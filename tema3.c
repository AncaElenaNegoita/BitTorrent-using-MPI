#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stddef.h>
#include "../../../../../usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h"

#define TRACKER_RANK 0
#define MAX_FILES 10
#define MAX_FILENAME 15
#define HASH_SIZE 33
#define MAX_CHUNKS 100

int numFiles = 0, clients = 0, c = 0;
int hasLeft = 0;
char segments[MAX_FILES][MAX_CHUNKS][HASH_SIZE];
int numSeg[MAX_FILES];
int hasFile[MAX_FILES][MAX_FILES] = {0};
int needsFiles[MAX_FILES] = {0};

typedef struct {
	int user;
	int file;
	char hashes[MAX_CHUNKS][HASH_SIZE];
	int sizeFile;
}torrentInfo;

void *download_thread_func(void *arg)
{
	int rank = *(int*) arg;
	int numberSegFile = 0, curr = 0;
	
	char hash[MAX_CHUNKS][HASH_SIZE];
	
	while (1) {
		// Send to the tracker the files it needs
		MPI_Send(needsFiles, MAX_FILES, MPI_INT, TRACKER_RANK, 5, MPI_COMM_WORLD);
		if (hasLeft == 0) {
			return NULL;
		}
		torrentInfo t[MAX_CHUNKS];
		
		int numEl, n = 0, m = 0;

		MPI_Datatype customtype, oldtypes[4];
		int blockcounts[4];
		MPI_Aint offsets[4];

		offsets[0] = offsetof(torrentInfo, user);
		oldtypes[0] = MPI_INT;
		blockcounts[0] = 1;

		offsets[1] = offsetof(torrentInfo, file);
		oldtypes[1] = MPI_INT;
		blockcounts[1] = 1;

		offsets[2] = offsetof(torrentInfo, hashes);
		oldtypes[2] = MPI_CHAR;
		blockcounts[2] = HASH_SIZE * MAX_CHUNKS;

		offsets[3] = offsetof(torrentInfo, sizeFile);
		oldtypes[3] = MPI_INT;
		blockcounts[3] = 1;

		MPI_Type_create_struct(4, blockcounts, offsets, oldtypes, &customtype);
		MPI_Type_commit(&customtype);

		// The number of elements in the array above
		MPI_Recv(&numEl, 1, MPI_INT, TRACKER_RANK, 1, MPI_COMM_WORLD, NULL);

		// Receive the list of users, files and their owned segments
		MPI_Recv(t, numEl, customtype, TRACKER_RANK, 0, MPI_COMM_WORLD, NULL);

		// The current wanted file to download
		int file = needsFiles[curr];
		// How many users have segments of this file
		for (n = 0; n < numEl; n++) {
			if (t[n].file != file) {
				break;
			}
		}

		int i = 0;
		while (i < 10) {
			// The remaining segments the file needs to download
			int neededSeg = t[m].sizeFile - numberSegFile;
			
			int minSeg;
			// Verify if the needed segmnents are less than 10 - i 
			if (neededSeg < 10 - i) {
				minSeg = neededSeg;
			} else {
				minSeg = 10 - i;
			}

			// If there are more clients having the wanted file than segments needed, take 1
			// segment from each one of them
			if (minSeg <= n - m) {
				for (int k = 0; k < minSeg; k++) {
					// Send to the seed/peer the hash and wait for confirmation
					char s[MAX_FILES];
					MPI_Send(t[k + m].hashes[numberSegFile], HASH_SIZE, MPI_CHAR, t[k + m].user, 6, MPI_COMM_WORLD);
					// Receive the confirmation for the hash
					MPI_Recv(s, MAX_FILES, MPI_CHAR, t[k + m].user, 3, MPI_COMM_WORLD, NULL);
					if (strcmp(s, "ok") == 0) {
						numSeg[file]++;
						strcpy(hash[numberSegFile], t[k + m].hashes[numberSegFile]);
						numberSegFile++;
					}
					i++;
				}
			} else {
				int j = m;
				while (n != j && i < 10) {
					int userSeg = minSeg / (n - m);
					// Store the hashes from each client, in order
					for (int k = 0; k < userSeg; k++) {
						char st[MAX_FILES];
						// Ask the client for confirmation of the hash
						MPI_Send(t[j].hashes[numberSegFile], HASH_SIZE, MPI_CHAR, t[j].user, 6, MPI_COMM_WORLD);
						// Receive the confirmation and add the segment
						MPI_Recv(st, MAX_FILES, MPI_CHAR, t[j].user, 3, MPI_COMM_WORLD, NULL);
						if (strcmp(st, "ok") == 0) {
							numSeg[file]++;
							strcpy(hash[numberSegFile], t[j].hashes[numberSegFile]);
							numberSegFile++;
						}
						i++;
					}
					j++;
				}
			}
			
			neededSeg = t[m].sizeFile - numberSegFile;
			// If the file completed downloading, print the result
			if (neededSeg == 0) {
				// Open the file where the hashes should be printed
				FILE* f;
				char outputFile[HASH_SIZE];
				sprintf(outputFile, "client%d_file%d", rank, file);
				f = fopen(outputFile, "w");
				if (f == NULL) {
					printf("Could not open file.\n");
					exit(-1);
				}
				
				// Print the file
				for (int j = 0; j < t[m].sizeFile; j++) {
					fprintf(f, "%s\n", hash[j]);
				}

				fclose(f);
				// Reset the variables and decrease the number of wanted files
				needsFiles[curr] = 0;
				numberSegFile = 0;
				curr++;
				hasLeft--;

				// If there aren't any files left for downloading, tell the tracker
				if (hasLeft == 0) {
					MPI_Type_free(&customtype);
					char s[MAX_CHUNKS] = "gata";
					MPI_Send(s, MAX_CHUNKS, MPI_CHAR, TRACKER_RANK, 2, MPI_COMM_WORLD);
					// Send all the remaining segments
					MPI_Send(numSeg, MAX_FILES, MPI_INT, TRACKER_RANK, 4, MPI_COMM_WORLD);
					return NULL;
				}
				// Shift to the next file needed
				m = n;
				file = needsFiles[curr];
				for (n = m + 1; n < numEl; n++) {
					if (t[n].file != file) {
						break;
					}
				}
				
			}
		}
		MPI_Type_free(&customtype);
		// Send to the tracker the new updated list of segments the client has
		char s[MAX_FILES] = "continue";
		MPI_Send(s, MAX_CHUNKS, MPI_CHAR, TRACKER_RANK, 2, MPI_COMM_WORLD);
		MPI_Send(numSeg, MAX_FILES, MPI_INT, TRACKER_RANK, 4, MPI_COMM_WORLD);
	}
	return NULL;
}

void *upload_thread_func(void *arg)
{
	int rank = *(int*) arg;
	while (1) {
		char hash[HASH_SIZE];
		MPI_Status status;
		// Receive the hash or the finished message
		MPI_Recv(hash, HASH_SIZE, MPI_CHAR, MPI_ANY_SOURCE, 6, MPI_COMM_WORLD, &status);
		int senderRank = status.MPI_SOURCE;
		if (senderRank == TRACKER_RANK) {
			return NULL;
		}
		
		// Confirm the received hash
		strcpy(hash, "ok");
		MPI_Send(hash, MAX_FILES, MPI_CHAR, senderRank, 3, MPI_COMM_WORLD);
	}
	return NULL;
}

void tracker(int numtasks, int rank) {
	c = numtasks;
	// Array that stores the user, the file it has and how many segments
	torrentInfo t[MAX_CHUNKS];
	// Array that stores the files a client needs
	int needed[MAX_FILES];
	// Array that stored the finished clients
	int finished[MAX_FILES];

	// Create the datatype to send the structure
	MPI_Datatype customtype, oldtypes[4];
	int blockcounts[4];
	MPI_Aint offsets[4];

	offsets[0] = offsetof(torrentInfo, user);
	oldtypes[0] = MPI_INT;
	blockcounts[0] = 1;

	offsets[1] = offsetof(torrentInfo, file);
	oldtypes[1] = MPI_INT;
	blockcounts[1] = 1;

	offsets[2] = offsetof(torrentInfo, hashes);
	oldtypes[2] = MPI_CHAR;
	blockcounts[2] = HASH_SIZE * MAX_CHUNKS;

	offsets[3] = offsetof(torrentInfo, sizeFile);
	oldtypes[3] = MPI_INT;
	blockcounts[3] = 1;

	MPI_Type_create_struct(4, blockcounts, offsets, oldtypes, &customtype);
	MPI_Type_commit(&customtype);

	while(1) {
		char s[MAX_CHUNKS];
		for (int us = 1; us < numtasks; us++) {
			if (finished[us] == 0) {
				int numEl = 0;
				// The request of a client to receive informations about the files it needs
				MPI_Recv(needed, MAX_FILES, MPI_INT, us, 5, MPI_COMM_WORLD, NULL);

				for (int i = 0; i < MAX_FILES; i++) {
					if (needed[i] != 0) {
						for (int j = 1; j < numtasks; j++) {
							// If the client has the file and it's not the requesting client, store its information
							if (hasFile[j][needed[i]] != 0 && j != us && hasFile[us][needed[i]] < hasFile[j][needed[i]]) {
								t[numEl].user = j;
								t[numEl].file = needed[i];
								t[numEl].sizeFile = numSeg[needed[i]];
								for (int l = 0; l < hasFile[j][needed[i]]; l++) {
									strcpy(t[numEl].hashes[l], segments[needed[i]][l]);
								}
								numEl++;
							}
						}
					}
				}
				// If the user initially didn't have any files, mark it as finished
				if (numEl == 0) {
					finished[us] = 1;
					c--;
				} else {
					// Send the needed data that the client neeeds in order to download a file
					MPI_Send(&numEl, 1, MPI_INT, us, 1, MPI_COMM_WORLD);
					MPI_Send(t, numEl, customtype, us, 0, MPI_COMM_WORLD);

					// Receive a message from the current client that he has received all the files it needed
					MPI_Recv(s, MAX_CHUNKS, MPI_CHAR, us, 2, MPI_COMM_WORLD, NULL);
					if (strcmp(s, "gata") == 0) {
						c--;
						finished[us] = 1;
						// Receive the segments that the curret client has
						int a[MAX_FILES];
						MPI_Recv(a, MAX_FILES, MPI_INT, us, 4, MPI_COMM_WORLD, NULL);
						for (int i = 1; i < MAX_FILES; i++) {
							hasFile[us][i] = a[i];
						}
						// The program finished
						if (c == 1) {
							char s[HASH_SIZE] = "gata";
							for (int i = 1; i < numtasks; i++) {
								MPI_Send(s, HASH_SIZE, MPI_CHAR, i, 6, MPI_COMM_WORLD);
							}
							MPI_Type_free(&customtype);
							return;
						}
					} else {
						// Receive the segments that the curret client has
						int a[MAX_FILES];
						MPI_Recv(a, MAX_FILES, MPI_INT, us, 4, MPI_COMM_WORLD, NULL);
						for (int i = 1; i < MAX_FILES; i++) {
							hasFile[us][i] = a[i];
						}
					}
				}
			}
		}
	}
}

void peer(int numtasks, int rank) {
	pthread_t download_thread;
	pthread_t upload_thread;
	void *status;
	int r;

	r = pthread_create(&download_thread, NULL, download_thread_func, (void *) &rank);
	if (r) {
		printf("Eroare la crearea thread-ului de download\n");
		exit(-1);
	}

	r = pthread_create(&upload_thread, NULL, upload_thread_func, (void *) &rank);
	if (r) {
		printf("Eroare la crearea thread-ului de upload\n");
		exit(-1);
	}

	r = pthread_join(download_thread, &status);
	if (r) {
		printf("Eroare la asteptarea thread-ului de download\n");
		exit(-1);
	}

	r = pthread_join(upload_thread, &status);
	if (r) {
		printf("Eroare la asteptarea thread-ului de upload\n");
		exit(-1);
	}
}
 
int main (int argc, char *argv[]) {
	int numtasks, rank;
 
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if (provided < MPI_THREAD_MULTIPLE) {
		fprintf(stderr, "MPI nu are suport pentru multi-threading\n");
		exit(-1);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	FILE *file;
	char inputFile[HASH_SIZE];
	char name[HASH_SIZE + 3];
	int number, n;

	if (rank == TRACKER_RANK) {
		// Receive all the hashes of the files from the clients
		for (int i = 1; i < numtasks; i++) {
			int nRecv, numRecv, kRecv;
			char nameRecv[HASH_SIZE + 3];
			MPI_Recv(&nRecv, 1, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
			for (int j = 0; j < nRecv; j++) {
				MPI_Recv(&numRecv, 1, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
				MPI_Recv(&kRecv, 1, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
				if (numSeg[kRecv] == 0) {
					numFiles++;
				}
				// Total number of segments a file has
				numSeg[kRecv] = numRecv;
				// Store the number of segments the client has of each file
				hasFile[i][kRecv] = numRecv;

				// Store the hashes, associating the row with the file
				for (int k = 0; k < numRecv; k++) {
					MPI_Recv(nameRecv, HASH_SIZE, MPI_CHAR, i, 0, MPI_COMM_WORLD, NULL);
					strcat(nameRecv, "\0");
					strcpy(segments[kRecv][k], nameRecv);
				}
			}
		}

		// Send to every process that they can start downloading
		char s[3] = "ok";
		for (int i = 1; i < numtasks; i++) {
			MPI_Send(s, 3, MPI_CHAR, i, 0, MPI_COMM_WORLD);
		}
		
		tracker(numtasks, rank);
	} else {
		// Open the file of the client
		sprintf(inputFile, "in%d.txt", rank);
		file = fopen(inputFile, "r");
		if (file == NULL) {
			printf("Could not open file.\n");
			exit(-1);
		}
		clients++;
		// Read each entry file to get the files a user has
		fscanf(file, "%d", &n);
		// Send the number of files the current client has
		MPI_Send(&n, 1, MPI_INT, TRACKER_RANK, 0, MPI_COMM_WORLD);
		for (int i = 0; i < n; i++) {
			int k = 0;
			fscanf(file, "%s %d", name, &number);

			// Send the number of segments from the current file
			MPI_Send(&number, 1, MPI_INT, TRACKER_RANK, 0, MPI_COMM_WORLD);
			while (!isdigit(name[k])) {
				k++;
			}
			// Convert the number of the file to int
			k =  atoi(&name[k]);
			// How many segments of the file the user has
			numSeg[k] = number;
			// Send the number of the file
			MPI_Send(&k, 1, MPI_INT, TRACKER_RANK, 0, MPI_COMM_WORLD);
			
			
			for (int j = 0; j < number; j++) {
				// Send all the owned hashes of the file to the tracker
				fscanf(file, "%s", name);
				MPI_Send(name, HASH_SIZE, MPI_CHAR, TRACKER_RANK, 0, MPI_COMM_WORLD);
			}
		}

		// Read the files a user needs
		fscanf(file, "%d", &hasLeft);
		for (int i = 0; i < hasLeft; i++) {
			fscanf(file, "%s", name);
			int k = 0;
			while (!isdigit(name[k])) {
				k++;
			}
			k =  atoi(&name[k]);
			needsFiles[i] =  k;
		}

		fclose(file);

		// Wait for the response of the tracker that it gathered all the information
		char s[3];
		MPI_Recv(s, 3, MPI_CHAR, TRACKER_RANK, 0, MPI_COMM_WORLD, NULL);

		// Start the program if the tracker received all the information it needed
		if (strcmp(s, "ok") == 0) {
			peer(numtasks, rank);
		}
	}

	// Wait for all processes to finish
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
