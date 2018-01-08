/*
	Author: Dosaru Daniel-Florin 331CA
	APD, Tema 3, 2017-2018

	Compile:
		mpic++ src.cpp -o filter -Wall
	Run:
		mpirun -np 12 filtru topologie1.in imagini.in statistica.out
		or
		mpirun -np 29 filtru topologie2.in imagini.in statistica.out
*/

#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>

#define LINE_MAX_LENGTH 1024
#define NR_TEST_TAG 0
#define ERROR_TAG 999
#define IMG_INFO_TAG 1001
#define SOBEL_TAG 1002
#define MEANR_TAG 1003
#define RESULT_TAG 1004
#define TERMINATION_TAG 1005

using namespace std;

// Allocate continously memory for images
int** my_alloc(int m, int n) {
	int *continue_space = (int *) malloc (m*n*sizeof(int));
	if (continue_space == NULL) {
		fprintf(stderr, "Out of memory at continue_space");
	}
	int **my_space = (int **) calloc (m, sizeof(int *));
	if (my_space == NULL) {
		fprintf(stderr, "Out of memory at my_space");
	}
	for (int i = 0; i < m; ++i) {
		my_space[i] = &(continue_space[n*i]);
	}
	return my_space;
}

int main(int argc, char * argv[]) {
	int rank;
	int nProcesses;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

	int * parent, *topo;
	parent = (int *) malloc(nProcesses*sizeof(int));
	if (parent == NULL) {
		fprintf(stderr, "Out of memory at parent");
	}

	topo = (int *)malloc(nProcesses * sizeof(int));
	if (topo == NULL) {
		fprintf(stderr, "Out of memory at topo");
	}

	// stat[i] = number of lines that had been processed by node i
	int * stat = (int *)malloc(nProcesses * sizeof(int));
	if (stat == NULL) {
		fprintf(stderr, "Out of memory at topo");
	}

	// initial values for topology, parent and stat array
	for (int i = 0; i < nProcesses; i++){
		topo[i] = 0;
		stat[i] = 0;
		parent[i] = -1;
	}

	const char sep[2] = " ";
	char *token;

	if (rank == 0) { // only the leader of the network executes this:
		//reading only the rank th line from topology file
		FILE *f_topo = fopen(argv[1], "r");
		int count = 0;
		char line[LINE_MAX_LENGTH];

		while (fgets(line, sizeof(line), f_topo) != NULL) {
			if (count == rank) {
				token = strtok(line, sep); // I already know that this is "rank: "  
				token = strtok(NULL, sep);
				while( token != NULL ) {
					int nr = atoi(token);
					topo[nr] = 1;
					token = strtok(NULL, sep);
				}

				break;
			} else {
				count++;
			}
		}

		fclose(f_topo);

		vector<int> kids;
		for (int i = 0; i < nProcesses; ++i) {
			if (topo[i] == 1) {
				kids.push_back(i);
			}
		}

		ifstream infile(argv[2]);
		ofstream outfile(argv[3]);
		int nr_tests = 0;
		infile >> nr_tests;
		
		// Sending the number of images that will be processed by the network
		for (int i = 0; i < (int)kids.size(); ++i) {
			MPI_Send(&nr_tests, 1, MPI_INT, kids[i], NR_TEST_TAG, MPI_COMM_WORLD);
		}

		for (int k = 0; k < nr_tests; k++) {
			printf("Photo number %d:\n", k+1);
			string type, img_in, img_out, inputLine = "";

			int i = 0, j = 0, H = 0, W = 0, maxValue, info_img[3];
			stringstream ss;
			infile >> type >> img_in >> img_out;
			ifstream in_f_img(img_in.c_str());
			FILE * out_f_img;

			out_f_img = fopen (img_out.c_str(),"w");

			getline(in_f_img, inputLine);
			if(inputLine.compare("P2") != 0)
				cerr << "Version error, this is not a image exported with GIMP" << endl;

			fprintf(out_f_img, "%s\n", inputLine.c_str());
			getline(in_f_img, inputLine); // comment
			fprintf(out_f_img, "%s\n", inputLine.c_str());
			ss << in_f_img.rdbuf();
			ss >> W >> H;
			fprintf(out_f_img, "%d %d\n", W, H);

			ss >> maxValue; // pure white of a pixel
			fprintf(out_f_img, "%d\n", maxValue);

			info_img[1] = W + 2; // number of columns of a line + adding 2 zeros
			info_img[2] = maxValue;

			for (int i = 0; i < (int)kids.size(); ++i) {
				if (i == (int)kids.size() - 1) {
					info_img[0] = H / kids.size() + 2 + H % kids.size(); // last kid
				} else {
					info_img[0] = H / kids.size() + 2;
				}
				MPI_Send(info_img, 3, MPI_INT, kids[i], IMG_INFO_TAG, MPI_COMM_WORLD);
			}

			int** m = my_alloc(H+2, W+2);

			//read pixel values
			printf("Reading image's number %d pixels... ", k+1);
			for(i = 1; i < H+1; ++i) {
				for (j = 1; j < W+1; ++j) {
					ss >> m[i][j];  
				}
			}
			printf("Done.\n");
			in_f_img.close();
			
			printf("Processing image number %d......... ", k+1);
			int tag = ERROR_TAG;
			if (type.compare("sobel") == 0) {
				tag = SOBEL_TAG;
			} else if (type.compare("mean_removal") == 0) {
				tag = MEANR_TAG;
			} else {
				fprintf(stderr, "Eu, %d, nu recunosc acest tip de filtru. ", rank);
			}

			if (kids.size() == 1) { // I have only one kid, I will send him the entire image
				MPI_Send(m[0], (H + 2) * (W + 2), MPI_INT, kids[0], tag, MPI_COMM_WORLD);
			} else {
				for (int i = 0; i < (int)kids.size(); i++) {
					if (i == (int)kids.size() - 1) { // last kid
						MPI_Send(m[i * ((int)H/kids.size())], 
								((int)H / kids.size() + 2 + H % kids.size()) * (W + 2), 
								MPI_INT, kids[i], tag, MPI_COMM_WORLD);
					} else {// other kids
						MPI_Send(m[i * ((int)H/kids.size())],
							((int)H / kids.size() + 2) * (W + 2), MPI_INT, kids[i], tag, MPI_COMM_WORLD);
					}
				}
			}

			free(m[0]);
			free(m);

			int** mNew = my_alloc(H, W);

			// Receiving and join the result(s) from my kids
			for (int i = 0; i < (int)kids.size(); i++) {
				if (i == (int)kids.size() - 1) { // last kid
					MPI_Recv(mNew[i * ((int)H/kids.size())], (((int)H / kids.size()) + H % kids.size())* W, 
						MPI_INT, kids[i], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				} else {
					MPI_Recv(mNew[i * ((int)H/kids.size())], ((int)H / kids.size()) * W, MPI_INT, 
						kids[i], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}
			printf("Done.\n");// Done processing message

			printf("Writing results into a new file...");
			for(i = 0; i < H; ++i) {
				for (j = 0; j < W; ++j) {
					fprintf(out_f_img, "%d\n", mNew[i][j]);
				}
			}
			printf(" Done.\n");
			
			free(mNew[0]);
			free(mNew);
			fclose(out_f_img);
		}

		infile.close();
		outfile.close();
		// Anounce that the work is done
		int am_terminat = 1;
		for (int i = 0; i < (int)kids.size(); ++i) {
			MPI_Send(&am_terminat, 1, MPI_INT, kids[i], TERMINATION_TAG, MPI_COMM_WORLD);
		}

		int * statRcv = (int *) malloc(nProcesses* sizeof(int));
		for (int i = 0; i < (int)kids.size(); ++i) {
			MPI_Recv(statRcv, nProcesses, MPI_INT, 
						kids[i], TERMINATION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int i = 0; i < nProcesses; ++i) {
				if (statRcv[i] > 0) {
					stat[i] = statRcv[i];
				}
			}
		}
		free(statRcv);
		FILE * f_stat = fopen(argv[3], "w");
		// write in statistica.out:
		for (int i = 0; i < nProcesses; ++i) {
			fprintf(f_stat, "%d: %d\n", i, stat[i]); 
		}
		fclose(f_stat);

	} else {
		int nr_neigh = 0;
		FILE *f_topo = fopen(argv[1], "r");
		int count = 0;
		char line[LINE_MAX_LENGTH];
		while (fgets(line, sizeof(line), f_topo) != NULL) {
			if (count == rank) {
				token = strtok(line, sep); // "rank: "  
				token = strtok(NULL, sep);

				while( token != NULL ) {
					int nr = atoi(token);
					topo[nr] = 1;
					nr_neigh++;
					token = strtok(NULL, sep);
				}
				break;
			} else {
				count++;
			}
		}
		fclose(f_topo);

		// primirea numarul de imagini care vor fi procesate
		int nr_tests = 0;
		MPI_Status my_status;
		MPI_Recv(&nr_tests, 1, MPI_INT, MPI_ANY_SOURCE, NR_TEST_TAG, MPI_COMM_WORLD, &my_status);
		if (parent[rank] == -1) {
			parent[rank] = my_status.MPI_SOURCE;
		} else {
			fprintf(stderr, "Me, %d th process, I already had a parent node \
					when I got this message from you, %d", rank, my_status.MPI_SOURCE);
		}

		vector<int> kids;
		for (int i = 0; i < nProcesses; ++i) {
			if (topo[i] == 1 && parent[rank] != i) {
				kids.push_back(i);
			}
		}

		for (int i = 0; i < (int)kids.size(); i++) {
			MPI_Send(&nr_tests, 1, MPI_INT, kids[i], NR_TEST_TAG, MPI_COMM_WORLD);
		}

		for (int k = 0; k < nr_tests; k++) {
			int info_img[3];
			MPI_Recv(info_img, 3, MPI_INT, parent[rank], IMG_INFO_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			int H = info_img[0] - 2;
			int W = info_img[1] - 2;
			int maxValue = info_img[2];

			if (kids.size() > 0) {
				if (kids.size() == 1) {
					info_img[0] = (int)H / kids.size() + 2;
					MPI_Send(info_img, 3, MPI_INT, kids[0], IMG_INFO_TAG, MPI_COMM_WORLD);
				} else {
					for (int i = 0; i < (int)kids.size(); i++) {
						if (i < (int)kids.size() - 1) {
							info_img[0] = (int)H / kids.size() + 2;
						} else {
							info_img[0] = (int) H / kids.size() + 2 + H % kids.size();
						}
						MPI_Send(info_img, 3, MPI_INT, kids[i], IMG_INFO_TAG, MPI_COMM_WORLD);
					}
				}
			}

			int** m = my_alloc(H+2, W+2);
			if (m == NULL) {
				fprintf(stderr, "Out of memory!\n");
			}
			int** mNew = my_alloc(H, W);
			if (mNew == NULL) {
				fprintf(stderr, "Out of memory!\n");
			}

			// receive pixels
			MPI_Status stat_filter;
			MPI_Recv(m[0], (H+2) * (W+2), 
				MPI_INT, parent[rank], MPI_ANY_TAG, MPI_COMM_WORLD, &stat_filter);
			int tag = stat_filter.MPI_TAG;

			if (kids.size() > 0) {
				if (kids.size() == 1) { // I have only one kid, I will send him the entire image
					MPI_Send(m[0], (H+2) * (W+2), 
						MPI_INT, kids[0], tag, MPI_COMM_WORLD);
					// Receiving the processed image
					MPI_Recv(mNew[0], H * W, 
						MPI_INT, kids[0], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Send(mNew[0], H * W, 
						MPI_INT, parent[rank], RESULT_TAG, MPI_COMM_WORLD);
				} else { // this node (with more than one child) executes this:
					for (int i = 0; i < (int)kids.size(); i++) {
						if (i == (int)kids.size() - 1) { // last child
							MPI_Send(m[i * ((int)H/kids.size())], (((int)H / kids.size()) + 2 + H % kids.size()) * (W+2), 
									MPI_INT, kids[i], tag, MPI_COMM_WORLD);
						} else {
							MPI_Send(m[i * ((int)H/kids.size())], (((int)H / kids.size()) + 2) * (W+2),
								MPI_INT, kids[i], tag, MPI_COMM_WORLD);
						}
					}

					// Receive and join the result(s) from my kids
					for (int i = 0; i <  (int)kids.size(); i++) {
						if (i ==  (int)kids.size() - 1) { // last child
							MPI_Recv(mNew[i * ((int)H/kids.size())], (((int)H / kids.size()) + H % kids.size()) * W, 
								MPI_INT, kids[i], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						} else {
							MPI_Recv(mNew[i * ((int)H/kids.size())], ((int)H / kids.size()) * W, 
								MPI_INT, kids[i], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
					}

					// Send the result to my parent
					MPI_Send(mNew[0], H * W, MPI_INT, parent[rank], RESULT_TAG, MPI_COMM_WORLD);
				}
			} else { // This node is a leaf, applying the filter:
				if (tag == SOBEL_TAG) {
					for(int i = 1; i < H+1; ++i) {
						for (int j = 1; j < W+1; ++j) {
							mNew[i-1][j-1] = m[i-1][j-1] - m[i-1][j+1] +
								2 * m[i][j-1] - 2 * m[i][j+1] +
								m[i+1][j-1] - m[i+1][j+1] + 127;
							if (mNew[i-1][j-1] < 0) {
								mNew[i-1][j-1] = 0;
							}

							if (mNew[i-1][j-1] > maxValue) {
								mNew[i-1][j-1] = maxValue;
							}
						}
					}
				} else if(tag == MEANR_TAG) {
					for(int i = 1; i < H+1; ++i) {
						for (int j = 1; j < W+1; ++j) {
							mNew[i-1][j-1] = -m[i-1][j-1] - m[i-1][j]- m[i-1][j+1] +
								- m[i][j-1]   + 9 * m[i][j]   - m[i][j+1] +
								- m[i+1][j-1] - m[i+1][j] - m[i+1][j+1];
							if (mNew[i-1][j-1] < 0) {
								mNew[i-1][j-1] = 0;
							}
							
							if (mNew[i-1][j-1] > maxValue) {
								mNew[i-1][j-1] = maxValue;
							}
						}
					}
				} else {
					fprintf(stderr, "Error, unrecognised filter tag (%d) received by %d", tag, rank);
				}
				stat[rank] += H;

				// Send the result to my parent
				MPI_Send(mNew[0], H * W, 
					MPI_INT, parent[rank], RESULT_TAG, MPI_COMM_WORLD);
			}
			free(mNew[0]);
			free(mNew);
			free(m[0]);
			free(m);
		}

		int finished = 0;
		MPI_Recv(&finished, 1, MPI_INT, parent[rank], TERMINATION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if (kids.size() == 0){ // leaf node
			MPI_Send(stat, nProcesses, MPI_INT, parent[rank], TERMINATION_TAG, MPI_COMM_WORLD);
		} else {
			int am_terminat = 1;
			for (int i = 0; i < (int)kids.size(); ++i) {
				MPI_Send(&am_terminat, 1, MPI_INT, kids[i], TERMINATION_TAG, MPI_COMM_WORLD);
			}
			int * statRcv = (int *) malloc(nProcesses* sizeof(int));
			for (int i = 0; i < (int)kids.size(); ++i) {
				MPI_Recv(statRcv, nProcesses, MPI_INT, 
							kids[i], TERMINATION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			printf("Done.\n");// Done processing message

			printf("Writing results into a new file...");
			for(i = 0; i < H; ++i) {
				for (j = 0; j < W; ++j) {
					fprintf(out_f_img, "%d\n", mNew[i][j]);
				}
			}
			printf(" Done.\n");
			
			free(mNew[0]);
			free(mNew);
			fclose(out_f_img);
		}

		infile.close();
		outfile.close();
		// Anounce that the work is done
		int am_terminat = 1;
		for (int i = 0; i < (int)kids.size(); ++i) {
			MPI_Send(&am_terminat, 1, MPI_INT, kids[i], TERMINATION_TAG, MPI_COMM_WORLD);
		}

		int * statRcv = (int *) malloc(nProcesses* sizeof(int));
		for (int i = 0; i < (int)kids.size(); ++i) {
			MPI_Recv(statRcv, nProcesses, MPI_INT, 
						kids[i], TERMINATION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int i = 0; i < nProcesses; ++i) {
				if (statRcv[i] > 0) {
					stat[i] = statRcv[i];
				}
			}
		}
		free(statRcv);
		FILE * f_stat = fopen(argv[3], "w");
		// write in statistica.out:
		for (int i = 0; i < nProcesses; ++i) {
			fprintf(f_stat, "%d: %d\n", i, stat[i]); 
		}
		fclose(f_stat);

	} else {
		int nr_neigh = 0;
				for (int i = 0; i < nProcesses; ++i) {
					if (statRcv[i] > 0) {
						stat[i] = statRcv[i];
					}
				}
			}
			free(statRcv);
			// Send the statistics to my parent
			MPI_Send(stat, nProcesses, MPI_INT, parent[rank], TERMINATION_TAG, MPI_COMM_WORLD);
		}
	}
	free(stat);
	free(topo);
	free(parent);

	printf("Bye from %i/%i\n", rank, nProcesses);
	MPI_Finalize();
	return 0;
}
