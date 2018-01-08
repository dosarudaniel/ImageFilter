/*
	Author: Dosaru Daniel-Florin 331CA
	APD, Tema 3, 2017-2018

	Compile:
		mpic++ src.cpp -o filter -Wall
	Run:
		mpirun -np 12 filter topologie1.in imagini.in statistica.out
		or
		mpirun -np 29 filter topologie2.in imagini.in statistica.out

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

	printf("%d ", 6*333/7);
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

	// stat[i] = cate linii a procesat nodul i
	int * stat = (int *)malloc(nProcesses * sizeof(int));
	if (stat == NULL) {
		fprintf(stderr, "Out of memory at topo");
	}

	// initial values for topologi array and parent array
	for (int i = 0; i < nProcesses; i++){
		topo[i] = 0;
		stat[i] = 0;
		parent[i] = -1;
	}

	const char sep[2] = " ";
	char *token;

	if (rank == 0) {
		// citesc prima linie din topologie.in <- copii nodului radacina.
		FILE *f_topo = fopen(argv[1], "r");
		int count = 0;
		char line[LINE_MAX_LENGTH];

		while (fgets(line, sizeof(line), f_topo) != NULL) {
			if (count == rank) {
				token = strtok(line, sep); // "rank: "  
				token = strtok(NULL, sep);
				/* walk through other tokens */
				while( token != NULL ) {
					int nr = atoi(token);
					topo[nr] = 1;
					token = strtok(NULL, sep);
					//printf( "%s_", token );
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

		for (int i = 0; i < (int)kids.size(); ++i) {
			MPI_Send(&nr_tests, 1, MPI_INT, kids[i], NR_TEST_TAG, MPI_COMM_WORLD);
		}

		for (int k = 0; k < nr_tests; k++) {
			string type, img_in, img_out, inputLine = "";

			int i = 0, j = 0, H = 0, W = 0, maxValue, info_img[3];
			stringstream ss;
			printf("aici1\n");
			infile >> type >> img_in >> img_out;
			ifstream in_f_img(img_in.c_str());
			FILE * out_f_img;

			out_f_img = fopen (img_out.c_str(),"w");

			getline(in_f_img, inputLine);
			if(inputLine.compare("P2") != 0)
				cerr << "Version error" << endl;

			printf("aici2\n");
			fprintf(out_f_img, "%s\n", inputLine.c_str());

			getline(in_f_img, inputLine); // comment
			fprintf(out_f_img, "%s\n", inputLine.c_str());

			ss << in_f_img.rdbuf();
			ss >> W >> H;
			fprintf(out_f_img, "%d %d\n", W, H);

			printf("aici3\n");
			ss >> maxValue;
			fprintf(out_f_img, "%d\n", maxValue);

			info_img[1] = W + 2; // nr de coloane ale unei linii + bordarea cu 0
			info_img[2] = maxValue;

			printf("aici4\n");
			for (int i = 0; i < (int)kids.size(); ++i) {
				if (i == (int)kids.size() - 1) {
					info_img[0] = H / kids.size() + 2 + H % kids.size(); // last kid
				} else {
					info_img[0] = H / kids.size() + 2;
				}
				MPI_Send(info_img, 3, MPI_INT, kids[i], IMG_INFO_TAG, MPI_COMM_WORLD);
			}

			printf("aici5\n");
			int** m = my_alloc(H+2, W+2);

			printf("aici6\n");
			//read pixel values
			for(i = 1; i < H+1; ++i) {
				for (j = 1; j < W+1; ++j) {
					ss >> m[i][j];  
				}
			}
			in_f_img.close();

			printf("aici7\n");
			int tag = ERROR_TAG;
			if (type.compare("sobel") == 0) {
				tag = SOBEL_TAG;
			} else if (type.compare("mean_removal") == 0) {
				tag = MEANR_TAG;
			} else {
				fprintf(stderr, "Eu, %d, nu recunosc acest tip de filtru. ", rank);
			}

			printf("aici8\n");
			if (kids.size() == 1) { // am un singur copil, ii trimit toată imaginea
				MPI_Send(m[0], (H + 2) * (W + 2), MPI_INT, kids[0], tag, MPI_COMM_WORLD);
			} else {
				for (int i = 0; i < (int)kids.size(); i++) {
					if (i == 0) { // primul copil
						MPI_Send(m[0], ((int)H / kids.size() + 2) * (W + 2), MPI_INT, kids[0], tag, MPI_COMM_WORLD);

					} else if (i == (int)kids.size() - 1) { // ultimul copil
						MPI_Send(m[i * ((int)H/kids.size())], 
								((int)H / kids.size() + 2 + H % kids.size()) * (W + 2), 
								MPI_INT, kids[i], tag, MPI_COMM_WORLD);
					} else {// ceilalti copii
						MPI_Send(m[i * ((int)H/kids.size())],
							((int)H / kids.size() + 2) * (W + 2), MPI_INT, kids[i], tag, MPI_COMM_WORLD);
					}
				}
			}

			printf("aici9\n");
			free(m[0]);
			free(m);

			int** mNew = my_alloc(H, W);

			printf("aici10\n");
			// primesc rezultatele de la copii si trebuie sa le unesc.
			for (int i = 0; i < (int)kids.size(); i++) {
				if (i == 0) { // primul copil
					MPI_Recv(mNew[0], ((int)H / kids.size()) * W, MPI_INT, 
						kids[i], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				} else if (i == (int)kids.size() - 1) { // ultimul copil
					MPI_Recv(mNew[i * ((int)H/kids.size())], (((int)H / kids.size()) + H % kids.size())* W, 
						MPI_INT, kids[i], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				} else {
					MPI_Recv(mNew[i * ((int)H/kids.size())], ((int)H / kids.size()) * W, MPI_INT, 
						kids[i], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}

			// apoi le scriu in fisier:
			printf("Acum scriu imaginea in fisier\n");
			for(i = 0; i < H; ++i) {
				for (j = 0; j < W; ++j) {
					fprintf(out_f_img, "%d\n", mNew[i][j]);
				}
			}

			free(mNew[0]);
			free(mNew);
			fclose(out_f_img);

			printf("next photo: %d\n", k+1);
		}

		infile.close();
		outfile.close();
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
					printf("%d: %d linii\n", i, stat[i]); 
				}
			}
		}
		free(statRcv);
		FILE * f_stat = fopen(argv[3], "w");
		// scrie in statistica.out datele:
		for (int i = 0; i < nProcesses; ++i) {
			fprintf(f_stat, "%d: %d\n", i, stat[i]); 
		}
		fclose(f_stat);

		// trimite msg cu tag de terminare.
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
			// printf("Eu , %d, am parintele %d si sunt pregatit sa procesez %d imagini\n", rank, parent[rank], nr_tests);
		} else {
			fprintf(stderr, "Eu, %d, aveam deja parinte când am primit de la %d.", rank, my_status.MPI_SOURCE);
		}

		vector<int> kids;
		for (int i = 0; i < nProcesses; ++i) {
			if (topo[i] == 1 && parent[rank] != i) { // trimit doar la copii.
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

			printf("Eu, %d, am %d copii şi voi primi, %d linii\n",rank, (int) kids.size(), H);

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

			//printf("%d, Am primit %d linii\n", rank, H+2);

			if (kids.size() > 0) {
				if (kids.size() == 1) { // am un singur copil, ii trimit toată imaginea
					MPI_Send(m[0], (H+2) * (W+2), 
						MPI_INT, kids[0], tag, MPI_COMM_WORLD);
					// primesc imaginea procesata
					printf("Eu, %d, primesc imaginea procesata de la singurul meu copil..\n", rank);
					MPI_Recv(mNew[0], H * W, 
						MPI_INT, kids[0], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					printf("Eu, %d, am primit imaginea procesata de la copil. trimit la parinte...\n", rank);
					MPI_Send(mNew[0], H * W, 
						MPI_INT, parent[rank], RESULT_TAG, MPI_COMM_WORLD);
				} else { // nod intermediar cu mai multi copii
					for (int i = 0; i < (int)kids.size(); i++) {
						if (i == (int)kids.size() - 1) { // ultimul copil

							MPI_Send(m[i * ((int)H/kids.size())], (((int)H / kids.size()) + 2 + H % kids.size()) * (W+2), 
									MPI_INT, kids[i], tag, MPI_COMM_WORLD);
						} else {
							MPI_Send(m[i * ((int)H/kids.size())], (((int)H / kids.size()) + 2) * (W+2),
								MPI_INT, kids[i], tag, MPI_COMM_WORLD);
						}
						printf("Eu, %d, trimit linii la %d, al %d lea copil\n",  rank, kids[i], i);
					}
					printf("Eu, %d, am impartit treaba la fii mei. \n", rank);

					int nr_elem = 0;

					// primesc rezultatele de la copii si trebuie sa le unesc.
					for (int i = 0; i <  (int)kids.size(); i++) {
						if (i ==  (int)kids.size() - 1) { // ultimul copil
							MPI_Recv(mNew[i * ((int)H/kids.size())], (((int)H / kids.size()) + H % kids.size()) * W, 
								MPI_INT, kids[i], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								nr_elem += (H / kids.size() + H % kids.size()) * W;
						} else {
							MPI_Recv(mNew[i * ((int)H/kids.size())], ((int)H / kids.size()) * W, 
								MPI_INT, kids[i], RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							nr_elem += ((int)H / kids.size()) * W;
						}
					}

					printf("Eu, %d, am primit rezultatele copiilor: %d linii, trimit la parinte %d\n", rank, nr_elem, H*W);
					// si trimit rezultatul la părinte
					MPI_Send(mNew[0], H * W, MPI_INT, parent[rank], RESULT_TAG, MPI_COMM_WORLD);

					printf("Eu, %d, am trimis rezultatul la %d \n", rank, parent[rank]);
				} // nod intermediar// nod intermediar

			} else { // frunza, aplic filtrul
				printf("eu, %d, Mă pregătesc pentru a procesa liniile\n", rank);
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
					fprintf(stderr, "Error, am primit un tag necunoscut %d", tag);
				}
				stat[rank] += H;

				printf("Eu, %d, trimit rezultatul la %d \n", rank, parent[rank]);
				// apoi trimit inapoi rezultatul..
				MPI_Send(mNew[0], H * W, 
					MPI_INT, parent[rank], RESULT_TAG, MPI_COMM_WORLD);
				printf("Eu, %d, am trimis rezultatul\n", rank);
			}

			printf("Eu, %d, aici sunt 1\n", rank);
			free(mNew[0]);
			printf("Eu, %d, aici sunt 2\n", rank);
			free(mNew);
			printf("Eu, %d, aici sunt 3\n", rank);
			free(m[0]); // try to delete this.
			printf("Eu, %d, aici sunt 4\n", rank);
			free(m);
			printf("Eu, %d, aici sunt 5\n", rank);
		}

		int am_terminat = 0;
		MPI_Recv(&am_terminat, 1, MPI_INT, parent[rank], TERMINATION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("Eu, %d, Am primit mesaj cu tag de terminare: %d\n", rank, am_terminat);

		if (kids.size() == 0){ // nod frunză {
			MPI_Send(stat, nProcesses, MPI_INT, parent[rank], TERMINATION_TAG, MPI_COMM_WORLD);
		} else { // nod intermediar
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
					//	printf("De la %d: %d: %d linii\n", rank, i, stat[i]); 
					}
				}
			}
			free(statRcv);
			// trimite statistica mai departe la parintele meu
			MPI_Send(stat, nProcesses, MPI_INT, parent[rank], TERMINATION_TAG, MPI_COMM_WORLD);
		}
	
	}


	free(stat);
	printf("Eu, %d, before free of topo\n", rank);
	free(topo);
	free(parent);

	printf("Bye from %i/%i\n", rank, nProcesses);
	MPI_Finalize();
	return 0;
}