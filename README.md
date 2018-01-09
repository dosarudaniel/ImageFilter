# ImageFilter
Image processing using MPI networks

Before compiling the source file make sure you have installed libopenmpi-dev,
openmpi-bin, openmpi-doc, openmpi-common packages:
For linux, ubuntu execute:
	sudo apt-get install libopenmpi-dev openmpi-bin openmpi-doc openmpi-common


#Compile:
	mpic++ src.cpp -o filtru -Wall
#Run:
	mpirun -np 12 filtru topologie1.in imagini.in statistica.out
	or
	mpirun -np 29 filtru topologie2.in imagini.in statistica.out

