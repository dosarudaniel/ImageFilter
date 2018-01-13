# ImageFilter
Image processing using MPI networks

Before compiling the source file make sure you have installed libopenmpi-dev,
openmpi-bin, openmpi-doc, openmpi-common packages:
For linux, ubuntu execute:
	`sudo apt-get install libopenmpi-dev openmpi-bin openmpi-doc openmpi-common`


### Compile:
	mpic++ src.cpp -o filter -Wall
### Run:
	mpirun -np 12 filter topology1.in images.in statistics.out
	or
	mpirun -np 29 filter topology2.in images.in statistics.out

