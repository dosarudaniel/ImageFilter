all:build

build:src.cpp
	mpic++ src.cpp -o filter -Wall
clean:
	rm filter