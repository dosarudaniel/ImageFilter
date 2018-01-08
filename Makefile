all:build

build:src.cpp
	mpic++ src.cpp -o filtru -Wall
clean:
	rm filtru