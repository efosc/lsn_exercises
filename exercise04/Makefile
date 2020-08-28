CC = g++
CFLAGS = -Wall -O3 --std=c++11


Mol_functions.o: Mol_functions.cpp Mol_functions.h
	$(CC) -c -o $@ $< $(CFLAGS)

%.exe : %.cpp 
	$(CC) $^ -o $@ $(CFLAGS)

cancella:
	rm config.0
	rm config.final
	rm old.0
	rm old.final

fcc0:
	cp config.fcc config.0 