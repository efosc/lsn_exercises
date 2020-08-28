CC = g++
CFLAGS = -Wall -O3 --std=c++11


Monte_Carlo_NVT_new.exe: Monte_Carlo_NVT_new.cpp random.o
	$(CC) $^ -o $@ $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

%.exe : %.cpp 
	$(CC) $^ -o $@ $(CFLAGS)

cancella:
	rm config.0
	rm config.final
	rm old.0
	rm old.final

fcc0:
	cp config.fcc config.0 