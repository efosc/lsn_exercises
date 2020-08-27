CC = g++
CFLAGS = -Wall -O3 --std=c++11


random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

funzioni9.o : funzioni9.cpp funzioni9.h
	$(CC) -c funzioni9.cpp -o funzioni9.o $(CFLAGS)

%.exe : %.o random.o funzioni9.o
	$(CC) $^ -o $@ $(CFLAGS)

%.o : %.cpp
	$(CC) -c $^ -o $@ $(CFLAGS)



