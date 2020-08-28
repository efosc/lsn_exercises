CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out

main2.exe : main2.o random.o posizione.o
	$(CC) random.o main2.o posizione.o -o main2.exe
main2.o : main2.cpp
	$(CC) -c main2.cpp -o main2.o $(CFLAGS)

