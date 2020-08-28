CC = g++
CFLAGS = -Wall -O3 --std=c++11

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
posizione.o : posizione.cpp posizione.h random.h
	$(CC) -c posizione.cpp -o posizione.o $(CFLAGS)

%.exe : %.cpp random.o
	$(CC) $^ -o $@ $(CFLAGS)

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp main.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

medie.exe: random.o
	$(CC) -o $@ $^ $(CFLAGS)
 
