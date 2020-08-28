CC = g++
CFLAGS = -Wall -O3 --std=c++11


random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out

%.exe : %.cpp random.o
	$(CC) $^ -o $@ $(CFLAGS)

