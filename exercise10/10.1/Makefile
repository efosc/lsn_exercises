CC = g++
CFLAGS = -Wall -O3 --std=c++11

#Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
#	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
#Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
#	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
#clean :
#	rm *.o Monte_Carlo_NVT.exe

funzioni10.o : funzioni10.cpp funzioni10.h
	$(CC) -c funzioni10.cpp -o funzioni10.o $(CFLAGS)

%.exe : %.o random.o funzioni10.o
	$(CC) $^ -o $@ $(CFLAGS)

%.o : %.cpp
	$(CC) -c $^ -o $@ $(CFLAGS)




