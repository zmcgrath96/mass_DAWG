# Variables 
CC = g++
CFLAGS = -Wall -g -std=c++11

# Executable
main: main.o MassDawg.o MassDawgNode.o 
	$(CC) $(CFLAGS) -o main main.o MassDawg.o MassDawgNode.o

# Object files
main.o: main.cpp MassDawg.hpp
	$(CC) $(CFLAGS) -c main.cpp

MassDawg.o: MassDawg.hpp MassDawgNode.hpp
	$(CC) $(CFLAGS) -c MassDawg.cpp

MassDawgNode.o: MassDawgNode.hpp 
	$(CC) $(CFLAGS) -c MassDawgNode.cpp

clean:
	rm main *.o