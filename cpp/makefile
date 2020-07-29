# Variables 
CC = g++
CFLAGS = -Wall -g -std=c++11

# Executable
main: main.o MassDawg.o MassDawgNode.o Edge.o
	$(CC) $(CFLAGS) -o main main.o MassDawg.o MassDawgNode.o Edge.o

# Object files
main.o: main.cpp MassDawg.hpp
	$(CC) $(CFLAGS) -c main.cpp

MassDawg.o: MassDawg.hpp MassDawgNode.hpp Edge.hpp
	$(CC) $(CFLAGS) -c MassDawg.cpp

MassDawgNode.o: MassDawgNode.hpp Edge.hpp
	$(CC) $(CFLAGS) -c MassDawgNode.cpp

Edge.o: Edge.hpp
	$(CC) $(CFLAGS) -c Edge.cpp

clean:
	rm main *.o