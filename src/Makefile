# Variables 
CC = g++
CFLAGS = -Wall -g -std=c++11

# Executable
main: main.o MassDawg.o MassDawgNode.o utils.o
	$(CC) $(CFLAGS) -o main main.o MassDawg.o MassDawgNode.o utils.o

test: test.o MassDawg.o MassDawgNode.o utils.o
	$(CC) $(CFLAGS) -o test test.o MassDawg.o MassDawgNode.o utils.o

# Object files
main.o: main.cpp MassDawg.hpp
	$(CC) $(CFLAGS) -c main.cpp

test.o: test.cpp MassDawg.hpp
	$(CC) $(CFLAGS) -c test.cpp

MassDawg.o: MassDawg.hpp MassDawgNode.hpp utils.hpp
	$(CC) $(CFLAGS) -c MassDawg.cpp 

MassDawgNode.o: MassDawgNode.hpp 
	$(CC) $(CFLAGS) -c MassDawgNode.cpp

utils.o: utils.hpp
	$(CC) $(CFLAGS) -c utils.cpp

clean:
	rm main test *.o