CC = g++
CFLAGS = -Wall -g -std=c++11

testmain: tests-main.o tests-MassDawgNode.o ../src/MassDawgNode.o
	${CC} ${CFLAGS} -o testmain tests-main.o tests-MassDawgNode.o ../src/MassDawgNode.o

tests-main.o: tests-main.cpp catch.hpp
	${CC} ${CFLAGS} -c tests-main.cpp 

tests-MassDawgNode.o: tests-MassDawgNode.cpp tests-main.cpp 
	${CC} ${CFLAGS} -c tests-MassDawgNode.cpp 

clean:
	rm testmain *.o