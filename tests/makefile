CC = g++
CFLAGS = -Wall -g -std=c++11

testmain: tests-main.o tests-MassDawgNode.o tests-MassDawg.o ../src/MassDawgNode.o ../src/MassDawg.o ../src/utils.o
	${CC} ${CFLAGS} -o testmain tests-main.o tests-MassDawgNode.o tests-MassDawg.o ../src/MassDawgNode.o ../src/MassDawg.o ../src/utils.o

tests-main.o: tests-main.cpp catch.hpp
	${CC} ${CFLAGS} -c tests-main.cpp 

tests-MassDawgNode.o: tests-MassDawgNode.cpp tests-main.cpp 
	${CC} ${CFLAGS} -c tests-MassDawgNode.cpp 

tests-MassDawg.o: tests-MassDawg.cpp tests-main.cpp 
	${CC} ${CFLAGS} -c tests-MassDawg.cpp

clean:
	rm testmain *.o