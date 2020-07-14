CC=gcc
CFLAGS=-I.
DEPS = MassPrefixTree.h
OBJ = main.o MassPrefixTree.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -r *.o
