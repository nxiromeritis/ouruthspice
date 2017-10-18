CC = gcc
CFLAGS = -g -Wall
OBJ = spicy.o cir_parser.o hashtable.o lists.o
EXECUTABLE = spicy

all: $(OBJ)
	$(CC) $(OBJ) -o $(EXECUTABLE)

test_hashtable: hashtable.o test.o
	$(CC) hashtable.o test.o -o test	

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@


.PHONY: clean

clean:
	rm -rvf *.o $(EXECUTABLE) test lists
