CC = gcc
CFLAGS = -g -Wall
OBJ = spicy.o cir_parser.o hashtable.o lists.o mna.o
EXECUTABLE = spicy
DFLAGS = -DCOLORS_ON


all: $(OBJ)
	$(CC) $(OBJ) -o $(EXECUTABLE)


colors_off: DFLAGS = 
colors_off: $(OBJ)
	@echo "\nBuild with colored output.."
	$(CC) $(OBJ) -o $(EXECUTABLE)

%.o: %.c %.h
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@


.PHONY: clean

clean:
	rm -rvf *.o $(EXECUTABLE) test lists
