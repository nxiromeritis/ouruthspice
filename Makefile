CC = gcc
CFLAGS = -g -Wall
OBJ = spicy.o cir_parser.o hashtable.o lists.o mna.o
EXECUTABLE = spicy
DFLAGS = -DCOLORS_ON

CLINK = -lgsl -lgslcblas

all: $(OBJ)
	$(CC) $(OBJ) -o $(EXECUTABLE) $(CLINK)


colors_off: DFLAGS = 
colors_off: $(OBJ)
	@echo "\nBuild with colored output.."
	$(CC) $(OBJ) -o $(EXECUTABLE) $(CLINK)

%.o: %.c %.h
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@


.PHONY: clean clean1 clean2

clean: 
	rm -rvf *.o $(EXECUTABLE) test lists
	rm -rvf *.txt *.png draw.sh

clean1:
	rm -rvf *.o $(EXECUTABLE) test lists

clean2:
	rm -rvf *.txt *.png draw.sh
