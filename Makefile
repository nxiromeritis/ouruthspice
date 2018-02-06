CC = gcc
CFLAGS = -g -Wall 
OBJ = build/spicy.o build/cir_parser/cir_parser.o build/hashtable/hashtable.o build/lists/lists.o build/mna/mna.o build/csparse/csparse.o
BFOLDERS = build/ build/cir_parser/ build/hashtable/ build/lists/ build/mna/ build/csparse/
EXECUTABLE = spicy
DFLAGS = -DCOLORS_ON

CLINK = -lgsl -lgslcblas -lm

all: $(OBJ)
	@mkdir -p build
	$(CC) $(OBJ) -o $(EXECUTABLE) $(CLINK)


colors_off: DFLAGS = 
colors_off: $(OBJ)
	@echo "\nBuild with colored output.."
	$(CC) $(OBJ) -o $(EXECUTABLE) $(CLINK)

build/%/%.o: src/%/%.c src/%/%.h
	@mkdir -p $(BFOLDERS)
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@
build/%.o: src/%.c src/%.h
	@mkdir -p $(BFOLDERS)
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@


.PHONY: clean clean_obj clean_outputs

clean: 
	rm -rvf $(OBJ) $(EXECUTABLE) $(BFOLDERS)
	rm -rvf *.txt *.out *.png draw.sh

clean_obj:
	rm -rvf $(OBJ) $(EXECUTABLE) $(BFOLDERS)

clean_outputs:
	rm -rvf *.txt *.out *.png draw.sh
