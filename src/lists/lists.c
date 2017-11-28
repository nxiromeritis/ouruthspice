#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cir_parser/cir_parser.h"
#include "../spicy.h"
#include "lists.h"
#include "../mna/mna.h"


list_head team1_list;
list_head team2_list;
sec_list_head sec_list;


// and array of strings that represent
// the commands that will be executed
char **command_list = NULL;
unsigned int command_list_len = 0;




// adds the given command to the command list
// NOTE1: Among many subsequent .DC commands, only the last is stored, no point solving the previous MNAs
// NOTE2: Each .PRINT or .PLOT command is refering the previous DC command
// NOTE3: Plot and print commands that have no DC command before them will be ignored DURING EXECUTION
// NOTE4: DC commands that have not print or plot commands after them will be ignored DURING EXECUTION
// NOTE5: Possible errors contained in command parameters will be unveiled during execution (interpreter logic)
void add_command_to_list(char *command) {
	static byte prev_command_is_dc = 0;
	int strcmp_res;

	char *token = NULL;
	const char delim[5] = " \r\t\n";
	unsigned int tok_count = 0;

	/*printf("COMMAND: %s\n", command);*/


	// if the parsed command is .OPTIONS execute the command during parsing
	strcmp_res = strncmp(command, ".OPTIONS", 8);
	if (strcmp_res == 0) {

		// bypass command name
		token = strtok(command, delim);

		token = strtok(NULL, delim);
		while (token != NULL) {
			tok_count++;

			if (strcmp(token, "SPD") == 0) {

				// we assume that the user knows what he does and he does not give
				// the same option multiple times
				if (solver_type == BI_CG_SOLVER)
					solver_type = CG_SOLVER;
				else if (solver_type  != CG_SOLVER)
					solver_type = CHOL_SOLVER;

			}
			else if (strcmp(token, "ITER") == 0) {

				// we assume that the user knows what he does and he does not give
				// the same option multiple times
				if (solver_type == CHOL_SOLVER)
					solver_type = CG_SOLVER;
				else if (solver_type != CG_SOLVER)
					solver_type = BI_CG_SOLVER;
			}
			else if (strncmp(token, "ITOL", 4) == 0) {
				// update itol
				parse_double(&itol, &token[5]);
			}
			else if (strncmp(token, "SPARSE", 6) == 0) {
				// sparse matrixes
				is_sparse = 1;
			}
			// else bypass argument (future arguments: ITER SPARSE)

			token = strtok(NULL, delim);
		}
		return;
	}


	strcmp_res = strncmp(command, ".DC ", 4);

	// if current command is .dc and the previous one is also dc..
	if ( (strcmp_res == 0) && (prev_command_is_dc == 1) ) {

		// the new dc command will overwrite the previous one
		free(command_list[command_list_len-1]);
	}
	else { // allocate one more command entry in the list

		command_list_len++;
		command_list = (char **) realloc(command_list, command_list_len*sizeof(char *));
		if (command_list == NULL) {
			printf("Error. Memory allocation problems. Exiting..\n");
			exit(EXIT_FAILURE);
		}

	}

	// store the command

	printf("\t\t\t\t\tcommand : %s\n",command);
	command_list[command_list_len-1] = strdup(command);
	if (command_list[command_list_len-1] == NULL) {
		printf("Error. Memory allocation problems. Exiting..\n");
		exit(EXIT_FAILURE);
	}

	// update the .dc flag
	/*prev_command_is_dc = !(strcmp_res); // correct but lacks readability*/
	if (strcmp_res == 0)
		prev_command_is_dc = 1;
	else
		prev_command_is_dc = 0;

}


// prints the command list
void print_command_list() {
	unsigned int i;

	printf(BLU "\n*COMMAND LIST*\n" NRM);

	for (i=0; i < command_list_len; i++)
		printf("#%u: %s\n", i, command_list[i]);

	printf(BLU "*END OF COMMAND LIST*\n\n" NRM);
}



// frees the command list
void free_command_list() {
	unsigned int i;

	for (i=0; i < command_list_len; i++) {
		free(command_list[i]);
	}
	free(command_list);
}




// Initialize the lists
void init_lists(){
	team1_list.size = 0;
	team2_list.size = 0;
	sec_list.size = 0;
	team1_list.list = NULL;
	team2_list.list = NULL;
	sec_list.list = NULL;
}

// Free the lists and their allocated fields
void free_lists(){
	unsigned long i;

	for(i=0; i < team1_list.size; i++) {
		free(team1_list.list[i].name);
	}
	free(team1_list.list);
	team1_list.size = 0;

	for(i=0; i < team2_list.size; i++) {
		free(team2_list.list[i].name);
	}
	free(team2_list.list);
	team2_list.size = 0;

	for(i=0; i < sec_list.size; i++) {
		free(sec_list.list[i].name);
		free(sec_list.list[i].model_name);
	}
	free(sec_list.list);
	sec_list.size = 0;
}

// Insert an element into one of the two lists
int insert_element(list_head *list, c_type type, char *name, element_h *node_plus, element_h *node_minus, double value){
	list_element *tmp;

	tmp = realloc(list->list, ((list->size + 1) * sizeof(list_element)));
	if (tmp == NULL)
		return -1;
	list->list = tmp;

	list->list[list->size].type = type;

	list->list[list->size].name = strdup(name);
	if (list->list[list->size].name == NULL)
		return -1;
	strtoupper(list->list[list->size].name);

	// nodes
	list->list[list->size].node_plus = node_plus;
	list->list[list->size].node_minus = node_minus;

	list->list[list->size].op_point_val = value;
	list->list[list->size].value = value;

	list->size++;
	printf("%lu\n", list->size);
	return 0;
}

// Insert a BJT Transistor to the sec_list
int insert_bjt(char *name, element_h *node_c, element_h *node_b, element_h *node_e, char *model_name){
	sec_list_element *tmp;

	tmp = realloc(sec_list.list, ((sec_list.size + 1) * sizeof(sec_list_element)));
	if (tmp == NULL)
		return -1;
	sec_list.list = tmp;

	sec_list.list[sec_list.size].type = Q;

	sec_list.list[sec_list.size].name = strdup(name);
	if (sec_list.list[sec_list.size].name == NULL)
		return -1;
	strtoupper(sec_list.list[sec_list.size].name);

	// nodes
	sec_list.list[sec_list.size].character.bjt.node_c = node_c;
	sec_list.list[sec_list.size].character.bjt.node_b = node_b;
	sec_list.list[sec_list.size].character.bjt.node_e = node_e;

	sec_list.list[sec_list.size].model_name = strdup(model_name);
	if (sec_list.list[sec_list.size].model_name == NULL)
		return -1;

	sec_list.size++;

	return 0;
}

// Inser a Diode into the sec_list
int insert_diode(char *name, element_h *node_plus, element_h *node_minus, char *model_name){
	sec_list_element *tmp;

	tmp = realloc(sec_list.list, ((sec_list.size + 1) * sizeof(sec_list_element)));
	if (tmp == NULL)
		return -1;
	sec_list.list = tmp;

	sec_list.list[sec_list.size].type = D;

	sec_list.list[sec_list.size].name = strdup(name);
	if (sec_list.list[sec_list.size].name == NULL)
		return -1;
	strtoupper(sec_list.list[sec_list.size].name);

	// nodes
	sec_list.list[sec_list.size].character.diode.node_plus = node_plus;
	sec_list.list[sec_list.size].character.diode.node_minus = node_minus;

	sec_list.list[sec_list.size].model_name = strdup(model_name);
	if (sec_list.list[sec_list.size].model_name == NULL)
		return -1;

	sec_list.size++;

	return 0;
}

// Insert a MOS Transistor into the sec_list
int insert_mos(char *name, element_h *node_d, element_h *node_g, element_h *node_s, element_h *node_b, long l, long w, char *model_name){
	sec_list_element *tmp;

	tmp = realloc(sec_list.list, ((sec_list.size + 1) * sizeof(sec_list_element)));
	if (tmp == NULL)
		return -1;
	sec_list.list = tmp;

	sec_list.list[sec_list.size].type = M;

	sec_list.list[sec_list.size].name = strdup(name);
	if (sec_list.list[sec_list.size].name == NULL)
		return -1;
	strtoupper(sec_list.list[sec_list.size].name);

	// nodes
	sec_list.list[sec_list.size].character.mos.node_d = node_d;
	sec_list.list[sec_list.size].character.mos.node_g = node_g;
	sec_list.list[sec_list.size].character.mos.node_s = node_s;
	sec_list.list[sec_list.size].character.mos.node_b = node_b;
	sec_list.list[sec_list.size].character.mos.l = l;
	sec_list.list[sec_list.size].character.mos.w = w;

	sec_list.list[sec_list.size].model_name = strdup(model_name);
	if (sec_list.list[sec_list.size].model_name == NULL)
		return -1;

	sec_list.size++;

	return 0;
}

// Print the elements of the team1_list
void print_list1 (){
	unsigned long i;

	printf("\n" BLU
	       "-----------------------------\n"
	       "-- List of Team 1 Elements --\n"
	       "-----------------------------\n"
	       "---->Size: " CYN "%lu\n" BLU
	       "-----------------------------\n"
	       "-----------------------------\n",
	       team1_list.size);

	for (i = 0; i < team1_list.size; i++){
		printf("------>Name: " CYN "%s\n" BLU
			   "-----------------------------\n", team1_list.list[i].name);
		switch(team1_list.list[i].type){
			case R:
				printf("-------->Type: " CYN "Resistor\n" BLU);
				break;
			case C:
				printf("-------->Type: " CYN "Capacitor\n" BLU);
				break;
			case I:
				printf("-------->Type: " CYN "Current Source\n" BLU);
				break;
			default:
				break;
		}
		printf("-------->Node (+): " CYN "%s\n" BLU
		       "-------->Node (-): " CYN "%s\n" BLU
		       "-------->Value: " CYN "%lf\n" BLU
		       "-----------------------------\n",
		       team1_list.list[i].node_plus->name,
			   team1_list.list[i].node_minus->name,
			   team1_list.list[i].value);
	}
	printf("------- End of List 1 -------\n"
	       "-----------------------------\n" NRM);
}

// Print the elements of the team2_list
void print_list2 (){
	unsigned long i;

	printf("\n" BLU
	       "-----------------------------\n"
	       "-- List of Team 2 Elements --\n"
	       "-----------------------------\n"
	       "---->Size: " CYN "%lu\n" BLU
	       "-----------------------------\n"
	       "-----------------------------\n",
	       team2_list.size);

	for (i = 0; i < team2_list.size; i++){
		printf("------>Name: " CYN "%s\n" BLU
			   "-----------------------------\n", team2_list.list[i].name);
		switch(team2_list.list[i].type){
			case R:
				printf("-------->Type: " CYN "Resistor [G2]\n" BLU);
				break;
			case L:
				printf("-------->Type: " CYN "Inductor\n" BLU);
				break;
			case C:
				printf("-------->Type: " CYN "Capacitor [G2]\n" BLU);
				break;
			case V:
				printf("-------->Type: " CYN "Voltage Source\n" BLU);
				break;
			case I:
				printf("-------->Type: " CYN "Current Source [G2]\n" BLU);
				break;
			default:
				break;
		}
		printf("-------->Node (+): " CYN "%s\n" BLU
		       "-------->Node (-): " CYN "%s\n" BLU
		       "-------->Value: " CYN "%lf\n" BLU
		       "-----------------------------\n",
		       team2_list.list[i].node_plus->name,
			   team2_list.list[i].node_minus->name,
			   team2_list.list[i].value);
	}
	printf("------- End of List 2 -------\n"
	       "-----------------------------\n" NRM);
}

// Print the elements of the sec_list
void print_sec_list (){
	unsigned long i;

	printf("\n" BLU
	       "-----------------------------\n"
	       "-- List of Unused Elements --\n"
	       "-----------------------------\n"
	       "---->Size: " CYN "%lu\n" BLU
	       "-----------------------------\n"
	       "-----------------------------\n",
	       sec_list.size);

	for (i = 0; i < sec_list.size; i++){
		printf("------>Name: " CYN "%s\n" BLU
			   "-----------------------------\n", sec_list.list[i].name);
		switch(sec_list.list[i].type){
			case D:
				printf("-------->Type: " CYN "Diode\n" BLU
				       "-------->Node (+): " CYN "%s\n" BLU
				       "-------->Node (-): " CYN "%s\n" BLU
				       "-------->Model name: " CYN "%s\n" BLU
				       "-----------------------------\n",
				       sec_list.list[i].character.diode.node_plus->name,
				       sec_list.list[i].character.diode.node_minus->name,
				       sec_list.list[i].model_name);
				break;
			case M:
				printf("-------->Type: " CYN "MOS Transistor\n" BLU
				       "-------->Node (D): " CYN "%s\n" BLU
				       "-------->Node (G): " CYN "%s\n" BLU
				       "-------->Node (S): " CYN "%s\n" BLU
				       "-------->Node (B): " CYN "%s\n" BLU
				       "-------->Model name: " CYN "%s\n" BLU
				       "-----------------------------\n",
				       sec_list.list[i].character.mos.node_d->name,
				       sec_list.list[i].character.mos.node_g->name,
				       sec_list.list[i].character.mos.node_s->name,
				       sec_list.list[i].character.mos.node_b->name,
				       sec_list.list[i].model_name);
				break;
			case Q:
				printf("-------->Type: " CYN "BJT Transistor\n" BLU
				       "-------->Node (C): " CYN "%s\n" BLU
				       "-------->Node (B): " CYN "%s\n" BLU
				       "-------->Node (E): " CYN "%s\n" BLU
				       "-------->Model name: " CYN "%s\n" BLU
				       "-----------------------------\n",
				       sec_list.list[i].character.bjt.node_c->name,
				       sec_list.list[i].character.bjt.node_b->name,
				       sec_list.list[i].character.bjt.node_e->name,
				       sec_list.list[i].model_name);
				break;
			default:
				break;
		}
	}
	printf("-------- End of List --------\n"
	       "-----------------------------\n" NRM);
}

#ifdef STANDALONE
int main(){
	init_lists();

	int i = insert_element(&team1_list, R, "1", 5e-3);
	if (!(i==0))
		printf("Shit!\n");

	print_list1();
	print_list2();
	print_sec_list();

	free_lists();

	return 0;
}
#endif
