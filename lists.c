#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spicy.h"
#include "lists.h"


list_head team1_list;
list_head team2_list;
sec_list_head sec_list;


void init_lists(){
	team1_list.size = 0;
	team2_list.size = 0;
	sec_list.size = 0;
	team1_list.list = NULL;
	team2_list.list = NULL;
	sec_list.list = NULL;
}

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

	for(i=0; i < team2_list.size; i++) {
		free(sec_list.list[i].name);
		free(sec_list.list[i].model_name);
	}
	free(sec_list.list);
	sec_list.size = 0;
}

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

	// nodes
	list->list[list->size].node_plus = node_plus;
	list->list[list->size].node_minus = node_minus;

	list->list[list->size].value = value;

	list->size++;
	printf("%lu\n", list->size);
	return 0;
}

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

	// nodes
	sec_list.list[sec_list.size].character.diode.node_plus = node_plus;
	sec_list.list[sec_list.size].character.diode.node_minus = node_minus;
	
	sec_list.list[sec_list.size].model_name = strdup(model_name);
	if (sec_list.list[sec_list.size].model_name == NULL)
		return -1;

	sec_list.size++;

	return 0;
}

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

void print_list1 (){
	printf("\n" BLU
	       "-----------------------------\n"
	       "-- List of Team 1 Elements --\n"
	       "-----------------------------\n"
	       "---->Size: " CYN "%lu\n" BLU
	       "-----------------------------\n"
	       "-----------------------------\n",
	       team1_list.size);

	for (unsigned long i = 0; i < team1_list.size; i++){
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

void print_list2 (){
	printf("\n" BLU
	       "-----------------------------\n"
	       "-- List of Team 2 Elements --\n"
	       "-----------------------------\n"
	       "---->Size: " CYN "%lu\n" BLU
	       "-----------------------------\n"
	       "-----------------------------\n",
	       team2_list.size);

	for (unsigned long i = 0; i < team2_list.size; i++){
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

void print_sec_list (){
	printf("\n" BLU
	       "-----------------------------\n"
	       "-- List of Unused Elements --\n"
	       "-----------------------------\n"
	       "---->Size: " CYN "%lu\n" BLU
	       "-----------------------------\n"
	       "-----------------------------\n",
	       sec_list.size);

	for (unsigned long i = 0; i < sec_list.size; i++){
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

/*int main(){*/
	/*init_lists();*/

	/*int i = insert_element(&team1_list, R, "1", 5e-3);*/
	/*if (!(i==0))*/
		/*printf("Shit!\n");*/

	/*print_list1();*/
	/*print_list2();*/
	/*print_sec_list();*/

	/*free_lists();*/

	/*return 0;*/
/*}*/
