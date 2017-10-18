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

int insert_bjt(char *name, char *model_name){
	sec_list_element *tmp;

	tmp = realloc(sec_list.list, ((sec_list.size + 1) * sizeof(sec_list_element)));
	if (tmp == NULL)
		return -1;
	sec_list.list = tmp;

	sec_list.list[sec_list.size].type = Q;

	// TODO strdup
	sec_list.list[sec_list.size].name = name;

	// nodes
	// TODO
	sec_list.list[sec_list.size].model_name = model_name;

	sec_list.size++;

	return 0;
}

int insert_diode(char *name, char *model_name){
	sec_list_element *tmp;

	tmp = realloc(sec_list.list, ((sec_list.size + 1) * sizeof(sec_list_element)));
	if (tmp == NULL)
		return -1;
	sec_list.list = tmp;

	sec_list.list[sec_list.size].type = D;

	// TODO strdup
	sec_list.list[sec_list.size].name = name;

	// nodes
	// TODO
	sec_list.list[sec_list.size].model_name = model_name;

	sec_list.size++;

	return 0;
}

int insert_mos(char *name, char *model_name){
	sec_list_element *tmp;

	tmp = realloc(sec_list.list, ((sec_list.size + 1) * sizeof(sec_list_element)));
	if (tmp == NULL)
		return -1;
	sec_list.list = tmp;

	sec_list.list[sec_list.size].type = M;

	// TODO strdup
	sec_list.list[sec_list.size].name = name;

	// nodes
	// TODO
	sec_list.list[sec_list.size].model_name = model_name;

	sec_list.size++;

	return 0;
}

void print_list1 (){
	printf("\n%s"
	       "-----------------------------\n"
	       "-- List of Team 1 Elements --\n"
	       "-----------------------------\n"
	       "---->Size: %s%lu%s\n"
	       "-----------------------------\n"
	       "-----------------------------\n",
	       BLU, CYN, team1_list.size, BLU);

	for (unsigned long i = 0; i < team1_list.size; i++){
		printf("------>Name: %s%s%s\n"
			   "-----------------------------\n", CYN, team1_list.list[i].name, BLU);
		switch(team1_list.list[i].type){
			case R:
				printf("-------->Type: %sResistor%s\n", CYN, BLU);
				break;
			case C:
				printf("-------->Type: %sCapacitor%s\n", CYN, BLU);
				break;
			case I:
				printf("-------->Type: %sCurrent Source%s\n", CYN, BLU);
				break;
			default:
				break;
		}
		printf("-------->Node (+): %s%s%s\n"
		       "-------->Node (-): %s%s%s\n"
		       "-------->Value: %s%lf%s\n"
		       "-----------------------------\n",
		       CYN, team1_list.list[i].node_plus->name, BLU,
			   CYN, team1_list.list[i].node_minus->name, BLU,
			   CYN, team1_list.list[i].value, BLU);
	}
	printf("------- End of List 1 -------\n"
	       "-----------------------------%s\n", NRM);
}

void print_list2 (){
	printf("\n%s"
	       "-----------------------------\n"
	       "-- List of Team 2 Elements --\n"
	       "-----------------------------\n"
	       "---->Size: %s%lu%s\n"
	       "-----------------------------\n"
	       "-----------------------------\n",
	       BLU, CYN, team2_list.size, BLU);

	for (unsigned long i = 0; i < team2_list.size; i++){
		printf("------>Name: %s%s%s\n"
			   "-----------------------------\n", CYN, team2_list.list[i].name, BLU);
		switch(team2_list.list[i].type){
			case L:
				printf("-------->Type: %sInductor%s\n", CYN, BLU);
				break;
			case V:
				printf("-------->Type: %sVoltage Source%s\n", CYN, BLU);
				break;
			default:
				break;
		}
		printf("-------->Node (+): %s%s%s\n"
		       "-------->Node (-): %s%s%s\n"
		       "-------->Value: %s%lf%s\n"
		       "-----------------------------\n",
		       CYN, team2_list.list[i].node_plus->name, BLU,
			   CYN, team2_list.list[i].node_minus->name, BLU,
			   CYN, team2_list.list[i].value, BLU);
	}
	printf("------- End of List 2 -------\n"
	       "-----------------------------%s\n", NRM);
}

void print_sec_list (){
	printf("\n%s"
	       "-----------------------------\n"
	       "-- List of Unused Elements --\n"
	       "-----------------------------\n"
	       "---->Size: %s%lu%s\n"
	       "-----------------------------\n"
	       "-----------------------------\n",
	       BLU, CYN, sec_list.size, BLU);

	for (unsigned long i = 0; i < sec_list.size; i++){
		printf("------>Name: %s%s%s\n"
			   "-----------------------------\n", CYN, sec_list.list[i].name, BLU);
		switch(sec_list.list[i].type){
			case D:
				printf("-------->Type: %sDiode%s\n"
				       "-------->Node (+): %s%s\n"
				       "-------->Node (-): %s%s\n"
				       "-------->Model name: %s%s%s\n"
				       "-----------------------------\n",
				       CYN, BLU, CYN, BLU, CYN, BLU, CYN, sec_list.list[i].model_name, BLU);
				break;
			case M:
				printf("-------->Type: %sMOS Transistor%s\n"
				       "-------->Node (D): %s%s\n"
				       "-------->Node (G): %s%s\n"
				       "-------->Node (S): %s%s\n"
				       "-------->Node (B): %s%s\n"
				       "-------->Model name: %s%s%s\n"
				       "-----------------------------\n",
				       CYN, BLU, CYN, BLU, CYN, BLU, CYN, BLU, CYN, BLU, CYN, sec_list.list[i].model_name, BLU);
				break;
			case Q:
				printf("-------->Type: %sBJT Transistor%s\n"
				       "-------->Node (C): %s%s\n"
				       "-------->Node (B): %s%s\n"
				       "-------->Node (E): %s%s\n"
				       "-------->Model name: %s%s%s\n"
				       "-----------------------------\n",
				       CYN, BLU, CYN, BLU, CYN, BLU, CYN, BLU, CYN, sec_list.list[i].model_name, BLU);
				break;
			default:
				break;
		}
	}
	printf("-------- End of List --------\n"
	       "-----------------------------%s\n", NRM);
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
