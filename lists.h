#ifndef _LISTS_
#define _LISTS_

#include "hashtable.h"

typedef struct RLCS_list_element {
	c_type type;
	char *name;

	element_h *node_plus;	// node +
	element_h *node_minus;	// node -

	double value;
}list_element;


typedef struct RLCS_list_head{
	unsigned long size;
	list_element *list;
}list_head;

typedef union{
	struct {
		element_h *node_plus;	// node +
		element_h *node_minus;	// node -
	}diode;
	struct {
		element_h *node_d;		// node D
		element_h *node_g;		// node G
		element_h *node_s;		// node S
		element_h *node_b;		// node B
		long l;
		long w;
	}mos;
	struct {
		element_h *node_c;		// node C
		element_h *node_b;		// node B
		element_h *node_e;		// node E
	}bjt;
}characteristics;

typedef struct DTran_list_element {
	c_type type;
	char *name;
	char *model_name;
	characteristics character;
}sec_list_element;

typedef struct DTran_list_head{
	unsigned long size;
	sec_list_element *list;
}sec_list_head;



extern list_head team1_list;
extern list_head team2_list;
extern sec_list_head sec_list;

void init_lists();
void free_lists();

int insert_element(list_head *list, c_type type, char *name, element_h *node_plus, element_h *node_minus, double value);

int insert_bjt(char *name, element_h *node_c, element_h *node_b, element_h *node_e, char *model_name);
int insert_diode(char *name, element_h *node_plus, element_h *node_minus, char *model_name);
int insert_mos(char *name, element_h *node_d, element_h *node_g, element_h *node_s, element_h *node_b, long l, long w, char *model_name);

void print_list1();
void print_list2();
void print_sec_list();

#endif
