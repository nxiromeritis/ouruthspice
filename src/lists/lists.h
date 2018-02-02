#ifndef _LISTS_
#define _LISTS_

#include "../hashtable/hashtable.h"


//typedef enum component_type {undefined, R, L, C, V, I, D, M, Q} c_type;
typedef enum component_type {undefined = 0 , R ='R' , L ='L', C ='C', V ='V', I ='I' , D ='D',M ='M', Q ='Q'} c_type;

#define TR_TYPE_NONE	0
#define TR_TYPE_SIN		1
#define TR_TYPE_EXP		2
#define TR_TYPE_PULSE	3
#define TR_TYPE_PWL		4


typedef struct ExpInfo {
	double i1;
	double i2;
	double td1;
	double tc1;
	double td2;
	double tc2;
} ExpInfoT;

typedef struct SinInfo {
	double i1;
	double ia;
	double fr;
	double td;
	double df;
	double ph;
} SinInfoT;

typedef struct PulseInfo {
	double i1;
	double i2;
	double td;
	double tr;
	double tf;
	double pw;
	double per;
} PulseInfoT;

typedef struct PwlInfo {
	double *times;
	double *values;
	unsigned int total_tuples;
} PwlInfoT;


typedef union tranSpec {
	void *data;
	struct ExpInfo *exp_data;
	struct SinInfo *sin_data;
	struct PulseInfo *pulse_data;
	struct PwlInfo *pwl_data;
} tranSpecT;



typedef struct RLCS_list_element {
	c_type type;
	char *name;

	element_h *node_plus;	// node +
	element_h *node_minus;	// node -

	double op_point_val;
	double value;

	// fields for transient analysis
	int tr_type;
	union tranSpec tran_spec;
}list_element;




typedef struct RLCS_list_head {
	unsigned long size;
	list_element *list;
} list_head;



typedef union {
	struct {
		element_h *node_plus;	// node +
		element_h *node_minus;	// node -
	} diode;
	struct {
		element_h *node_d;		// node D
		element_h *node_g;		// node G
		element_h *node_s;		// node S
		element_h *node_b;		// node B
		long l;
		long w;
	} mos;
	struct {
		element_h *node_c;		// node C
		element_h *node_b;		// node B
		element_h *node_e;		// node E
	} bjt;
} characteristics;



typedef struct DTran_list_element {
	c_type type;
	char *name;
	char *model_name;
	characteristics character;
} sec_list_element;

typedef struct DTran_list_head {
	unsigned long size;
	sec_list_element *list;
} sec_list_head;



extern list_head team1_list;
extern list_head team2_list;
extern sec_list_head sec_list;

extern char **command_list;
extern unsigned int command_list_len;

void parse_plot(char *command);
extern void add_command_to_list(char *command);
extern void free_command_list();
extern void print_command_list();

extern void init_lists();
extern void free_lists();

extern int insert_element(list_head *list, c_type type, char *name, element_h *node_plus, element_h *node_minus, double value, int tr_type, void *tran_spec_data);

extern int insert_bjt(char *name, element_h *node_c, element_h *node_b, element_h *node_e, char *model_name);
extern int insert_diode(char *name, element_h *node_plus, element_h *node_minus, char *model_name);
extern int insert_mos(char *name, element_h *node_d, element_h *node_g, element_h *node_s, element_h *node_b, long l, long w, char *model_name);

extern void test_tran_spec();
extern void print_list1();
extern void print_list2();
extern void print_sec_list();

#endif
