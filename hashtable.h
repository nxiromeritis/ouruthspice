#ifndef __HASHTABLE_H__
#define __HASHTABLE_H__


typedef struct element_h{
	char *name;	// name
	unsigned long id;
	// TODO: add lists of pointers to elements
	// according to the connection type

	double val;
	struct element_h *next;

}element_h;


typedef struct hashtable_t{
	int size;
	int capacity;
	struct element_h **table;

}hashtable_t;


extern hashtable_t *HashTable;
extern element_h **id_to_node;
extern unsigned long total_ids;

extern void ht_init(unsigned long size);
extern unsigned long hs_function(unsigned long size_hs, char *name);
extern element_h *newElement(char *name,unsigned long id);
extern element_h *ht_put(char *name, unsigned long id);
extern element_h * ht_get(char *name);
extern void printHastable();
extern void freeHashTable();

extern void add_id_to_list();
extern void free_id_list();
extern void print_id_list();

#endif
