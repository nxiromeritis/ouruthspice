#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "../cir_parser/cir_parser.h"
#include "hashtable.h"
#include "../spicy.h"


hashtable_t *HashTable = NULL;
element_h **id_to_node = NULL;
unsigned long total_ids = 0;


void add_id_to_list(element_h *node, unsigned long id) {
	total_ids = id+1;

	id_to_node = (element_h **)realloc(id_to_node, sizeof(element_h *)*total_ids);
	if (id_to_node == NULL) {
		printf("Error. Memory allocation problems. Exiting..\n");
		exit(EXIT_FAILURE);
	}

	id_to_node[total_ids-1] = node;
}


void free_id_list() {
	free(id_to_node);
	id_to_node = NULL;
	total_ids = 0;
}


void print_id_list() {
	unsigned long i;

	printf("\n\n***\nPrinting ID array\n");
	for (i = 0; i < total_ids; i++) {
		printf("id(" GRN "%lu" NRM ")-> node(" RED "%s " GRN "%lu" NRM ")\n",
				i, id_to_node[i]->name, id_to_node[i]->id);
	}
	printf("***\n");
}


void ht_init(unsigned long size){

	int i;


	if(size <= 0)
		return;


	HashTable = (hashtable_t *)malloc(sizeof(hashtable_t));
	if(HashTable==NULL){
		 perror("Error malloc: ");
	}

	HashTable->table = (element_h **)malloc(size*sizeof(element_h *));
	if(HashTable->table==NULL){
		 perror("Error malloc: ");
	}

	for(i=0 ; i < size; i++){
		HashTable->table[i]=NULL;
	}

	/*
		prosthiki kai allwn stoixeiwn

	*/

	HashTable->size = size;
	HashTable->capacity = 0;
}


unsigned long hs_function(unsigned long size_hs, char *name ){

	unsigned long val=0;
	int i = 0;
	int size = (int)strlen(name);

	while(i < size){
		val += (val << 8) + name[ i ];
		i++;
	}

	return val%size_hs;
}


element_h *newElement(char *name, unsigned long id){

	element_h *element;

	element = (element_h *)malloc(sizeof(element_h));
	if(element == NULL){
		perror("Error malloc: ");
	}

	element->name=strdup(name);
	if(element->name == NULL){
		perror("Error strdup: ");
	}

	element->id = id;
	element->next = NULL;

	return element;
}



element_h *ht_put(char *name, unsigned long id) {

	element_h *curr= NULL;
	element_h *prev= NULL;
	unsigned long index = 0;

	strtoupper(name);
	index = hs_function(HashTable->size, name);


	if(HashTable->table[index] == NULL){
		HashTable->table[index]=newElement(name,id);
		HashTable->capacity++;

		// add the id to the id array
		add_id_to_list(HashTable->table[index], id);
		return HashTable->table[index];

	}else{
		curr = HashTable->table[index];
		while((curr!=NULL) && (strcmp(curr->name,name)!=0)){

			prev=curr;
			curr=curr->next;
		}


		if(curr == NULL){
			prev->next=newElement(name,id);

			// add the id to the id array
			add_id_to_list(prev->next, id);
			return prev->next;
		}
	}

	return NULL;
}


// searches for name_element.
// - getElement is NULL if name not found
//	 otherwise it contains the pointer of the node
element_h *ht_get(char *name) {

	unsigned long index = 0;
	element_h *curr=NULL;

	strtoupper(name);
	index = hs_function(HashTable->size,name);

	curr = HashTable->table[index];
	while((curr!=NULL) && (strcmp(curr->name,name)!=0)){
			curr = curr->next;
	}

	return curr;
}



void printHastable(){
	int i;
	element_h *curr= NULL;


	printf("\t**************HashTable**************** \n");

	printf("\tSize  : %d\n",HashTable->size);
	for(i=0;i < HashTable->size;i++){

		curr = HashTable->table[i];

		printf(" %d ",i);
		while(curr != NULL){
			printf("-> (%s%s%s , %s%lu%s)",RED,curr->name, NRM, GRN, curr->id, NRM);

			curr=curr->next;
		}
		printf("\n");

	}
}

void freeHashTable(){

	int i;
	element_h *curr = NULL;
	element_h *prev = NULL;

	for(i = 0; i < HashTable->size ; i++){
		curr = HashTable->table[i];
		while(curr!= NULL){
			prev=curr;
			curr = curr->next;

			free(prev->name);
			prev->name = NULL;
			free(prev);
			prev = NULL;

		}

	}

	free_id_list();

	free(HashTable->table);
	HashTable->table = NULL;

	free(HashTable);
	HashTable = NULL;

}
