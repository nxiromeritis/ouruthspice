#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "hashtable.h"
#include "spicy.h"


int main ( int argc,char *argv[]){

    int ret;
    element_h *getElement;

  	if(argc <2){
  		printf("Wrong number of arguments ");
  		exit(0);
  	}



  	printf("\t\tINSERTING ELEMENTS\n");

	ht_init( (unsigned long)atol(argv[1]));

	ht_put("name1", 15,&ret);
	if(ret==1){
		printf("one element inserted\n");
	}
	ht_put("name2", 1,&ret);
	ht_put("name3", 2,&ret);
	ht_put("name4", 3 ,&ret);
	ht_put("name5", 4,&ret);
	ht_put("name6", 5,&ret);;
	ht_put("name7", 6,&ret);
	ht_put("name8", 7,&ret);
	ht_put("name9", 8,&ret);
	ht_put("name10", 9,&ret);
	ht_put("name11", 10,&ret);
	ht_put("name12", 11,&ret);
	ht_put("name13", 12,&ret);
	ht_put("name14", 13,&ret);
	ht_put("name15", 14,&ret);
	/*ht_put( hashtable, "name16", "rtuuuutt",&ret);*/
	/*ht_put( hashtable, "name10", "rtytwrwer",&ret);*/
	/*ht_put( hashtable, "name11", "vbcvbmvbm",&ret);*/
	/*ht_put( hashtable, "name12", "werwerw",&ret);*/
	/*ht_put( hashtable, "name13", "bcvbcvbcv",&ret);*/
	/*ht_put( hashtable, "name14", "pinsdfsky",&ret);*/
	/*ht_put( hashtable, "name15", "werwe",&ret);*/
	/*ht_put( hashtable, "name16", "rtuuuutt",&ret);*/
	/*ht_put( hashtable, "name1", "pavlos",&ret);*/
	/*ht_put( hashtable, "name2", "nikos",&ret);*/
	/*ht_put( hashtable, "name3", "george",&ret);*/
	/*ht_put( hashtable, "name4", "thanos",&ret);*/
	/*ht_put( hashtable, "name5", "sifou",&ret);*/

	printf("\t\tGETTING ELEMENTS\n");

	ht_get("name6",&getElement);
	if(getElement!=NULL){
		printf("(%s%s%s , %s%lu%s)\n",RED,getElement->name,NRM,GRN,getElement->id,NRM);

	}

	ht_get("name3",&getElement);
	if(getElement!=NULL){
		printf("(%s%s%s , %s%lu%s)\n",RED,getElement->name,NRM,GRN,getElement->id,NRM);
	}

	ht_get("name16",&getElement);
	if(getElement!=NULL){
		printf("(%s%s%s , %s%lu%s)\n",RED,getElement->name,NRM,GRN,getElement->id,NRM);
	}
	ht_get("name17",&getElement);
	if(getElement!=NULL){
		printf("(%s%s%s , %s%lu%s)\n",RED,getElement->name,NRM,GRN,getElement->id,NRM);
	}


	printHastable();

	printf("\t\t Deleting HashTable\n");
	freeHashTable();
}















