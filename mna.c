#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "spicy.h"
#include "hashtable.h"
#include "lists.h"
#include "mna.h"

double *mna_array = NULL;
double *mna_vector = NULL;
unsigned long mna_dimension_size = 0;



void init_MNA_system() {

	mna_dimension_size = team2_list.size + total_ids - 1;
	printf("MNA Dimension size: %lu (%lu + %lu)\n", mna_dimension_size, total_ids-1, team2_list.size);

	// mna array dimensions: ((n-1) + m2)x((n-1) + m2)
	mna_array = (double *)calloc(mna_dimension_size * mna_dimension_size, sizeof(double));
	if (mna_array == NULL) {
		printf("Error. Memory allocation problems. Exiting..\n");
		exit(EXIT_FAILURE);
	}

	// mna vector dimension: ((n-1) + m2)
	mna_vector = (double *)calloc(mna_dimension_size, sizeof(double));
	if (mna_vector == NULL) {
		printf("Error. Memory allocation problems. Exiting..\n");
		exit(EXIT_FAILURE);
	}

	printf("first array element: %lf\n", mna_array[0]);

}


void fill_MNA_system() {
	unsigned long i;
	unsigned long node_plus_idx;
	unsigned long node_minus_idx;
	double component_value;

	// iterate through the Group1 List and initialise the MNA system
	for (i=0; i < team1_list.size; i++) {


		// TODO TODO TODO NOW
		// handle GROUNDING.. node_minus_idx causes undeflow... add if statements to bypass GND


		node_plus_idx = team1_list.list[i].node_plus->id - 1;
		node_minus_idx = team1_list.list[i].node_minus->id - 1;
		printf("node_plus_idx: %lu\nnode_minus_idx: %lu\n", node_plus_idx, node_minus_idx);

		component_value = team1_list.list[i].value;

		switch(team1_list.list[i].type) {
			case R:
				// array[<+>][<+>] -> +gk
				mna_array[node_plus_idx * mna_dimension_size + node_plus_idx] += 1/component_value;

				// array[<->][<->] -> +gk
				mna_array[node_minus_idx * mna_dimension_size + node_minus_idx] += 1/component_value;

				// array[<+>][<->] -> -gk
				mna_array[node_plus_idx * mna_dimension_size + node_minus_idx] -= 1/component_value;

				// array[<->][<+>] -> -gk
				mna_array[node_minus_idx * mna_dimension_size + node_plus_idx] -= 1/component_value;

				break;
			case I:
				// vector[<->] -> +sk
				mna_vector[node_minus_idx] += component_value;

				// vector[<+>] -> -sk
				mna_vector[node_plus_idx] -= component_value;

				break;
			case C:		// ignored at DC analysis
				break;
			default:
				printf("Unknown type (%d) in list1\n", team1_list.list[i].type);
				exit(EXIT_FAILURE);

		}
	}


	// iterate though the Group2 List and initialise the MNA system
	for (i=0; i < team2_list.size; i++) {

		node_plus_idx = team1_list.list[i].node_plus->id - 1;
		node_minus_idx = team1_list.list[i].node_minus->id - 1;
		component_value = team1_list.list[i].value;

		switch(team2_list.list[i].type) {
			case V:
				// vector[k] -> +sk
				// ... where k = (total_ids-1+i) = (n-1+i)
				mna_vector[(total_ids-1+i)] += component_value;

				// it is intended to enter the case L code
			case L:

				// TODO no need to increment as the ith entry is accessed only be one component
				// array[k][<+>] -> +1
				mna_array[(total_ids-1+i)*mna_dimension_size + node_plus_idx] += 1;

				// array[k][<->] -> -1
				mna_array[(total_ids-1+i)*mna_dimension_size + node_minus_idx] -= 1;

				// array[<+>][k] -> +1
				mna_array[(node_plus_idx)*mna_dimension_size + (total_ids-1+i)] += 1;

				// array[<->][k] -> -1
				mna_array[(node_minus_idx)*mna_dimension_size + (total_ids-1+i)] -= 1;

				// vector[k] -> 0
				// ... where k = (total_ids-1+i) = (n-1+i)
				break;

			// these cases are here because of the optional .spice I R C field G2
			// ...not yet implemented in our MNA
			case I:
			case R:
			case C:
				break;
			default:
				printf("Unknown type (%d) in list2\n", team2_list.list[i].type);
				exit(EXIT_FAILURE);

		}
	}

}


void free_MNA_system() {
	free(mna_array);
	free(mna_vector);
}
