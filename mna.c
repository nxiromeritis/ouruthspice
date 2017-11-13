#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

#include "cir_parser.h"
#include "spicy.h"
#include "hashtable.h"
#include "lists.h"
#include "mna.h"

double *mna_array = NULL;
double *mna_vector = NULL;
unsigned long mna_dimension_size = 0;

gsl_matrix_view gsl_mna_array;
gsl_vector_view gsl_mna_vector;
gsl_permutation *gsl_p = NULL;

byte solver_type = LU_SOLVER;
double itol = ITOL_DEFAULT;

void decomp_lu_MNA() {
	int s;

	gsl_mna_array = gsl_matrix_view_array(mna_array, mna_dimension_size, mna_dimension_size);
	gsl_mna_vector = gsl_vector_view_array(mna_vector, mna_dimension_size);
	gsl_p = gsl_permutation_alloc(mna_dimension_size);

	// user's responsibility is this fails
	gsl_linalg_LU_decomp(&gsl_mna_array.matrix, gsl_p, &s);

}

void decomp_cholesky_MNA() {

	gsl_mna_array = gsl_matrix_view_array(mna_array, mna_dimension_size, mna_dimension_size);
	gsl_mna_vector = gsl_vector_view_array(mna_vector, mna_dimension_size);

	// user's responsibility if this fails
	gsl_linalg_cholesky_decomp(&gsl_mna_array.matrix);
}


// dont forget to free the permutation after the last call of this function
void solve_lu_MNA() {
	unsigned long i;
	gsl_vector *gsl_x = NULL;

	gsl_x = gsl_vector_alloc(mna_dimension_size);

	gsl_linalg_LU_solve(&gsl_mna_array.matrix, gsl_p, &gsl_mna_vector.vector, gsl_x);

	// id_to_nodes has the GND located at idx 0
	// although MNA array and vector ignore GND and the first node starts from 0
	for (i=1; i < total_ids; i++) {
		id_to_node[i]->val = gsl_vector_get(gsl_x, i-1);
	}

	gsl_vector_free(gsl_x);
}


void solve_cholesky_MNA() {
	unsigned long i;
	gsl_vector *gsl_x = NULL;

	gsl_x = gsl_vector_alloc(mna_dimension_size);

	gsl_linalg_cholesky_solve(&gsl_mna_array.matrix, &gsl_mna_vector.vector, gsl_x);

	// id_to_nodes has the GND located at idx 0
	// although MNA array and vector ignore GND and the first node starts from 0
	for (i=1; i < total_ids; i++) {
		id_to_node[i]->val = gsl_vector_get(gsl_x, i-1);
	}

	gsl_vector_free(gsl_x);
}


// executes the command_list (command .OPTIONS is excluded from the list as it is executed during the parsing phase)
void execute_commands() {
	unsigned long i,k;
	const char delim[5] = " \r\t\n";
	char *token = NULL;

	// variables used for DC command
	FILE *fp_draw = NULL;
	FILE *node_fp = NULL;
	char *var_name = NULL;
	char *node_name = NULL;
	char *filename = NULL;
	double j = 0;
	double start = 0;
	double end = 0;
	double jump = 0;
	byte var_found = 0; // 0 if not found, 1 if found in list1, 2 if found in list2
	element_h *node = NULL;
	unsigned long idx1;
	unsigned long idx2;
	list_element *var = NULL;

	// TODO -> do with malloc
	//char tmp_name[128];



	for (i = 0; i < command_list_len; i++) {
		if (strncmp(command_list[i], ".DC ", 4) == 0) {

			if( var_name != NULL){

				start = 0;
				end = 0;
				jump = 0;
				free(var_name);
				var_name = NULL;
			}


			// Command name
			token = strtok(command_list[i], delim);
			if (token == NULL) {
				printf(RED "Error" NRM ": Not enough arguments (%s)\n Bypassing..\n", command_list[i]);
				continue;
			}
			// nothing to store here


			// v/i source name
			token = strtok(NULL, delim);
			if (token == NULL) {
				printf(RED "Error" NRM ": Not enough arguments (%s)\n Bypassing..\n", command_list[i]);
				continue;
			}
			// before storing var name check if it exist inside the list
			var_found = 0;


			if( toupper(token[0]) == 'I'){
				for (k=0; k < team1_list.size; k++) {

					// token + 1 becase we bypass the first character that refers to the component type
					if ((strcmp(token+1, team1_list.list[k].name) == 0 ) && (team1_list.list[k].type == 'I')) {
						var_found = 1;
						var = &team1_list.list[k];
						idx1 = var->node_plus->id - 1;
						idx2 = var->node_minus->id -1;
						break;
					}
				}
			}

			if (var_found == 0) {

				if(toupper(token[0]) == 'V'){

					for (k=0; k < team2_list.size; k++) {
						if ((strcmp(token+1, team2_list.list[k].name) == 0 ) && (team2_list.list[k].type == 'V')) {
							var_found = 2;
							var = &team2_list.list[k];
							idx1 = k + total_ids - 1;
							break;
						}
					}
				}
			}
			// not found in any of those lists
			if (var_found == 0) {
				printf(RED "Error" NRM ": .DC variable not found\n Bypassing..\n");
				continue;
			}
			var_name = strdup(token);
			if (var_name == NULL) {
				printf("Error. Memory allocation problems. Exiting..\n");
				exit(EXIT_FAILURE);
			}



			// source start value
			token = strtok(NULL, delim);
			if (token == NULL) {
				printf(RED "Error" NRM ": Not enough arguments (%s)\n Bypassing..\n", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}
			if (parse_double(&start, token) == 0) {
				printf(RED "Error" NRM ": Invalid argument value (%s)\n Bypassing..\n", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}



			// source end value
			token = strtok(NULL, delim);
			if (token == NULL) {
				printf(RED "Error" NRM ": Not enough arguments (%s)\n Bypassing..", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}
			if (parse_double(&end, token) == 0) {
				printf(RED "Error" NRM ": Invalid argument value (%s)\n Bypassing..", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}



			// source increment step
			token = strtok(NULL, delim);
			if (token == NULL) {
				printf(RED "Error" NRM ": Not enough arguments (%s)\n Bypassing..", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}
			if (parse_double(&jump, token) == 0) {
				printf(RED "Error" NRM ": Invalid argument value (%s)\n Bypassing..", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}


			// possible extra false arguments
			token = strtok(NULL, delim);
			if (token != NULL) {
				printf(RED "Error" NRM ": Command contains extra false arguments (%s)\n Bypassing", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}

			// the iteration values are set up
			// the next plot command will solve the system
			// and write the results to the file

			// no need to free var_name here. it is needed during the plot command
		}
		if ((strncmp(command_list[i], ".PRINT ", 7) == 0) || (strncmp(command_list[i], ".PLOT ", 6) == 0)) {
			// TODO: multiple nodes in one line of .print or .plot command
			// Make the parser execute the options command and create a list of dc-print commands

			if(var_found == 0 ){
				continue;
			}

			// bypass command name
			token = strtok(command_list[i], delim);
			if (token == NULL) {
				printf(RED "Error" NRM ": Not enough arguments (%s)\n Bypassing..\n", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}

			// store node name
			token = strtok(NULL, delim);
			printf("token: %s\n", token);
			if (token == NULL) {
				printf(RED "Error" NRM ": Not enough arguments (%s)\n Bypassing..\n", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}
			if ((toupper(token[0]) != 'V') || (token[1] != '(') || (token[strlen(token)-1] != ')') ) {
				printf(RED "Error" NRM ": Invalid argument value (%s)\n Bypassing..\n", token);
				free(var_name);
				var_name = NULL;
				continue;
			}
			// +1 for '\0', -1 for 'v', -1 for '(' and -1 for ')'
			node_name = malloc( (strlen(token) - 2)*sizeof(char));
			if (node_name == NULL) {
				printf("Error. Memory allocation problems. Exiting..\n");
				exit(EXIT_FAILURE);
			}
			snprintf(node_name, strlen(token)-2, "%s", &token[2]);
			node = ht_get(node_name);
			if (node == NULL) {
				printf(RED "Error" NRM ": Node not found (%s)\n Bypassing\n", token);
				free(var_name);
				free(node_name);
				var_name = NULL;
				node_name = NULL;
				continue;
			}

			/***********/
			node_name = realloc(node_name,(strlen(token)+1)*sizeof(char));
			if (node_name == NULL) {
				printf("Error. Memory allocation problems. Exiting..\n");
				exit(EXIT_FAILURE);
			}
			snprintf(node_name, strlen(token)+1, "%s", &token[0]);
			/*********/

			// there should be no more arguments
			token = strtok(NULL, delim);
			if (token != NULL) {
				printf(RED "Error" NRM ": Command contains extra false arguments (%s)\n Bypassing\n", command_list[i]);
				free(var_name);
				free(node_name);
				var_name = NULL;
				node_name = NULL;
				continue;
			}


			// +9 -> strlen("_DC_") + strlen(".txt") + 1 // 1 for '\0'
			filename = (char *) malloc((strlen(var_name) + strlen(node_name) + 9)*sizeof(char));
			if (filename == NULL) {
				printf("Error. Memory allocation problems. Exiting..\n");
				exit(EXIT_FAILURE);
			}
			sprintf(filename, "%s_DC_%s.txt", node_name, var_name);

			node_fp = fopen(filename, "w");
			if (node_fp == NULL) {
				perror("fopen");
				exit(EXIT_FAILURE);
			}


			/*  Script for plots */

			fp_draw = fopen("draw.sh", "a");

			if (fp_draw == NULL) {
				perror("fopen");
				exit(EXIT_FAILURE);
			}

			/***************/

			if (var_found == 1) {


				for (j=start; j < end + 0.000000001; j = j + jump) {

					if( (idx1+1) != 0 ){
						mna_vector[idx1] += var->value;	// eliminate old value from b vector
						mna_vector[idx1] -= j;			// add new value to b vector
					}
					if( (idx2+1) != 0 ){
						mna_vector[idx2] -= var->value;	// eliminate old value from b vector
						mna_vector[idx2] += j;			// add new value to b vector
					}

					var->value = j;
					if (solver_type == LU_SOLVER) {
						solve_lu_MNA();
						fprintf(node_fp, "%lf\t\t%lf\n", j, node->val);
					}
					else {
						solve_cholesky_MNA();
						fprintf(node_fp, "%lf\t\t%lf\n", j, node->val);
					}

				}


				if( (idx1+1) != 0 ){
					mna_vector[idx1] += var->value;	// eliminate old value from b vector
					mna_vector[idx1] -= var->op_point_val;			// add new value to b vector
				}
				if( (idx2+1) != 0 ){
					mna_vector[idx2] -= var->value;	// eliminate old value from b vector
					mna_vector[idx2] += var->op_point_val;			// add new value to b vector
				}
				var->value = var->op_point_val;

			}
			else {  // it is guaranteed that var_found == 2
				for (j=start; j < end + 0.00000001; j = j + jump) {

					// position in list is unique
					mna_vector[idx1] = j;
					var->value = j;

					if (solver_type == LU_SOLVER) {
						solve_lu_MNA();
						fprintf(node_fp, "%lf\t\t%lf\n", j, node->val);
					}
					else {
						solve_cholesky_MNA();
						fprintf(node_fp, "%lf\t\t%lf\n", j, node->val);
					}
				}

				mna_vector[idx1] = var->op_point_val;
				var->value = var->op_point_val;
			}

			fprintf(fp_draw, "gnuplot -e \"set terminal png size 1024, 1024;");
			fprintf(fp_draw, "set output \\\"%s_DC_%s.png\\\";",node_name,var_name);
			fprintf(fp_draw, "plot \\\"%s\\\" using 1:2 with linespoints;\"\n", filename);

			// redirect sterr to stdout and redirect stdout to /dev/null to avoid viewing xdg-open warnings
			fprintf(fp_draw, "xdg-open \"%s_DC_%s.png\" > /dev/null 2>&1\n",node_name,var_name);


			fclose(node_fp);
			fclose(fp_draw);
			free(node_name);
			free(filename);
			node_name = NULL;
			filename = NULL;
		}

	}

	system("bash draw.sh");

	if(var_name != NULL){

		start = 0;
		end = 0;
		jump = 0;
		free(var_name);
		var_name = NULL;
	}

}


void init_MNA_system() {

	mna_dimension_size = team2_list.size + total_ids - 1;

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


}


// iterate through the node hash table and dump the nodes into a file
void dump_MNA_nodes() {
	FILE *fp = NULL;
	int i;

	fp = fopen("nodes_op_point_all.txt", "w");
	if (fp == NULL) {
		perror("fopen");
		exit(EXIT_FAILURE);
	}

	for (i=1; i < total_ids; i++) {
		fprintf(fp, "%s\t\t%lf\n", id_to_node[i]->name, id_to_node[i]->val);
	}

	fclose(fp);
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

		component_value = team1_list.list[i].value;

		switch(team1_list.list[i].type) {
			case R:
				if ((node_plus_idx + 1) == 0){
					// array[<->][<->] -> +gk
					mna_array[node_minus_idx * mna_dimension_size + node_minus_idx] += 1/component_value;
				}else if ((node_minus_idx + 1) == 0){
					// array[<+>][<+>] -> +gk
					mna_array[node_plus_idx * mna_dimension_size + node_plus_idx] += 1/component_value;
				}else{
					// array[<+>][<+>] -> +gk
					mna_array[node_plus_idx * mna_dimension_size + node_plus_idx] += 1/component_value;

					// array[<->][<->] -> +gk
					mna_array[node_minus_idx * mna_dimension_size + node_minus_idx] += 1/component_value;

					// array[<+>][<->] -> -gk
					mna_array[node_plus_idx * mna_dimension_size + node_minus_idx] -= 1/component_value;

					// array[<->][<+>] -> -gk
					mna_array[node_minus_idx * mna_dimension_size + node_plus_idx] -= 1/component_value;
				}
				break;
			case I:
				// underflow handling
				if ((node_minus_idx + 1) != 0){
					// vector[<->] -> +sk
					mna_vector[node_minus_idx] += component_value;
				}

				// underflow handling
				if ((node_plus_idx + 1) != 0){
					// vector[<+>] -> -sk
					mna_vector[node_plus_idx] -= component_value;
				}
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

		node_plus_idx = team2_list.list[i].node_plus->id - 1;
		node_minus_idx = team2_list.list[i].node_minus->id - 1;
		component_value = team2_list.list[i].value;

		switch(team2_list.list[i].type) {
			case V:
				// vector[k] -> +sk
				// ... where k = (total_ids-1+i) = (n-1+i)
				mna_vector[(total_ids-1+i)] += component_value;

				// it is intended to enter the case L code
			case L:
				// TODO no need to increment as the ith entry is accessed only be one component
				if ((node_minus_idx + 1) != 0){
					// array[k][<->] -> -1
					mna_array[(total_ids-1+i)*mna_dimension_size + node_minus_idx] -= 1;

					// array[<->][k] -> -1
					mna_array[(node_minus_idx)*mna_dimension_size + (total_ids-1+i)] -= 1;
				}

				if ((node_plus_idx + 1) != 0){
					// array[k][<+>] -> +1
					mna_array[(total_ids-1+i)*mna_dimension_size + node_plus_idx] += 1;

					// array[<+>][k] -> +1
					mna_array[(node_plus_idx)*mna_dimension_size + (total_ids-1+i)] += 1;
				}
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


void print_MNA_system(){
	unsigned long i, j;

	printf("\n\n");

	for (i = 0; i < mna_dimension_size; i++){
		if (i < (total_ids - 1)){
			for (j = 0; j < (total_ids -1); j++){
				printf(RED "%.2lf " NRM, mna_array[(i * mna_dimension_size) + j]);
			}

			for (j = (total_ids - 1); j < mna_dimension_size; j++){
				printf(GRN "%.2lf " NRM, mna_array[(i * mna_dimension_size) + j]);
			}
			putchar('\n');
		}else{
			for (j = 0; j < (total_ids -1); j++){
				printf(GRN "%.2lf " NRM, mna_array[(i * mna_dimension_size) + j]);
			}
			for (j = (total_ids - 1); j < mna_dimension_size; j++){
				printf(YEL "%.2lf " NRM, mna_array[(i * mna_dimension_size) + j]);
			}
			putchar('\n');
		}
	}

	printf("\n\n");

	for (i = 0; i < (total_ids - 1); i++){
		printf(BLU "%.4lf\n" NRM, mna_vector[i]);
	}
	for (i = (total_ids - 1); i < mna_dimension_size; i++){
		printf(GRN "%.4lf\n" NRM, mna_vector[i]);
	}

	printf("\n");
}
