#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

#include "spicy.h"
#include "lists.h"
#include "hashtable.h"
#include "cir_parser.h"
#include "mna.h"


int main(int argc, char *argv[]) {
	char gnd_name[2] = "0";
	char filename[BUF_MAX];
	unsigned long components_num;

	if (argc != 2) {
		printf("Error. Invalid number of arguments..\n");
		printf("Use: %s <filename>\n", argv[0]);
		return 1;
	}

	if (strlen(argv[1]) > BUF_MAX-1) {
		printf("Error. Filename too long. Exiting..\n");
		return 1;
	}


	strcpy(filename, argv[1]);

	components_num = get_components_num(filename);
	ht_init(components_num >> 1);

	// add grounding into the hash table (we want it to be reserved)
	ht_put(gnd_name, 0);

	init_lists();

	printf("Total number of components: %lu\n\n", components_num);
	parse_cir(filename);

	printHastable();
	print_id_list();
	//print_list1();
	//print_list2();
	//print_sec_list();

	// print the variables that were set by reading Options
	printf("\nSOLVER: %s\n", ((solver_type==0)?"lu_decomp":
						   ((solver_type == 1)?"cholesky_decomp":
						   ((solver_type == 2)?"cg_solver":
						   ((solver_type == 3)?"bi_cg_solver":"unknown_solver")))));
	printf("ITOL: %e\n", itol);

	init_MNA_system();
	fill_MNA_system();
	print_MNA_system();

	switch(solver_type) {
		case LU_SOLVER:
			decomp_lu_MNA();
			break;
		case CHOL_SOLVER:
			decomp_cholesky_MNA();
			break;
		// iterative solving method. No need to decompose
		case CG_SOLVER:
		case BI_CG_SOLVER:
			break;
		default:
			printf(RED "Error uknown solver type specified..\n" NRM);
			exit(EXIT_FAILURE);
	}



	print_command_list();
	execute_commands();


	// TODO: Write another switch case and call possible iterative functions
	if (solver_type == LU_SOLVER) {
		solve_lu_MNA();
		dump_MNA_nodes();
	}
	else {	// cholesky
		solve_cholesky_MNA();
		dump_MNA_nodes();
	}


	if (solver_type == LU_SOLVER) {
		gsl_permutation_free(gsl_p);
	}
	free_MNA_system();
	freeHashTable();
	free_lists();
	free_command_list();

	return 0;
}
