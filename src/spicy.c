#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

#include "spicy.h"
#include "lists/lists.h"
#include "hashtable/hashtable.h"
#include "cir_parser/cir_parser.h"
#include "mna/mna.h"


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

	init_list_trans();
	print_list_trans();

	//print_list1();
	//print_list2();
	//print_sec_list();

	// print the variables that were set by reading Options
	printf("\nSOLVER: %s\n", ((solver_type==0)?"lu_decomp":
						   ((solver_type == 1)?"cholesky_decomp":
						   ((solver_type == 2)?"cg_solver":
						   ((solver_type == 3)?"bi_cg_solver":"unknown_solver")))));
	printf("ITOL: %e\n", itol);
	printf("%sSPARSE\n", is_sparse?"":"NOT ");
	printf("%sTRANSIENT ANALYSIS\n", is_trans?"":"NO ");
	if (is_trans) {
		printf("TRANSIENT_METHOD: %s\n", (tr_method == TRAPEZOIDAL)?"TRAPEZOIDAL":"BACKWARD_EULER");
	}



	if (is_sparse) {
		init_triplet();
		print_sparse_matrix(triplet_A);
		create_compressed_column();
		print_sparse_matrix(compr_col_A);
	}
	else {
		init_MNA_system();
		fill_MNA_system();
		print_MNA_array();
		print_MNA_vector();
	}




	// TODO ADD if is_sparse == 1 and call new solvers/functions
	switch(solver_type) {
		case LU_SOLVER:
			gsl_x_vector = gsl_vector_alloc(mna_dimension_size);
			gsl_p = gsl_permutation_alloc(mna_dimension_size);
			decomp_lu();
			break;
		case CHOL_SOLVER:
			gsl_x_vector = gsl_vector_alloc(mna_dimension_size);
			decomp_cholesky();
			break;
		// iterative solving method. No need to decompose
		case CG_SOLVER:
		case BI_CG_SOLVER:
			gsl_x_vector = gsl_vector_calloc(mna_dimension_size);
			initialise_iter_methods();
			break;
		default:
			printf(RED "Error uknown solver type specified..\n" NRM);
			exit(EXIT_FAILURE);
	}


	// this function call will generate a unique file for each component
	// that has a transient spec for gnuplot plotting
	test_tran_spec(); //debug
	switch(solver_type) {
		case LU_SOLVER:
			solve_lu();

			//gsl_permutation_free(gsl_p);
			//gsl_vector_free(gsl_x_vector);
			break;
		case CHOL_SOLVER:
			solve_cholesky();

			//gsl_vector_free(gsl_x_vector);
			break;
		case CG_SOLVER:

			solve_CG_iter_method();
			free_gsl_vectors();
			break;
		case BI_CG_SOLVER:

			solve_BI_CG_iter_method();
			/*free_gsl_vectors();*/
			break;
		default:
			break;
	}

	default_X_vector_copy = gsl_vector_alloc(mna_dimension_size);
	default_mna_vector_copy = (double *)malloc(mna_dimension_size*sizeof(double));
	gsl_vector_memcpy(default_X_vector_copy,gsl_x_vector);
	memcpy(default_mna_vector_copy,mna_vector,mna_dimension_size*sizeof(double));


	dump_MNA_nodes();

	print_command_list();
	execute_commands();

	gsl_vector_free(default_X_vector_copy);
	free(default_mna_vector_copy);
	gsl_vector_free(gsl_x_vector);

	if (solver_type == LU_SOLVER)
		gsl_permutation_free(gsl_p);


	free_gsl_vectors();
	free_MNA_system();
	freeHashTable();
	free_lists();
	free_command_list();
	free(Trans_list.list);
	free(Trans_list.k);


	if (css_S)
		cs_sfree(css_S);
	if (csn_N)
		cs_nfree(csn_N);

	if(compr_col_A)
		cs_spfree(compr_col_A);
	if(p_vector)
		free(p_vector);
	if(q_vector)
		free(q_vector);


	return 0;
}
