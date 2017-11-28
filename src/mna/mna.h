#ifndef _MNA_H_
#define _MNA_H_

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include "../csparse/csparse.h"

extern double *mna_array;
extern double *mna_vector;
extern double *default_mna_vector_copy;
extern unsigned long mna_dimension_size;

#define MAX_ITERATIONS 	100	//Placeholder
#define ITOL_DEFAULT	10e-6
#define EPS_DEFAULT 	10e-16
#define LU_SOLVER		0
#define CHOL_SOLVER		1
#define CG_SOLVER		2
#define BI_CG_SOLVER	3

extern double itol;

extern byte solver_type;
extern byte is_sparse;

extern cs *triplet_A;
extern cs *compr_col_A;
extern css *css_S;
extern csn *csn_N;

extern gsl_matrix_view gsl_mna_array;
extern gsl_vector_view gsl_mna_vector;
extern gsl_vector *gsl_M_array;
extern gsl_vector *gsl_x_vector;
extern gsl_vector *gsl_z_vector;
extern gsl_vector *gsl_r_vector;
extern gsl_vector *gsl_p_vector;
extern gsl_vector *gsl_q_vector;
extern gsl_vector *gsl_zT_vector;
extern gsl_vector *gsl_rT_vector;
extern gsl_vector *gsl_pT_vector;
extern gsl_vector *gsl_qT_vector;
extern gsl_permutation *gsl_p;

// functions for sparse matrixes
extern void init_triplet();
extern void create_compressed_column();
extern void print_sparse_matrix(cs *A);

extern void init_MNA_system();
extern void fill_MNA_system();
extern void free_MNA_system();
extern void print_MNA_array();
extern void print_MNA_vector();
extern void execute_commands();


extern void dump_MNA_nodes();

extern void decomp_lu();
extern void decomp_cholesky();
extern void initialise_iter_methods();
extern void solve_lu();
extern void solve_cholesky();
extern void solve_CG_iter_method();
extern void solve_precond();
extern void solve_q();
extern void free_gsl_vectors();
extern void solve_BI_CG_iter_method();
extern void Transpose_solve_precond();
extern void Transpose_solve_q();

#endif
