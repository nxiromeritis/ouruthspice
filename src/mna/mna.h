#ifndef _MNA_H_
#define _MNA_H_

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include "../csparse/csparse.h"
#include "../lists/lists.h"
#include <math.h>

extern double *mna_array;
extern double *mna_vector;
extern double *default_mna_vector_copy;
extern gsl_vector *default_X_vector_copy;
extern double *old_mna_vector;
extern gsl_vector *gsl_old_x_vector;
extern double *B_vector;
extern double factor;
extern unsigned long mna_dimension_size;

extern gsl_matrix_complex *gsl_complex_mna_array;
extern gsl_vector_complex *gsl_complex_x_vector;
extern gsl_vector_complex *gsl_complex_b_vector;

extern double freq;
#define OMEGA (2*M_PI*freq)
#define PHASOR_REAL(m,f) (m*cos(f))
#define PHASOR_IMG(m,f) (m*sin(f))

extern int sweep;
extern int ac_points;
extern double start_freq;
extern double end_freq;

extern double *G_array;
extern double *C_array;
extern double *C_AC_array;

#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))

#define MIN_ITER		100
#define MAX_ITER		100
#define ITOL_DEFAULT	10e-6
#define EPS_DEFAULT		10e-16
#define LU_SOLVER		0
#define CHOL_SOLVER		1
#define CG_SOLVER		2
#define BI_CG_SOLVER	3
#define TRAPEZOIDAL		0
#define BACKWARD_EULER	1
#define DC_PLOT			0
#define TRAN_PLOT		1
#define AC_PLOT			2 // TODO

extern double itol;

extern byte solver_type;
extern byte tr_method;
extern byte is_sparse;
extern byte is_trans;
extern byte is_ac;

extern int plot_type;

extern double timestep;
extern double end_time;

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
extern double *p_vector;
extern double *q_vector;

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

extern double get_exp_val(ExpInfoT *data, double t);
extern double get_sin_val(SinInfoT *data, double t);
extern double get_pulse_val(PulseInfoT *data, double t);
extern double get_pwl_val(PwlInfoT *data, double t);

void create_trans_MNA_array();
void reset_MNA_array();
void print_C_array();
void print_G_array();

void init_AC_MNA_system();
void fill_AC_MNA_array();
void free_AC_MNA_array();
void print_complex_mna();
#endif
