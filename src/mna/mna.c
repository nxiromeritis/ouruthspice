#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <math.h>

#include "../cir_parser/cir_parser.h"
#include "../spicy.h"
#include "../hashtable/hashtable.h"
#include "mna.h"

// variables regarding the MNA system
double *mna_array = NULL;
double *mna_vector = NULL;
double *default_mna_vector_copy = NULL;
unsigned long mna_dimension_size = 0;

// variable regarding the plot state
int plot_type = DC_PLOT;

// variables used for the complete MNA system
double *G_array = NULL;
double *C_array = NULL;

// timeste and rnd time of transient analysys
double timestep = 0.0;
double end_time = 0.0;

// variables used with sparse matrixes
cs *triplet_A = NULL;
cs *compr_col_A = NULL;
css *css_S = NULL;
csn *csn_N = NULL;

gsl_matrix_view gsl_mna_array;
gsl_vector_view gsl_mna_vector;
gsl_vector *gsl_M_array = NULL;
gsl_vector *gsl_x_vector = NULL;
gsl_vector *gsl_z_vector = NULL;
gsl_vector *gsl_r_vector = NULL;
gsl_vector *gsl_p_vector = NULL;
gsl_vector *gsl_q_vector = NULL;
gsl_vector *gsl_zT_vector = NULL;
gsl_vector *gsl_rT_vector = NULL;
gsl_vector *gsl_pT_vector = NULL;
gsl_vector *gsl_qT_vector = NULL;
gsl_permutation *gsl_p = NULL;
double *p_vector = NULL;
double *q_vector = NULL;

byte solver_type = LU_SOLVER;
byte tr_method = TRAPEZOIDAL;
byte is_sparse = 0;
byte is_trans = 0;

double itol = ITOL_DEFAULT;

void decomp_lu() {
	int s;

	if (is_sparse) {
		css_S = cs_sqr(2, compr_col_A, 0);
		csn_N = cs_lu(compr_col_A, css_S, 1);
		cs_spfree(compr_col_A);
		compr_col_A = NULL;
	}
	else {
		gsl_x_vector = gsl_vector_alloc(mna_dimension_size);

		gsl_mna_array = gsl_matrix_view_array(mna_array, mna_dimension_size, mna_dimension_size);
		gsl_mna_vector = gsl_vector_view_array(mna_vector, mna_dimension_size);
		gsl_p = gsl_permutation_alloc(mna_dimension_size);

		// user's responsibility is this fails
		gsl_linalg_LU_decomp(&gsl_mna_array.matrix, gsl_p, &s);
	}

}

void decomp_cholesky() {

	if (is_sparse) {
		css_S = cs_schol(1, compr_col_A);
		csn_N = cs_chol(compr_col_A, css_S);
		cs_spfree(compr_col_A);
		compr_col_A = NULL;
	}
	else {
		gsl_x_vector = gsl_vector_alloc(mna_dimension_size);

		gsl_mna_array = gsl_matrix_view_array(mna_array, mna_dimension_size, mna_dimension_size);
		gsl_mna_vector = gsl_vector_view_array(mna_vector, mna_dimension_size);

		// user's responsibility if this fails
		gsl_linalg_cholesky_decomp(&gsl_mna_array.matrix);
	}
}


void initialise_iter_methods() {
	unsigned long i,p;
	int k;


	if(is_sparse){

		p_vector = (double *) calloc(mna_dimension_size,sizeof(double));
		if (p_vector == NULL) {
			printf("Error. Memory allocation problems. Exiting..\n");
			exit(EXIT_FAILURE);
		}

		// this is necessary as a helping/temporary vector.
		// the same helper vector can be used in solve qT aswell
		q_vector = (double *) calloc(mna_dimension_size,sizeof(double));
		if (q_vector == NULL) {
			printf("Error. Memory allocation problems. Exiting..\n");
			exit(EXIT_FAILURE);
		}

		gsl_mna_vector = gsl_vector_view_array(mna_vector, mna_dimension_size);

		gsl_M_array = gsl_vector_alloc(mna_dimension_size);
		gsl_vector_set_all(gsl_M_array,1);

		for(k = 0; k < compr_col_A->n; k++){

			for(p = compr_col_A ->p[k]; p < compr_col_A->p[k+1];p++){
				if( k==compr_col_A->i[p]){
					gsl_vector_set(gsl_M_array,k,compr_col_A->x[p]);
				}
			}
		}

		for(k = 0; k<mna_dimension_size;k++){
			printf("M_array[%d] = %lf\n",k,gsl_vector_get(gsl_M_array,k));
		}

	}
	else {

		gsl_mna_array = gsl_matrix_view_array(mna_array, mna_dimension_size, mna_dimension_size);
		gsl_mna_vector = gsl_vector_view_array(mna_vector, mna_dimension_size);

		gsl_M_array = gsl_vector_alloc(mna_dimension_size);
		gsl_vector_set_all(gsl_M_array,1);
		for (i = 0; i < mna_dimension_size; i++){

			// TODO add if is_sparse
			if(  mna_array[(i * mna_dimension_size) + i] != 0 )
				gsl_vector_set(gsl_M_array,i,mna_array[(i * mna_dimension_size) + i]);

		}
		for(k = 0; k<mna_dimension_size;k++){
			printf("M_array[%d] = %lf\n",k,gsl_vector_get(gsl_M_array,k));
		}
	}


	gsl_x_vector = gsl_vector_calloc(mna_dimension_size);
	gsl_z_vector = gsl_vector_calloc(mna_dimension_size);
	gsl_r_vector = gsl_vector_calloc(mna_dimension_size);
	gsl_p_vector = gsl_vector_calloc(mna_dimension_size);
	gsl_q_vector = gsl_vector_calloc(mna_dimension_size);

	if(solver_type == BI_CG_SOLVER){
		gsl_zT_vector = gsl_vector_calloc(mna_dimension_size);
		gsl_rT_vector = gsl_vector_calloc(mna_dimension_size);
		gsl_pT_vector = gsl_vector_calloc(mna_dimension_size);
		gsl_qT_vector = gsl_vector_calloc(mna_dimension_size);
	}
}


// dont forget to free the permutation after the last call of this function
void solve_lu() {
	unsigned long i;
	double *x;



	if (is_sparse) {

		x = (double *) malloc(sizeof(double)*mna_dimension_size);
		if (x == NULL) {
			printf("Error. Memory allocation problems. Exiting..\n");
			exit(EXIT_FAILURE);
		}

		cs_ipvec(csn_N->pinv, mna_vector, x, mna_dimension_size);
		cs_lsolve(csn_N->L, x);
		cs_usolve(csn_N->U, x);
		cs_ipvec(css_S->q, x, mna_vector, mna_dimension_size);

		// mna_vector will contain the solution
		// copy mna_vector to id_to_node here
		for (i=1; i < total_ids; i++) {
			id_to_node[i]->val = mna_vector[i-1];
		}

		free(x);
	}
	else {
		gsl_linalg_LU_solve(&gsl_mna_array.matrix, gsl_p, &gsl_mna_vector.vector, gsl_x_vector);

		// id_to_nodes has the GND located at idx 0
		// although MNA array and vector ignore GND and the first node starts from 0
		for (i=1; i < total_ids; i++) {
			id_to_node[i]->val = gsl_vector_get(gsl_x_vector, i-1);
		}
	}


}


void solve_cholesky() {
	unsigned long i;
	double *x;


	if (is_sparse) {

		x = (double *) malloc(sizeof(double)*mna_dimension_size);
		if (x == NULL) {
			printf("Error. Memory allocation problems. Exiting..\n");
			exit(EXIT_FAILURE);
		}

		cs_ipvec(css_S->pinv, mna_vector, x, mna_dimension_size);
		cs_lsolve(csn_N->L, x);
		cs_ltsolve(csn_N->L, x);
		cs_pvec(css_S->pinv, x, mna_vector, mna_dimension_size);

		// mna_vector will contain the solution
		// copy mna_vector to id_to_node here
		for (i=1; i < total_ids; i++) {
			id_to_node[i]->val = mna_vector[i-1];
		}

		free(x);
	}
	else {
		gsl_linalg_cholesky_solve(&gsl_mna_array.matrix, &gsl_mna_vector.vector, gsl_x_vector);

		// id_to_nodes has the GND located at idx 0
		// although MNA array and vector ignore GND and the first node starts from 0
		for (i=1; i < total_ids; i++) {
			id_to_node[i]->val = gsl_vector_get(gsl_x_vector, i-1);
		}
	}

}



void solve_CG_iter_method() {
	unsigned long i;
	unsigned int iter;
	double normR = 0.0, normB = 0.0;
	double alpha = 0.0, beta = 0.0, tmp = 0.0;
	double rho = 0.0, rho1 = 0.0;

	gsl_vector_set_zero(gsl_x_vector);

	//r = b
	gsl_vector_memcpy(gsl_r_vector, &gsl_mna_vector.vector);

	//r = r - Ax
	// our current initial x is zero so r=r
	/*gsl_blas_dgemv(CblasNoTrans, -1.0, &gsl_mna_array.matrix, gsl_x_vector, 1.0, gsl_r_vector);*/

	for (iter = 0; iter < MAX(MIN_ITER, mna_dimension_size); iter++) {
		normR = gsl_blas_dnrm2(gsl_r_vector);
		normB = gsl_blas_dnrm2(&gsl_mna_vector.vector);
		if (normB == 0.0)
			normB = 1.0;

		if ((normR / normB) <= itol)
			break;

		solve_precond();

		gsl_blas_ddot(gsl_r_vector, gsl_z_vector, &rho);

		if (iter == 0) {
			gsl_vector_memcpy(gsl_p_vector, gsl_z_vector);
		}
		else {
			beta = rho / rho1;
			gsl_vector_scale(gsl_p_vector, beta);
			gsl_vector_add(gsl_p_vector, gsl_z_vector);
		}

		rho1 = rho;

		solve_q();

		gsl_blas_ddot(gsl_p_vector, gsl_q_vector, &tmp);
		alpha = rho / tmp;


		gsl_blas_daxpy(alpha, gsl_p_vector, gsl_x_vector);
		gsl_blas_daxpy((0.0 - alpha), gsl_q_vector, gsl_r_vector);
	}

	for (i=1; i < total_ids; i++) {
		id_to_node[i]->val = gsl_vector_get(gsl_x_vector, i-1);
	}


}


void solve_BI_CG_iter_method() {
	unsigned long i;
	unsigned int iter;
	double normR = 0.0, normB = 0.0;
	double alpha = 0.0, beta = 0.0, omega = 0.0;
	double rho = 0.0, rho1 = 0.0;



	gsl_vector_set_zero(gsl_x_vector);
	//r = b
	gsl_vector_memcpy(gsl_r_vector, &gsl_mna_vector.vector);

	//r = r - Ax
	// our current initial x is zero so r=r
	/*gsl_blas_dgemv(CblasNoTrans, -1.0, &gsl_mna_array.matrix, gsl_x_vector, 1.0, gsl_r_vector);*/

	//rT = r
	gsl_vector_memcpy(gsl_rT_vector, gsl_r_vector);

	for (iter = 0; iter < MAX(MIN_ITER, mna_dimension_size) ; iter++) {
		normR = gsl_blas_dnrm2(gsl_r_vector);
		normB = gsl_blas_dnrm2(&gsl_mna_vector.vector);
		if (normB == 0.0)
			normB = 1.0;

		if ((normR / normB) <= itol)
			break;

		solve_precond();	                              // M*z = r
		Transpose_solve_precond();						  // M(T)*zT = rT

		gsl_blas_ddot(gsl_rT_vector, gsl_z_vector, &rho); //rho = rT . z

		if(fabs(rho) < EPS_DEFAULT){
			printf(RED" 1) i : %d , Bi-CG failed\n"NRM,iter);
			exit(EXIT_FAILURE);
		}

		if (iter == 0) {
			gsl_vector_memcpy(gsl_p_vector, gsl_z_vector);  // p = z
			gsl_vector_memcpy(gsl_pT_vector, gsl_zT_vector);  // pT = zT
		}
		else {
			beta = rho / rho1;
			gsl_vector_scale(gsl_p_vector, beta);        // p = p*beta
			gsl_vector_add(gsl_p_vector, gsl_z_vector);	 // p = p +z

			gsl_vector_scale(gsl_pT_vector, beta);		 // pT = pT*beta
			gsl_vector_add(gsl_pT_vector, gsl_zT_vector); // pT = pT +zT
		}

		rho1 = rho;

		solve_q();                               // q = A . p

		Transpose_solve_q();					// qT = A(T) . pT

		gsl_blas_ddot(gsl_pT_vector, gsl_q_vector, &omega);    // omega = pT . q

		/*printf("OMEGA = %lf\n", omega);*/
		if(fabs(omega) < EPS_DEFAULT){
			printf(RED" 2) i : %d , Bi-CG failed\n"NRM,iter);
			exit(EXIT_FAILURE);
		}

		alpha = rho/omega;

		gsl_blas_daxpy(alpha, gsl_p_vector, gsl_x_vector);			 // x = x + alpha*p
		gsl_blas_daxpy((0.0 - alpha), gsl_q_vector, gsl_r_vector);   // r = r -alpha*q
		gsl_blas_daxpy((0.0 - alpha), gsl_qT_vector, gsl_rT_vector); // rT = rT -alpha*qT
	}

	for (i=1; i < total_ids; i++) {
		id_to_node[i]->val = gsl_vector_get(gsl_x_vector, i-1);
	}


}

void solve_precond() {
	gsl_vector_memcpy(gsl_z_vector,gsl_r_vector);
	gsl_vector_div(gsl_z_vector, gsl_M_array);
}

// q = A.p
void solve_q() {

	unsigned long i;

	if(is_sparse) {

		for( i=0; i<mna_dimension_size; i++){
			p_vector[i] = gsl_vector_get(gsl_p_vector, i);
		}

		// q must be zero so that q = A.p + q = A.p
		cs_gaxpy(compr_col_A,p_vector,q_vector);

		for( i=0; i<mna_dimension_size; i++){
			gsl_vector_set(gsl_q_vector,i,q_vector[i]);
		}

		// restore q
		memset(q_vector, 0, mna_dimension_size*sizeof(double));

	}
	else {
		gsl_blas_dgemv(CblasNoTrans, 1.0, &gsl_mna_array.matrix, gsl_p_vector, 0.0, gsl_q_vector);
	}


	// DEBUG
	/*printf("solve_q (%s)\n", is_sparse?"sparse":"not sparse");*/
	/*for (i=0; i < mna_dimension_size; i++) {*/
		/*printf("q[%lu] = %lf\n", i, gsl_vector_get(gsl_q_vector, i));*/
	/*}*/
}

void Transpose_solve_precond() {
	gsl_vector_memcpy(gsl_zT_vector,gsl_rT_vector);
	gsl_vector_div(gsl_zT_vector, gsl_M_array);        //M_transpose = M
}


void Transpose_solve_q(){

	int i,p;

	if(is_sparse) {

		for( i=0; i<mna_dimension_size; i++){
			p_vector[i] = gsl_vector_get(gsl_pT_vector, i);
		}

		// q must be zero so that q = A.p + q = A.p


		for(i = 0; i < compr_col_A->n; i++){

			for(p = compr_col_A ->p[i]; p < compr_col_A->p[i+1];p++){

				q_vector[i]=q_vector[i]+compr_col_A->x[p]*p_vector[compr_col_A->i[p]];
			}
		}

		for( i=0; i<mna_dimension_size; i++){
			gsl_vector_set(gsl_qT_vector,i,q_vector[i]);
		}

		// restore q
		memset(q_vector, 0, mna_dimension_size*sizeof(double));

	}
	else {
		gsl_blas_dgemv(CblasTrans, 1.0, &gsl_mna_array.matrix, gsl_pT_vector, 0.0, gsl_qT_vector);


	}
}


// get the value of the exp transient function at time t
double get_exp_val(ExpInfoT *data, double t) {
	double res;

	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}
	else if (t < data->td1) {
		res = data->i1;
	}
	else if (t < data->td2) {
		res = data->i1 + (data->i2 - data->i1)*(1.0 - exp(-1*(t - data->td1)/data->tc1));
	}
	else {
		res = data->i1 + (data->i2 - data->i1)*(exp(-1*(t- data->td2)/data->tc2) - exp(-1*(t - data->td1)/data->tc1));
	}
	return res;
}


// get the value of the sin tansient function at time t
double get_sin_val(SinInfoT *data, double t) {
	double res;

	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}
	else if (t < data->td) {
		res = data->i1 + data->ia*sin(2*M_PI/360.0);
	}
	else {
		res = data->i1 + data->ia * sin(2*M_PI * data->fr*(t - data->td) + 2*M_PI * data->ph/360.0)* \
			  exp(-1*(t - data->td)*data->df);
	}
	return res;
}


// get the value of the pulse transient function at time t
double get_pulse_val(PulseInfoT *data, double t) {
	double res;
	double y2, y1;
	double x1;

	// shift the time to take the corresponding value from the first period of the pulse
	// (only when time is bigger than td)

	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}
	else if (t < data->td) {
		res = data->i1;
	}
	else {
		t = t - (int)((t-data->td)/data->per)*data->per;

		if (t < data->td + data->tr) {
			// i1->i2 linearly
			/*x2 = data->td + data->tr;*/
			x1 = data->td;
			y2 = data->i2;
			y1 = data->i1;
			/*res = ((y2-y1)/(x2-x1)) * (t - x1) + y1;*/
			res = ((y2 - y1)/data->tr)* (t - x1) + y1;
		}
		else if (t < data->td + data->tr + data->pw) {
			res = data->i2;
		}
		else if (t < data->td + data->tr + data->pw + data->tf) {
			// i2->i1 linearly
			/*x2 = data->td + data->tr + data->pw + data->tf;*/
			x1 = data->td + data->tr + data->pw;
			y2 = data->i1;
			y1 = data->i2;
			res = ((y2 - y1)/data->tf) * (t - x1) + y1;
		}
		else {
			res = data->i1;
		}
	}

	return res;
}



// get the value of pwl transient function at time t
double get_pwl_val(PwlInfoT *data, double t) {
	double res;
	double y2, y1;
	double x2, x1;
	int idx;

	// negative time -> return zero
	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;;
	}

	// only one tuple is equivalent to a line parallel to x axis (or a constant DC value)
	if (data->total_tuples == 1)
		return data->values[0];

	// t is smaller than the first tuple time (return first value)
	if (t < data->times[0])
		return data->values[0];

	// t is bigger that the last tuple time (return the last value)
	if (t > data->times[data->total_tuples-1])
		return data->values[data->total_tuples-1];


	// find the very FIRST tuple with a time bigger than t
	idx = 0;
	while (1) {
		idx++;
		if (data->times[idx] > t)
			break;
	}

	x2 = data->times[idx];
	x1 = data->times[idx-1];
	y2 = data->values[idx];
	y1 = data->values[idx-1];
	res = ((y2-y1)/(x2-x1))*(t - x1) + y1;

	/*printf("%lf < %lf\n",t,  data->times[idx]);*/
	/*printf("\t(x1,x2),(y1,y2)(res) = (%lf, %lf),(%lf, %lf)(%lf)\n", x1, x2, y1, y2,res);*/

	return res;
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


	// this is a global variable that indicates the length list that contains
	// the commands to be executed. Therefore in this case there are no commands
	if (command_list_len == 0)
		return;


	for (i = 0; i < command_list_len; i++) {
		if (strncmp(command_list[i], ".DC ", 4) == 0) {

			// in this case variable var_name is aready allocated by a previous
			// command and must be freed and re-initialised
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

			// the following plot commands will plot DC graphs
			plot_type = DC_PLOT;

			// the iteration values are set up
			// the next plot command will solve the system
			// and write the results to the file

			// no need to free var_name here. it is needed during the plot command
		}
		if (strncmp(command_list[i], ".TRAN ", 6) == 0) {

			// bypass command name
			token = strtok(command_list[i], delim);
			if (token == NULL) {
				printf(RED "Error" NRM ": Not enough arguments (%s)\n Bypassing..\n", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}

			// time step
			if (parse_double(&timestep, token) == 0) {
				printf(RED "Error" NRM ": Invalid argument value (%s)\n Bypassing..\n", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}

			token = strtok(command_list[i], delim);
			if (token == NULL) {
				printf(RED "Error" NRM ": Not enough arguments (%s)\n Bypassing..\n", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}

			// end time
			if (parse_double(&end_time, token) == 0) {
				printf(RED "Error" NRM ": Invalid argument value (%s)\n Bypassing..\n", command_list[i]);
				free(var_name);
				var_name = NULL;
				continue;
			}

			// the following plot commands will plot TRAN graphs
			plot_type = TRAN_PLOT;
		}
		if ((plot_type == DC_PLOT) &&  
			((strncmp(command_list[i], ".PRINT ", 7) == 0) || (strncmp(command_list[i], ".PLOT ", 6) == 0))) {
			// TODO: multiple nodes in one line of .print or .plot command
			// Make the parser execute the options command and create a list of dc-print commands


			// I or V was not found in the previous DC command  (TODO what should we do with TRAN?)
			// TRANSIENT ANALYSIS?
			if(var_found == 0 ) {
				printf("Bypassing PLOT command. No DC command before or unknown DC source\n");
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

			// check if the node is written correctly in command (syntax check)
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

			// checks were successsful. store node name into a variable
			// +1 for '\0', -1 for 'v', -1 for '(' and -1 for ')'
			node_name = malloc( (strlen(token) - 2)*sizeof(char));
			if (node_name == NULL) {
				printf("Error. Memory allocation problems. Exiting..\n");
				exit(EXIT_FAILURE);
			}
			snprintf(node_name, strlen(token)-2, "%s", &token[2]);


			// search for the node in the hashtable
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
			// now that we have successfully located the node in the hash table and aquired
			// its pointer we are rereading the name including the parantheses. This is used only for
			// file name generation
			node_name = realloc(node_name,(strlen(token)+1)*sizeof(char));
			if (node_name == NULL) {
				printf("Error. Memory allocation problems. Exiting..\n");
				exit(EXIT_FAILURE);
			}
			snprintf(node_name, strlen(token)+1, "%s", &token[0]);
			/***********/


			// there should be no more arguments (syntax check)
			token = strtok(NULL, delim);
			if (token != NULL) {
				printf(RED "Error" NRM ": Command contains extra false arguments (%s)\n Bypassing\n", command_list[i]);
				free(var_name);
				free(node_name);
				var_name = NULL;
				node_name = NULL;
				continue;
			}


			// reminder: variable var_name is the I or V that changes value during the DC or TRAN analysis
			// strlen(node_name) +9 = strlen(node_name) + strlen("_DC_") + strlen(".txt") + 1 for '\0'
			filename = (char *) malloc((strlen(var_name) + strlen(node_name) + 9)*sizeof(char));
			if (filename == NULL) {
				printf("Error. Memory allocation problems. Exiting..\n");
				exit(EXIT_FAILURE);
			}

			// TODO add TRAN instead of DC if we are plotting for TRAN analysis
			sprintf(filename, "%s_DC_%s.txt", node_name, var_name);

			node_fp = fopen(filename, "w");
			if (node_fp == NULL) {
				perror("fopen");
				exit(EXIT_FAILURE);
			}


			/* GNUPLOT script */
			fp_draw = fopen("draw.sh", "a");

			if (fp_draw == NULL) {
				perror("fopen");
				exit(EXIT_FAILURE);
			}


			// TODO handle transient (TRAN)  analysis after this point

			// at this point it is guaranteed that var_found will be either 1 or 2
			// var_found == 1 -> I variations
			// var_found == 2 -> V variations
			if (var_found == 1) {

				for (j=start; j < end + 0.000000001; j = j + jump) {

					if ( (is_sparse) && ((solver_type == 0)||( solver_type ==1)) ) {
						// restore default b vector values
						memcpy(mna_vector, default_mna_vector_copy, mna_dimension_size*sizeof(double));
						if( (idx1+1) != 0 ) {
							mna_vector[idx1] += var->op_point_val;
							mna_vector[idx1] -= j;			// add new value to b vector
						}
						if( (idx2+1) != 0 ) {
							mna_vector[idx2] -= var->op_point_val;
							mna_vector[idx2] += j;			// add new value to b vector
						}
					}
					else {
						if( (idx1+1) != 0 ){
							mna_vector[idx1] += var->value;	// eliminate old value from b vector
							mna_vector[idx1] -= j;			// add new value to b vector
						}
						if( (idx2+1) != 0 ){
							mna_vector[idx2] -= var->value;	// eliminate old value from b vector
							mna_vector[idx2] += j;			// add new value to b vector
						}
					}

					var->value = j;

					switch(solver_type) {
						case LU_SOLVER:
							solve_lu();
							break;
						case CHOL_SOLVER:
							solve_cholesky();
							break;
						case CG_SOLVER:
							solve_CG_iter_method();
							break;
						case BI_CG_SOLVER:
							solve_BI_CG_iter_method();
							break;
						default:
							break;
					}
					/*fprintf(node_fp, "%lf\t\t%lf\n", j, node->val);*/
					fprintf(node_fp, "%lf\t\t%e\n", j, node->val);
				}

				if ( (is_sparse) && ((solver_type == 0)||( solver_type ==1)) ){
					// restore default b vector values
					memcpy(mna_vector, default_mna_vector_copy, mna_dimension_size*sizeof(double));
				}
				else {
					if( (idx1+1) != 0 ){
						mna_vector[idx1] += var->value;	// eliminate old value from b vector
						mna_vector[idx1] -= var->op_point_val;			// add new value to b vector
					}
					if( (idx2+1) != 0 ){
						mna_vector[idx2] -= var->value;	// eliminate old value from b vector
						mna_vector[idx2] += var->op_point_val;			// add new value to b vector
					}
				}

				var->value = var->op_point_val;

			}
			else {  // it is guaranteed that var_found == 2
				for (j=start; j < end + 0.00000001; j = j + jump) {

					// position in list is unique
					if ( (is_sparse) && ((solver_type == 0)||( solver_type ==1)) ){
						// restore default b vector values
						memcpy(mna_vector, default_mna_vector_copy, mna_dimension_size*sizeof(double));
					}

					mna_vector[idx1] = j;
					var->value = j;

					switch(solver_type) {
						case LU_SOLVER:
							solve_lu();
							break;
						case CHOL_SOLVER:
							solve_cholesky();
							break;
						case CG_SOLVER:
							solve_CG_iter_method();
							break;
						case BI_CG_SOLVER:
							solve_BI_CG_iter_method();
							break;
						default:
							break;
					}
					/*fprintf(node_fp, "%lf\t\t%lf\n", j, node->val);*/
					fprintf(node_fp, "%lf\t\t%e\n", j, node->val);
				}

				if ( (is_sparse) && ((solver_type == 0)||( solver_type ==1)) ){
					// restore default b vector values
					memcpy(mna_vector, default_mna_vector_copy, mna_dimension_size*sizeof(double));
				}
				else
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

		if ((plot_type == TRAN_PLOT) &&  
			((strncmp(command_list[i], ".PRINT ", 7) == 0) || (strncmp(command_list[i], ".PLOT ", 6) == 0))) {
			
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

void create_trans_MNA_array() {
	gsl_matrix_view gsl_mna = gsl_matrix_view_array(mna_array, mna_dimension_size, mna_dimension_size);
	gsl_matrix_view gsl_G_array = gsl_matrix_view_array(G_array, mna_dimension_size, mna_dimension_size);
	gsl_matrix_view gsl_C_array = gsl_matrix_view_array(C_array, mna_dimension_size, mna_dimension_size);

	memset(mna_array, 0, (mna_dimension_size * mna_dimension_size * sizeof(double)));
	gsl_matrix_add(&gsl_mna.matrix, &gsl_C_array.matrix);
	if (tr_method == BACKWARD_EULER) {
		gsl_matrix_scale (&gsl_mna.matrix, (1 / timestep));
	} else { // tr_method ==TRAPEZOIDAL
		gsl_matrix_scale (&gsl_mna.matrix, (2 / timestep));
	}

	gsl_matrix_add(&gsl_mna.matrix, &gsl_G_array.matrix);
}

void reset_MNA_array() {
	memcpy(mna_array, G_array, (mna_dimension_size * mna_dimension_size * sizeof(double)));
}

void print_sparse_matrix(cs *A) {
	int i;
	int p_bound;

	if (A == NULL)
		return;

	printf("%s:\n", (A->nz == -1)?"Compressed Column Form":"Triplet Form");

	printf("\tnzmax: %d\n"
		   "\tm    : %d\n"
		   "\tn    : %d\n"
		   "\tnz   : %d\n",
		   A->nzmax, A->m, A->n, A->nz);


	if (A->nz == -1)
		p_bound = A->n + 1;
	else
		p_bound = A->nzmax;

	printf("\ti = { ");
	for (i = 0; i < A->nzmax; i++) {
		printf("%d%s ", A->i[i], (i == (A->nz)-1)?"":",");
	}
	printf("}\n");

	printf("\tp = { ");
	for (i = 0; i < p_bound; i++) {
		printf("%d%s ", A->p[i], (i == p_bound-1)?"":",");
	}
	printf("}\n");

	printf("\tx = { ");
	for (i = 0; i < A->nzmax; i++) {
		printf("%.2lf%s ", A->x[i], (i == (A->nz)-1)?"":",");
	}
	printf("}\n\n");

}


void create_compressed_column() {
	//printf("\n\n\n\n\n\n\\n\n\n\n\n");
	compr_col_A = cs_compress(triplet_A);
	cs_dupl(compr_col_A);
	cs_spfree(triplet_A);
	triplet_A = NULL;
}




void init_triplet() {
	unsigned long i;
	unsigned long node_plus_idx;
	unsigned long node_minus_idx;
	double component_value;
	int nz = 0;
	int ret;

	// this global is used inside cs_spalloc
	mna_dimension_size = team2_list.size + total_ids - 1;

	mna_vector = (double *) calloc(mna_dimension_size, sizeof(double));
	if (mna_vector == NULL) {
		printf("Error. Memory allocation problems. Exiting..\n");
		exit(EXIT_FAILURE);
	}

	default_mna_vector_copy = (double *) calloc(mna_dimension_size, sizeof(double));
	if (default_mna_vector_copy == NULL) {
		printf("Error. Memory allocation problems. Exiting..\n");
		exit(EXIT_FAILURE);
	}

	// create triplet with initial nz set to 1
	triplet_A = cs_spalloc(mna_dimension_size, mna_dimension_size, 1, 1, 1);
	if (triplet_A == NULL) {
		printf("Error. Memory allocation problems. Exiting..\n");
		exit(EXIT_FAILURE);
	}


	for (i=0; i < team1_list.size; i++) {

		node_plus_idx = team1_list.list[i].node_plus->id - 1;
		node_minus_idx = team1_list.list[i].node_minus->id - 1;
		component_value = team1_list.list[i].value;

		switch(team1_list.list[i].type) {
			case R:
				if ( ((node_plus_idx + 1) == 0) && ((node_minus_idx + 1 != 0)) ) {
					nz++;

					// do sprealloc
					ret = cs_sprealloc(triplet_A, nz);
					if (ret == 0) {
						printf("Error. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}

					triplet_A->i[nz-1] = node_minus_idx;
					triplet_A->p[nz-1] = node_minus_idx;
					triplet_A->x[nz-1] = 1/component_value;

					// array[<->][<->] -> +gk
					/*mna_array[node_minus_idx * mna_dimension_size + node_minus_idx] += 1/component_value;*/
				} else if ( ((node_minus_idx + 1) == 0) && ((node_plus_idx + 1) != 0) ) {
					nz++;

					// do sprealloc
					ret = cs_sprealloc(triplet_A, nz);
					if (ret == 0) {
						printf("Error. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}

					triplet_A->i[nz-1] = node_plus_idx;
					triplet_A->p[nz-1] = node_plus_idx;
					triplet_A->x[nz-1] = 1/component_value;

					// array[<+>][<+>] -> +gk
					/*mna_array[node_plus_idx * mna_dimension_size + node_plus_idx] += 1/component_value;*/
				} else if ( ((node_plus_idx +1) != 0) && ((node_minus_idx + 1) != 0) ) {
					nz += 4;

					// do sprealloc
					ret = cs_sprealloc(triplet_A, nz);
					if (ret == 0) {
						printf("Error. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}

					triplet_A->i[nz-4] = node_plus_idx;
					triplet_A->p[nz-4] = node_plus_idx;
					triplet_A->x[nz-4] = 1/component_value;
					// array[<+>][<+>] -> +gk
					/*mna_array[node_plus_idx * mna_dimension_size + node_plus_idx] += 1/component_value;*/

					triplet_A->i[nz-3] = node_minus_idx;
					triplet_A->p[nz-3] = node_minus_idx;
					triplet_A->x[nz-3] = 1/component_value;
					// array[<->][<->] -> +gk
					/*mna_array[node_minus_idx * mna_dimension_size + node_minus_idx] += 1/component_value;*/

					triplet_A->i[nz-2] = node_plus_idx;
					triplet_A->p[nz-2] = node_minus_idx;
					triplet_A->x[nz-2] = -1/component_value;
					// array[<+>][<->] -> -gk
					/*mna_array[node_plus_idx * mna_dimension_size + node_minus_idx] -= 1/component_value;*/

					triplet_A->i[nz-1] = node_minus_idx;
					triplet_A->p[nz-1] = node_plus_idx;
					triplet_A->x[nz-1] = -1/component_value;
					// array[<->][<+>] -> -gk
					/*mna_array[node_minus_idx * mna_dimension_size + node_plus_idx] -= 1/component_value;*/
				}
				// else do nothing... R nodes are both connected to GND (R 0 0 <num>)
				// TODO fugure out what to do if R has 0 value.
				// (does it have to be in list grp2 and be considered a V source? [like L]
				// or completely ignored??)
				break;
			case I:
				// underflow handling
				if ((node_minus_idx + 1) != 0){
					// vector[<->] -> +sk
					default_mna_vector_copy[node_minus_idx] += component_value;
					mna_vector[node_minus_idx] += component_value;
				}

				// underflow handling
				if ((node_plus_idx + 1) != 0){
					// vector[<+>] -> -sk
					default_mna_vector_copy[node_plus_idx] -= component_value;
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
				default_mna_vector_copy[(total_ids-1+i)] += component_value;
				mna_vector[(total_ids-1+i)] += component_value;

				// NOTE: NO BREAK HERE!!!!!
				// it is intended to enter the case L code
			case L:
				// TODO no need to increment as the ith entry is accessed only be one component
				if ((node_minus_idx + 1) != 0){
					nz += 2;

					// do sprealloc
					ret = cs_sprealloc(triplet_A, nz);
					if (ret == 0) {
						printf("Error. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}

					triplet_A->i[nz-2] = total_ids-1+i;
					triplet_A->p[nz-2] = node_minus_idx;
					triplet_A->x[nz-2] = -1;
					// array[k][<->] -> -1
					/*mna_array[(total_ids-1+i)*mna_dimension_size + node_minus_idx] -= 1;*/

					triplet_A->i[nz-1] = node_minus_idx;
					triplet_A->p[nz-1] = total_ids-1+i;
					triplet_A->x[nz-1] = -1;
					// array[<->][k] -> -1
					/*mna_array[(node_minus_idx)*mna_dimension_size + (total_ids-1+i)] -= 1;*/
				}

				if ((node_plus_idx + 1) != 0){
					nz += 2;

					// do sprealloc
					ret = cs_sprealloc(triplet_A, nz);
					if (ret == 0) {
						printf("Error. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}

					triplet_A->i[nz-2] = total_ids-1+i;
					triplet_A->p[nz-2] = node_plus_idx;
					triplet_A->x[nz-2] = 1;
					// array[k][<+>] -> +1
					/*mna_array[(total_ids-1+i)*mna_dimension_size + node_plus_idx] += 1;*/

					triplet_A->i[nz-1] = node_plus_idx;
					triplet_A->p[nz-1] = total_ids-1+i;
					triplet_A->x[nz-1] = 1;
					// array[<+>][k] -> +1
					/*mna_array[(node_plus_idx)*mna_dimension_size + (total_ids-1+i)] += 1;*/
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

	printf("Before Assignment: (nz, nzmax) = (%d, %d)\n", triplet_A->nz, triplet_A->nzmax);
	triplet_A->nz = nz;
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

	if (is_trans) {
		// G array dimensions: ((n-1) + m2)x((n-1) + m2)
		G_array = (double *)calloc(mna_dimension_size * mna_dimension_size, sizeof(double));
		if (G_array == NULL) {
			printf("Error. Memory allocation problems. Exiting..\n");
			exit(EXIT_FAILURE);
		}

		// C array dimensions: ((n-1) + m2)x((n-1) + m2)
		C_array = (double *)calloc(mna_dimension_size * mna_dimension_size, sizeof(double));
		if (C_array == NULL) {
			printf("Error. Memory allocation problems. Exiting..\n");
			exit(EXIT_FAILURE);
		}
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

	// write grounding with name G, and value 0
	// changing the groudning name in main leads to errors
	fprintf(fp, "G\t\t0.00000e+00\n");
	for (i=1; i < total_ids; i++) {
		/*fprintf(fp, "%s\t\t%lf\n", id_to_node[i]->name, id_to_node[i]->val);*/
		fprintf(fp, "%s\t\t%.5e\n", id_to_node[i]->name, id_to_node[i]->val);
	}

	fclose(fp);
}



void fill_MNA_system() {
	unsigned long i;
	unsigned long node_plus_idx;
	unsigned long node_minus_idx;
	double component_value;
	void *dest;

	// iterate through the Group1 List and initialise the MNA system
	for (i=0; i < team1_list.size; i++) {


		// TODO TODO TODO NOW
		// handle GROUNDING.. node_minus_idx causes undeflow... add if statements to bypass GND


		node_plus_idx = team1_list.list[i].node_plus->id - 1;
		node_minus_idx = team1_list.list[i].node_minus->id - 1;

		component_value = team1_list.list[i].value;

		switch(team1_list.list[i].type) {
			case R:
				if ( ((node_plus_idx + 1) == 0) && ((node_minus_idx + 1 != 0)) ) {
					// array[<->][<->] -> +gk
					mna_array[node_minus_idx * mna_dimension_size + node_minus_idx] += 1/component_value;
				} else if ( ((node_minus_idx + 1) == 0) && ((node_plus_idx + 1) != 0) ) {
					// array[<+>][<+>] -> +gk
					mna_array[node_plus_idx * mna_dimension_size + node_plus_idx] += 1/component_value;
				} else if ( ((node_plus_idx +1) != 0) && ((node_minus_idx + 1) != 0) ) {
					// array[<+>][<+>] -> +gk
					mna_array[node_plus_idx * mna_dimension_size + node_plus_idx] += 1/component_value;

					// array[<->][<->] -> +gk
					mna_array[node_minus_idx * mna_dimension_size + node_minus_idx] += 1/component_value;

					// array[<+>][<->] -> -gk
					mna_array[node_plus_idx * mna_dimension_size + node_minus_idx] -= 1/component_value;

					// array[<->][<+>] -> -gk
					mna_array[node_minus_idx * mna_dimension_size + node_plus_idx] -= 1/component_value;
				}
				// else do nothing... R nodes are both connected to GND
				// TODO fugure out what to do if R has 0 value.
				// (does it have to be in list grp2 and be considered a V source? [like L]
				// or completely ignored??)
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
			case C:
				// ignored at DC analysis
				if (is_trans) {
					if ( ((node_plus_idx + 1) == 0) && ((node_minus_idx + 1 != 0)) ) {
						// C_array[<->][<->] -> +Ck
						C_array[node_minus_idx * mna_dimension_size + node_minus_idx] += component_value;
					} else if ( ((node_minus_idx + 1) == 0) && ((node_plus_idx + 1) != 0) ) {
						// C_array[<+>][<+>] -> +Ck
						C_array[node_plus_idx * mna_dimension_size + node_plus_idx] += component_value;
					} else if ( ((node_plus_idx +1) != 0) && ((node_minus_idx + 1) != 0) ) {
						// C_array[<+>][<+>] -> +Ck
						C_array[node_plus_idx * mna_dimension_size + node_plus_idx] += component_value;

						// C_array[<->][<->] -> +Ck
						C_array[node_minus_idx * mna_dimension_size + node_minus_idx] += component_value;

						// C_array[<+>][<->] -> -Ck
						C_array[node_plus_idx * mna_dimension_size + node_minus_idx] -= component_value;

						// C_array[<->][<+>] -> -Ck
						C_array[node_minus_idx * mna_dimension_size + node_plus_idx] -= component_value;
					}
				}
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

				if (is_trans) {
					// C_array[k][k] -> Lk
					C_array[(total_ids-1+i)*mna_dimension_size + (total_ids-1+i)] += component_value;
				}
				// vector[k] -> 0
				// ... where k = (total_ids-1+i) = (n-1+i)
				break;

			// these cases are here because of the optional .spice I R C field G2
			// described into Najm`s book: "CIRCUIT SIMULATION" 
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

	if (is_trans) {
		dest = memcpy(G_array, mna_array, ((mna_dimension_size * mna_dimension_size) * sizeof(double)));
		if (dest != G_array) {
			printf("Unable to copy MNA array to G array\n");
			exit(EXIT_FAILURE);
		}
	}

}

void free_gsl_vectors(){


	gsl_vector_free(gsl_M_array);
	gsl_M_array = NULL;
	gsl_vector_free(gsl_x_vector);
	gsl_x_vector = NULL;
	gsl_vector_free(gsl_z_vector);
	gsl_z_vector = NULL;
	gsl_vector_free(gsl_r_vector);
	gsl_r_vector = NULL;
	gsl_vector_free(gsl_p_vector);
	gsl_p_vector = NULL;
	gsl_vector_free(gsl_q_vector);
	gsl_q_vector = NULL;

	if(solver_type == BI_CG_SOLVER){
		gsl_vector_free(gsl_zT_vector);
		gsl_zT_vector = NULL;
		gsl_vector_free(gsl_rT_vector);
		gsl_rT_vector = NULL;
		gsl_vector_free(gsl_pT_vector);
		gsl_pT_vector = NULL;
		gsl_vector_free(gsl_qT_vector);
		gsl_qT_vector = NULL;
	}
}


void free_MNA_system() {
	free(mna_array);
	mna_array = NULL;
	free(mna_vector);
	mna_vector = NULL;

	if (default_mna_vector_copy)
		free(default_mna_vector_copy);

	if (is_trans) {
		free(G_array);
		G_array = NULL;
		free(C_array);
		C_array = NULL;
	}
}



void print_MNA_vector() {
	unsigned long i;

	for (i = 0; i < (total_ids - 1); i++){
		printf(BLU "%.4lf\n" NRM, mna_vector[i]);
	}
	for (i = (total_ids - 1); i < mna_dimension_size; i++){
		printf(GRN "%.4lf\n" NRM, mna_vector[i]);
	}

	printf("\n");
}



void print_MNA_array(){
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
}
