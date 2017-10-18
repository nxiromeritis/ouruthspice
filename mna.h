#ifndef _MNA_H_
#define _MNA_H_

extern double *mna_array;
extern double *mna_vector;
extern unsigned long mna_dimension_size;

extern void init_MNA_system();
extern void fill_MNA_system();
extern void free_MNA_system();
#endif
