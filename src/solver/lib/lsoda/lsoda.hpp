
#ifndef SOLVER_LIB_LSODA_LSODA_HPP_
#define SOLVER_LIB_LSODA_LSODA_HPP_

extern "C" {
void lsoda(void (*)(double t, double *y, double *ydot, void *data), int neq,
           double *y, double *t, double tout, int itol, double *rtol,
           double *atol, int itask, int *istate, int iopt, int jt, int iwork1,
           int iwork2, int iwork5, int iwork6, int iwork7, int iwork8,
           int iwork9, double rwork1, double rwork5, double rwork6,
           double rwork7, void *_data);
}

#endif
