/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
    for (int i = 0; i < *la; i++) {
        AB[(*kv) + 1 + i*(*lab)] = 2.0;
        if (i != 0) AB[(*kv) + i*(*lab)] = -1.0;
        if (i != (*la) - 1) AB[(*kv) + 2 + i*(*lab)] = -1.0;
    }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  for (int i = 0; i < *la; i++) AB[(*kv) + 1 + i*(*lab)] = 1.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  RHS[0] = *BC0;
  for (int i=1; i < *la - 1; i++){
    RHS[i] = 0.0;
  }
  RHS[*la - 1] = *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  for (int i=0; i < *la; i++) EX_SOL[i] = *BC0 + X[i]*(*BC1 - *BC0);
}  

void set_grid_points_1D(double* x, int* la){
  double h = 1.0/(1.0*(*la + 1));
  for (int i=0; i < *la; i++){
    x[i] = (i+1) * h;
  }
}

// ici implémenter erreur cblas_dnrm2(la, x, 1)
double relative_forward_error(double* x, double* y, int* la){
  return 0;
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}