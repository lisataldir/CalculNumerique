/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int tmp;
  for (int i = 0; i < *la; i++) {
      tmp = i*(*lab);
      /* diagonale */
      AB[(*kv) + 1 + tmp] = 2.0;
      /* sous-diagonale */
      if (i != 0) AB[(*kv) + tmp] = -1.0;
      /* sur-diagonale */
      if (i != (*la) - 1) AB[(*kv) + 2 + tmp] = -1.0;
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

double relative_forward_error(double* x, double* y, int* la){
  double* A = (double*)malloc(sizeof(double)* (*la));
  int lab = 1;
  int kv = 0;

  set_GB_operator_colMajor_poisson1D_Id(A, &lab, la, &kv);
  A[0] = 1.0;

  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, 0, 0, -1.0, A, lab, x, 1, 1.0, y, 1);
  
  double norm_res = cblas_dnrm2(*la, y, 1); 
  double norm_x = cblas_dnrm2(*la, x, 1);       
  free(A);

  return norm_res/norm_x; 
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){

  AB[(*kl) + 2] = AB[(*kl) + 2]/AB[(*kl) + 1];
  ipiv[0] = 1;

  int tmp;

  for (int i = 1; i < *la; i++) {

    tmp = i*(*lab);

    AB[(*kl) + 1 + tmp] = AB[(*kl) + 1 + tmp] - AB[(*kl) + tmp]*AB[(*kl) + 2 + tmp - *lab] ;

    if (i != (*la) - 1) AB[(*kl) + 2 + tmp] = AB[(*kl) + 2 + tmp] / AB[(*kl) + 1 + tmp];

    ipiv[i] = i + 1;
  }
  *info = 0;
  return *info;
}

// On suppose que AB est en col major
void set_CSR_operator_poisson1D(double* AB, int *lab, int *la, double* AA, int* JA, double* IA){

  int iter = 0;
  int k = 0;
  int tmp;

  for (int i = 0; i < (*la); i++) {

    tmp = 0;

    for (int j = 0; j < (*lab); j++){

      if(AB[i*(*lab) + j] != 0){

        AA[iter] = AB[i*(*lab) + j];
        JA[iter] = j;

        if (tmp = 0){

          IA[k] = AB[i*(*lab) + j];
          k++;
          tmp = 1;
          printf("k : %d\n", k);

        }
        printf("iter : %d\n", iter);
        iter++;
      }
    }  
  }
  printf("ok\n");
}

void set_CSC_operator_poisson1D(double* AB, int *lab, int *la, double* AA, int* JA, double* IA){
  
  int iter = 0;
  int k = 0;
  int tmp;

  for (int j = 0; j < (*lab); j++) {

    tmp = 0;

    for (int i = 0; i < (*la); i++){

      if(AB[i*(*lab) + j] != 0){

        AA[iter] = AB[i*(*lab) + j];
        JA[iter] = i;

        if (tmp = 0){

          IA[k] = AB[i*(*lab) + j];
          k++;
          tmp = 1;

        }

        iter++;
      }
    }  
  }
  printf("ok\n");
}