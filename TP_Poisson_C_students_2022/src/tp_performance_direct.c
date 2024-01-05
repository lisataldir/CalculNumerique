#include "lib_poisson1D.h"

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  double* relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS=1;
  T0=-5.0;
  T1=5.0;

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

	int nbpoints = 200;

	relres = (double *) malloc(sizeof(double)*200);
  for (int n=3; n < nbpoints + 3; n++){
    la=n-2;

    RHS=(double *) malloc(sizeof(double)*la);
    EX_SOL=(double *) malloc(sizeof(double)*la);
    X=(double *) malloc(sizeof(double)*la);

  
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  	AB = (double *) malloc(sizeof(double)*lab*la);

  	set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

  	ipiv = (int *) calloc(la, sizeof(int));

  	/* LU Factorization */
  	if (IMPLEM == TRF) {
    	dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  	}

  	/* LU for tridiagonal matrix  (can replace dgbtrf_) */
  	if (IMPLEM == TRI) {
    	dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  	}

  	if (IMPLEM == TRI || IMPLEM == TRF){
    	/* Solution (Triangular) */
    	if (info==0){
      	dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    	}else{
      	printf("\n INFO = %d\n",info);
    	}
  	}

  	/* It can also be solved with dgbsv */
  	if (IMPLEM == SV) {
    	dgbsv_(&lab, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
  	}

  	/* Relative forward error */
  	relres[la-1] = relative_forward_error(RHS, EX_SOL, &la);

  	free(RHS);
  	free(EX_SOL);
  	free(X);
  	free(AB);
  }

	write_vec(relres, &nbpoints, "erreur.dat");
	free(relres);
}