#include "lib_poisson1D.h"
#include <time.h>

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

  double* time;

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

  nbpoints = 1000;

  time = (double *) malloc(sizeof(double)*nbpoints);
  for (int n=6; n < nbpoints + 6; n++){
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

	/* Mesure de performance */
	clock_t begin, end;

  	/* LU Factorization */
  	if (IMPLEM == TRF) {
		begin = (double)clock();
    	dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  	}

  	/* LU for tridiagonal matrix  (can replace dgbtrf_) */
  	if (IMPLEM == TRI) {
		begin = (double)clock();
    	dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  	}

  	if (IMPLEM == TRI || IMPLEM == TRF){
    	/* Solution (Triangular) */
    	if (info==0){
      		dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
			end = (double)clock();
			time[n-6] = (end-begin)*100000000.0/CLOCKS_PER_SEC;
    	}
  	}

  	/* It can also be solved with dgbsv */
  	if (IMPLEM == SV) {
		begin = (double)clock();
    	dgbsv_(&lab, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
		end = (double)clock();
		time[n-6] = (double)(end-begin)*1000.0/CLOCKS_PER_SEC;
  	}

  	free(RHS);
  	free(EX_SOL);
  	free(X);
  	free(AB);
  }

  char file[50];
  sprintf(file, "time%d.dat", IMPLEM);

  write_vec(time, &nbpoints, file);
  free(time);
}