/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
    for (int k = 0; k < (*la); k++){
        double thetak = (k+1)*M_PI/((*la)+1);
        eigval[k] = -2.0*cos(thetak) + 2.0;
    }
}

// Comme cos décroissant sur [0, pi], -2*cos + 2 croissant sur [0, pi]
// donc plus grande vp = -2*cos(n*pi/(n+1)) + 2
// et plus petite vp = -2*cos(pi/(n+1)) + 2 où n = *la

double eigmax_poisson1D(int *la){
  return -2.0*cos((*la)*M_PI/(*la+1)) + 2.0;
}

double eigmin_poisson1D(int *la){
  return -2.0*cos(M_PI/(*la+1)) + 2.0;
}

double richardson_alpha_opt(int *la){
  return 2.0/(eigmax_poisson1D(la) + eigmin_poisson1D(la));
}

// Formule : x = x + alpha*(b-A*x)
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    double* b = (double*)malloc(sizeof(double)*(*la));
    double normb = 1.0/cblas_dnrm2(*la, RHS, 1);
    cblas_dcopy(*la, RHS, 1, b, 1);

    // b = beta*b + alpha*A*x, ici on calcule b - A*x
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, b, 1);

    // Calcul du résidu
    resvec[*nbite] = cblas_dnrm2(*la, b, 1)*normb;

    while(*nbite < (*maxit) && resvec[*nbite] > (*tol)){
      // y = alpha*x + y, ici x vaut b - A*x et a été stocké dans b
      cblas_daxpy(*la, *alpha_rich, b, 1, X, 1);
      cblas_dcopy(*la, RHS, 1, b, 1);

      // étape suivante
      cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1, AB, *lab, X, 1, 1, b, 1);

      *nbite = *nbite + 1;
      resvec[*nbite] = cblas_dnrm2(*la, b, 1)*normb;
    }
    printf("Nombre d'itérations : %d\n", *nbite);
    free(b);
}

// MB = D où D est la diagonale de AB
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  for(int i=0; i < *la; i++){
    MB[(*kv) + 1 + i*(*lab)] = AB[(*kv) + 1 + i*(*lab)];
  }
}

// MB = D-E où D est la diagonale de AB et E la sous diagonale
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  for(int i=0; i < *la; i++){
    MB[(*kv) + i*(*lab)] = AB[(*kv) + i*(*lab)];
    MB[(*kv) + 1 + i*(*lab)] = AB[(*kv) + 1 + i*(*lab)];
  }
}

// Généralisation de richardson_alpha : x = x + M^{-1}*(b-A*x)
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {

    int info = 1;
    int NRHS = 1;
    int* ipiv = (int*)malloc(sizeof(int)*(*la));
    int ku_mb = 0;

    double* b = (double*)malloc(sizeof(double)*(*la));
    double normb = 1.0/cblas_dnrm2(*la, RHS, 1);
    // b = RHS 
    cblas_dcopy(*la, RHS, 1, b, 1);

    // b = b - A*x
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, b, 1);

    // Calcul du résidu
    resvec[*nbite] = cblas_dnrm2(*la, b, 1)*normb;

    while(*nbite < (*maxit) && resvec[*nbite] > (*tol)){
      // b = M^{-1}*b
      dgbsv_(la, kl, &ku_mb, &NRHS, MB, lab, ipiv, b, la, &info);
      // x = b + x
      cblas_daxpy(*la, 1, b, 1, X, 1);
      cblas_dcopy(*la, RHS, 1, b, 1);

      // étape suivante
      cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, b, 1);
      *nbite = *nbite + 1;
      resvec[*nbite] = cblas_dnrm2(*la, b, 1)*normb;
    }
    printf("Nombre d'itérations : %d\n", *nbite);
    free(b);
    free(ipiv);
}