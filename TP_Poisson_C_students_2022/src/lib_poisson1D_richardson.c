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

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){

    /*
    double normb = 1.0;
    resvec[*nbite] = 1/normb;
    while(*nbite < (*maxit) && resvec[*nbite] > (*tol)){
        double* tmp = (double *) malloc(sizeof(double)*(*la));
        *tmp = *X + (dgbmv_("N", la, lab, kl, ku, alpha_rich, AB, lab, X, 1, -1, RHS, 1));
        *nbite++;
        resvec[*nbite] = relative_forward_error(tmp, X, la);
        *X = *tmp;
    }
    */
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}