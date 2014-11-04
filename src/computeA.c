#include <float.h> //DBL_EPSILON
#include <R_ext/Lapack.h>
#include <R.h>
#include <Rinternals.h>
//#include <R_ext/Applic.h>
//#include <R_ext/BLAS.h>
//#include <R_ext/RS.h> //R definitions for 'extending' R, registering functions,...

#define PRINTF Rprintf
#define MAX(A,B)    ((A) > (B) ? (A) : (B))
#define MIN(A,B)    ((A) < (B) ? (A) : (B))



SEXP compute_A(SEXP _c, SEXP _d, SEXP _b, SEXP _f, SEXP _dil_r, SEXP _sigma_sq, SEXP _support, SEXP _xvar, SEXP _yvar){
     int n=length(_xvar);
     int K=length(_support);
     SEXP _ans=PROTECT(allocMatrix(REALSXP, n, K));
     double c=*REAL(_c), d=*REAL(_d), b=*REAL(_b), f=*REAL(_f);
     double dil_r=*REAL(_dil_r), sigma_sq=*REAL(_sigma_sq);
     double *support=REAL(_support), *xvar=REAL(_xvar), *yvar=REAL(_yvar), *ans=REAL(_ans);

//    PRINTF("%f %f %f %f \n",c,d,b,f);
//    PRINTF("%f %f \n",dil_r,sigma_sq);
//    PRINTF("%i %i \n",n,K);
     
    int i,k;
    double tmp;
    for(i = 0;i < n;i++){
        for(k = 0; k < K; k++){
            tmp = pow(xvar[i] - support[k], 2) + pow(yvar[i] - log(c+(d-c)/pow(1+(pow((d-c)/(exp(support[k])-c), 1/f)-1)*pow(dil_r,b), f)), 2);
            ans[i+k*n] = 1/sigma_sq * exp(-tmp/sigma_sq/2);
        }
    }
    
    UNPROTECT(1);
    return _ans;
}

SEXP compute_four_pl_prc (SEXP _c,SEXP _d,SEXP _b,SEXP _f, SEXP _xx, SEXP _dil_r) {
     int n=length(_xx);
     SEXP _ans=PROTECT(allocVector(REALSXP, n));
     double c=*REAL(_c), d=*REAL(_d), b=*REAL(_b), f=*REAL(_f);
     double dil_r=*REAL(_dil_r);
     double *xx=REAL(_xx), *ans=REAL(_ans);
     
     int i;
    for(i = 0;i < n;i++){
         ans[i]=log(c+(d-c)/pow(1+(pow((d-c)/(exp(xx[i])-c), 1/f)-1)*pow(dil_r,b), f));
    }
    
    UNPROTECT(1);
    return _ans;
}
