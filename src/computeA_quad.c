#include <float.h> //DBL_EPSILON
#include <R_ext/Lapack.h>
#include <R.h>
#include <Rinternals.h>
//#include <R_ext/Applic.h>
//#include <R_ext/BLAS.h>
//#include <R_ext/RS.h> //R definitions for 'extending' R, registering functions,...

#define PRINTF Rprintf
#define	MAX(A,B)	((A) > (B) ? (A) : (B))
#define	MIN(A,B)	((A) < (B) ? (A) : (B))



SEXP compute_A_quad(SEXP _a, SEXP _b, SEXP _c, SEXP _sigma_sq, SEXP _support, SEXP _xvar, SEXP _yvar){
     int n=length(_xvar);
     int K=length(_support);
     SEXP _ans=PROTECT(allocMatrix(REALSXP, n, K));
     double a=*REAL(_a), b=*REAL(_b), c=*REAL(_c);
     double sigma_sq=*REAL(_sigma_sq);
     double *support=REAL(_support), *xvar=REAL(_xvar), *yvar=REAL(_yvar), *ans=REAL(_ans);

//    PRINTF("%f %f %f %f \n",c,d,b,f);
//    PRINTF("%f %f \n",dil_r,sigma_sq);
//    PRINTF("%i %i \n",n,K);
     
	int i,k;
	double tmp;
    for(i = 0;i < n;i++){
        for(k = 0; k < K; k++){
            tmp = pow(xvar[i] - support[k], 2) + pow(yvar[i] - (a*support[k]*support[k]+b*support[k]+c), 2);
			ans[i+k*n] = 1/sigma_sq * exp(-tmp/sigma_sq/2);
		}
	}
	
    UNPROTECT(1);
    return _ans;
}

