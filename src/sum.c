#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>



SEXP out(SEXP x, SEXP y) {
    R_len_t i, j, nx = length(x), ny = length(y);
    double tmp, *rx = REAL(x), *ry = REAL(y), *rans;
    SEXP ans;

    PROTECT(ans = allocMatrix(REALSXP, nx, ny));
    rans = REAL(ans);
    for(i = 0; i < nx; i++) {
        tmp = rx[i];
        for(j = 0; j < ny; j++)
            rans[i + nx*j] = tmp * ry[j];
    }
    UNPROTECT(1);
    return(ans);
}


SEXP absallpairs(SEXP x, SEXP y) {
    int i, j, nx = length(x), ny = length(y);
    double tmp, *rx = REAL(x), *ry = REAL(y), *rans;
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, nx));
    rans = REAL(ans);
    
    for(i = 0; i < nx; i++) {
        tmp = rx[i];
        rans[i] = 0;
        for(j = 0; j < ny; j++) {
            rans[i] = rans[i] +  fabs( tmp - ry[j] );
        }
        
    }
    UNPROTECT(1);
    return(ans);
}


SEXP sumabsallpairs(SEXP x, SEXP y) {
    int i, j, nx = length(x), ny = length(y);
    double tmp, *rx = REAL(x), *ry = REAL(y), *rans;
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    rans = REAL(ans);
    rans[0] = 0;
    
    for(i = 0; i < nx; i++) {
        tmp = rx[i];
        for(j = 0; j < ny; j++) {
            rans[0] +=  fabs( tmp - ry[j] );
        }
        
    }
    UNPROTECT(1);
    return(ans);
}




SEXP absallsortedpairs(SEXP x, SEXP y) {
    // x and y are sorted vectors.

    int i, k, nx = length(x);
    double tval,  *rx = REAL(x), *ry = REAL(y), *rans;
    int curR, preR;
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, nx));
  
    double *upper = (double*)malloc( sizeof(double)*(nx) ); // vector
    double *lower = (double*)malloc( sizeof(double)*(nx) ); // vector
    
    rans = REAL(ans);
    
    curR =1;
    lower[0] =0;
    upper[0] =0;
    for (i = 0; i < nx; i++) {
        if (rx[i] < ry[1-1] ) {
          curR += 1;
          lower[1-1] += ry[1-1] - rx[i];
        } else {
          upper[1-1] += rx[i] - ry[1-1];
        }
    }
    rans[0] = lower[0] + upper[0];
    
    for (k = 1; k < nx; k++) {
      tval = 0;
      preR = curR;
      
      if (curR <= nx & ry[k] >= rx[curR-1]  ) {
        while (curR <= nx & ry[k] >= rx[curR-1]  ) {
           tval = tval + (ry[k] - rx[curR-1] );
           curR = curR + 1;
        }
      }
      
      upper[k] = upper[k-1] - (nx-preR+1)*(ry[k]-ry[k-1]) + tval;
      lower[k] = lower[k-1] +    (preR-1)*(ry[k]-ry[k-1]) + tval;

        
    rans[k] = lower[k] + upper[k];
    }
    
    free(upper);
    free(lower);
    
    UNPROTECT(1);
    return(ans);
}



SEXP sumabsallsortedpairs(SEXP x, SEXP y) {
    // x and y are sorted vectors.

    int i, k, nx = length(x);
    double tval,  *rx = REAL(x), *ry = REAL(y), *rans;
    int curR, preR;
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, nx));
  
    rans = REAL(ans);
    
    curR =1;
    rans[0] =0;
    for (i = 0; i < nx; i++) {
        if (rx[i] < ry[0] ) {
          curR += 1;
          rans[0] += ry[0] - rx[i];
        } else {
          rans[0] += rx[i] - ry[0];
        }
    }
    
    for (k = 1; k < nx; k++) {
      tval = 0;
      preR = curR;
      
      if (curR <= nx & ry[k] >= rx[curR-1]  ) {
        while (curR <= nx & ry[k] >= rx[curR-1]  ) {
           tval = tval + (ry[k] - rx[curR-1] );
           curR = curR + 1;
        }
      }
       
      rans[k] = rans[k-1] + 2*tval + (2*preR-2-nx)*(ry[k]-ry[k-1]);
    }
    
    UNPROTECT(1);
    return(ans);
}









SEXP allpairsdiffmax(SEXP x) {
    // x is G x n matrix
    
    SEXP Rdim;
    Rdim = getAttrib(x, R_DimSymbol);
    //I = INTEGER(Rdim)[0];
    //J = INTEGER(Rdim)[1];

    int curR, i, j, k, n= INTEGER(Rdim)[1], G = INTEGER(Rdim)[0];
    double *rx = REAL(x), *rans, curVal, maxVal;
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, n*(n-1)/2 ));
    rans = REAL(ans);
    curR = 0;

    for (i = 0; i < n; i++) {
        for (j = i+1; j < n; j++) {
        // i&j
            maxVal = fabs(rx[i*G + 0] - rx[j*G + 0 ]);
            for (k = 1; k < G; k++) {
                curVal = fabs(rx[i*G + k] - rx[j*G + k]);
                if (curVal > maxVal) maxVal = curVal;
            }
            rans[curR]=maxVal;
            curR += 1;
        }

    }
    
  
    UNPROTECT(1);
    return(ans);
}






SEXP elementsix(SEXP x) {
    SEXP Rdim;
    Rdim = getAttrib(x, R_DimSymbol);
    //I = INTEGER(Rdim)[0];
    //J = INTEGER(Rdim)[1];
    double *rx = REAL(x), *rans;
    SEXP ans;
    
   
    PROTECT(ans = allocVector(REALSXP, 1));
    rans = REAL(ans);
    
    rans[0] = rx[5];
    
    UNPROTECT(1);
    return ans;
}


SEXP dimC(SEXP x) {
    SEXP Rdim;
    Rdim = getAttrib(x, R_DimSymbol);
    //I = INTEGER(Rdim)[0];
    //J = INTEGER(Rdim)[1];
    
    
    return Rdim;
}







SEXP createMat(SEXP x) {
    SEXP Rdim;
    Rdim = getAttrib(x, R_DimSymbol);
    //I = INTEGER(Rdim)[0];
    //J = INTEGER(Rdim)[1];
 
    SEXP ans = PROTECT(allocMatrix(REALSXP, INTEGER(Rdim)[0], INTEGER(Rdim)[1]));
    UNPROTECT(1);
    return ans;
}

SEXP createMat2(SEXP x) {
    SEXP Rdim;
    Rdim = getAttrib(x, R_DimSymbol);
    //I = INTEGER(Rdim)[0];
    //J = INTEGER(Rdim)[1];
    
    SEXP ans = PROTECT(allocMatrix(REALSXP, INTEGER(Rdim)[0], INTEGER(Rdim)[1]));
    double *pans = REAL(ans);
    for (int i = 0; i < length(ans); i++ ) pans[i] = 0;
    
    UNPROTECT(1);
    return ans;
}




