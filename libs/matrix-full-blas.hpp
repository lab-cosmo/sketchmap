/* Hooks for BLAS DGEMM with FMatrix
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/

#ifndef __MF_BLAS_H
#define __MF_BLAS_H 1

#include "tbdefs.hpp"
#include "matrix-full.hpp"
#include <complex>
namespace tbblas {
extern "C" {
    typedef std::complex<double> complex;
#define dgemm dgemm_
#define zgemm zgemm_
    int dgemm(char *transa, char* transb, int* m, 
              int *n, int *k, double* alpha, 
              double * a, int *lda, double* b, 
              int *ldb, double* beta, double*c,
              int* ldc);
    int zgemm(char *transa, char* transb, int* m, 
              int *n, int *k, complex* alpha, 
              complex * a, int *lda, complex* b, 
              int *ldb, complex* beta, complex*c,
              int* ldc);
}
}
/*
#ifdef _INTEL
    #include "mkl_cblas.h"
#else
    #include "cblas.h"
#endif
*/
namespace toolbox {
    //specializations of full matrix functions to use BLAS
    template <> void mult<double,double>(const FMatrix<double>& b, const FMatrix<double>& c, FMatrix<double>& a);
    template <> void mult<std::complex<double>,std::complex<double> >
            (const FMatrix<std::complex<double> >& b, const FMatrix<std::complex<double> >& c, FMatrix<std::complex<double> >& a);
    
}; //ends namespace toolbox

#endif // ends __MF_BLAS_H
