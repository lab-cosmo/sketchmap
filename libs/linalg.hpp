#ifndef __LINALG_HPP
#define __LINALG_HPP
#include <complex>
#include "matrix-full.hpp"

namespace tblapack{
    typedef std::complex<double> complex;
    extern "C" {
#define dgesv dgesv_
#define dgesvd dgesvd_
#define dsyev dsyev_
#define dsyevx dsyevx_
#define zgeev zgeev_
#define dgeev dgeev_
#define dgetrf dgetrf_
#define dgetri dgetri_

    int  dgesv (int* n, 
        int* nrhs, double* a, 
        int* lda,  int *ipiv,  
        double* b,  int* ldb,  int* info);
    int dsyev(char *jobz, char *uplo, int *n, double *a,
                int *lda, double *w, double *work, int *lwork, 
                int *info);
    int dsyevx(char *jobz, char* range, char *uplo, int *n, 
               double *a, int *lda,  double* vl, double* vu, int* il, int* iu,
               double* abtol, int *nfound,  double *w, double *z, int *ldz, 
               double *work, int *lwork, int*iwork, int *ifail, int *info);
    int dgeev(char *jobvl, char *jobvr, int *n, double *a,
              int *lda, double *wr, double *wi, double* vl, int* ldvl,
              double* vr, int*ldvr, double *work, int *lwork, int *info);
    int zgeev(char *jobvl, char *jobvr, int *n, complex *a,
              int *lda, complex *w, complex* vl, int* ldvl,
              complex* vr, int*ldvr, complex *work, int *lwork, 
              double *rwork, int *info);
    int dgesvd( char * jobu, char * jobvt, int *nrows, int *ncols, double *AT, 
                int *n3, double *S, double *UT, int *n1, double *V, int *n2, 
                double *work, int *lwork, int *info );
    int dgetrf(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
            
    int dgetri(int *n, double *a, int *lda, int *ipiv, double *work, int* lwork, int *info);
    }
}


namespace toolbox {
    void LinearSystem(const FMatrix<double>& A, const std::valarray<double>& b, std::valarray<double>& x); 
    enum LAEigenMode {LAEMAll, LAEMValue, LAEMIndex };
    void EigenSolverSym(const FMatrix<double>& A, FMatrix<double>& Q, std::valarray<double>& q, LAEigenMode emode=LAEMAll, double vmin=0.0, double vmax=1.0);
    void EigenSolver(const FMatrix<tblapack::complex>& A, FMatrix<tblapack::complex>& RQ, FMatrix<tblapack::complex>& LQ, std::valarray<tblapack::complex>& q);
    void EigenDecomposition(const FMatrix<tblapack::complex>& A, FMatrix<tblapack::complex>& O, FMatrix<tblapack::complex>& O1, std::valarray<tblapack::complex>& q);
    void EigenSolver(const FMatrix<double>& A, FMatrix<tblapack::complex>& RQ, FMatrix<tblapack::complex>& LQ, std::valarray<tblapack::complex>& q);
    void EigenDecomposition(const FMatrix<double>& A, FMatrix<tblapack::complex>& O, FMatrix<tblapack::complex>& O1, std::valarray<tblapack::complex>& q);
    void SVDecomposition(const FMatrix<double>& A, FMatrix<double>& U, FMatrix<double>& VT, std::valarray<double>& s);
    
    void MatrixInverse(const FMatrix<double>& A, FMatrix<double>& IA);
    void PseudoInverse(const FMatrix<double>& A, FMatrix<double>& IA);
    void MatrixFunctionSym(const FMatrix<double>& A,  double (*f) (double), FMatrix<double>& FA);
    void MatrixFunction(const FMatrix<double>& A,  tblapack::complex (*f) (tblapack::complex), FMatrix<double>& FA);
    void Cholesky(const FMatrix<double>& MMT, FMatrix<double>& M);
    void StabCholesky(const FMatrix<double>& MMT, FMatrix<double>& M);
} //ends namespace toolbox
#endif //ends #ifndef __LINALG_HPP
