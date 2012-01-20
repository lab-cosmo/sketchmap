/* Hooks for BLAS DGEMM for the full-matrix template.
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/
   
#include "matrix-full-blas.hpp"
using namespace tbblas;
namespace toolbox{
template <> void mult<double,double>(const FMatrix<double>& b, const FMatrix<double>& c, FMatrix<double>& a)
{
    //B and C are exchanged, because of BLAS using fortran row/col convention
#ifdef BENCHMARK
    TBBenchmarks["full_mult_dgemm"].n_calls++; TBBenchmarks["full_mult_dgemm"].timer.start(); 
#endif
#ifdef DEBUG
    if (c.wr!=b.wc) ERROR("Incompatible matrix size")
#endif
    a.resize(b.wr,c.wc);
    
    //because of BLAS using fortran row/col convention we need to make A^T=C^T B^T
    int m=c.wc, n=b.wr, k=c.wr;
    int lda=m, ldb=k, ldc=m;
    char trans='N';
    double alpha=1., beta=0.;
    dgemm(&trans, &trans, &m, &n, &k,
          &alpha, &(const_cast<FMatrix<double>&>(c).data[0]), &lda,
          &(const_cast<FMatrix<double>&>(b).data[0]), &ldb, &beta, 
          &(a.data[0]), &ldc);

    
#ifdef BENCHMARK
    TBBenchmarks["full_mult_dgemm"].timer.stop(); TBBenchmarks["full_mult_dgemm"].tot_time+=TBBenchmarks["full_mult_dgemm"].timer;
#endif
}

template <> void 
mult<std::complex<double>,std::complex<double> >
        (const FMatrix<std::complex<double> >& b, const FMatrix<std::complex<double> >& c, FMatrix<std::complex<double> >& a)
{
#ifdef BENCHMARK
    TBBenchmarks["full_mult_zgemm"].n_calls++; TBBenchmarks["full_mult_zgemm"].timer.start(); 
#endif
#ifdef DEBUG
    if (c.wr!=b.wc) ERROR("Incompatible matrix size")
#endif
    a.resize(b.wr,c.wc);
    
    //because of BLAS using fortran row/col convention we need to make A^T=C^T B^T
    int m=c.wc, n=b.wr, k=c.wr;
    int lda=m, ldb=k, ldc=m;
    char trans='N';
    complex alpha=complex(1.,0.), beta=complex(0.,0.);
    zgemm(&trans, &trans, &m, &n, &k,
         &alpha, &(const_cast<FMatrix<complex>&>(c).data[0]), &lda,
         &(const_cast<FMatrix<complex>&>(b).data[0]), &ldb, &beta, 
         &(a.data[0]), &ldc);

    
    
    /*
    int lda=a.wr, ldb=b.wr, ldc=c.wr;
    std::complex<double> one(1.,0.), zero(0.,0.);

    cblas_zgemm(CblasRowMajor, 
                CblasNoTrans, CblasNoTrans, b.wr, c.wc, b.wc,
                &one, &(const_cast<FMatrix<std::complex<double> >&>(b).data[0]), ldb, 
                        &(const_cast<FMatrix<std::complex<double> >&>(c).data[0]), ldc, &zero, & (a.data[0]), lda);

    std::cerr<<"\n";*/
#ifdef BENCHMARK
    TBBenchmarks["full_mult_zgemm"].timer.stop(); TBBenchmarks["full_mult_zgemm"].tot_time+=TBBenchmarks["full_mult_zgemm"].timer;
#endif
}
}; //ends namespace toolbox
