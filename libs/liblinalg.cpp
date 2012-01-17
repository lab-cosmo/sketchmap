#include "linalg.hpp"
using namespace tblapack;

namespace toolbox {
void LinearSystem(const FMatrix<double>& A, const std::valarray<double>& b, std::valarray<double>& x)
{
    int n=b.size(), nrhs=1, info=0;
#ifdef DEBUG
    if (A.rows()!=n || A.cols()!=n) ERROR("Matrix and vector size mismatch.");
#endif
    x.resize(n); x=b;
    std::valarray<int> ipiv(n);
    FMatrix<double> AT;
    transpose(A, AT);  //transpose because of c/fortan row convention
    //double *dbl=&(const_cast<FMatrix<double> &>(A)(0,0));
    int ilap=tblapack::dgesv(&n,&nrhs,&AT(0,0),&n,&ipiv[0],&x[0],&n,&info);
    if (info!=0) {std::cerr<<"LAPACK ERROR "<<info<<"\n"; ERROR("Error in dgesv call");}
}
        
//solver for symmetric matrix!!!
void EigenSolverSym(const FMatrix<double>& A, FMatrix<double>& Q, std::valarray<double>& q, LAEigenMode emode, double vmin, double vmax)
{
    unsigned long n=A.rows();
    if (A.cols()!=n) ERROR("Argument should be a square matrix\n");
    
    if (emode==LAEMAll) 
    {
        char JOBZ='V', UPLO='U';
        int N=n, LDA=n, LWORK=100*n, INFO;
        std::valarray<double> WORK(LWORK);
        toolbox::FMatrix<double> lA(A);
        q.resize(n);
        tblapack::dsyev(&JOBZ,&UPLO,&N,&(lA(0,0)),&LDA,&(q[0]),&WORK[0],&LWORK,&INFO);
        if (INFO!=0) ERROR("DSYEV returned an error");
        transpose(lA,Q);
    }
    else
    {
        char JOBZ='V', RANGE=(emode==LAEMValue?'V':'I'), UPLO='U';
        int N=n, LDA=n, INFO, LDZ=n, LWORK=100*n, NFOUND;
        double VL=vmin, VU=vmax, ABTOL=0.0; int IL=((int)vmin)+1, IU=((int)vmax)+1;
        std::valarray<int> IFAIL(n), IWORK(5*n);
        std::valarray<double> WORK(LWORK), lq(n);
        toolbox::FMatrix<double> lA(A), Z(n,n);
        
        tblapack::dsyevx(&JOBZ,&RANGE,&UPLO,&N,
                          &(lA(0,0)),&LDA,&VL,&VU,&IL,&IU,
                            &ABTOL, &NFOUND, &(lq[0]), &Z(0,0),
                           &LDZ, &WORK[0],&LWORK,&IWORK[0],
                           &IFAIL[0],&INFO);
        if (INFO!=0) ERROR("DSYEV returned an error: failing vector was "<<IFAIL);
        Q.resize(n,NFOUND); q.resize(NFOUND);
        for (unsigned long i=0; i<NFOUND; ++i)
        { q[i]=lq[i];    for (unsigned long j=0; j<n; ++j) Q(j,i)=Z(i,j); }
        
    }
}


void StabCholesky(const toolbox::FMatrix<double>& MMt, toolbox::FMatrix<double>& M)
{
    //"stabilized" cholesky which can handle blocks of zeroes and zeroes out negative "eigenvalues"
    unsigned long n=MMt.rows();
    if (MMt.cols()!=n) ERROR("Argument should be a square matrix\n");
    
    FMatrix<double> L(n,n), D(n,n);
    D*=0.; L*=0.;
    unsigned long i,j,k;
    for(i=0; i<n; ++i)
    {
        L(i,i)=1.;
        for (j=0; j<i; j++)
        {
            L(i,j)=MMt(i,j);
            for (k=0; k<j; ++k) L(i,j)-=L(i,k)*L(j,k)*D(k,k);
            if (D(j,j)!=0.) L(i,j)/=D(j,j); else L(i,j)=0.;
        }
        D(i,i)=MMt(i,i);
        for (k=0; k<i; ++k) D(i,i)-=L(i,k)*L(i,k)*D(k,k);
    }
    //FMatrix<double> LT(n,n), R(n,n);
    //transpose(L,LT);
    //mult(L,D,R); mult(R,LT,M); M-=MMt;
    
    for(i=0; i<n; ++i) D(i,i)=(D(i,i)>0.?sqrt(D(i,i)):0.);
    mult(L,D,M);
}

void Cholesky(const toolbox::FMatrix<double>& MMt, toolbox::FMatrix<double>& M)
{
    unsigned long n=MMt.rows();
    if (MMt.cols()!=n) ERROR("Argument should be a square matrix\n");
    
    M.resize(n,n); M*=0.;
    M(0,0)=sqrt(MMt(0,0));
    unsigned long i,j,k;
    for(i=0; i<n; ++i)
    {
        for (j=0; j<i; j++)
        {
            M(i,j)=MMt(i,j);
            for (k=0; k<j; ++k) M(i,j)-=M(i,k)*M(j,k);
            M(i,j)*=1./M(j,j);
        }
        M(i,i)=MMt(i,i);
        for (k=0; k<i; ++k) M(i,i)-=M(i,k)*M(i,k);
        M(i,i)=sqrt(M(i,i));
    }
}

void MatrixFunctionSym(const FMatrix<double>& A,  double (*f) (double), FMatrix<double>& FA)
{
    unsigned long n=A.rows();
    if (A.cols()!=n) ERROR("Argument should be a square matrix\n");
    FMatrix<double> W, WT; std::valarray<double> v;
    EigenSolverSym(A,W,v);
    transpose(W,WT);
    for(int i=0; i<n; ++i) 
    {
        v[i]=f(v[i]);
        for(int j=0; j<n; ++j) WT(i,j)*=v[i];
    }
    mult(W,WT,FA);
}

void MatrixFunction(const FMatrix<double>& A,  tblapack::complex (*f) (tblapack::complex), FMatrix<double>& FA)
{
    unsigned long n=A.rows();
    if (A.cols()!=n) ERROR("Argument should be a square matrix\n");
    FMatrix<tblapack::complex> W, W1, CA; std::valarray<complex> v;
    EigenDecomposition(A,W,W1,v);
    
    for(int i=0; i<n; ++i) 
    {
        v[i]=f(v[i]);
        for(int j=0; j<n; ++j) W1(i,j)*=v[i];
    }
    mult(W,W1,CA);
    FA.resize(n,n);
    for(int i=0; i<n; ++i) 
        for(int j=0; j<n; ++j) 
            FA(i,j)=CA(i,j).real();
}

void MatrixInverse(const FMatrix<double>& A, FMatrix<double>& IA)
{
    int n=A.rows();
    if (A.cols()!=n) ERROR("Argument should be a square matrix\n");
    FMatrix<double> lA; transpose(A,lA);  //accounts for storage of FORTRAN in col-major
    std::valarray<int> ipiv(n); 
    int info;
    tblapack::dgetrf(&n,&n,&lA(0,0),&n,&ipiv[0],&info);
    if (info!=0) ERROR("Error "<<info<<" in call to dgetrf.");
    
    int lwork=-1; std::valarray<double> work(n);
    tblapack::dgetri(&n,&lA(0,0),&n,&ipiv[0],&work[0],&lwork,&info);
    if (info!=0) ERROR("Error "<<info<<" in call to dgetri.");
    
    lwork=(int)work[0]; work.resize(lwork);
    tblapack::dgetri(&n,&lA(0,0),&n,&ipiv[0],&work[0],&lwork,&info);
    if (info!=0) ERROR("Error "<<info<<" in call to dgetri.");
    
    transpose(lA,IA); //once again FORT -> C conversion
}

void EigenDecomposition(const FMatrix<tblapack::complex>& A, FMatrix<tblapack::complex>& O, FMatrix<tblapack::complex>& O1, std::valarray<tblapack::complex>& q)
{
    
    int n=A.rows();
    if (A.cols()!=n) ERROR("Argument should be a square matrix\n");
    //finds O and D=diag(q) such that A=O D O1
    EigenSolver(A, O, O1, q);
    //almost there: we only need to renormalize O (which contains now right eigv) and O1 (which contains left eigv)
    std::valarray<complex> w(n); w=0.;
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<n; ++j) w[i]+=O1(i,j)*O(j,i);
        w[i]=complex(1.,0.)/sqrt(w[i]);
    }
    for (int i=0; i<n; ++i) for (int j=0; j<n; ++j) { O(j,i)*=w[i]; O1(i,j)*=w[i]; }
}

void SVDecomposition(const FMatrix<double>& A, FMatrix<double>& U, FMatrix<double>& VT, std::valarray<double>& s)
{
    FMatrix<double> AT; transpose(A,AT);
    
  // Create column and row information on the matrix
    int nsv, nrows, ncols, info; nrows=A.rows(); ncols=A.cols();
    if(nrows>ncols){nsv=ncols;}else{nsv=nrows;}

  // Create some containers for stuff from single value decomposition
    s.resize(nsv); 
    FMatrix<double> UT(nrows,nrows);
    FMatrix<double> V(ncols,ncols);

    char jobu='A', jobvt='A';     // These decide the mode of operation of lapack

  // This optimizes the size of the work array used in lapack singular value decomposition
    int lwork=-1; std::valarray<double> work(1);
    dgesvd( &jobu, &jobvt, &nrows, &ncols, &AT(0,0), &nrows, &s[0], &UT(0,0), &nrows, &V(0,0), &ncols, &work[0], &lwork, &info );
    if(info!=0) ERROR("Return "<<info<<" in work optimization call to dgesvd");

  // Retrieve correct sizes for work and rellocate
    lwork=(int) work[0]; work.resize(lwork);

  // This does the singular value decomposition
    dgesvd( &jobu, &jobvt, &nrows, &ncols, &AT(0,0), &nrows, &s[0], &UT(0,0), &nrows, &V(0,0), &ncols, &work[0], &lwork, &info );
    if(info!=0) ERROR("Return "<<info<<" in call to dgesvd");

    transpose(UT,U); transpose(V,VT);
}
        
void PseudoInverse(const FMatrix<double>& A, FMatrix<double>& IA)
{
    FMatrix<double> U, VT; std::valarray<double> s;
    SVDecomposition(A,U,VT,s);
    
    unsigned long nsv=s.size();
    // Compute the tolerance on the singular values ( machine epsilon * nsv * maximum singular value )
    double tol=s.max();
    tol*=nsv*1.11E-16; //machine precision
    
    // Get the inverses of the singlular values
    FMatrix<double> IS( A.rows() , A.cols(), 0.0); 
    
    for(unsigned long i=0;i<nsv;++i){ if( s[i]>tol ){ IS(i,i)=1./s[i]; }else{ IS(i,i)=0.0; } }

  // And now compute the psedoinverse
    FMatrix<double> tmp; //, V, UT;
    //transpose(VT,V); transpose(U,UT);    
    mult(U,IS,tmp); mult(tmp,VT,IS);  transpose(IS,IA);
    
    //CHECK
    //mult(A,IA,tmp); mult(tmp,A,IS); IS-=A; std::cerr<<"PSEUDO ERRR " <<normfrob(IS)<<"\n";
}

//solver for generic matrix!!!
/**********************************************************
 let D be the diagonal matrix of the elements of q. then
 A RQ = RQ D
 LQ A = D LQ
************************************************************/
void EigenSolver(const FMatrix<complex>& A, FMatrix<complex>& RQ, FMatrix<complex>& LQ, std::valarray<complex>& q)
{
    
    int n=A.rows();
    if (A.cols()!=n) ERROR("Argument should be a square matrix\n");
    
    char jobvl='V', jobvr='V';
    int lda=n, ldvr=n, ldvl=n;
    std::valarray<complex> work(n);
    FMatrix<complex> la;
    transpose(A,la);
    std::valarray<double> rwork(2*n);
    FMatrix<complex> vr(n,n), vl(n,n);
    q.resize(n);
    
    int info, lwork=-1;
    //here we only get the size of the work array
    tblapack::zgeev(&jobvl, &jobvr, &n, &(la(0,0)),
                     &lda, &q[0], &(vl(0,0)), &ldvl,
                                    &(vr(0,0)), &ldvr, &work[0], &lwork, 
                                      &rwork[0], &info);
    if (info!=0) 
        ERROR("Error in cgeev: error code: "<<info);
    lwork=(int) real(work[0]);  work.resize(lwork);
    
    tblapack::zgeev(&jobvl, &jobvr, &n, &la(0,0),
                     &lda, &q[0], &vl(0,0), &ldvl,
                                      &vr(0,0), &ldvr, &work[0], &lwork, 
                                          &rwork[0], &info);
    if (info!=0) 
        ERROR("Error in cgeev: error code: "<<info);
    transpose(vr, RQ);
    LQ=vl; toolbox::map(LQ,conj);
}

//solver for generic matrix!!!
/**********************************************************
 let D be the diagonal matrix of the elements of q. then
 A RQ = RQ D
 LQ A = D LQ
************************************************************/
void EigenSolver(const FMatrix<double>& A, FMatrix<complex>& RQ, FMatrix<complex>& LQ, std::valarray<complex>& q)
{
    //std::cerr<<"real eigensolver\n";
    int n=A.rows();
    if (A.cols()!=n) ERROR("Argument should be a square matrix\n");
    
    char jobvl='V', jobvr='V';
    int lda=n, ldvr=n, ldvl=n;
    std::valarray<double> work(n), wr(n), wi(n);
    FMatrix<double> la;
    transpose(A,la);
    FMatrix<double> vr(n,n), vl(n,n);
    
    int info, lwork=-1;
    //here we only get the size of the work array
    tblapack::dgeev(&jobvl, &jobvr, &n, &(la(0,0)),
                     &lda, &wr[0], &wi[0], &(vl(0,0)), &ldvl,
                    &(vr(0,0)), &ldvr, &work[0], &lwork, &info);
    if (info!=0) 
        ERROR("Error in cgeev: error code: "<<info);
    lwork=(int) work[0];  work.resize(lwork);
    tblapack::dgeev(&jobvl, &jobvr, &n, &(la(0,0)),
                    &lda, &wr[0], &wi[0], &(vl(0,0)), &ldvl,
                    &(vr(0,0)), &ldvr, &work[0], &lwork, &info);
    if (info!=0) 
        ERROR("Error in cgeev: error code: "<<info);
    //now, lines of vr contains eigenvectors, with splitted real/imag part. we must put them together
    q.resize(n); RQ.resize(n,n); RQ*=0.; LQ=RQ;
    for (unsigned long i=0; i<n; ++i)
    {
        q[i]=complex(wr[i],wi[i]);
        if(wi[i]==0.)
        {
            //real eigenvalue!
            for(unsigned long j=0; j<n; ++j) 
            { RQ(j,i)=complex(vr(i,j),0.); LQ(i,j)=complex(vl(i,j),0.); }
        }
        else
        {
            q[i+1]=complex(wr[i],-wi[i]);
            for(unsigned long j=0; j<n; ++j) 
            { RQ(j,i)=complex(vr(i,j),vr(i+1,j)); RQ(j,i+1)=complex(vr(i,j),-vr(i+1,j));  
              LQ(i,j)=complex(vl(i,j),-vl(i+1,j)); LQ(i+1,j)=complex(vl(i,j),vl(i+1,j)); }
            ++i;
        }
    }
}

void EigenDecomposition(const FMatrix<double>& A, FMatrix<tblapack::complex>& O, FMatrix<tblapack::complex>& O1, std::valarray<tblapack::complex>& q)
{
    
    int n=A.rows();
    if (A.cols()!=n) ERROR("Argument should be a square matrix\n");
    //finds O and D=diag(q) such that A=O D O1
    EigenSolver(A, O, O1, q);
    //almost there: we only need to renormalize O (which contains now right eigv) and O1 (which contains left eigv)
    std::valarray<complex> w(n); w=0.;
    
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<n; ++j) w[i]+=O1(i,j)*O(j,i);
        w[i]=complex(1.,0.)/w[i];
    }
    for (int i=0; i<n; ++i) for (int j=0; j<n; ++j) { O1(i,j)*=w[i]; }
}

}; //ends namespace toolbox

