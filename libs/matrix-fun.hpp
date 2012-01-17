#ifndef __MATRIX_FUN_H
#define __MATRIX_FUN_H

#include "matrices.hpp"
#include <valarray>

namespace toolbox{
/**********************************************************************
 a collection of functions to compute matrix polynomials with different
 methods. these are conceived to work on any matrix type defining the 
 same functions as the CrsMatrix, so are made up as templates.
 for what concerns extra properties such as autotruncation, these are
 inherited from R. 
***********************************************************************/
template <class MC, class U>
void poly_standard(const MC& M, const std::valarray<std::valarray<U>* >& b, std::valarray<MC*>& R )
{
#ifdef DEBUG
    if (M.cols()!=M.rows()) ERROR("Input matrix must be square!");
#endif

    unsigned long ncmp=b.size();
#ifdef DEBUG
    if (R.size()!=ncmp) ERROR("Return array must be initialized to as many valid pointers to matrices as the size of coefficients list");
#endif

    typename MC::index_type n=M.rows(), q;

    for (unsigned long icmp=0; icmp<ncmp; ++icmp)
    {
        q=b[icmp]->size()-1;

        if (R[icmp]->rows()!=M.rows() || R[icmp]->cols()!=M.cols()) R[icmp]->resize(n,n);
        MC tmp; MatrixOptions<MC> rops; R[icmp]->getops(rops); tmp.setops(rops);
        (*R[icmp])*=0.;
        if (q==0) 
        {
            (*R[icmp])+=(*b[icmp])[0];
            return;
        }
        else
        {
            incr(*R[icmp],M,(*b[icmp])[q]);
            (*R[icmp])+=(*b[icmp])[q-1];
        }
    
    
        for (long k=q-2; k>=0; --k)
        {
            //avoids unnecessary copying
            if ((q-k)%2==0)
            {
                mult(*R[icmp],M,tmp);
                tmp+=(*b[icmp])[k];
            }
            else
            {
                mult(tmp,M,*R[icmp]);
                (*R[icmp])+=(*b[icmp])[k];
            }
        }
        if (q%2==0) (*R[icmp])=tmp;
    }
}

template <class MC, class U>
void poly_standard(const MC& M, const std::valarray<U>& b, MC& R )
{
    const std::valarray<std::valarray<U>* > rb((const_cast<std::valarray<U>* >(&b)),1);
    std::valarray< MC* >rR(&R,1);
    poly_standard(M, rb, rR);
}

//computes many polynomials doing the matrix powers just once. The properties of the first matrix are used 
//to set autotruncation and similar in the results.
template <class MC, class U>
void poly_standard_nonhorner(const MC& M, const std::valarray<std::valarray<U>* >& b, std::valarray<MC*>& R )
{
#ifdef DEBUG
    if (M.cols()!=M.rows()) ERROR("Input matrix must be square!")
#endif
    unsigned long ncmp=b.size();
#ifdef DEBUG
    if (R.size()!=ncmp) ERROR("Return array must be initialized to as many valid pointers to matrices as the size of coefficients list");
#endif
    
    typename MC::index_type n=M.rows(), q;
#ifdef DEBUG
    unsigned long q0=b[0]->size()-1;
#endif
    
    MC tmp, tmp2;
    MatrixOptions<MC> rops; (*R[0]).getops(rops); tmp.setops(rops); tmp2.setops(rops);
    tmp.resize(n,n);
    for (unsigned long icmp=0; icmp<ncmp; ++icmp)
    {
        q=b[icmp]->size()-1;
#ifdef DEBUG
        if (q!=q0) ERROR("All the polynomials must be of the same degree");
#endif
        
        (*R[icmp]).resize(n,n);
        if (q==0) {  (*R[icmp])+=(*b[icmp])[0];  continue;  }
        
        incr((*R[icmp]),M,(*b[icmp])[1]);
        (*R[icmp])+=(*b[icmp])[0];
    }
 
    if (q<2) return;
        
    incr(tmp,M);

    for (unsigned int i=2; i<q; ++i)
    {
        //we avoid to copy from tmp to tmp2, and we upgrade the results at once
        if (i%2==0)
        {
            mult(tmp,M,tmp2); tmp.resize(0,0); 
            for (unsigned long icmp=0; icmp<ncmp; ++icmp)
                incr((*R[icmp]),tmp2,(*b[icmp])[i]);
        }
        else
        {
            mult(tmp2,M,tmp); tmp2.resize(0,0); 
            for (unsigned long icmp=0; icmp<ncmp; ++icmp)
                incr((*R[icmp]),tmp,(*b[icmp])[i]);
        }
    }
}

template <class MC, class U>
void poly_standard_nonhorner(const MC& M, const std::valarray<U>& b, MC& R )
{
    const std::valarray<std::valarray<U>* > rb((const_cast<std::valarray<U>* >(&b)),1);
    std::valarray< MC* >rR(&R,1);
    poly_standard_nonhorner(M, rb, rR);
}

/******************************************************************
  matrix plynomial computation using ~2 sqrt(q) matrix multiplications,
  but requiring to hold O(sqrt(q)) indermediate matrices
******************************************************************/
template <class MC, class U>
void poly_vanloan1(const MC& M, const std::valarray<std::valarray<U>* >& b, std::valarray<MC*>& R )
{
#ifdef BENCHMARK
    unsigned long nm=0;
    TBBenchmarks["mf_polyvanl1"].n_calls++; TBBenchmarks["mf_polyvanl1"].timer.start(); 
#endif
#ifdef DEBUG
    if (M.cols()!=M.rows()) ERROR("Input matrix must be square!");
#endif
    unsigned long ncmp=b.size();
#ifdef DEBUG
    if (R.size()!=ncmp) ERROR("Return array must be initialized to as many valid pointers to matrices as the size of coefficients list");
#endif
    
    typename MC::index_type n=M.rows(), q=b[0]->size()-1;
    
    unsigned long s=(unsigned long) sqrt(1. *q), r=q/s, sr=s*r;
    
    MC tmp; 
    MatrixOptions<MC> rops; (*R[0]).getops(rops); tmp.setops(rops); tmp.resize(n,n);
    std::valarray<MC>Y(tmp,s+1);
    //prepares the array with powers up to s
    Y[0]+=(typename MC::data_type) 1.;
    incr(Y[1],M);
    for (unsigned long i=2;i<=s;++i) 
    { 
        mult(Y[i-1],M,Y[i]); 
        std::perr<<"#POWER "<<i<<" el. per row: "<<Y[i].size()/n<<"  norm: "<<norminf(Y[i])<<"\n";
    }
#ifdef BENCHMARK
    nm+=s-1;
#endif
    
    for (unsigned long icmp=0; icmp<ncmp; ++icmp)
    {
#ifdef DEBUG
        if (b[icmp]->size()-1!=q) ERROR("All the polynomials must be of the same degree");
#endif
        (*R[icmp]).resize(n,n);
        if (q-sr>= 1)
        {
            incr(*R[icmp],M,(*b[icmp])[sr+1]); 
            incr(*R[icmp],Y[0],(*b[icmp])[sr]);  //ok, no need to store the identity matrix. this should be fixed sooner or later
            for (unsigned long i=2; i<=q-sr; ++i) incr(*R[icmp],Y[i],(*b[icmp])[sr+i]);
        }
        else
        {
            *R[icmp]=Y[0]; *R[icmp]*=(*b[icmp])[sr];
        }
        for (unsigned long k=1; k<=r; ++k)
        {
            
            mult(Y[s],*R[icmp],tmp);
            for (unsigned long i=0; i<s; ++i)
                incr(tmp,Y[i],(*b[icmp])[i+s*(r-k)]);
            *R[icmp]=tmp;
            std::perr<<"# POLY STEP "<<k<<" el. per row: "<<R[icmp]->size()/n<<"  norm: "<<norminf(*R[icmp])<<" coeff: "<< (*b[icmp])[k]<<"\n";
        }
#ifdef BENCHMARK
        nm+=r+1;
#endif
    }
#ifdef BENCHMARK
    TBBenchmarks["mf_polyvanl1"].timer.stop(); TBBenchmarks["mf_polyvanl1"].tot_time+=TBBenchmarks["mf_polyvanl1"].timer;
    TBBenchmarks["mf_polyvanl1"].tot_prop["mm_mult"]+=nm;  //n. of multiplications
#endif
}

template <class MC, class U>
void poly_vanloan1(const MC& M, const std::valarray<U>& b, MC& R )
{
    const std::valarray<std::valarray<U>* > rb((const_cast<std::valarray<U>* >(&b)),1);
    std::valarray< MC* >rR(&R,1);
    poly_vanloan1(M, rb, rR);
}

/******************************************************************
  matrix polynomial computation using ~2 sqrt(q) matrix multiplications,
  but requiring to hold O(sqrt(q)) indermediate matrices. should be
  SLIGHTLY faster than vanloan 1...
******************************************************************/
template <class MC, class U>
void poly_liang(const MC& M, const std::valarray< std::valarray<U> * >& b, std::valarray<MC*>& R )
{
#ifdef BENCHMARK
    unsigned long nm=0;
    TBBenchmarks["mf_polyliang"].n_calls++; TBBenchmarks["mf_polyliang"].timer.start(); 
#endif
#ifdef DEBUG
    if (M.cols()!=M.rows()) ERROR("Input matrix must be square!");
#endif
    
    unsigned long ncmp=b.size();
#ifdef DEBUG
    if (R.size()!=ncmp) ERROR("Return array must be initialized to as many valid pointers to matrices as the size of coefficients list");
#endif

    typename MC::index_type n=M.rows();
    unsigned long q=b[0]->size()-1;
    unsigned long s=(unsigned long) floor(sqrt(1. + q)), r=(unsigned long) (1+floor(1.*q/sqrt(1.+q)));
    if (q>r*s) ++s;
    
    MC tmp; 
    MatrixOptions<MC> rops; (*R[0]).getops(rops); tmp.setops(rops); tmp.resize(n,n);
    std::valarray<MC> Y(tmp,r+1); 
    incr(Y[0], 1.);
    //fills up the array of matrices up to r 
    if (r>0) incr(Y[1],M);
    for (unsigned long i=2;i<=r;++i) mult(Y[i-1],M,Y[i]); 
#ifdef BENCHMARK
    nm+=r-1;
#endif
    
    unsigned long sr=(s-1)*r;
    
    for (unsigned long icmp=0; icmp<ncmp; ++icmp)
    {
#ifdef DEBUG
        if (b[icmp]->size()-1!=q) ERROR("All the polynomials must be of the same degree");
#endif
        (*R[icmp]).resize(n,n);
    
        for (unsigned long i=1; i<=q-sr; ++i) incr(*R[icmp],Y[i],(*b[icmp])[sr+i]); 
    
        for (long k=s-2; k>=0; --k)
        {
            mult(*R[icmp],Y[r],tmp); 
#ifdef BENCHMARK
            nm++;
#endif
            for (unsigned long i=1; i<=r; ++i)
                incr(tmp,Y[i],(*b[icmp])[i+r*k]);
            (*R[icmp])=tmp;
        }
        incr(*R[icmp],Y[0],(*b[icmp])[0]);
    }
#ifdef BENCHMARK
    TBBenchmarks["mf_polyliang"].timer.stop(); TBBenchmarks["mf_polyliang"].tot_time+=TBBenchmarks["mf_polyliang"].timer;
    TBBenchmarks["mf_polyliang"].tot_prop["mm_mult"]+=nm;  //n. of multiplications
#endif
}

template <class MC, class U>
void poly_liang(const MC& M, const std::valarray<U>& b, MC& R )
{
    const std::valarray<std::valarray<U>* > rb((const_cast<std::valarray<U>* >(&b)),1);
    std::valarray< MC* >rR(&R,1);
    poly_liang(M, rb, rR);
}

#define __POLY_BREAKEVEN 10
template <class MC, class U>
void poly(const MC& M, const std::valarray<U>& b, MC& R )
{
    if (b.size()<= __POLY_BREAKEVEN) poly_standard(M,b,R);
    else poly_vanloan1(M,b,R);
}

/******************************************************************
  matrix exponential through scale and square trotter-like method
  plus Taylor expansion. The order of the Taylor expansion and of
  the squaring is chosen with an extimate according to 
  "19 dubious ways" paper
******************************************************************/
#define __TB_DQ_ACCU 1e-2
template <class MC>
void exp(const MC& M, MC& R, const double eps=__TB_STD_EPS)
{
#ifdef BENCHMARK
    TBBenchmarks["mf_exp"].n_calls++; TBBenchmarks["mf_exp"].timer.start(); 
#endif
    /*************************************************************
      we are going to find exp[M] so that if instead of the exact 
      exponential we get exp[M+E], |E|/|M|<eps, with the least
      number of matrix multiplications possible
    *************************************************************/
    
    double dq, de, nm; 
    nm=norminf(M);
    //even a simplified expression to find the total number of
    //multiplies requires a trascendental equation. well, we solve it somehow
    //initial guess
    if (nm<=0.)  //exponential of the null matrix!
    {
        R.resize(M.rows(),M.cols());
        R*=0.;
        R+=1.;
        return;
    }
                
    double l2nm=dq=log(nm)/log(2.), le=log(eps)/log(2.);
    if (dq<1.) dq=1.;
    
    //newton method on the logarithm of the error for the solution
    double odq=dq+2*__TB_DQ_ACCU, u, m;
    while (fabs(dq-odq)>__TB_DQ_ACCU)
    {
        odq=dq; 
        u=sqrt(12.+l2nm*(l2nm-2*(dq+2))+dq*(4+dq));
        de=1./16.*(36-3*l2nm*l2nm+dq*(4-3*dq)-2*u*(2+dq)+u*u+
                2*l2nm*(u+3*dq-2))-log(2-l2nm+dq+u)/log(2.);
        m=-1./(4*u)*(4/log(2.) +8+ l2nm*l2nm-2*u+dq*(4+dq+u)-l2nm*(4+2*dq+u));
        dq=odq+(le-de)/m;
    }
    //this defines a minimum n. of operations for the exponential: one scaling of first order taylor!
    if (dq<2.) dq=2.;
    double dj=(0.25*(2+l2nm+3*dq-sqrt(12.+l2nm*(l2nm-2.*(dq+2.))+dq*(4.+dq))));
    if (dj<1) dj=1;
    double dk=dq-dj;

    //ally we find the best combination of k (order of taylor) and j(scale-square) given q.
    //we round up both to be on the safe side
    long k, j;
    j=(unsigned long) ceil(dj); k=(unsigned long) ceil(dk);

#ifdef DEBUG
    std::perr<<" * TEST MATRIX EXPONENTIAL: for |M|= "<<nm<<"  eventually we chose j= "<<j<<"  k= "<<k<<"\n";
#endif
    
    double s=1.;  for (unsigned long i=1; i<=j;++i) s*=2.;    //scaling factor
        
    std::valarray<double> taycf(k+1);
    taycf[0]=1.;  for (unsigned long i=1; i<=k;++i) taycf[i]=taycf[i-1]/i;   //coefficients for taylor expansion
    
    R.resize(0,0); MC tmp=R;
    poly(M*(1/s), taycf, R); 
    
    for (unsigned long i=1; i<=j;++i) 
    {
        if (i%2==1)  { mult(R,R,tmp); R.resize(0,0);}
        else { mult(tmp,tmp,R); tmp.resize(0,0);}
    }
    if (j%2==1) R=tmp;
#ifdef BENCHMARK
    TBBenchmarks["mf_exp"].timer.stop(); TBBenchmarks["mf_exp"].tot_time+=TBBenchmarks["mf_exp"].timer;
    TBBenchmarks["mf_exp"].tot_prop["mm_mult"]+=j+k-1;  //n. of multiplications
#endif
}

/*****************************************************************************
  compute iteratively a matrix's inverse by using a first order newton method.
  if a nonempty matrix is passed as the return value, this is used as initial
  guess, otherwise a guaranteed-to-converge-with-infinite-precision guess is
  generated. we use the norm of the residual |1-MR|_inf as a measure for 
  convergence,  but there are two issues: first, the method is diverging
  when the spectral radius of (1-MR) is greater than one, but we need a way
  of extimating this reliably. most matrix norms can get larger than one even 
  if the s.r. is smaller, so we need something which is a LOWER bound to 
  the s.r. we don't need a good extimate, just something we can use to bailout
  when we are diverging. abs(Tr A)/n is a good choice, cheap to compute as well
  after all, we hope that users are not so silly as trowing in guesses which are 
  not expected to converge, so we can waste a couple of iterations before 
  detecting divergency.
  convergency is computed based on inf-norm of the residual. one can use different
  combinations of targets, thresholds, ecc.
*******************************************************************************/
template <template<class U> class MC, class U, unsigned long N>
bool inverse_newton(const MC<U>& M, MC<U>& R, 
                    IterOptions<double,N> iops=IterOptions<double,N>(__TB_STD_MAXSTEP, __TB_STD_EPS, 0., ichk_default)
                    )
{
#ifdef BENCHMARK
    unsigned long nm=0;
    TBBenchmarks["mf_inv_newt"].n_calls++; TBBenchmarks["mf_inv_newt"].timer.start(); 
#endif
    
    bool fguess=false;
    std::cerr<<R.rows()<<"\n";
    //we can pass a guess for the inverse in G, if we pass an empty matrix, we start from scratch
    if (R.rows()!=0)
    {
        fguess=true;
    }
    else
    {
        //"guaranteed convergency" initial guess
        //std::perr<<"guessing safe\n";
        transpose(M,R); //beware! transpose is broke for PCRS....
        typedef U(*fmap_type) (const U&);
        double t=norminf(R)*norminf(M);
        toolbox::map<U>(R,fmap_type(conj));
        scale(R,1./t);
    } 
    typename MC<U>::index_type n=M.rows();
    
    MC<U> tmp=R; tmp.resize(0,0); MC<U> tmp2=tmp;
    double errn;
    
    mult(R,M,tmp);
    while (!iops)
    {
        
#ifdef BENCHMARK
        nm+=2;
#endif
        mult(tmp,R,tmp2); 
        scale(R,-2.);
        incr(R, tmp2);
        neg(R);  //here we have the new guess
        
        mult(R,M,tmp);  //by computing this here, we get the error 'almost' for free
        tmp-=1.;  //for (typename MC::index_type i=0; i<n; ++i) tmp(i,i)-=1.;
        errn=norminf(tmp); 
        //check for divergence
        if (abs(trace(tmp))/n >1.) 
        {
            //it is not conceptually right to fallback on 'safe guess'.
            //we must return false, then caller may recall w/o guess
/*            if (!fguess) ERROR("Matrix inversion is diverging even with 'safe guess'. M should be VERY ill-conditioned");
            R.resize(0,0);
#ifdef BENCHMARK
            TBBenchmarks["mf_inv_newt"].timer.stop(); TBBenchmarks["mf_inv_newt"].tot_time+=TBBenchmarks["mf_inv_newt"].timer;
            TBBenchmarks["mf_inv_newt"].tot_prop["mm_mult"]+=nm;  //n. of multiplications
#endif
            inverse_newton(M,R,iops); //try to invert from scratch
            */
            return false;
        }
        for (unsigned long i=0; i<N; ++i)
            iops.setval(errn,i);
#ifdef DEBUG
        std::cerr<<"error: "<<errn<<"  nR: " <<norm1(R)<<"\n";
        std::perr<<"el. per row: "<< R.size()/n<<"  error: "<<errn<<" TRACE: "<<trace(R)<<"\n";
#endif
        //recover the shift
        tmp+=1.; //for (typename MC::index_type i=0; i<n; ++i) tmp(i,i)+=1.;
    }
    if (iops.maxstep_reached()) { return false; } //ERROR("Newton inversion reached max number of steps without converging within desired accuracy\n");
#ifdef BENCHMARK
    TBBenchmarks["mf_inv_newt"].timer.stop(); TBBenchmarks["mf_inv_newt"].tot_time+=TBBenchmarks["mf_inv_newt"].timer;
    TBBenchmarks["mf_inv_newt"].tot_prop["mm_mult"]+=nm;  //n. of multiplications
#endif
    return true;
}

/*****************************************************************
                  CHEBYSHEV POLYNOMIALS ROUTINES
******************************************************************/
template <class MC, class U>
void chebyshev_poly(const MC& M, const std::valarray<std::valarray<U>* >& b, std::valarray<MC*>& R )
{
#ifdef DEBUG
    if (M.cols()!=M.rows()) ERROR("Input matrix must be square!");
#endif

    unsigned long ncmp=b.size();
#ifdef DEBUG
    if (R.size()!=ncmp) ERROR("Return array must be initialized to as many valid pointers to matrices as the size of coefficients list");
#endif
                
    typename MC::index_type n=M.rows(), q;
#ifdef DEBUG
    unsigned long q0=b[0]->size()-1;
#endif
    
    MC tmp, T, Told;
    MatrixOptions<MC> rops; (*R[0]).getops(rops); T.setops(rops); Told.setops(rops); tmp.setops(rops);
    tmp.resize(n,n); 
    
    for (unsigned long icmp=0; icmp<ncmp; ++icmp)
    {
        q=b[icmp]->size()-1;
#ifdef DEBUG
        if (q!=q0) ERROR("All the polynomials must be of the same degree");
#endif
        
        (*R[icmp]).resize(n,n);
        if (q==0) {  (*R[icmp])+=(*b[icmp])[0];  continue;  }
        
        incr((*R[icmp]),M,(*b[icmp])[1]);    // R=b[0] T_0(M) + b[1] T_1(M)
        (*R[icmp])+=(*b[icmp])[0];
    }
 
    if (q<2) return;
    
    //sets up T_0 and T_1
    T.resize(n,n); T*=0.; incr(T,M); Told.resize(n,n); Told*=0.;
    Told+=1.;
    
    for (unsigned int i=2; i<=q; ++i)
    {
        //makes up T[i]
        mult(T,M,tmp); 
        tmp*=2.; 
        incr(tmp,Told,-1.);
        Told=T; T=tmp;
        //we increment the results
        for (unsigned long icmp=0; icmp<ncmp; ++icmp)
            incr((*R[icmp]),T,(*b[icmp])[i]);
        
        std::perr<<"# CHEB POLY STEP "<<i<<" el. per row: "<<R[0]->size()/n<<"  norm: "<<norminf(*R[0])<<
                "  trace: "<<trace(*R[0])<<"\n";
    }
}

template <class MC, class U>
void chebyshev_poly(const MC& M, const std::valarray<U>& b, MC& R )
{
    const std::valarray<std::valarray<U>* > rb((const_cast<std::valarray<U>* >(&b)),1);
    std::valarray< MC* >rR(&R,1);
    chebyshev_poly(M, rb, rR);
}

template <class MC, class U>
void chebyshev_fastpoly(const MC& M, const std::valarray< std::valarray<U> * >& b, std::valarray<MC*>& R )
{
#ifdef BENCHMARK
    unsigned long nm=0;
    TBBenchmarks["mf_fastcheb"].n_calls++; TBBenchmarks["mf_fastcheb"].timer.start(); 
#endif
#ifdef DEBUG
    if (M.cols()!=M.rows()) ERROR("Input matrix must be square!");
#endif
    
    unsigned long ncmp=b.size();
#ifdef DEBUG
    if (R.size()!=ncmp) ERROR("Return array must be initialized to as many valid pointers to matrices as the size of coefficients list");
#endif

    
    typename MC::index_type N=M.rows();
    unsigned long q=b[0]->size()-1;
    unsigned long s=(unsigned long) floor(sqrt(1. + q)), r=(unsigned long) (1+floor(1.*q/sqrt(1.+q)));
    if (q>r*s) ++s; unsigned long ns=q-(s-1)*r;
    unsigned long n,m;
    //here comes a painful part, we must compute the cheb coefficients. the original paper is
    //complete crap, so full of misprints and inaccuracies I wasted one day redoing all the stuff
    //from scratch. the procedure down here at least is working.
    std::valarray<std::valarray<double> > c1(std::valarray<double>(0.,r+1),s), c2(c1);
    //A matrix can be initialized here.
    std::valarray<std::valarray<double> > A(std::valarray<double>(0.,s),s);
    A[0][0]=1.;
    if (s>1) A[1][1]=1.; if (s>2) { A[2][2]=2.; A[0][2]=-1.; }
    if (s>3) for(n=3; n<s; ++n) for(m=0; m<s; ++m) A[m][n]=(m>0?2.*A[m-1][n-1]:0.)-(n==m || m>n-2?0.:A[m][n-2]); 
    
    MC tmp;
    MatrixOptions<MC> rops; (*R[0]).getops(rops); tmp.resize(N,N); tmp.setops(rops);
    
    std::valarray<MC> T(tmp,r+1);
    //sets up T_0, T_1 and the cheb polynomials up to r
    T[0]+=1.; if (r>0) T[1]+=M;
    
    T[1].getops(rops);
    for (unsigned int i=2; i<=r; ++i)
    {
        //makes up T[i]
        mult(T[i-1],M,T[i]); T[i]*=2.; incr(T[i],T[i-2],-1.);
    }
    
#ifdef BENCHMARK
    nm+=r-1;
#endif
    
    //unsigned long sr=(s-1)*r;
    for (unsigned long icmp=0; icmp<ncmp; ++icmp)
    {
#ifdef DEBUG
        if (b[icmp]->size()-1!=q) ERROR("All the polynomials must be of the same degree");
#endif
        //here we do compute the c coefficients from the b. once again, trust this and not liang03
        std::valarray<U>& c=*b[icmp];
        for (n=0; n<s; ++n) c1[n]=c2[n]=0.;  //init, zeroing out
        for(unsigned long i=1; i<=ns; ++i) c1[s-1][i]=2.*c[(s-1)*r+i];
        for(n=s-2; n>0; --n) c1[n][r]=2.*c[(n+1)*r]-(n<s-2?c1[n+2][r]:0.);
        c1[0][r]=c[r]-(s>2?0.5*c1[2][r]:0.);
        for (n=s-2;n>0;--n) for(unsigned long i=r-1; i>0;--i) c1[n][i]=2.*c[n*r+i]-c1[n+1][r-i];
        for(long i=0; i<(long) r; ++i) c1[0][i]=c[i]-0.5*c1[1][r-i];
        
        for(n=0;n<s;++n)
            for (unsigned long i=0; i<=r; ++i)
                for (m=n;m<s;++m)
                    c2[n][i]+=c1[m][i]*A[n][m];
        
        (*R[icmp]).resize(N,N);  (*R[icmp])*=0.;
        for (unsigned long i=1; i<=ns; ++i) incr(*R[icmp],T[i],c2[s-1][i]); 
    
        for (long k=s-2; k>=0; --k)
        {
            mult(*R[icmp],T[r],tmp); 
#ifdef BENCHMARK
            nm++;
#endif
            for (unsigned long i=1; i<=r; ++i)
                incr(tmp,T[i],c2[k][i]);
            (*R[icmp])=tmp;
            std::perr<<"# CHEB FASTPOLY STEP "<<k<<" el. per row: "<<R[icmp]->size()/N<<"  norm: "<<norminf(*R[icmp])<<"\n";
        }
        incr(*R[icmp],c2[0][0]);
    }
#ifdef BENCHMARK
    TBBenchmarks["mf_fastcheb"].timer.stop(); TBBenchmarks["mf_fastcheb"].tot_time+=TBBenchmarks["mf_fastcheb"].timer;
    TBBenchmarks["mf_fastcheb"].tot_prop["mm_mult"]+=nm;  //n. of multiplications
#endif
}

template <class MC, class U>
void chebyshev_fastpoly(const MC& M, const std::valarray<U>& b, MC& R )
{
    const std::valarray<std::valarray<U>* > rb((const_cast<std::valarray<U>* >(&b)),1);
    std::valarray< MC* >rR(&R,1);
    chebyshev_fastpoly(M, rb, rR);
}




/******************************************************************
  matrix polynomial computation using ~2 sqrt(2 q) matrix multiplications,
  but requiring to hold only O(1) indermediate matrices. This
  is much tougher to implement properly, and to keep sparse and
  easy to parallelize!!! so here I just list thoughts about how
  it could be implemented:
  1) sparse-matrix vector multiply is needed. probably the most
  long-term approach would be to implement a sparse vector class
  in the first place
  2) the result could be written as if rowwise and then transposed
  3) also the necessary parallelism can be addressed at the 
  sparse vector class level!!!
******************************************************************/
/*template <class MC, class U>
void poly_vanloan2(const MC& M, const std::valarray<U>& b, MC& R )
{
#ifdef DEBUG
    if (M.cols()!=M.rows()) ERROR("Input matrix must be square!")
#endif
    }*/
} //ends namespace toolbox
#endif // ends ifdef __MATRIX_FUN_H
