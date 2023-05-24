/* Full-matrix library
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/

#ifndef __MATRIX_FULL_H
#define __MATRIX_FULL_H

#include "tbdefs.hpp"
#include "matrices.hpp"
#include <valarray>
#include <string>
#include <sstream>
#include <iomanip>

namespace toolbox {
/****************************************************
 * A class for full-matrices. This also contains    *
 * all the relevant functions which should be imple-*
 * mented by a matrix class.                        *
 ****************************************************/
template <class U>
        class MatrixOptions<FMatrix<U> >
{
    public:
        double atthresh;
        bool atrel;
        ATNormMode atnorm;
        bool atnodiag;    //if it is set to true, we don't include the diagonal into the computation of the norms.
    
MatrixOptions(double natthresh=0., bool natrel=false, bool natnodiag=false, ATNormMode natnorm=at_norminf) : 
        atthresh(natthresh), atrel(natrel), atnorm(natnorm), atnodiag(natnodiag) {}

template <class T> 
        MatrixOptions<FMatrix<U> >
        (MatrixOptions<FMatrix<T> > const& nops)
{ atthresh=nops.atthresh; atrel=nops.atrel; atnodiag=nops.atnodiag; atnorm=nops.atnorm; }
};


    
template <class U>
class FMatrix
{
    template <class T> friend class CoordMatrix;
    template <class T> friend class CrsMatrix;
    template <class T> friend class FMatrix;
    friend class IField<FMatrix<U> >;
    friend class IField<const FMatrix<U> >;
    template <class T> friend bool operator == (const FMatrix<T>& A, const FMatrix<T>& B);
    template <class T, class V> friend void scale(FMatrix<T>& a, const V& b);
    template <class T> friend void neg(FMatrix<T>& a);
    template <class T> friend void map(FMatrix<T>& a, T (*f)(const T));
    template <class T> friend void map(FMatrix<T>& a, T (*f)(const T&));
    template <class T, class V> friend void copymap(const FMatrix<T>& b, FMatrix<V>&a, V (*f)(const T));
    template <class T, class V> friend void copymap(const FMatrix<T>& b, FMatrix<V>&a, V (*f)(const T&));
    template <class T, class V> friend void transpose(const FMatrix<T>& a, FMatrix<V>& b);
    template <class T, class V> friend void incr(FMatrix<T>& a, const FMatrix<V>& b);
    template <class T, class V> friend void incr(FMatrix<T>& a, const V& b);
    template <class T, class V> friend void incr(FMatrix<T>& a, const std::valarray<V>& b);    
    template <class Z, class V, class T> friend void incr(FMatrix<Z>& a, const FMatrix<V>& b, const T& s);
    template <class T, class V> friend void add(const FMatrix<T>& b, const FMatrix<V>& c, FMatrix<T>& a);
    template <class T, class V> friend void mult(const FMatrix<T>& b, const FMatrix<V>& c, FMatrix<T>& a);
    friend double norm1<>(const FMatrix<U>& a);
    friend double norminf<>(const FMatrix<U>& a);
    friend double normfrob<>(const FMatrix<U>& a);
    friend U trace<>(const FMatrix<U>& a);
    
    friend std::ostream& operator<< <> (std::ostream& ostr, const FMatrix& cm);
    
public:
    typedef unsigned long index_type;
    typedef U data_type;
    
private:
    index_type wc, wr;
    std::valarray<U> data;
    MatrixOptions<FMatrix> pops;
    
public:
    //matrix properties
    index_type rows() const { return wr; }
    index_type cols() const { return wc; }
    index_type size() const { return wc*wr; }
        
    //extra properties transfer functions
    inline void do_at()
    {
        double thf;
        if (pops.atthresh==0.) return;
        
        if (!pops.atrel) thf=1.;
        else 
        {
            switch(pops.atnorm)
            {
                case at_normone:
                    thf=norm1(*this); break;
                case at_norminf:
                    thf=norminf(*this); break;
                case at_normfrob:
                    thf=normfrob(*this); break;
            }
            if (pops.atshift<thf) thf-=pops.atshift; else thf=0.; //shift if required, but if shift brings to negative, don't truncate
        }
        trunc(pops.atthresh*thf);
    }

    MatrixOptions<FMatrix> getops() const  { return pops; }
    void getops(MatrixOptions<FMatrix>& rops) const  { rops=pops; }
    void setops(const MatrixOptions<FMatrix>& rops)
    { 
#ifdef DEBUG
        if (rops.atthresh<0.) ERROR("Negative truncation threshold")
#endif
        
                    pops=rops;
    }
    
    //element access functions
    inline bool exists(const index_type& i, const index_type& j) const
    { 
#ifdef DEBUG
        if (i>=wr || j>=wc) ERROR("Element index out of bounds")
#endif        
        return true; 
    }

    inline void remove(const index_type& i, const index_type& j) 
    {
#ifdef DEBUG
        if (i>=wr || j>=wc) ERROR("Element index out of bounds")
#endif
        (*this).access(i,j)=(U) 0.; 
    }

    inline void insert(const index_type& i, const index_type& j, const data_type& val) 
    { 
#ifdef DEBUG
        if (i>=wr || j>=wc) ERROR("Element index out of bounds")
#endif        
        (*this).access(i,j)=val; 
    }
    
    inline data_type& access(const index_type& i, const index_type& j) 
    { 
#ifdef DEBUG
        if (i>=wr || j>=wc) ERROR("Element index out of bounds")
#endif
        return data[i*wc+j]; 
    }
    inline const data_type& access(const index_type& i, const index_type& j) const
    { 
#ifdef DEBUG
        if (i>=wr || j>=wc) ERROR("Element index out of bounds")
#endif
        return data[i*wc+j]; 
    }
    
    inline std::slice_array<U> all()
    {
        return data[std::slice(0,data.size(),1)];
    }
    
    inline std::slice_array<U> row(const index_type& i)
    {
        return data[std::slice(i*wc,wc,1)];
    }
    
    inline std::slice_array<U> col(const index_type& j)
    {
        return data[std::slice(j,wr,wc)];
    }

    inline std::slice_array<U> diag()
    {
        return data[std::slice(0,wr<wc?wr:wc,wc+1)];
    }
            
    void resize(const index_type& r, const index_type& c)
    { wr=r; wc=c; data.resize(r*c, (U) 0.); }
    
    //constructors
    FMatrix(const index_type& r=0, const index_type& c=0, const data_type& v=(data_type) 0.) : wc(c), wr(r),  data(v, r*c) {}
    FMatrix(const FMatrix<U>& a) :  wc(a.wc), wr(a.wr), data(a.data), pops(a.pops) {}
    template <class T> 
    FMatrix(const FMatrix<T>& a) : data(a.data.size()), wc(a.wc), wr(a.wr), pops(a.pops) 
    { for(index_type i=0; i<data.size(); ++i) data[i]=(data_type) a.data[i]; }
    
    FMatrix(const index_type& r, const index_type& c,
              const std::valarray<std::valarray<index_type> >& iind, 
              const std::valarray<std::valarray<data_type> >& idata
             ) :  wc(c), wr(r), pops()
    {
        resize(r,c); 
        
#ifdef DEBUG
        if (wr!=iind.size() || wr!=idata.size()) ERROR("Data arrays does not match number of rows") 
        for(index_type i=0; i<wr; ++i) 
            if (iind[i].size()!=idata[i].size()) ERROR("Data and indices arrays mismatch")
#endif
        for(index_type i=0; i<wr; ++i) 
        { 
            for (index_type j=0;j<idata[i].size(); ++j) access(i,iind[i][j])=idata[i][j];
        }
    }
    
    //copy operator
    FMatrix<U>& operator=(const FMatrix<U>& a)
    { if (this==&a) return *this; resize(a.wr, a.wc); data=a.data; pops=a.pops; return *this;}
    
    template <class T> 
    FMatrix<U>& operator=(const FMatrix<T>& a)
    { 
        if (this==(FMatrix<U>*) &a) return *this; 
        resize(a.wr, a.wc); 
        pops=a.pops;
        for(index_type i=0; i<data.size(); ++i) data[i]=(data_type) a.data[i]; 
        return *this;
    }
    
    //conversions from other matrix types
    template <class T> FMatrix(const std::valarray<std::valarray<T> >& s);
    template <class T> FMatrix(const CrsMatrix<T>& s);
    template <class T> FMatrix(const CoordMatrix<T>& s);
    
    operator std::string() 
    {
        std::ostringstream ostr;
        
        ostr.setf(std::ios::scientific);
        ostr.setf(std::ios::right);
        ostr.precision(5);
        
        index_type k=0;
        for (index_type i=0; i<wr; ++i)
        {
            for (index_type j=0; j<wc; ++j)
            { ostr<<std::setw(12)<<data[k]<<" "; ++k; }
            ostr<<"\n";
        }
        
        return ostr.str();
    }
    
    //matrix operations
    inline void trunc(const double& thresh)
    {
#ifdef DEBUG
        if (thresh<0.) ERROR("Negative truncation threshold")
#endif
        
        for (index_type k=0; k<data.size(); ++k) if (abs(data[k])<=thresh) data[k]=(U) 0.;
    }
    
    //overloaded operators
    inline U& operator() (const unsigned long& i, const unsigned long& j) 
    { if (!exists(i,j))  insert(i,j,(U) 0.);  return access(i,j); }
    
    inline U operator() (const unsigned long& i, const unsigned long& j)  const
    { if (exists(i,j)) return access(i,j); else return 0;}
       
    template <class V> inline FMatrix<U>& operator += (const V& s)
    { incr(*this, s); return *this; }
    
    template <class V> inline FMatrix<U>& operator -= (const V& s)
    { incr(*this, -s); return *this; }
    
    template <class V> inline FMatrix<U>& operator *= (const V& s)
    { scale(*this, s); return *this; }
    
    inline const FMatrix<U> operator -() 
    {  FMatrix<U> r(*this); neg(r); return r; }
    
    template <class V> inline FMatrix<U>& operator += (const FMatrix<V>& a)
    { incr(*this,a); return *this; }
    
    template <class V> inline FMatrix<U>& operator -= (const FMatrix<V>& a)
    { neg(*this); incr(*this,a); neg(*this); return *this; }  //this way we avoid creating a temporary
    
    template <class V> inline FMatrix<U>& operator *= (const FMatrix<V>& a) //we cannot avoid a temporary!
    { FMatrix<U> r; r=*this; mult(r,a,*this); return *this; }
    
}; //ends class CrsMatrix

template <class V, class T> inline const FMatrix<V> operator + (const FMatrix<V>& a, const T& b)
{ FMatrix<V> r(a); incr(r,b); return r; }
    
template <class V, class T> inline const FMatrix<V> operator + (const FMatrix<V>& a, const FMatrix<T> b)
{ FMatrix<V> r(a); incr(r,b); return r; }

template <class V, class T> inline const FMatrix<V> operator * (const FMatrix<V>& a, const T& b)
{ FMatrix<V> r(a); scale(r,b); return r; }

template <class V, class T> inline const FMatrix<V> operator * (const T& b, const FMatrix<V>& a)
{ FMatrix<V> r(a); scale(r,b); return r; }

template <class V, class T> inline const FMatrix<V> operator * (const FMatrix<V>& a, const FMatrix<T> b)
{ 
    FMatrix<V> r; double att; MatrixOptions<FMatrix<V> > mops; 
    a.getops(mops); r.setops(mops);  mult(a,b,r); return r; 
}

template <class U> bool operator == (const FMatrix<U>& A, const FMatrix<U>& B) 
{ 
	if (A.wc != B.wc || A.wr != B.wr) return false;
	for (unsigned long i=0; i<A.data.size(); ++i)
	  if (A.data[i]!=B.data[i]) return false;
	return true;
}

template <class U, class V> 
void incr(FMatrix<U>& a, const FMatrix<V>& b)
{
#ifdef BENCHMARK
    TBBenchmarks["full_incr_mtx"].n_calls++; TBBenchmarks["full_incr_mtx"].timer.start(); 
#endif
#ifdef DEBUG
    if (a.wr!=b.wr || a.wc!=b.wc) ERROR("Incompatible matrix size")
#endif 
    for (unsigned long k=0; k<a.data.size(); ++k)
        a.data[k]+=b.data[k];
#ifdef BENCHMARK
    TBBenchmarks["full_incr_mtx"].timer.stop(); TBBenchmarks["full_incr_mtx"].tot_time+=TBBenchmarks["full_incr_mtx"].timer;
#endif
}

template <class U, class V> void incr(FMatrix<U>& a, const V& b)
{
#ifdef BENCHMARK
        TBBenchmarks["full_incr_scal"].n_calls++; TBBenchmarks["full_incr_scal"].timer.start(); 
#endif
#ifdef DEBUG
    if (a.wr!=a.wc) ERROR("Adding a diagonal part to a non-square matrix is nonsense.");
#endif 
    typename FMatrix<U>::index_type k, dk=a.wr+1;
    for (k=0;k<a.data.size();k+=dk) { a.data[k]+=b; }
#ifdef BENCHMARK
    TBBenchmarks["full_incr_scal"].timer.stop(); TBBenchmarks["full_incr_scal"].tot_time+=TBBenchmarks["full_incr_scal"].timer;
#endif
}

template <class U, class V> void incr(FMatrix<U>& a, const std::valarray<V>& b)
{
#ifdef DEBUG
    if (a.wr!=a.wc) ERROR("Adding a diagonal part to a non-square matrix is nonsense.")
    if (a.wr!=b.size()) ERROR("Incompatible matrix and vector size.")
#endif 
    typename FMatrix<U>::index_type k=0, i, dk=a.wr+1;
    for (i=0;i<a.wr;++i) { a.data[k]+=b[i]; k+=dk; }
}

template <class U, class V, class T>
void incr(FMatrix<U>& a, const FMatrix<V>& b, const T& s) 
{
#ifdef BENCHMARK
    TBBenchmarks["full_incr_mtx_s"].n_calls++; TBBenchmarks["full_incr_mtx_s"].timer.start(); 
#endif
#ifdef DEBUG
    if (a.wr!=b.wr || a.wc!=b.wc) ERROR("Incompatible matrix size")
#endif 
    typename FMatrix<U>::index_type k;
    for (k=0; k<a.data.size();++k) { a.data[k]+=b.data[k]*s; }
#ifdef BENCHMARK
    TBBenchmarks["full_incr_mtx_s"].timer.stop(); TBBenchmarks["full_incr_mtx_s"].tot_time+=TBBenchmarks["full_incr_mtx_s"].timer;
#endif
}
    
template <class U, class V> void scale(FMatrix<U>& a, const V& b) { a.data*=b; }
template <class T> void neg(FMatrix<T>& a) { a.data=-a.data; }
template <class T> void map(FMatrix<T>& a, T (*f)(const T)) { for (typename FMatrix<T>::index_type k=0; k<a.data.size();++k) { a.data[k]=f(a.data[k]); } }
template <class T> void map(FMatrix<T>& a, T (*f)(const T&)) { for (typename FMatrix<T>::index_type k=0; k<a.data.size();++k) { a.data[k]=f(a.data[k]); } }
template <class U, class V> void transpose(const FMatrix<U>& a, FMatrix<V>& b)
{
    //yeah, we should do it blockwise for cache locality. consider this just a stub
    b.resize(a.wc,a.wr);
    for (typename FMatrix<U>::index_type i=0; i<a.wr;++i)
        for (typename FMatrix<U>::index_type j=0; j<a.wc;++j)  
        { b(j,i)=a(i,j); }
}

template <class U> double norm1(const FMatrix<U>& a)
{
    std::valarray<double> maxc((double) 0., a.wc);
    typename FMatrix<U>::index_type k=0;
    for (typename FMatrix<U>::index_type i=0; i<a.wr; ++i)
        for (typename FMatrix<U>::index_type j=0; j<a.wc; ++j)
            maxc[j]+=abs(a.data[k++]);
    
    return maxc.max();
}
template <class U> double norminf(const FMatrix<U>& a)
{    
#ifdef BENCHMARK
    TBBenchmarks["full_norminf"].n_calls++; TBBenchmarks["full_norminf"].timer.start(); 
#endif
    double max=0, rtot;
    typename FMatrix<U>::index_type k=0;
    for (typename FMatrix<U>::index_type i=0; i<a.wr; ++i)
    {
        rtot=0;
        for (typename FMatrix<U>::index_type j=0; j<a.wc; ++j) rtot+=abs(a.data[k++]);
        if (rtot>max) max=rtot;
    }
    
#ifdef BENCHMARK
    TBBenchmarks["full_norminf"].timer.stop(); TBBenchmarks["full_norminf"].tot_time+=TBBenchmarks["full_norminf"].timer;
#endif
    return max;
}
template <class U> double normfrob(const FMatrix<U>& a)
{ 
    double nf=0,av; 
    for (typename FMatrix<U>::index_type k=0; k<a.data.size();++k) nf+=(av=abs(a.data[k]))*av;
    return sqrt(nf);
}

template <class U>  U trace(const FMatrix<U>& a)
{
#ifdef BENCHMARK
    TBBenchmarks["full_trace"].n_calls++; TBBenchmarks["full_trace"].timer.start(); 
#endif
#ifdef DEBUG
        if (a.wr!=a.wc) ERROR("Making the trace of a non-square matrix is nonsense.")
#endif 
    U rtr=0;
    typename FMatrix<U>::index_type k, dk=a.wr+1;
    for (k=0;k<a.data.size();k+=dk) { rtr+=a.data[k]; }
#ifdef BENCHMARK
    TBBenchmarks["full_trace"].timer.stop(); TBBenchmarks["full_trace"].tot_time+=TBBenchmarks["full_trace"].timer;
#endif
    return rtr;
}

template <class T, class V> 
void copymap(const FMatrix<T>& b, FMatrix<V>&a, V (*f)(const T))
{
    a.resize(b.wr, b.wc);
    for (typename FMatrix<V>::index_type k=0;k<a.data.size();++k) { a.data[k]=f(b.data[k]); }
}

template <class T, class V> 
void copymap(const FMatrix<T>& b, FMatrix<V>&a, V (*f)(const T&))
{
    a.resize(b.wr, b.wc);
    for (typename FMatrix<V>::index_type k=0;k<a.data.size();++k) { a.data[k]=f(b.data[k]); }
}


template <class U, class V> void add(const FMatrix<U>& b, const FMatrix<V>& c, FMatrix<U>& a)
{
#ifdef BENCHMARK
        TBBenchmarks["full_add"].n_calls++; TBBenchmarks["full_add"].timer.start(); 
#endif
#ifdef DEBUG
    if (c.wr!=b.wr || c.wc!=b.wc) ERROR("Incompatible matrix size")
#endif
    a.resize(b.wr,b.wc);
    a.data=b.data+c.data;
#ifdef BENCHMARK
    TBBenchmarks["full_add"].timer.stop(); TBBenchmarks["full_add"].tot_time+=TBBenchmarks["full_add"].timer;
#endif
}

template <class U, class V> void mult(const FMatrix<U>& b, const FMatrix<V>& c, FMatrix<U>& a)
{
#ifdef BENCHMARK
    TBBenchmarks["full_mult"].n_calls++; TBBenchmarks["full_mult"].timer.start(); 
#endif
#ifdef DEBUG
    if (c.wr!=b.wc) ERROR("Incompatible matrix size")
#endif
    U tel;
    a.resize(b.wr,c.wc);
    //std::cerr<<"MMMult: ";
    typename FMatrix<U>::index_type i,j,k,kb,kc;
    for (i=0; i<a.wr;++i)
    {
       // std::cerr<<".";
    for (j=0; j<a.wc;++j)
    {
        tel=(U) 0;
        kb=i*b.wc; kc=j;
        for (k=0;k<b.wc;++k) { tel+=b.data[kb++]*c.data[kc]; kc+=c.wc; }
        a.data[i*a.wc+j]=tel;
    }
    }
#ifdef BENCHMARK
    TBBenchmarks["full_mult"].timer.stop(); TBBenchmarks["full_mult"].tot_time+=TBBenchmarks["full_mult"].timer;
#endif
}
} //ends namespace toolbox
#endif //ends ifdef __MATRIX_FULL_H
