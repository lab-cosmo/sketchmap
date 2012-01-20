/* Rudimentary tensor class
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/

#ifndef __TENSOR_FULL_H
#define __TENSOR_FULL_H

#include "tbdefs.hpp"
#include "matrices.hpp"
#include <valarray>
#include <string>
#include <sstream>
#include <iomanip>

namespace toolbox {
/****************************************************
 * A class for full tensors. This also contains    *
 * all the relevant functions which should be imple-*
 * mented by a matrix class.                        *
 ****************************************************/

template <class U, unsigned long order=1u>
class FTensor
{
public:
    typedef unsigned long index_type;
    typedef U data_type;
    typedef fixarray<index_type, order> point_type;

private:
    point_type wdim; 
    std::valarray<U> data;

public:
    //matrix properties
    index_type dim( index_type i ) const 
    { 
#ifdef DEBUG
        if (i>=order) ERROR("Dimension index out of bounds")
#endif
        return wdim[i]; 
    }
    void dims(point_type& d) const { d=wdim; }
    index_type size() const { index_type sz=wdim[0]; for (index_type i=1; i<order; i++) sz*=wdim[i]; return sz; }
    
    inline data_type& operator()(const point_type& p) 
    { 
#ifdef DEBUG
        for (index_type i=0; i<order; ++i) if (p[i]>=wdim[i]) ERROR("Element index "<<i<<" out of bounds");
#endif
        index_type pos=p[0]; for (index_type i=1; i<order; ++i) pos=pos*wdim[i]+p[i];
        return data[pos]; 
    }
    
    inline const data_type& operator()(const point_type& p) const
    { 
#ifdef DEBUG
        for (index_type i=0; i<order; ++i) if (p[i]>=wdim[i]) ERROR("Element index "<<i<<" out of bounds");
#endif
        index_type pos=p[0]; for (index_type i=1; i<order; ++i) pos=pos*wdim[i]+p[i];
        return data[pos]; 
    }
    
    //variadic access
    inline data_type& operator()(const index_type first, ...)
    {
        index_type pos=first; 
        if (order<=1) return data[pos];
        
        va_list ap; va_start(ap, first);
        for (unsigned long i=1; i<order; ++i) pos=pos*wdim[i]+va_arg(ap, index_type);
        va_end(ap);
        return data[pos];
    }
    
    inline const data_type& operator()(const index_type first, ...) const 
    {
        index_type pos=first; 
        
        if (order<=1) return data[pos];
        va_list ap; va_start(ap, first);
        for (unsigned long i=1; i<order; ++i) pos=pos*wdim[i]+va_arg(ap, index_type);
        va_end(ap);
        return data[pos];
    }
    
    
    void resize(const point_type& ndim)
    { wdim=ndim; data.resize((*this).size(), (U) 0.); }
    
    //copy operator
    FTensor(const FTensor<U,order>& a) { resize(a.wdim); data=a.data; }
    
    FTensor<U,order>& operator=(const FTensor<U,order>& a)
    { if (this==&a) return *this; resize(a.wdim); data=a.data; return *this;}
    
    template <class T> 
    FTensor<U,order>& operator=(const FTensor<T,order>& a)
    { 
        if (this==(FTensor<U>*) &a) return *this; 
        resize(a.wdim); 
        for(index_type i=0; i<data.size(); ++i) data[i]=(data_type) a.data[i]; 
        return *this;
    }
    
    FTensor(
            const point_type& ndim=(fixarray<index_type, order >)(0), 
            const data_type& v=(data_type) 0.) : wdim(ndim) {  data.resize(size()); if (size()>0) data=v; }

    /*
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
    */
}; //ends class CrsMatrix

/*
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
}*/
} //ends namespace toolbox
#endif //ends ifdef __TENSOR_FULL_H
