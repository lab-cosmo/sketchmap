/* Compressed row storage matrix format
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/


#ifndef __MATRIX_CRS_H
#define __MATRIX_CRS_H

#include "tbdefs.hpp"
#include "matrices.hpp"
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

#define __CRS_STORAGE_CHUNK 65536
#define __CRS_STORAGE_FREE  131072
namespace toolbox {
/****************************************************
 * A class for CRS-matrices.                        *
 ****************************************************/
template <class U>
class MatrixOptions<CrsMatrix<U> >
{
public:
    double atthresh;
    bool atrel;
    ATNormMode atnorm;
    bool atnodiag;    //if it is set to true, we don't include the diagonal into the computation of the norms.
        
    MatrixOptions(double natthresh=0., bool natrel=false, bool natnodiag=false, ATNormMode natnorm=at_norminf) : 
            atthresh(natthresh), atrel(natrel), atnorm(natnorm), atnodiag(natnodiag) {}

    template <class T> 
    MatrixOptions<CrsMatrix<U> >
    (MatrixOptions<CrsMatrix<T> > const& nops)
    { atthresh=nops.atthresh; atrel=nops.atrel; atnodiag=nops.atnodiag; atnorm=nops.atnorm; }
};

template <class U>
class CrsMatrix
{
    template <class T> friend class CoordMatrix;
    template <class T> friend class CrsMatrix;
    template <class T> friend class PCrsMatrix;
    template <class T> friend class FMatrix;
    
    template <class T, class V> friend void scale(CrsMatrix<T>& a, const V& b);
    template <class T> friend void neg(CrsMatrix<T>& a);
    template <class T> friend void map(CrsMatrix<T>& a, T (*f)(const T));
    template <class T> friend void map(CrsMatrix<T>& a, T (*f)(const T&));
    template <class T, class V> friend void copymap(const CrsMatrix<T>& b, CrsMatrix<V>&a, V (*f)(const T));
    template <class T, class V> friend void copymap(const CrsMatrix<T>& b, CrsMatrix<V>&a, V (*f)(const T&));
    template <class T, class V> friend void transpose(const CrsMatrix<T>& a, CrsMatrix<V>& b);
    template <class T, class V> friend void incr(CrsMatrix<T>& a, const CrsMatrix<V>& b);
    template <class T, class V> friend void incr(CrsMatrix<T>& a, const V& b);
    template <class T, class V> friend void incr(CrsMatrix<T>& a, const std::valarray<V>& b);    
    template <class Z, class V, class T> friend void incr(CrsMatrix<Z>& a, const CrsMatrix<V>& b, const T& s);
    template <class T, class V> friend void add(const CrsMatrix<T>& b, const CrsMatrix<V>& c, CrsMatrix<T>& a);
    template <class T, class V> friend void mult(const CrsMatrix<T>& b, const CrsMatrix<V>& c, CrsMatrix<T>& a);
    template <class T, class V> friend void mult(const CrsMatrix<T>& a, const std::valarray<V>& x, std::valarray<T>& y);
    template <class T, class V> friend void Tmult(const CrsMatrix<T>& a, const std::valarray<V>& x, std::valarray<T>& y);
    friend double norm1<>(const CrsMatrix<U>& a, const bool nodiag);
    friend double norminf<>(const CrsMatrix<U>& a, const bool nodiag);
    friend double normfrob<>(const CrsMatrix<U>& a, const bool nodiag);
    friend U trace<>(const CrsMatrix<U>& a);
    
    friend class toolbox::IField<CrsMatrix>;
    friend class toolbox::IField<const CrsMatrix>;
    friend std::ostream& operator<< <> (std::ostream& ostr, const CrsMatrix& cm);
    
    //pcrs matrix functions which should be able to access private CrsMatrix fields
    template <class T, class V> friend void incr(PCrsMatrix<T>& a, const V& b);
    friend U trace<>(const PCrsMatrix<U>& a);
    friend double normfrob<>(const PCrsMatrix<U>& a, const bool nodiag);
    friend double norminf<>(const PCrsMatrix<U>& a, const bool nodiag);
    template <class T, class V> friend void mult(const PCrsMatrix<T>& b, const PCrsMatrix<V>& c, PCrsMatrix<T>& a);
    template <class T, class V> friend void mult(const CrsMatrix<T>& a, const std::valarray<V>& x, std::valarray<T>& y);
    template <class T, class V> friend void Tmult(const CrsMatrix<T>& a, const std::valarray<V>& x, std::valarray<T>& y);
    template <class T, class V> friend void incr(PCrsMatrix<T>& a, const PCrsMatrix<V>& b);
    template <class Z, class V, class T> friend void incr(PCrsMatrix<Z>& a, const PCrsMatrix<V>& b, const T& s);
    template <class T, class V> friend void transpose(const PCrsMatrix<T>& a, PCrsMatrix<V>& b);
public:
    typedef unsigned long index_type;
    typedef U data_type;
    
private:
    index_type wc, wr;
    //! these variables are to help with making parallel CRS matrix operations consistent with varying n. of procs
    index_type RBASE; //!this is to implement CRSMatrix as a sub-block of another matrix. makes easy finding diagonal el.
    double LASTNORM;
    
    MatrixOptions<CrsMatrix<U> > pops;

    std::valarray<data_type> values;
    std::valarray<index_type> indices;
    std::valarray<index_type> rpoints;
    
    //resizes the storage non-destructively
    void presize(const index_type& nsz) 
    {
        int myrank=0;
#ifdef TB_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
#ifdef PCFS_MDBG
        std::ofstream cm;
        cm.open((std::string("cresize.")+int2str(myrank)).c_str(),std::ios_base::app);
        cm<<"presize: "<<indices.size()<<","<<rpoints[wr]<<">>"<<nsz<<std::endl;
#endif
        index_type sz=indices.size();
        if (rpoints[wr]<sz) sz=rpoints[wr];
        
#ifdef PCFS_MDBG
        cm<<"making temporaries"<<std::endl;
#endif
        std::valarray<data_type> dtmp(values);
        std::valarray<index_type> itmp(indices);
#ifdef PCFS_MDBG
        cm<<"resizing"<<std::endl;
#endif
        values.resize(nsz, (data_type) 0.);
        indices.resize(nsz, (index_type) 0.);
        
        if (nsz<sz) sz=nsz; //if we shrink we must not get out of the borders
#ifdef PCFS_MDBG
        cm<<"copying "<<sz<<std::endl;
#endif
        for (index_type i=0; i<sz; ++i) { values[i]=dtmp[i]; indices[i]=itmp[i]; }
#ifdef PCFS_MDBG
        cm<<"success"<<std::endl<<std::endl;
        cm.close();
#endif
    }
    
    
    inline index_type pinsert(const index_type& i, const index_type& j, const data_type& val, const index_type& prek=0) 
    { 
        //std::<< "P INSERTING "<<i<<","<<j<< "\n";
        #ifdef DEBUG
            if (i>=wr || j>=wc) ERROR("Element index out of bounds")
        #endif
        if (indices.size()<=rpoints[wr]) //we must grow the storage allocation
            presize(rpoints[wr]+__CRS_STORAGE_CHUNK);
        
        index_type k;
        if (rpoints[i]>prek) k=rpoints[i]; else k=prek;
        for (; j>indices[k] && k<rpoints[i+1]; ++k); 
        
        //we check for double insertions, and behave gracefully
        //the second check makes sure that the row is not empty, so that we are 
        //not fiddling with the next one
        if (j==indices[k] && k!=rpoints[i+1]) {values[k]=val; return k; }
        
        //std::cerr<<rpoints[wr]<<":"<<k<<">>"<<indices[k]<<"\n";
        //shifts values and indices
        for (index_type h=rpoints[wr]; h>k; --h)
        { indices[h]=indices[h-1]; values[h]=values[h-1]; }
        
        //correct pointers
        for (index_type h=i+1; h<=wr; ++h) ++rpoints[h];
        indices[k]=j; values[k]=val;
        return k;
    } 
    
    inline index_type pexists(const index_type& i, const index_type& j, bool& fex) 
    { 
        #ifdef DEBUG
            if (i>=wr || j>=wc) ERROR("Element index out of bounds")
        #endif
        
        //bisection search
        index_type ka=rpoints[i], kc=rpoints[i+1];
         
        if (ka==kc) { fex=false; return ka; }
        index_type kb=(ka+kc)/2, jb=indices[kb];
        while (ka<kb)
        {
            if (jb==j) {fex=true; return kb;}
            if (jb<j) ka=kb;
            else kc=kb;
            kb=(ka+kc)/2;
            jb=indices[kb];
        }

        fex=(jb==j); 
        return kb;
        
        /*
        index_type k;
        for (k=rpoints[i]; k<rpoints[i+1] && j>indices[k]; ++k);
        if (k>=rpoints[i+1]) { fex=false; return k; }
        if (indices[k]==j) { fex=true;  return k;}
        fex=false; return k; 
        */
    }
    

public:
    //sanity check
    bool sanity() const
    {
        //performs a sanity check over the matrix structure
        if (values.size()!=indices.size()) 
            std::cerr<<"indices and values mismatch\n";
        for (unsigned long i=0; i<wr;++i)
        {
            if(rpoints[i+1]>indices.size())
                std::cerr<<"row pointer overflow in row "<<i<<"\n";
            for (unsigned long k=rpoints[i]; k<rpoints[i+1]; ++k)
            {
                if (indices[k]>wc) 
                    std::cerr<<"col. index overflow on row "<<i<<" pos "<<k<<"\n";
                if (k>rpoints[i] && indices[k]<=indices[k-1]) 
                    std::cerr<<"index error on row "<<i<<" pos "<<k<<"\n";
            }
        }   
        return true;
    }

    /*******************************************************************************
      if a matrix is supposed to be constantly thresholded (i.e. all the elements
      below a given threshold are ripped off), often big savings can be made if 
      these elements are not inserted in the first place, when computing certain
      operations. different models are available; beware that setting autotresholding
      does not guarantee that a matrix is thresholded at any time. for example, elements
      insertion does not launch thresholding, but it can be enforced by do_at().
      anyway, autotresholding is implemented so that the result of a given operation
      is the same as the one obtained without thresholding, and applying a treshold
      on the untuncated result.
     *******************************************************************************/
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
                    thf=norm1(*this,pops.atnodiag); break;
                case at_norminf:
                    thf=norminf(*this,pops.atnodiag); break;
                case at_normfrob:
                    thf=normfrob(*this,pops.atnodiag); break;
            }
        }
        trunc(pops.atthresh*thf);
    }

    MatrixOptions<CrsMatrix<U> > getops() const  { return pops; }
    void getops(MatrixOptions<CrsMatrix<U> >& rops) const  { rops=pops; }
    void setops(const MatrixOptions<CrsMatrix<U> >& rops)
    { 
        #ifdef DEBUG
            if (rops.atthresh<0.) ERROR("Negative truncation threshold")
        #endif
        
        pops=rops;
    }
    
    //matrix properties
    index_type rows() const { return wr; }
    index_type cols() const { return wc; }    
    index_type size() const { return rpoints[wr]; }
        
    //element access functions
    inline bool exists(const index_type& i, const index_type& j) 
    { 
        #ifdef DEBUG
            if (i>=wr || j>=wc) ERROR("Element index out of bounds")
        #endif
        bool fex;
        pexists(i,j,fex); return fex;
    }
    
    inline void remove(const index_type& i, const index_type& j) 
    { 
        #ifdef DEBUG
            if (i>=wr || j>=wc) ERROR("Element index out of bounds")
        #endif
        for (index_type k=rpoints[i]; k<rpoints[i+1]; ++k) 
            if(indices[k]==j) 
            {   
                //shifts values and indices
                for (index_type h=k; h<rpoints[wr]-1; ++h)
                {
                    indices[h]=indices[h+1];
                    values[h]=values[h+1];
                }
                //correct pointers
                for (index_type h=i+1; h<=wr; ++h) --rpoints[h];
                
                //if needed, free up some memory
                if (indices.size()-rpoints[wr]>__CRS_STORAGE_FREE) presize(rpoints[wr]);
                return;
            }
        
        #ifdef DEBUG
            ERROR("Trying to remove a non-existing element")
        #endif
    }
    
    inline void insert(const index_type& i, const index_type& j, const data_type& val) {pinsert(i,j,val);} 
    
    inline data_type& access(const index_type& i, const index_type& j) 
    { 
#ifdef DEBUG
        if (i>=wr || j>=wc) ERROR("Element index out of bounds")
#endif
        bool fex; index_type kv;
        kv=pexists(i,j,fex);
        if (!fex) return values[pinsert(i,j,(data_type) 0.,kv)]; 
        else return values[kv];
    }

    void resize(const index_type& r, const index_type& c)
    { 
        wr=r; wc=c; 
        values.resize(__CRS_STORAGE_CHUNK, (data_type) 0.); 
        indices.resize(__CRS_STORAGE_CHUNK, (index_type) 0);
        rpoints.resize(r+1, (index_type) 0);
    }
    
    void optimize()
    {
        presize(size());
    }
    //constructors
    class jdpair{  public: index_type j; data_type d; 
        bool operator < (const jdpair &rhs) { return j< rhs.j; }
    };
    CrsMatrix(const index_type& r=0, const index_type& c=0) : wc(c), wr(r), RBASE(0), LASTNORM(0.), pops() { resize(r,c); }
    //from arrays of data. temporary solution until other sparse matrix formats are implemented
    CrsMatrix(const index_type& r, const index_type& c,
              const std::valarray<std::valarray<index_type> >& iind, 
              const std::valarray<std::valarray<data_type> >& idata
             ) :  wc(c), wr(r), RBASE(0), LASTNORM(0.), pops()
    { 
        resize(r,c); 
        
#ifdef DEBUG
        if (wr!=iind.size() || wr!=idata.size()) ERROR("Data arrays does not match number of rows") 
        for(index_type i=0; i<wr; ++i) 
            if (iind[i].size()!=idata[i].size()) ERROR("Data and indices arrays mismatch")
#endif
        index_type dtot=0;
        for(index_type i=0; i<wr; ++i) { rpoints[i]=dtot; dtot+=iind[i].size(); }
        rpoints[wr]=dtot;
        presize(dtot);
        index_type j,k=0;
        for(index_type i=0; i<wr; ++i) 
        { 
            j=0;
            //indices must be sorted! 
            std::valarray<jdpair> trow(iind[i].size()); 
            for(index_type h=0; h<iind[i].size(); ++h) 
            { trow[h].j=iind[i][h]; trow[h].d=idata[i][h]; }
            if (trow.size()>1) heapsort(trow);
            
            for (;k<rpoints[i+1]; ++k)
            { values[k]=trow[j].d; indices[k]=trow[j++].j; }
        }
    }
    /*
    CrsMatrix(const CrsMatrix<U>& a) 
    { 
        pops=a.pops;
        resize(a.wr,a.wc); 
        for(index_type i=0; i<=wr; ++i) rpoints[i]=a.rpoints[i];
        values.resize(rpoints[wr]); indices.resize(rpoints[wr]);
        for(index_type i=0; i<rpoints[wr]; ++i) { values[i]=a.values[i]; indices[i]=a.indices[i]; }
    }
    */
    template <class T> 
    CrsMatrix(const CrsMatrix<T>& a)
    { 
        RBASE=a.RBASE;
        LASTNORM=a.LASTNORM;
        pops=a.pops;
        resize(a.wr,a.wc); 
        for(index_type i=0; i<=wr; ++i) rpoints[i]=a.rpoints[i];
        values.resize(rpoints[wr]); indices.resize(rpoints[wr]);    
        for(index_type i=0; i<rpoints[wr]; ++i) { values[i]=(data_type)a.values[i]; indices[i]=a.indices[i]; }
    }
        
    //copy operator
    CrsMatrix<U>& operator = (const CrsMatrix<U>& a)
    { 
        if (this==&a) return *this;
        RBASE=a.RBASE;
        LASTNORM=a.LASTNORM;
        pops=a.pops;
        resize(a.wr,a.wc);         
        rpoints=a.rpoints;
        values.resize(rpoints[wr]); indices.resize(rpoints[wr]);
        for(index_type i=0; i<rpoints[wr]; ++i) { values[i]=a.values[i]; indices[i]=a.indices[i]; }
        return *this;
    }

    template <class T> 
    CrsMatrix<U>& operator = (const CrsMatrix<T>& a)
    { 
        if (this==(CrsMatrix<U>*) &a) return *this;
        RBASE=a.RBASE;
        LASTNORM=a.LASTNORM;
        resize(a.wr,a.wc); 
        pops=a.pops;
        rpoints=a.rpoints;
        values.resize(rpoints[wr]); indices.resize(rpoints[wr]);
        for(index_type i=0; i<rpoints[wr]; ++i) { values[i]=(data_type) a.values[i]; indices[i]=a.indices[i]; }
        return *this;
    }
    
    //comparison operator
    template<class T> bool operator == (const CrsMatrix<T>& a)
    {
        //this checks if the two matrices have the same values, NOT if the various members are 
        //strictly identical, i.e. the two matrices are equal even though they have different buffer sizes
        //or different auto-truncation thresholds
        if (wr!=a.wr || wc!=a.wc || a.rpoints[wr]!=rpoints[wr]) return false;
        for (index_type i=0 ; i<wr; ++i) if (a.rpoints[i]!=rpoints[i]) return false;
        for (index_type i=0 ; i<rpoints[wr]; ++i) if (a.indices[i]!=indices[i]) return false;
        for (index_type i=0 ; i<rpoints[wr]; ++i) if ((U) a.values[i]!=values[i]) return false;
        return true;
    }
    
    //conversion operators. these are just declarations, definitions are in the file matrix-conv.hpp
    //which must be included separately

    template<class T> CrsMatrix(const CoordMatrix<T>& s);
    template<class T> CrsMatrix(const FMatrix<T>& s);
    template<class T> CrsMatrix(const PCrsMatrix<T>& s);
    
    operator std::string() 
    {
        std::ostringstream ostr;
        
        ostr.setf(std::ios::scientific);
        ostr.setf(std::ios::right);
        ostr.precision(5);
        
        for (index_type i=0; i<wr; ++i)
        {
            for (index_type k=rpoints[i]; k<rpoints[i+1]; ++k)
            {
                ostr<<std::setw(5)<<k<<" "
                    <<std::setw(5)<<i<<" "<<std::setw(5)<<indices[k]<<" "
                    <<std::setw(12)<<values[k]<<"\n";
            }
        }
        return ostr.str();
    }
    
    //matrix operations
    inline void trunc(const double& thresh)
    {
        #ifdef DEBUG
            if (thresh<0.) ERROR("Negative truncation threshold")
        #endif

        index_type k=0, h=0, sh=0;
        for (index_type i=0; i<wr; ++i)
        {
            rpoints[i]-=sh;
            
            for (; k<rpoints[i+1]; ++k)
            {
                if (abs(values[k])<=thresh) {++sh; continue;}
                
                values[h]=values[k];
                indices[h]=indices[k];
                ++h;
            }
        }
        
        rpoints[wr]-=sh;
        if (sh>__CRS_STORAGE_FREE) presize(rpoints[wr]);
    }
    
    //matrix functions
    
    //overloaded operators, for syntactic sugar
    inline U& operator() (const unsigned long& i, const unsigned long& j) 
    {  return access(i,j); }
    
    template <class V> inline CrsMatrix<U>& operator += (const V& s)
    { incr(*this, s); return *this; }
    
    template <class V> inline CrsMatrix<U>& operator -= (const V& s)
    { incr(*this, -s); return *this; }
    
    template <class V> inline CrsMatrix<U>& operator *= (const V& s)
    { scale(*this, s); return *this; }
    
    inline const CrsMatrix<U> operator -() 
    {  CrsMatrix<U> r(*this); neg(r); return r; }
    
    template <class V> inline CrsMatrix<U>& operator += (const CrsMatrix<V>& a)
    { incr(*this,a); return *this; }
    
    template <class V> inline CrsMatrix<U>& operator -= (const CrsMatrix<V>& a)
    { neg(*this); incr(*this,a); neg(*this); return *this; }  //this way we avoid creating a temporary
    
    template <class V> inline CrsMatrix<U>& operator *= (const CrsMatrix<V>& a) //we cannot avoid a temporary!
    { CrsMatrix<U> r; r=*this; mult(r,a,*this); return *this; }
    
}; //ends class CrsMatrix

template <class V, class T> inline const CrsMatrix<V> operator + (const CrsMatrix<V>& a, const T& b)
{ CrsMatrix<V> r(a); incr(r,b); return r; }
    
template <class V, class T> inline const CrsMatrix<V> operator + (const CrsMatrix<V>& a, const CrsMatrix<T> b)
{ CrsMatrix<V> r(a); incr(r,b); return r; }

template <class V, class T> inline const CrsMatrix<V> operator * (const CrsMatrix<V>& a, const T& b)
{ CrsMatrix<V> r(a); scale(r,b); return r; }


template <class V, class T> inline const CrsMatrix<V> operator * (const T& b, const CrsMatrix<V>& a)
{ CrsMatrix<V> r(a); scale(r,b); return r; }

template <class V, class T> inline const CrsMatrix<V> operator * (const CrsMatrix<V>& a, const CrsMatrix<T> b)
{ 
    CrsMatrix<V> r; double att; MatrixOptions<CrsMatrix<V> > mops; 
    a.getops(mops); r.setops(mops);  mult(a,b,r); return r; 
}
//matrix operations

template <class U, class V>
void incr(CrsMatrix<U>& a, const CrsMatrix<V>& b) 
{
    //adds in place, this is more expensive if the two matrices have different
    //sparsity patterns, as some merging is necessary
#ifdef BENCHMARK
    TBBenchmarks["crs_incr_mtx"].n_calls++; TBBenchmarks["crs_incr_mtx"].timer.start(); 
    TBBenchmarks["crs_incr_mtx"].tot_avg["ops_sparsity"]+=(b.rpoints[b.wr]*1./b.wr+a.rpoints[a.wr]*1./a.wr)/2;  //mean number of el. per row of operands
#endif
#ifdef DEBUG
    if (a.wr!=b.wr || a.wc!=b.wc) ERROR("Incompatible matrix size")
#endif

    std::valarray<typename CrsMatrix<U>::index_type> a_insert((typename CrsMatrix<U>::index_type) 0, b.size());
    std::valarray<typename CrsMatrix<V>::index_type> b_insert((typename CrsMatrix<V>::index_type) 0, b.size());
    typename CrsMatrix<U>::index_type ki=0, ka=0, ma, mb;
    typename CrsMatrix<V>::index_type kb=0;
    double nguess=0., thresh=0., tr, ta;
    if (a.pops.atthresh<=0.)
    {
        //no truncation
        thresh=0.;
        ma=a.rpoints[0]; mb=b.rpoints[0];
        for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
        {
            ka=ma; kb=mb;
            a.rpoints[i]+=ki;
            ma=a.rpoints[i+1]; 
            mb=b.rpoints[i+1];
            for (; ka<ma; ++ka) 
            {   
                while (kb<mb && b.indices[kb]<a.indices[ka])
                { a_insert[ki]=ka; b_insert[ki]=kb; ++kb; ++ki; }
                if (kb==mb) break;
                if (b.indices[kb]==a.indices[ka]) {a.values[ka]+=(U) b.values[kb]; ++kb;}
            }
            while (kb<mb) { a_insert[ki]=ka; b_insert[ki]=kb; ++kb; ++ki; }
        }
    }
    else
    {
        if(!a.pops.atrel)
        {
            thresh=a.pops.atthresh;
            //absolute truncation
            for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
            {
                ka=a.rpoints[i]; kb=b.rpoints[i];
                a.rpoints[i]+=ki;
                for (; ka<a.rpoints[i+1]; ++ka) 
                {   
                    while (kb<b.rpoints[i+1] && b.indices[kb]<a.indices[ka])
                    { if (abs((U) b.values[kb])>thresh) { a_insert[ki]=ka; b_insert[ki]=kb; ++ki; } ++kb;  }
                    if (kb==b.rpoints[i+1]) continue;
                    if (b.indices[kb]==a.indices[ka]) {  a.values[ka]+=(U) b.values[kb]; ++kb;}
                }
                
                while (kb<b.rpoints[i+1]) { if (abs((U) b.values[kb])>thresh) { a_insert[ki]=ka; b_insert[ki]=kb; ++ki; } ++kb; }
            }
        }
        else
        {
            //relative truncation
            switch (a.pops.atnorm)
            {
            case at_norminf:
                nguess=thresh=0.;  //we could take as well |norm1(a)-norm1(b)|, but calculating that would cost so.... who knows.
                for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
                {
                    ka=a.rpoints[i]; kb=b.rpoints[i];
                    a.rpoints[i]+=ki; tr=0;
                    for (; ka<a.rpoints[i+1]; ++ka) 
                    {                   
                        while (kb<b.rpoints[i+1] && b.indices[kb]<a.indices[ka])
                        { 
                            ta=abs((U) b.values[kb]); if(!a.pops.atnodiag || b.indices[kb]!=i+a.RBASE) tr+=ta;
                            if (ta>thresh) { a_insert[ki]=ka; b_insert[ki]=kb; ++ki;} ++kb; 
                        }
                        if (kb==b.rpoints[i+1]) continue;
                        //if the element A is already present, it is a waste of time to remove it here.
                        if (b.indices[kb]==a.indices[ka]) {a.values[ka]+=(U) b.values[kb]; ++kb;}
                        if (!a.pops.atnodiag || a.indices[ka]!=i+a.RBASE) tr+=abs(a.values[ka]);
                    }
                    while (kb<b.rpoints[i+1]) { 
                        ta=abs((U) b.values[kb]);  if (!a.pops.atnodiag || b.indices[kb]!=i+a.RBASE) tr+=ta;
                        if (ta>thresh) { a_insert[ki]=ka; b_insert[ki]=kb; ++ki; } ++kb; }
                        if (nguess<tr) nguess=tr; thresh=nguess*a.pops.atthresh;
                }
                break;
            default:
                ERROR("Unsupported truncation model.");
            }
        }
    }
    
    a.LASTNORM=thresh;
    if (ki==0) 
    {
        a.trunc(thresh);
        return; // we are soooo lucky
    }
    typename CrsMatrix<U>::index_type j=a.rpoints[a.wr]-1;
    a.rpoints[a.wr]+=ki;
    //now we must insert the marked points. 
    //first we resize a
    
    if (a.indices.size()<a.rpoints[a.wr]||a.indices.size()>a.rpoints[a.wr]+__CRS_STORAGE_FREE) a.presize(a.rpoints[a.wr]);
    ki--;
    
    typename CrsMatrix<U>::index_type k=a.rpoints[a.wr]-1;
    while (ki+1>0)
    {
        if (j+1==a_insert[ki]) {a.values[k]=(U) b.values[b_insert[ki]]; a.indices[k]=b.indices[b_insert[ki]]; --ki; }
        else { a.values[k]=a.values[j]; a.indices[k]=a.indices[j]; --j; }
        --k;
    }
    while(k>j) { a.values[k]=a.values[j]; a.indices[k]=a.indices[j]; --j; --k; }
    
    a.trunc(thresh);
#ifdef BENCHMARK
    TBBenchmarks["crs_incr_mtx"].timer.stop(); TBBenchmarks["crs_incr_mtx"].tot_time+=TBBenchmarks["crs_incr_mtx"].timer;
    TBBenchmarks["crs_incr_mtx"].tot_avg["res_sparsity"]+=(a.rpoints[a.wr]*1./a.wr);  //mean number of el. per row of result
    if  ((a.rpoints[a.wr]*1./a.wr)>TBBenchmarks["crs_incr_mtx"].tot_prop["min_res_sparsity"]) TBBenchmarks["crs_incr_mtx"].tot_prop["min_res_sparsity"]=(a.rpoints[a.wr]*1./a.wr);
#endif
} 

template <class U, class V, class T>
void incr(CrsMatrix<U>& a, const CrsMatrix<V>& b, const T& s) 
{
    //adds in place, this is more expensive if the two matrices have different
    //sparsity patterns, as some merging is necessary
    
#ifdef DEBUG
    if (a.wr!=b.wr || a.wc!=b.wc) ERROR("Incompatible matrix size")
#endif

    std::valarray<typename CrsMatrix<U>::index_type> a_insert((typename CrsMatrix<U>::index_type) 0, b.size());
    std::valarray<typename CrsMatrix<V>::index_type> b_insert((typename CrsMatrix<V>::index_type) 0, b.size());
    typename CrsMatrix<U>::index_type ki=0, ka=0;
    typename CrsMatrix<V>::index_type kb=0;
    double nguess=0., thresh=0., tr, ta;
    if (a.pops.atthresh<=0.)
    {
        //no truncation
        thresh=0.;
        for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
        {
            ka=a.rpoints[i]; kb=b.rpoints[i];
            a.rpoints[i]+=ki;
            for (; ka<a.rpoints[i+1]; ++ka) 
            {   
                while (kb<b.rpoints[i+1] && b.indices[kb]<a.indices[ka])
                { a_insert[ki]=ka; b_insert[ki]=kb; ++kb; ++ki; }
                if (kb==b.rpoints[i+1]) continue;
                if (b.indices[kb]==a.indices[ka]) {a.values[ka]+=((U) (b.values[kb]*s)); ++kb;}
            }
            while (kb<b.rpoints[i+1]) { a_insert[ki]=ka; b_insert[ki]=kb; ++kb; ++ki; }
        }
    }
    else
    {
        if(! a.pops.atrel)
        {
            thresh=a.pops.atthresh;
            //absolute truncation
            for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
            {
                ka=a.rpoints[i]; kb=b.rpoints[i];
                a.rpoints[i]+=ki;
                for (; ka<a.rpoints[i+1]; ++ka) 
                {   
                    while (kb<b.rpoints[i+1] && b.indices[kb]<a.indices[ka])
                    { if (abs((U) (b.values[kb]*s))>thresh) { a_insert[ki]=ka; b_insert[ki]=kb; ++ki; } ++kb;  }
                    if (kb==b.rpoints[i+1]) continue;
                    if (b.indices[kb]==a.indices[ka]) {a.values[ka]+=(U) (b.values[kb]*s); ++kb;}
                }
                while (kb<b.rpoints[i+1]) { if (abs((U) (b.values[kb]*s))>thresh) { a_insert[ki]=ka; b_insert[ki]=kb; ++ki; } ++kb; }
            }
        }
        else
        {
            switch (a.pops.atnorm)
            {
            case at_norminf:
                nguess=thresh=0.; 
                for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
                {
                    ka=a.rpoints[i]; kb=b.rpoints[i];
                    a.rpoints[i]+=ki; tr=0;
                    for (; ka<a.rpoints[i+1]; ++ka) 
                    {
                        while (kb<b.rpoints[i+1] && b.indices[kb]<a.indices[ka])
                        {
                            ta=abs((U) (b.values[kb]*s)); if (!a.pops.atnodiag || b.indices[kb]!=i+a.RBASE) tr+=ta;
                            if (ta>thresh) { a_insert[ki]=ka; b_insert[ki]=kb; ++ki;} ++kb; 
                        }
                        if (kb==b.rpoints[i+1]) continue;
                    //if the element A is already present, it is a waste of time to remove it here, but we must update the trace
                        if (b.indices[kb]==a.indices[ka]) {a.values[ka]+=(U) (b.values[kb]*s); ++kb;}
                        if (!a.pops.atnodiag || a.indices[ka]!=i+a.RBASE) tr+=abs(a.values[ka]);
                    }
                    while (kb<b.rpoints[i+1]) { 
                        ta=abs((U) (b.values[kb]*s)); if (!a.pops.atnodiag || b.indices[kb]!=i+a.RBASE) tr+=ta;
                        if (ta>thresh) { a_insert[ki]=ka; b_insert[ki]=kb; ++ki; } ++kb; 
                    }
                    //also implement shifting of norm
                    if (nguess<tr) nguess=tr; thresh=nguess*a.pops.atthresh;
                }
                break;
            default:
                ERROR("Unsupported truncation model.");
            }
        }
    }
    a.LASTNORM=thresh;
    if (ki==0) 
    {
        a.trunc(thresh);
        return; // we are soooo lucky
    }
    typename CrsMatrix<U>::index_type j=a.rpoints[a.wr]-1;
    a.rpoints[a.wr]+=ki;
    //now we must insert the marked points. 
    //first we resize a
    if (a.indices.size()<a.rpoints[a.wr]||a.indices.size()>a.rpoints[a.wr]+__CRS_STORAGE_FREE) a.presize(a.rpoints[a.wr]);
    ki--;
    //then we insert all the new elements in a single sweep. SINCE WE USE UNSIGNED TYPES WE MUST USE BAILOUT CONDITIONS BASED ON UNDERFLOW
    for (typename CrsMatrix<U>::index_type k=a.rpoints[a.wr]-1; k+1>0; --k)
    {
        if (ki+1>0 && (j+1==0 || j+1==a_insert[ki])) {a.values[k]=((U) (b.values[b_insert[ki]]*s)); a.indices[k]=b.indices[b_insert[ki]]; --ki;}
        else { a.values[k]=a.values[j]; a.indices[k]=a.indices[j]; --j; }
    }
    a.trunc(thresh);
} 


template <class U, class V>
void add(const CrsMatrix<U>& b, const CrsMatrix<V>& c, CrsMatrix<U>& a) 
{
        
#ifdef DEBUG
    if (c.wr!=b.wr || c.wc!=b.wc) ERROR("Incompatible matrix size")
#endif
    
    a.resize(b.wr, c.wc);
    a.RBASE=b.RBASE;
    
    typename CrsMatrix<U>::index_type max_sz;
    max_sz=b.rpoints[b.wr]+c.rpoints[b.wr];
    if (max_sz>b.wr*b.wc) max_sz=b.wr*b.wc;
    a.values.resize(max_sz, (typename CrsMatrix<U>::data_type) 0);
    a.indices.resize(max_sz, (typename CrsMatrix<U>::index_type) 0);
    
    typename CrsMatrix<U>::index_type ka=0, kb=0;
    typename CrsMatrix<V>::index_type kc=0;
    typename CrsMatrix<U>::index_type i=0;
    
    double nguess=0., thresh=0., tr, ta;
    
    if (a.pops.atthresh<=0.)
    {
        //no truncation
        thresh=0.;
        while (i<a.wr)
        {
            if (kb==b.rpoints[i+1]) 
            {
                //exausts row of matrix c
                while (kc<c.rpoints[i+1]) 
                { a.indices[ka]=c.indices[kc]; a.values[ka]=(U) c.values[kc]; ++ka; ++kc; }
                a.rpoints[i+1]=ka; ++i;
            }
            else if (kc==c.rpoints[i+1])
            {
                //exausts row of matrix b
                while (kb<b.rpoints[i+1]) 
                { a.indices[ka]=b.indices[kb]; a.values[ka]=b.values[kb]; ++ka; ++kb; }
                a.rpoints[i+1]=ka; ++i;
            }
            else 
            {
                if (b.indices[kb]<c.indices[kc])
                { a.indices[ka]=b.indices[kb]; a.values[ka]=b.values[kb]; ++ka; ++kb; }
                else if (b.indices[kb]>c.indices[kc])
                { a.indices[ka]=c.indices[kc]; a.values[ka]=(U) c.values[kc]; ++ka; ++kc; }
                else
                { a.indices[ka]=b.indices[kb]; a.values[ka]=(U)c.values[kc]+b.values[kb]; ++ka; ++kb; ++kc; }
            }
        }
        a.trunc(0.); //get rid of unnecessary points
    }
    else
    {
        if(! a.pops.atrel)
        {
            thresh=a.pops.atthresh;
            //absolute truncation
            
            while (i<a.wr)
            {
                if (kb==b.rpoints[i+1]) 
                {
                //exausts row of matrix c
                    while (kc<c.rpoints[i+1]) 
                    { if (abs(c.values[kc])>thresh) { a.indices[ka]=c.indices[kc]; a.values[ka]=(U) c.values[kc]; ++ka;} ++kc; }
                    a.rpoints[i+1]=ka; ++i;
                }
                else if (kc==c.rpoints[i+1])
                {
                //exausts row of matrix b
                    while (kb<b.rpoints[i+1]) 
                    { if (abs(b.values[kb])>thresh) { a.indices[ka]=b.indices[kb]; a.values[ka]=b.values[kb]; ++ka; } ++kb; }
                    a.rpoints[i+1]=ka; ++i;
                }
                else 
                {
                    if (b.indices[kb]<c.indices[kc])
                    { if (abs(b.values[kb])>thresh) { a.indices[ka]=b.indices[kb]; a.values[ka]=b.values[kb]; ++ka; } ++kb; }
                    else if (b.indices[kb]>c.indices[kc])
                    { if (abs(c.values[kc])>thresh) { a.indices[ka]=c.indices[kc]; a.values[ka]=(U) c.values[kc]; ++ka; } ++kc; }
                    else
                    { if (abs(c.values[kc]+b.values[kb])>thresh) 
                    { a.indices[ka]=b.indices[kb]; a.values[ka]=(U)c.values[kc]+b.values[kb]; ++ka; } ++kb; ++kc; }
                }
            }
            a.trunc(thresh); //get rid of unnecessary points
        }
        else
        {
            switch (a.pops.atnorm)
            {
            case at_norminf:
                nguess=thresh=0.; tr=0.;
                while (i<a.wr)
                {
                    if (kb==b.rpoints[i+1]) 
                    {
                //exausts row of matrix c
                        while (kc<c.rpoints[i+1]) 
                        { 
                            ta=abs((U) c.values[kc]);
                            if(!a.pops.atnodiag || c.indices[kc]!=i+a.RBASE) tr+=ta;
                            if (ta>thresh) { a.indices[ka]=c.indices[kc]; a.values[ka]=(U) c.values[kc]; ++ka;} 
                            ++kc; 
                        }
                        a.rpoints[i+1]=ka; ++i; 
                        if (tr>nguess) nguess=tr; tr=0; thresh=nguess*a.pops.atthresh;
                    }
                    else if (kc==c.rpoints[i+1])
                    {
                //exausts row of matrix b
                        while (kb<b.rpoints[i+1]) 
                        { 
                            (ta=abs(b.values[kb])); if(!a.pops.atnodiag || b.indices[kb]!=i+a.RBASE) tr+=ta;
                            if (ta>thresh) { a.indices[ka]=b.indices[kb]; a.values[ka]=b.values[kb]; ++ka; } 
                            ++kb; 
                        }
                        a.rpoints[i+1]=ka; ++i;
                        if (tr>nguess) nguess=tr; tr=0.; thresh=nguess*a.pops.atthresh;
                    }
                    else 
                    {
                        if (b.indices[kb]<c.indices[kc])
                        { 
                            ta=abs(b.values[kb]); 
                            if(!a.pops.atnodiag || b.indices[kb]!=i+a.RBASE) tr+=ta;
                            if (ta>thresh) { a.indices[ka]=b.indices[kb]; a.values[ka]=b.values[kb]; ++ka; } ++kb; 
                        }
                        else if (b.indices[kb]>c.indices[kc])
                        { 
                            ta=abs((U) c.values[kc]); 
                            if(!a.pops.atnodiag || c.indices[kc]!=i+a.RBASE) tr+=ta;
                            if (ta>thresh) { a.indices[ka]=c.indices[kc]; a.values[ka]=(U) c.values[kc]; ++ka; } ++kc; }
                        else
                        { 
                            ta=abs((U) c.values[kc]+b.values[kb]);
                            if(!a.pops.atnodiag || b.indices[kb]!=i+a.RBASE) tr+=ta;
                            if (ta>thresh) 
                            { a.indices[ka]=b.indices[kb]; a.values[ka]=(U)c.values[kc]+b.values[kb]; ++ka; } ++kb; ++kc; }
                    }
                }
                if (tr>nguess) nguess=tr; thresh=nguess*a.pops.atthresh;
                a.trunc(thresh); //get rid of unnecessary points
                
               break;
            default:
                ERROR("Unsupported truncation model.");
            }
        }
    }
}

//computes A=B.C
template <class U, class V>
inline U bmerge(U *is1, V* vs1, const U l1, U *is2, V *vs2, const U l2, U *id, V *vd)
{
   U p1=0, p2=0, pd=0, i;
    while (p1<l1 && p2<l2)
    {
        if (is1[p1]==is2[p2])
        { id[pd]=is1[p1]; vd[pd]=vs1[p1]+vs2[p2]; ++p1; ++p2;}
        else if (is2[p2]<is1[p1])
        { id[pd]=is2[p2]; vd[pd]=vs2[p2]; ++p2;}
        else
        { id[pd]=is1[p1]; vd[pd]=vs1[p1]; ++p1;}
        ++pd;
    }
    for (i=p1;i<l1;++i) {id[pd]=is1[i]; vd[pd]=vs1[i]; ++pd;}
    for (i=p2;i<l2;++i) {id[pd]=is2[i]; vd[pd]=vs2[i]; ++pd;}

    return pd;
}

template <class T, class V> void mult(const CrsMatrix<T>& a, const std::valarray<V>& x, std::valarray<T>& y)
{
    //computes x= a^T y
#ifdef DEBUG
    if (a.wc!=x.size()) ERROR("Incompatible matrix and vector sizes");
#endif
    y.resize(a.wr); y=0.0;
    typename CrsMatrix<T>::index_type ir, k=a.rpoints[0];
    for (ir=0; ir<a.wr; ++ir)
    {
        while (k<a.rpoints[ir+1]) y[ir]+=x[a.indices[k]]*a.values[k++];
    }
}

template <class T, class V> void Tmult(const CrsMatrix<T>& a, const std::valarray<V>& x, std::valarray<T>& y)
{
    //computes y= a^T x
#ifdef DEBUG
    if (a.wr!=x.size()) ERROR("Incompatible matrix and vector sizes");
#endif
    y.resize(a.wc); y=0.0;
    typename CrsMatrix<T>::index_type ir, k=a.rpoints[0];
    for (ir=0; ir<a.wr; ++ir)
    {
        while (k<a.rpoints[ir+1]) y[a.indices[k]]+=x[ir]*a.values[k++];
    }
}

template <class U, class V>
void mult(const CrsMatrix<U>& b, const CrsMatrix<V>& c, CrsMatrix<U>& a) 
{
    /*!MDBG
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    std::string efname;
    efname=std::string("mult.")+int2str(myrank);
    std::ofstream efstream;
    efstream.open(efname.c_str(),std::ios_base::out);
    */
    
#ifdef BENCHMARK
    TBBenchmarks["crs_mult"].n_calls++; TBBenchmarks["crs_mult"].timer.start(); 
    TBBenchmarks["crs_mult"].tot_avg["ops_sparsity"]+=(b.rpoints[b.wr]*1./b.wr+c.rpoints[c.wr]*1./b.wr)/2;  //mean number of el. per row of operands
#endif
#ifdef DEBUG
    if (c.wr!=b.wc) ERROR("Incompatible matrix size");
#endif
    //resize a
    a.resize(b.wr,c.wc);
    a.RBASE=b.RBASE;
    
    
    //we must take the max n. of el. per _row_ of b an the max. el. per _col_ of c, which is slightly harder
    typename CrsMatrix<U>::index_type max_er_b=0, ib, jb, kb;
    typename CrsMatrix<V>::index_type max_er_c=0, ic, jc, kc;
    
    for (ib=0; ib<b.wr; ++ib) if (b.rpoints[ib+1]-b.rpoints[ib]>max_er_b) max_er_b=b.rpoints[ib+1]-b.rpoints[ib];
    for (ic=0; ic<c.wr; ++ic) if (c.rpoints[ic+1]-c.rpoints[ic]>max_er_c) max_er_c=c.rpoints[ic+1]-c.rpoints[ic];
    double nguess=0, thresh=0., tr, ta;
    //!efstream<<"MULT STARTED"<<std::endl;
    
    if ((max_er_b*max_er_c)<a.wc) // N*M^2 algorithm 
    {
        //!efstream<<"* NM2 algorithm chosen"<<std::endl;
        std::valarray<U> bc_val(max_er_b*max_er_c), a_val(max_er_b*max_er_c);
        std::valarray<typename CrsMatrix<U>::index_type> 
            bc_ind(max_er_b*max_er_c), a_ind(max_er_b*max_er_c),
            bc_len(max_er_b), a_len(max_er_b);
        typename CrsMatrix<U>::index_type kbc=0;
        
        //loop across all the rows in B, generating contributions to the same row in A
        for(ib=0; ib<b.wr; ++ib) 
        {
            //stride through nonzero elements in B[ib], computing a temporary array
            //which contains all the terms B[ib,j]*C[j,k] which will collapse into A[ib,k]
            kb=0;
            for (jb=b.rpoints[ib]; jb<b.rpoints[ib+1]; ++jb)
            {
                U bv=b.values[jb];
                kc=b.indices[jb];
                kbc=kb*max_er_c;
                bc_len[kb]=c.rpoints[kc+1]-c.rpoints[kc];

                for (jc=c.rpoints[kc]; jc<c.rpoints[kc+1]; ++jc)
                {
                    bc_val[kbc]=bv*c.values[jc];
                    bc_ind[kbc++]=c.indices[jc];
                }
                ++kb;
            }

            //since we are wise guys, the indices are already sorted, and we "just"
            //need to merge them into A[ib]
            //kb now contains the number of rows into bc temporary
            if (kb==0) 
            {
                //row [ib] is empty. we just update the row pointers and go ahead
                a.rpoints[ib+1]=a.rpoints[ib];
                continue; 
            }
            
            //this part is copied from previous work, so we have to duplicate a few variables 
            //to match the names
            typedef typename CrsMatrix<U>::index_type _int;
            _int block=1; _int kshift, k;
	  
            _int *is1, *is2, *id;
            U *vs1, *vs2, *vd;
            _int *i1=&bc_ind[0], *i2=&a_ind[0], *i3; 
            U *v1=&bc_val[0], *v2=&a_val[0], *v3; 
            _int *l1=&bc_len[0], *l2=&a_len[0], *l3;  
            _int *ls1, *ls2, *ld;
	  
            _int lx=0; 
            _int *ix=NULL; U *vx;	
            if (kb>1) while (block<=kb)
            {
                block<<=1;	
                for (k=0; k+block<=kb; k+=block)
                {
                    kshift=k*max_er_c;
                    is1=i1+kshift;
                    is2=is1+(block>>1)*max_er_c;
                    id=i2+kshift;

                    vs1=v1+kshift;
                    vs2=vs1+(block>>1)*max_er_c;
                    vd=v2+kshift;

                    ls1=l1+k; 
                    ls2=ls1+(block>>1); 
                    ld=l2+k;
                    *ld=bmerge(is1,vs1,*ls1,is2,vs2,*ls2,id,vd);
                }
                if (k<kb)
                {
                    kshift=k*max_er_c;
                    is1=i1+kshift;
                    if ( ix!=NULL && is1!=ix) //if a new block is lagging behind, we merge it with the previous one
                    {								
                        is2=ix;
                        id=i2+kshift;
                        vs1=v1+kshift;
                        vs2=vx;
                        vd=v2+kshift;
                        ls1=l1+k;
                        ld=l2+k;
		      
                        *ld=lx=bmerge(is1,vs1,*ls1,is2,vs2,lx,id,vd);
                        //set the new position of the extra block
                        ix=is1; vx=vs1;
                    }
                    else
                    {
                        ix=is1;
                        id=i2+kshift;
                        vx=vs1=v1+kshift;
                        vd=v2+kshift;
                        ls1=l1+k;
                        ld=l2+k;
                        *ld=lx=*ls1;									    
                        _int *wi=is1+lx;
                        while (is1<wi) { *id=*is1;  ++is1; ++id; }
                        U *wd=vs1+lx;
                        while (vs1<wd) {*vd=*vs1;  ++vs1; ++vd; }
                    }
                }	
	      
                if (ix!=NULL) {	ix=i2+(ix-i1); vx=v2+(vx-v1);} //extra blocks pointers must be "converted" to new work area
                i3=i1; v3=v1; l3=l1;
                i1=i2; v1=v2; l1=l2;
                i2=i3; v2=v3; l2=l3;
            }
            kc=(*l1);	 
            
            //we remove negligeable elements
            jc=0;
            if (a.pops.atthresh<=0.)
            { for (ic=0; ic<kc; ++ic) if (v1[ic]!=(U) 0.) {v1[jc]=v1[ic]; i1[jc++]=i1[ic]; } } // no trunc  
            else
            { 
                if (!a.pops.atrel)
                { for (ic=0; ic<kc; ++ic) if (abs(v1[ic])>a.pops.atthresh) {v1[jc]=v1[ic]; i1[jc++]=i1[ic];} } //absolute 
                else switch (a.pops.atnorm)
                {
                case at_norminf:
                    tr=0;
                    for (ic=0; ic<kc; ++ic) { ta=abs(v1[ic]); if(!a.pops.atnodiag || i1[ic]!=ib+a.RBASE) tr+=ta; if (ta>thresh) {v1[jc]=v1[ic]; i1[jc++]=i1[ic];} }
                    if(tr>nguess) nguess=tr; thresh=nguess*a.pops.atthresh;
                    break;
                default:
                    ERROR("Unsupported truncation mode");
                }
            }
            
            //now put the row into A: we must update the "total elements" value as well, otherwise the resize would be screwed
            a.rpoints[a.wr]=a.rpoints[ib+1]=a.rpoints[ib]+jc;
            
            if (a.rpoints[ib+1]>a.indices.size())
            {
                //we must resize the storage for A!!
                //we do so extimating the size of the rows so far, with an extra buffer
                //we convert to double to avoid overflows
                a.presize((unsigned long) ceil((a.rpoints[ib+1]*1.*a.wr)/(ib+1.))+__CRS_STORAGE_CHUNK);
            }
            
            
            kbc=a.rpoints[ib];

            for (ic=0; ic<jc; ++ic)
            { a.values[kbc]=v1[ic]; a.indices[kbc++]=i1[ic]; }
                            
        } //ends loop on rows
    } //ENDS N*M^2 ALGO
    else //M*N^2 ALGO (which is also waaay easier!)
    {
        //!efstream<<"* MN2 algorithm chosen"<<std::endl;
        std::valarray<U> a_val((U)0., a.wc);
        std::valarray<typename CrsMatrix<U>::index_type> a_ind((typename CrsMatrix<U>::index_type) 0 ,a.wc);
        
        for (ib=0; ib<b.wr;++ib)
        {
            //!efstream<<"+ row "<<ib<<" of "<<b.wr<<std::endl;
            for (unsigned long ja=0; ja<a_val.size(); ++ja) a_val[ja]=(U) 0.;
            /*!performs a sanity check, we must understand whether there are REALLY uninitialized values....
            for (jb=b.rpoints[ib]; jb< b.rpoints[ib+1]; ++jb)
            {
                if (jb>=b.indices.size() ) std::cerr<<" INDEX OUT OF BOUNDS (b1)\n";
                if (jb>=b.values.size() ) std::cerr<<" INDEX OUT OF BOUNDS (b2)\n";
                if (b.indices[jb]>=c.rpoints.size() ) std::cerr<<" INDEX OUT OF BOUNDS (b3)\n";
                if (abs(b.values[jb])>=1e100) std::cerr<<" WEIRD VALUE (b1)\n";
                kc=b.indices[jb];
                for (jc=c.rpoints[kc]; jc<c.rpoints[kc+1]; ++jc)
                {
                    if (jc>=c.indices.size() ) std::cerr<<" INDEX OUT OF BOUNDS (c1)\n";
                    if (jc>=c.values.size() ) std::cerr<<" INDEX OUT OF BOUNDS (c2)\n";
                    if (c.indices[jc]>=a_val.size() ) std::cerr<<" INDEX OUT OF BOUNDS (c3)\n";
                    if (abs(c.values[jc])>=1e100) std::cerr<<" WEIRD VALUE (c1)\n";
                    if (abs(a_val[c.indices[jc]])>=1e100) std::cerr<<" WEIRD VALUE (c2)\n";
                }
            }*/
            
            for (jb=b.rpoints[ib]; jb< b.rpoints[ib+1]; ++jb)
            {
                U bv=b.values[jb];
                kc=b.indices[jb];
                for (jc=c.rpoints[kc]; jc<c.rpoints[kc+1]; ++jc)
                {a_val[c.indices[jc]]+=(U)c.values[jc]*bv; }
            }
            //shrinks the list dropping empty or negligeable elements
            kc=0;
            if (a.pops.atthresh<=0.)
            { for (jc=0; jc<a.wc; ++jc) if (a_val[jc] != (U) 0.) {a_ind[kc]=jc; a_val[kc]=a_val[jc]; kc++; } } // no trunc  
            else
            { 
                if (!a.pops.atrel)
                { for (jc=0; jc<a.wc; ++jc) if (a_val[jc]!=(U) 0. && abs(a_val[jc])>a.pops.atthresh) { a_ind[kc]=jc; a_val[kc++]=a_val[jc]; } }
                else switch (a.pops.atnorm)
                {
                    case at_norminf:
                        tr=0;
                        for (jc=0; jc<a.wc; ++jc) 
                        { 
                            if (a_val[jc]!=(U) 0.) 
                            { 
                                (ta=abs(a_val[jc])); if (!a.pops.atnodiag || jc!=ib+a.RBASE) tr+=ta;
                                if (ta>thresh) { a_ind[kc]=jc; a_val[kc++]=a_val[jc]; } 
                            }
                        }
                        if(tr>nguess) nguess=tr; thresh=nguess*a.pops.atthresh;
                        break;
                    default:
                        ERROR("Unsupported truncation mode");
                }
            }
            //now put the row into A: we must update the "total elements" value as well, otherwise the resize would be screwed
            a.rpoints[a.wr]=a.rpoints[ib+1]=a.rpoints[ib]+kc;
            
            if (a.rpoints[ib+1]>a.indices.size())
            {
                unsigned long newsize=(a.rpoints[ib+1]/(ib+1)+1)*a.wr+__CRS_STORAGE_CHUNK;
                if (newsize< a.rpoints[ib+1]) newsize=a.rpoints[ib+1]+__CRS_STORAGE_CHUNK;
                
                //we must resize the storage for A!!
                //we do so extimating the size of the rows so far, with an extra buffer
                a.presize((a.rpoints[ib+1]/(ib+1)+1)*a.wr+__CRS_STORAGE_CHUNK);
            }
            kc=0;
            for (jc=a.rpoints[ib]; jc<a.rpoints[ib+1]; ++jc)
            { a.values[jc]=a_val[kc]; a.indices[jc]=a_ind[kc++]; }
        }
    
    } //ends N^2 algo
    //!efstream<<"NOW DOING FINAL TRUNCATION"<<std::endl;
    if (a.pops.atthresh>0.) {a.LASTNORM=thresh; a.trunc(a.pops.atrel?(thresh):a.pops.atthresh);}
    //!efstream<<"SUCCESSFULLY FINISHED"<<std::endl;
#ifdef BENCHMARK
    TBBenchmarks["crs_mult"].timer.stop(); TBBenchmarks["crs_mult"].tot_time+=TBBenchmarks["crs_mult"].timer;
    TBBenchmarks["crs_mult"].tot_avg["res_sparsity"]+=(a.rpoints[a.wr]*1./a.wr);  //mean number of el. per row of result
    if  ((a.rpoints[a.wr]*1./a.wr)>TBBenchmarks["crs_mult"].tot_prop["min_res_sparsity"]) TBBenchmarks["crs_mult"].tot_prop["min_res_sparsity"]=(a.rpoints[a.wr]*1./a.wr);
#endif
}

template <class U, class V> 
void transpose(const CrsMatrix<U>& a, CrsMatrix<V>& b)
{
    b.resize(a.wc,a.wr);
    b.presize(a.rpoints[a.wr]);
    
    //first, we build the index
    for (typename CrsMatrix<U>::index_type k=0; k<a.rpoints[a.wr]; ++k)
    {  b.rpoints[a.indices[k]+1]++;  }
    
    for (typename CrsMatrix<V>::index_type i=1; i<=b.wr; ++i)
        b.rpoints[i]+=b.rpoints[i-1];

    //then we populate it, using the index as a guideline
    typename CrsMatrix<U>::index_type ka=0;
    for (typename CrsMatrix<U>::index_type ia=0; ia<a.wr; ++ia)
    {  
        
        while (ka<a.rpoints[ia+1])
        {
            typename CrsMatrix<V>::index_type kb=b.rpoints[a.indices[ka]];
            b.indices[kb]=ia;
            b.values[kb]=a.values[ka];
            ++b.rpoints[a.indices[ka]];
            ++ka; 
        }
    }
    
    //we must recorrect the index
    for (typename CrsMatrix<V>::index_type i=b.wr; i>=1; --i)
        b.rpoints[i]=b.rpoints[i-1];
    b.rpoints[0]=0;
}

template <class U, class V> 
void scale(CrsMatrix<U>& a, const V& b)
{
    for (typename CrsMatrix<U>::index_type k=0; k<a.rpoints[a.wr]; ++k)
        a.values[k]*=b;
}

template <class U>
void neg(CrsMatrix<U>& a)
{
    for (typename CrsMatrix<U>::index_type k=0; k<a.rpoints[a.wr]; ++k)
        a.values[k]=-a.values[k];
}

template <class U>
void map(CrsMatrix<U>& a, U (*f) (const U)) 
{ 
    for (typename CrsMatrix<U>::index_type k=0; k<a.rpoints[a.wr]; ++k)
        a.values[k]=f(a.values[k]);
}

template <class U>
void map(CrsMatrix<U>& a, U (*f) (const U&)) 
{ 
    for (typename CrsMatrix<U>::index_type k=0; k<a.rpoints[a.wr]; ++k)
        a.values[k]=f(a.values[k]);
}

template <class T, class U> 
void copymap(const CrsMatrix<T>& b, CrsMatrix<U>&a, U (*f)(const T))
{
    if ((void*)&b== (void *)&a) 
    {
#ifdef DEBUG
        ERROR("Self-copying is not allowed")
#endif
        return;
    }
    //we don't copy matrix extra properties, this is a value copy, not an assignment operator
    a.resize(b.wr,b.wc);
    a.rpoints=b.rpoints;
    a.values.resize(a.rpoints[a.wr]); a.indices.resize(a.rpoints[a.wr]);
    for(typename CrsMatrix<U>::index_type i=0; i<a.rpoints[a.wr]; ++i) { a.values[i]=f(b.values[i]); a.indices[i]=b.indices[i]; }
}

template <class T, class U> 
void copymap(const CrsMatrix<T>& b, CrsMatrix<U>&a, U (*f)(const T&))
{
    if ((void*)&b== (void *)&a) 
    {
#ifdef DEBUG
        ERROR("Self-copying is not allowed")
#endif
                return;
    }
    //we don't copy matrix extra properties, this is a value copy, not an assignment operator
    a.resize(b.wr,b.wc);
    a.rpoints=b.rpoints;
    a.values.resize(a.rpoints[a.wr]); a.indices.resize(a.rpoints[a.wr]);
    for(typename CrsMatrix<U>::index_type i=0; i<a.rpoints[a.wr]; ++i) { a.values[i]=f(b.values[i]); a.indices[i]=b.indices[i]; }
}

template <class U, class V> 
void incr(CrsMatrix<U>& a, const V& b)
{
#ifdef DEBUG
    if (a.wr!=a.wc) ERROR("Adding a diagonal part to a non-square matrix is nonsense.")
#endif
    typename CrsMatrix<U>::index_type ka, kb, kc, j, ki=0;
    std::valarray<typename CrsMatrix<U>::index_type> iins(a.wr), kins(a.wr);
    U ub=(U) b;
    //choose between bisection and sequential search depending on how many average elements per row are present
    /*sequential search. might be better for VERY sparse matrices, but in general.... it isn't
        for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
        {
            ka=a.rpoints[i]; a.rpoints[i]+=ki;
            for (; ka<a.rpoints[i+1] && a.indices[ka]<i; ++ka);
            if (ka==a.rpoints[i+1]) { iins[ki]=i;  kins[ki++]=ka; continue;}
            if (ka<a.rpoints[i+1] && i==a.indices[ka]) a.values[ka]+=b;
            else { iins[ki]=i;  kins[ki++]=ka;}
        }
    */
    //do bisection search
    for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
    {
        ka=a.rpoints[i]; kc=a.rpoints[i+1];
        a.rpoints[i]+=ki;
        if (kc==ka)
        { iins[ki]=i;  kins[ki++]=ka; }
        else
        {
            //bisection search for the diagonal element
            kb=(kc+ka)/2; j=a.indices[kb];
            while (ka<kb) 
            {
                if (j==i) break; 
                if (j>i) kc=kb;
                else ka=kb;
                kb=(kc+ka)/2;
                j=a.indices[kb];
            }
            
            if (i==j) a.values[kb]+=b; else { iins[ki]=i; kins[ki++]=(i<j?kb:kb+1);}
        }
    }
    if (ki>0) //otherwise we're very lucky
    {    //create space if needed
        if (a.rpoints[a.wr]+ki > a.indices.size()) a.presize(a.rpoints[a.wr]+ki);
        
        kb=a.rpoints[a.wr]-1; 
        a.rpoints[a.wr]+=ki;
        --ki;

        ka=a.rpoints[a.wr]-1;
        while (ki+1>0)
        {
            if (kb+1==kins[ki]) {a.values[ka]=ub; a.indices[ka]=iins[ki]; --ki; }
            else { a.values[ka]=a.values[kb]; a.indices[ka]=a.indices[kb]; --kb; }
            --ka;
        }
        while(ka>kb) { a.values[ka]=a.values[kb]; a.indices[ka]=a.indices[kb]; --ka; --kb; }
    }
    //a.do_at(); //we don't want to do this, it screws up several things...
}

template <class U, class V> 
void incr(CrsMatrix<U>& a, const std::valarray<V>& b)
{
#ifdef DEBUG
    if (a.wr!=a.wc) ERROR("Adding a diagonal part to a non-square matrix is nonsense.")
    if (a.wr!=b.size()) ERROR("Size mismatch")
#endif
    
    //we look for the diagonal elements by bisection
    typename CrsMatrix<U>::index_type ka, kb, kc, j, ki=0;
    std::valarray<typename CrsMatrix<U>::index_type> iins(a.wr), kins(a.wr);
    for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
    {
        ka=a.rpoints[i]; kc=a.rpoints[i+1];
        a.rpoints[i]+=ki;
        if (kc==ka)
        { iins[ki]=i;  kins[ki++]=ka; }
        else
        {
            //bisection search for the diagonal element
            kb=(kc+ka)/2; j=a.indices[kb];
            while (ka<kb) 
            {
                if (j==i) break; 
                if (j>i) kc=kb;
                else ka=kb;
                kb=(kc+ka)/2;
                j=a.indices[kb];
            }
            
            if (i==j) a.values[kb]+=b[i]; else { iins[ki]=i; kins[ki++]=(i<j?kb:kb+1);}
        }
    }
    if (ki>0) //otherwise we're very lucky
    {    //create space if needed
        if (a.rpoints[a.wr]+ki > a.indices.size()) a.presize(a.rpoints[a.wr]+ki);
        
        kb=a.rpoints[a.wr]-1; 
        a.rpoints[a.wr]+=ki;
        --ki;

        ka=a.rpoints[a.wr]-1;
        while (ki+1>0)
        {
            if (kb+1==kins[ki]) {a.values[ka]=b[iins[ki]]; a.indices[ka]=iins[ki]; --ki; }
            else { a.values[ka]=a.values[kb]; a.indices[ka]=a.indices[kb]; --kb; }
            --ka;
        }
        while(ka>kb) { a.values[ka]=a.values[kb]; a.indices[ka]=a.indices[kb]; --ka; --kb; }
    //    hrt.stop(); std::cerr<<"Phase2 "<<hrt*1e-6<<" s\n"; hrt.start();
    }
   // a.do_at();
}

template <class U> 
U trace(const CrsMatrix<U>& a)
{
#ifdef BENCHMARK
    TBBenchmarks["crs_trace"].n_calls++; TBBenchmarks["crs_trace"].timer.start(); 
#endif
    U tr=(U) 0.;
    //we look for the diagonal elements by bisection
    typename CrsMatrix<U>::index_type ka, kb, kc, j;
    for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
    {
        
        ka=a.rpoints[i]; kc=a.rpoints[i+1];
        if (kc==ka) continue;
        
        kb=(kc+ka)/2; j=a.indices[kb];
        while (ka<kb)  
        {
            if (j==i) break;
            if (j>i) kc=kb;
            else ka=kb; 
            kb=(kc+ka)/2; 
            j=a.indices[kb];
        }
        
        if (j==i+a.RBASE) tr+=a.values[kb];
    }
#ifdef BENCHMARK
    TBBenchmarks["crs_trace"].timer.stop(); TBBenchmarks["crs_trace"].tot_time+=TBBenchmarks["crs_trace"].timer;
#endif
    return tr;
}


template <class U> 
double norm1(const CrsMatrix<U>& a, const bool nodiag) 
{
    std::valarray<double> maxc((double) 0., a.wc);
    if (nodiag)
        for (typename CrsMatrix<U>::index_type k=0; k<a.rpoints[a.wr]; ++k)
            if (a.indices[k]!=k+a.RBASE) maxc[a.indices[k]]+=abs(a.values[k]);
    else
        for (typename CrsMatrix<U>::index_type k=0; k<a.rpoints[a.wr]; ++k)
            maxc[a.indices[k]]+=abs(a.values[k]);
    
    return maxc.max();
}

template <class U> 
double norminf(const CrsMatrix<U>& a, const bool nodiag)
{
    double ti, maxc=0;
    typename CrsMatrix<U>::index_type k=0;
    if (nodiag)
        for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
        {
            ti=0;
            for (; k<a.rpoints[i+1]; ++k) if(a.indices[k]!=i+a.RBASE) ti+=abs(a.values[k]);
            if (ti>maxc) maxc=ti;
        }
    else
        for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
        {
            ti=0;
            for (; k<a.rpoints[i+1]; ++k) ti+=abs(a.values[k]);
            if (ti>maxc) maxc=ti;
        }
    //std::cerr<<"NORM IS "<<maxc<<" NODIAG: "<<nodiag<<" RBASE "<<a.RBASE<<"\n";
    return maxc;
}
    
template <class U> 
double normfrob(const CrsMatrix<U>& a, const bool nodiag)
{
    double ak=0,fn=0;
    typename CrsMatrix<U>::index_type k=0;
    if (nodiag)
        for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
        {
            for (; k<a.rpoints[i+1]; ++k) if (a.indices[k]!=i+a.RBASE) {ak=abs(a.values[k]); fn+=ak*ak; }
        }
    else 
    for (; k<a.rpoints[a.wr]; ++k)
        { ak=abs(a.values[k]); fn+=ak*ak; }

    return sqrt(fn);
}

/*//upper bound for the spectral radius according to
template <class U> 
double srad_ub(const CrsMatrix<U>& a)
{
#ifdef DEBUG
    if (a.wr!=a.wc) ERROR("This spectral radius upper bound is not defined for rectangular matrices.");
#endif
    std::valarray<double> C((double) 0., a.wc);
    std::valarray<double> R((double) 0., a.wr);
    double fn2, av2;
    
    typename CrsMatrix<U>::index_type i=0;
    for (typename CrsMatrix<U>::index_type k=0; k<a.rpoints[a.wr]; ++k)
    {
        if (k==a.rpoints[i+1]) ++i;
        av2=abs(a.values[k]); av2*=av2;
        if (i!=a.indices[k]) 
        {   R[i]+=av2;  C[a.indices[k]]+=av2; }
        fn2+=av2;
    }
    for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i) R[i]=abs(sqrt(R[i])-sqrt(C[i]));

    std::cerr<<"RMAX " <<R.max()<<"\n";
    double rho2=R.max(); rho2*=-rho2; rho2+=fn2;
    std::cerr <<sqrt(rho2) <<"extimate 1";
    
    double tra=abs(trace(a));
    
    return sqrt((1.-1./a.wr)*(rho2-tra*tra/a.wr))+tra/a.wr;
}
*/
} //ends namespace toolbox
#endif //ends ifdef __MATRIX_CRS_H
