#ifdef TB_MPI //we don't include if MPI is not active
#ifndef __MATRIX_PCRS_H
#define __MATRIX_PCRS_H
//!
#include <fstream>
#include "tbdefs.hpp"
#include "matrix-crs.hpp"

/*******************************************************************
 Simple parallel CRS matrix class. Every node store a few lines of 
 the whole matrix, in CRS format. We try to make the parallelism 
 as transparent as possible to the end user of the class.
********************************************************************/
namespace toolbox {
//this is ok as long as we don't specify particular options for the parallel version
template <class U>
class MatrixOptions<PCrsMatrix<U> > : public MatrixOptions<CrsMatrix<U> > 
{
public:
    MatrixOptions(double natthresh=0., bool natrel=false, bool natnodiag=false, ATNormMode natnorm=at_norminf) : 
        MatrixOptions<CrsMatrix<U> >(natthresh,natrel,natnodiag,natnorm) {}
    
    MatrixOptions(const MatrixOptions<CrsMatrix<U> >& cmm) : MatrixOptions<CrsMatrix<U> >(cmm) {}
    
    template <class T> MatrixOptions<PCrsMatrix<U> >(MatrixOptions<PCrsMatrix<T> > const& nops) :
    MatrixOptions<CrsMatrix<U> >(nops.atthresh,nops.atrel,nops.atnodiag,nops.atnorm) {}
};

template <class U>
void mpi_sum_function(void* a, void* b, int* len, MPI_Datatype* dt)
{  for (int i=0; i<*len; ++i) ((U*)b)[i]+=((U*)a)[i];  }

template <class U>
class PCrsMatrix {
    template <class T> friend class CoordMatrix;
    template <class T> friend class CrsMatrix;
    template <class T> friend class PCrsMatrix;
    template <class T> friend class FMatrix;
    template <class T, class V> friend void scale(PCrsMatrix<T>& a, const V& b);
    template <class T> friend void neg(PCrsMatrix<T>& a);
    template <class T> friend void map(PCrsMatrix<T>& a, T (*f)(const T));
    template <class T> friend void map(PCrsMatrix<T>& a, T (*f)(const T&));
    template <class T, class V> friend void copymap(const PCrsMatrix<T>& b, PCrsMatrix<V>&a, V (*f)(const T));
    template <class T, class V> friend void copymap(const PCrsMatrix<T>& b, PCrsMatrix<V>&a, V (*f)(const T&));
    template <class T, class V> friend void transpose(const PCrsMatrix<T>& a, PCrsMatrix<V>& b);
    template <class T, class V> friend void incr(PCrsMatrix<T>& a, const PCrsMatrix<V>& b);
    template <class T, class V> friend void incr(PCrsMatrix<T>& a, const V& b);
    template <class T, class V> friend void incr(PCrsMatrix<T>& a, const std::valarray<V>& b);    
    template <class Z, class V, class T> friend void incr(PCrsMatrix<Z>& a, const PCrsMatrix<V>& b, const T& s);
    template <class T, class V> friend void add(const PCrsMatrix<T>& b, const PCrsMatrix<V>& c, PCrsMatrix<T>& a);
    template <class T, class V> friend void mult(const PCrsMatrix<T>& b, const PCrsMatrix<V>& c, PCrsMatrix<T>& a);
    friend double norm1<>(const PCrsMatrix<U>& a, const bool nodiag);
    friend double norminf<>(const PCrsMatrix<U>& a, const bool nodiag);
    friend double normfrob<>(const PCrsMatrix<U>& a, const bool nodiag);
    friend U trace<>(const PCrsMatrix<U>& a);
    template <class T, class V> friend bool samenodes(const PCrsMatrix<T>& a, const PCrsMatrix<V>& b);
    
    friend class toolbox::IField<PCrsMatrix>;
    friend class toolbox::IField<const PCrsMatrix>;
    friend std::ostream& operator<< <> (std::ostream& ostr, const PCrsMatrix& cm);
    
public:
    typedef unsigned long index_type;
    typedef U data_type;
    
private:
    
    MPI_Datatype mpi_data, mpi_index;
    MPI_Op mpi_sum;
    
    //data are stored in a local CRS matrix. we just duplicate by linking upstream
    MPI_Comm mycomm;
    int myrank, mysize;

    std::valarray<index_type> nroots;
    CrsMatrix<U> pmat;
    index_type wc, wr;
    MatrixOptions<PCrsMatrix<U> > pops;
    U dummy;

    void mpi_init(const MPI_Comm ncomm)
    {
        mycomm=ncomm;
        MPI_Comm_rank(mycomm, &myrank);
        MPI_Comm_size(mycomm, &mysize);

       //define an MPI type for elements of type U and index 
        MPI_Type_contiguous(sizeof(data_type),MPI_CHAR, &mpi_data);
        MPI_Type_contiguous(sizeof(index_type),MPI_CHAR,&mpi_index);
        MPI_Type_commit(&mpi_data); MPI_Type_commit(&mpi_index);
        MPI_Op_create(&mpi_sum_function<U>,true,&mpi_sum);
    }
    
public:
    
    void dump(const std::string& fbase) const
    {
        std::string fname=fbase+std::string(".")+int2str(myrank);
        //dump one processor at a time to avoid flooding the filesystem
        for (int in=0; in<mysize; ++in)
        {
            if (in==myrank)
            {
                std::ofstream fmout(fname.c_str());        
                fmout << wr <<" "<< wc<<std::endl;
                fmout << pmat;
                fmout.close();
            }
            MPI_Barrier(mycomm);
        }
    }
    void read(const std::string& fbase)
    {
        //reads one processor at a time to avoid flooding the filesystem
        std::string fname=fbase+std::string(".")+int2str(myrank);
        for (int in=0; in<mysize; ++in)
        {
            if (in==myrank)
            {
                std::ifstream fmin(fname.c_str());
                fmin >> wr >> wc;
                resize(wr,wc);
                unsigned long mbase=pmat.RBASE;
                fmin>> pmat;
                std::cerr<<" READ "<<pmat.RBASE<<" saved " <<mbase<<" node "<<myrank<<std::endl;
                trunc(0.);
                fmin.close();
            }
            MPI_Barrier(mycomm);
        }
    }
    
    PCrsMatrix(const index_type& nwr=0u, const index_type& nwc=0u,
              const MPI_Comm ncomm=MPI_COMM_WORLD)
   {
       
       mpi_init(ncomm);
       resize(nwr, nwc);
       MPI_Barrier(mycomm); //we wait for all the nodes
   }
   
   ~PCrsMatrix()
   {
       MPI_Type_free(&mpi_data);
       MPI_Type_free(&mpi_index);
       MPI_Op_free(&mpi_sum);
       MPI_Barrier(mycomm); //we wait for all the nodes
   }
   
   void resize(const index_type& r, const index_type& c)
   { 
        wr=r; wc=c; 
        //stupid load-balancing: we put the same number of rows on each node; 
        //if there is a remainder, it is put spread among the last nodes.
        nroots.resize(mysize+1);
        index_type kr=0;
        for (int i=0; i<mysize; ++i)
        { nroots[i]=kr; kr+=r/mysize; if (i>=mysize-((long)r)%mysize) kr++; }
        nroots[mysize]=kr;
        pmat.resize(nroots[myrank+1]-nroots[myrank],c);
        pmat.RBASE=nroots[myrank];
        MPI_Barrier(mycomm); //we wait for all the nodes
   }
   
   void optimize()  { pmat.optimize(); }
   //copy operator
   PCrsMatrix<U>& operator = (const PCrsMatrix<U>& a)
   { 
        if (this==&a) return *this;
        mpi_init(a.mycomm);
        resize(a.wr,a.wc);
        pmat=a.pmat;
        setops(a.pops);
        MPI_Barrier(mycomm); //we wait for all the nodes to have copied their part
        return *this;
   }

   template <class T> 
   PCrsMatrix<U>& operator = (const PCrsMatrix<T>& a)
   { 
       if (this==(PCrsMatrix<U>*) &a) return *this;
       mpi_init(a.mycomm);
       resize(a.wr,a.wc);
       pmat=a.pmat;
       setops(a.pops); 
       MPI_Barrier(mycomm); //we wait for all the nodes
       return *this;
   }
   
   PCrsMatrix(const PCrsMatrix& a) : pops(a.pops)
   {
       mpi_init(a.mycomm);
       resize(a.wr, a.wc);
       pmat=a.pmat;
       MPI_Barrier(mycomm); //we wait for all the nodes
   }
   
   template<class T> PCrsMatrix(const CrsMatrix<T>& s, const MPI_Comm ncomm=MPI_COMM_WORLD);
   
   //matrix properties
   index_type rows() const { return wr; }
   index_type cols() const { return wc; }    
   //we want the global size. this means we first have to do some MPI things
   index_type size() const { 
       index_type lsize=pmat.size(), gsize;
       MPI_Allreduce(&lsize, &gsize, 1, MPI_UNSIGNED_LONG, MPI_SUM, mycomm);
       return gsize; }
   
    MatrixOptions<PCrsMatrix<U> > getops() const  { return pops; }
    void getops(MatrixOptions<PCrsMatrix<U> >& rops) const  { rops=pops; }
    void setops(const MatrixOptions<PCrsMatrix<U> >& rops)
    { 
#ifdef DEBUG
        if (rops.atthresh<0.) ERROR("Negative truncation threshold")
#endif
        pops=rops;
        pmat.setops(rops);
        MPI_Barrier(mycomm); //we wait for all the nodes
    }
    
    inline void do_at() { pmat.do_at(); }
    
    /* if we try to access an element on another node, we return a ref. 
    to a dummy element, which is always initialized to zero and which 
    can be written without harm, basically doing nothing to the "good"
    data stored in the local matrix. */
    inline data_type& access(const index_type& i, const index_type& j) 
    { 
#ifdef DEBUG
        if (i>=wr || j>=wc) ERROR("Element index out of bounds");
#endif
        bool fex; index_type kv;
        if (i<nroots[myrank] || i>=nroots[myrank+1])
        {
#ifdef DEBUG
        //!    WARNING("Unable to retrive r/w reference to a nonlocal element");
#endif
            dummy=0; 
            //MPI_Barrier(mycomm); //we wait for all the nodes
            return dummy;
        }
        else
        {
            unsigned long bi=i-nroots[myrank];
            kv=pmat.pexists(bi,j,fex);
            if (!fex) kv=pmat.pinsert(bi,j,(data_type) 0.,kv);
            //MPI_Barrier(mycomm); //we wait for all the nodes
            return pmat.values[kv];
        }
    }
    
    inline void trunc(const double& thresh)
    {
        pmat.trunc(thresh);
        MPI_Barrier(mycomm);
    }
    
    //overloaded operators, for syntactic sugar
    inline U& operator() (const unsigned long& i, const unsigned long& j) 
    {  return access(i,j); }
    
    template <class V> inline PCrsMatrix<U>& operator += (const V& s)
    { incr(*this, s); return *this; }
    
    template <class V> inline PCrsMatrix<U>& operator -= (const V& s)
    { incr(*this, -s); return *this; }
    
    template <class V> inline PCrsMatrix<U>& operator *= (const V& s)
    { scale(*this, s); return *this; }
    
    inline const PCrsMatrix<U> operator -() 
    {  PCrsMatrix<U> r(*this); neg(r); return r; }
    
    template <class V> inline PCrsMatrix<U>& operator += (const PCrsMatrix<V>& a)
    { incr(*this,a); return *this; }
    
    template <class V> inline PCrsMatrix<U>& operator -= (const PCrsMatrix<V>& a)
    { neg(*this); incr(*this,a); neg(*this); return *this; }  //this way we avoid creating a temporary
    
    template <class V> inline PCrsMatrix<U>& operator *= (const PCrsMatrix<V>& a) //we cannot avoid a temporary!
    { PCrsMatrix<U> r; r=*this; mult(r,a,*this); return *this; }
};

template <class V, class T> inline const PCrsMatrix<V> operator + (const PCrsMatrix<V>& a, const T& b)
{ PCrsMatrix<V> r(a); incr(r,b); return r; }
    
template <class V, class T> inline const PCrsMatrix<V> operator + (const PCrsMatrix<V>& a, const PCrsMatrix<T> b)
{ PCrsMatrix<V> r(a); incr(r,b); return r; }

template <class V, class T> inline const PCrsMatrix<V> operator * (const PCrsMatrix<V>& a, const T& b)
{ PCrsMatrix<V> r(a); scale(r,b); return r; }


template <class V, class T> inline const PCrsMatrix<V> operator * (const T& b, const PCrsMatrix<V>& a)
{ PCrsMatrix<V> r(a); scale(r,b); return r; }

template <class V, class T> inline const PCrsMatrix<V> operator * (const PCrsMatrix<V>& a, const PCrsMatrix<T> b)
{ 
    PCrsMatrix<V> r; double att; MatrixOptions<PCrsMatrix<V> > mops; 
    a.getops(mops); r.setops(mops);  mult(a,b,r); return r; 
}

//most of the operations are straightforward, as they are 100% local!!!!
//we really don't have to bother of being in MPI!
template <class U, class V> void scale(PCrsMatrix<U>& a, const V& b)
{ 
    scale(a.pmat, b); 
    MPI_Barrier(a.mycomm); //we wait for all the nodes
}

template <class T> void neg(PCrsMatrix<T>& a)
{ 
    neg(a.pmat); 
    MPI_Barrier(a.mycomm); //we wait for all the nodes
}

template <class T> void map(PCrsMatrix<T>& a, T (*f)(const T))
{ 
    map(a.pmat, f);
    MPI_Barrier(a.mycomm); //we wait for all the nodes
}
template <class T> void map(PCrsMatrix<T>& a, T (*f)(const T&))
{ 
    toolbox::map(a.pmat, f);
    MPI_Barrier(a.mycomm); //we wait for all the nodes
}

template <class T, class U>
bool samenodes(const PCrsMatrix<T>& a,const PCrsMatrix<U>&b)
{
    unsigned long nsz;
    if (a.mycomm!=b.mycomm) return false;
    if (a.nroots.size()!=(nsz=b.nroots.size())) return false;
    else
        for (unsigned long i=0; i<nsz; ++i) 
            if (a.nroots[i]!=b.nroots[i]) return false;
    return true;
}

template <class T, class U> void copymap(const PCrsMatrix<T>& b, PCrsMatrix<U>&a, U (*f)(const T))
{ 
    if (!samenodes(a,b)) ERROR("Unable to copy from matrices with different node maps.");
    map(b.pmat, a.pmat, f); 
    MPI_Barrier(a.mycomm); //we wait for all the nodes
}
template <class T, class U> void copymap(const PCrsMatrix<T>& b, PCrsMatrix<U>&a, U (*f)(const T&))
{ 
    if (!samenodes(a,b)) ERROR("Unable to copy from matrices with different node maps.");
    copymap(b.pmat, a.pmat, f); 
    MPI_Barrier(a.mycomm); //we wait for all the nodes
}

template <class U, class V> void incr(PCrsMatrix<U>& a, const PCrsMatrix<V>& b)
{ 
    if (!samenodes(a,b)) ERROR("Unable to sum matrices with different node maps.");
    incr(a.pmat, b.pmat); 
    
    MPI_Barrier(a.mycomm); 
    //!we need the TRUE norm, which is the max from what's scattered across processors
    double gthr; 
    MPI_Allreduce(&a.pmat.LASTNORM, &gthr, 1, MPI_DOUBLE, MPI_MAX, a.mycomm); 
    if (a.pmat.pops.atthresh>0.) {a.pmat.trunc(a.pmat.pops.atrel?(gthr):a.pmat.pops.atthresh);}
    MPI_Barrier(a.mycomm); //we wait for all the nodes
}

//we must take care not to add to elements not within our local block
//some code replication is therefore necessary.
template <class U, class V> 
void incr(PCrsMatrix<U>& pa, const V& b)
{
    CrsMatrix<U>&a=pa.pmat;
#ifdef DEBUG
    if (pa.wr!=pa.wc) ERROR("Adding a diagonal part to a non-square matrix is nonsense.")
#endif

    typename CrsMatrix<U>::index_type ka, kb, kc, pi, j, ki=0;
    std::valarray<typename CrsMatrix<U>::index_type> iins(a.wr), kins(a.wr);
    U ub=(U) b;
    //do bisection search
    unsigned long mynrows=pa.nroots[pa.myrank+1]-pa.nroots[pa.myrank];
    for (typename CrsMatrix<U>::index_type i=0; i<mynrows; ++i)
    {
        
        pi=pa.nroots[pa.myrank]+i;
        ka=a.rpoints[i]; kc=a.rpoints[i+1];
        a.rpoints[i]+=ki;
        if (kc==ka)
        { iins[ki]=pi;  kins[ki++]=ka; }
        else
        {
            //bisection search for the diagonal element
            kb=(kc+ka)/2; j=a.indices[kb];
            while (ka<kb) 
            {
                if (j==pi) break; 
                if (j>pi) kc=kb;
                else ka=kb;
                kb=(kc+ka)/2;
                j=a.indices[kb];
            }
            
            if (pi==j) a.values[kb]+=b; else { iins[ki]=pi; kins[ki++]=(pi<j?kb:kb+1);}
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
    //a.do_at();
    MPI_Barrier(pa.mycomm);
}

template <class U, class V> void incr(PCrsMatrix<U>& a, const std::valarray<V>& b)
{ 
    ERROR("Reached a stub!");
    incr(a.pmat, b.pmat); 
    //! THIS MUST BE MADE UP WITH CODE REPLICATION, AS WE WANT ONLY TO AFFECT THE RELEVANT PART
}
template <class U, class V, class T> void incr(PCrsMatrix<U>& a, const PCrsMatrix<V>& b, const T& s)
{ 
    if (!samenodes(a,b)) ERROR("Unable to sum matrices with different node maps.");
    incr(a.pmat, b.pmat, s); 
    
    MPI_Barrier(a.mycomm); 
    
    //!we need the TRUE norm, which is the max from what's scattered across processors
    double gthr; 
    MPI_Allreduce(&a.pmat.LASTNORM, &gthr, 1, MPI_DOUBLE, MPI_MAX, a.mycomm); 
    if (a.pmat.pops.atthresh>0.) {a.pmat.trunc(a.pmat.pops.atrel?(gthr):a.pmat.pops.atthresh);}
    
    MPI_Barrier(a.mycomm); //we wait for all the nodes
}
template <class U, class V> void add(const PCrsMatrix<U>& b, const PCrsMatrix<V>& c, PCrsMatrix<U>& a)
{
    a.resize(b.wr,b.wc);
    if (!samenodes(b,c)) ERROR("Unable to sum matrices with different node maps.");
    add(b.pmat, c.pmat, a.pmat); 
    MPI_Barrier(a.mycomm); //we wait for all the nodes
}

template <class U> U trace(const PCrsMatrix<U>& pa)
{
    //compute the local part of the trace
    const CrsMatrix<U>& a=pa.pmat;
    U tr=(U) 0.;
    //we look for the diagonal elements by bisection
    typename CrsMatrix<U>::index_type ka, kb, kc, j;
    typename CrsMatrix<U>::index_type prbase=pa.nroots[pa.myrank];
    for (typename CrsMatrix<U>::index_type i=0; i<a.wr; ++i)
    {
        typename CrsMatrix<U>::index_type pi=i+prbase;
        ka=a.rpoints[i]; kc=a.rpoints[i+1];
        if (kc==ka) continue;
        
        kb=(kc+ka)/2; j=a.indices[kb];
        while (ka<kb)  
        {
            if (j==pi) break;
            if (j>pi) kc=kb;
            else ka=kb; 
            kb=(kc+ka)/2; 
            j=a.indices[kb];
        }
        
        if (j==pi) tr+=a.values[kb];
    }
    
    U gtr;
    MPI_Allreduce(&tr, &gtr, 1, pa.mpi_data, pa.mpi_sum, pa.mycomm);  //this is also blocking
    return gtr;
}

//#define PCRS_MDBG 1
//computes A=B.C
template <class U, class V> void mult(const PCrsMatrix<U>& pb, const PCrsMatrix<V>& pc, PCrsMatrix<U>& pa) 
{
#ifdef BENCHMARK
    TBBenchmarks["pcrs_mult"].tot_avg["ops_sparsity"]+=(pb.size()*1./pb.wr+pc.size()*1./pb.wr)/2.;  //mean number of el. per row of operands
    TBBenchmarks["pcrs_mult"].n_calls++; TBBenchmarks["pcrs_mult"].timer.start(); 
    HRTimer hrt; 
    hrt.start();
#endif

    /* 
      this is a bit more complicated. it seems that the best is to combine the on-node and 
      off-node elements into a temporary matrix¸ which we'll then combine to do the multiply
    */
    //!
    int mpirv;
    static int icall=0;
#ifdef PCFS_MDBG
    //pb.dump("pmult-a");
    //pc.dump("pmult-b");
    MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
    std::string efname;
    efname=std::string("pmult.")+int2str(pc.myrank);
    std::ofstream efstream;
    efstream.open(efname.c_str(),std::ios_base::out | std::ios_base::app);
    efstream<<"PMATRIX mult CALLED "<< icall++ <<std::endl;
#endif

    const CrsMatrix<U>& b=pb.pmat;
    const CrsMatrix<V>& c=pc.pmat;
    //b.sanity();
    //c.sanity();
    CrsMatrix<U>& a=pa.pmat;
    CrsMatrix<U> wc;
    
#ifdef DEBUG
    if (pc.wr!=pb.wc) ERROR("Incompatible matrix size.")
#endif
    //resize a
    pa.resize(pb.wr,pc.wc);
    wc.resize(pc.wr,pc.wc);
    
    if (!samenodes(pa,pb)) ERROR("First matrix and result should have same node maps.");
    if (pa.mycomm!=pc.mycomm) ERROR("Unable to multiply matrices with different communicators.");

    unsigned long myrowbase=pc.nroots[pc.myrank];
    //first, we mark which lines are needed, by checking all the elements in b
    std::valarray<unsigned char> fneeded((unsigned char) 0, pc.wr); 
    for(unsigned long k=0; k<b.rpoints[b.wr]; ++k) fneeded[b.indices[k]]=1;
    
    //now we make up some vectors to query 
    //the required lines from each processor
    //s_*** arrays hold the send data, i.e. data concerning communication ME -> OTHERS
    //r_*** arrays hold the recv data, i.e. data concerning communication ME <- OTHERS
    
    //first we want to tell everybody how many lines we need from them. 
    //we make up an array holding this information, then briadcast it
    std::valarray<int> r_nlines(pc.mysize); //this holds the number of lines needed from node [i]
    std::valarray<int> r_disp(pc.mysize+1); //incremental displacement: will be a 'pointer' to the line info
    r_nlines==0; r_disp[0]=0;
    for (int in=0; in<pc.mysize; ++in)
    {
        for (unsigned long k=pc.nroots[in]; k<pc.nroots[in+1]; ++k) r_nlines[in]+=fneeded[k];
        r_disp[in+1]=r_disp[in]+r_nlines[in];
    }
    
    //then we prepare an array which contains the line numbers we need.
    //r_disp[in] will point to the beginning of the slice belonging to node in
    std::valarray<unsigned long> r_lines;
    r_lines.resize(r_disp[pc.mysize]);

    unsigned long i=0;
    for (int in=0; in<pc.mysize; ++in)
    {
        for (unsigned long k=pc.nroots[in]; k<pc.nroots[in+1]; ++k) 
            if (fneeded[k]!=0u) { r_lines[i]=k; i++; }
    }

    //we prepare the space to store the requests coming from the other nodes
    std::valarray<int> s_nlines(pc.mysize), s_disp(pc.mysize+1);
    //we do an MPI_Alltoall to know how many lines every node wants to have
    mpirv=MPI_Alltoall(&r_nlines[0], 1, MPI_INTEGER, &s_nlines[0], 1, MPI_INTEGER, pc.mycomm);
#ifdef PCFS_MDBG
        if (mpirv!=0) efstream<<"ERROR "<<mpirv<<" on node "<<pc.myrank<<" on MPI_Alltoall"<<std::endl;
#endif 

    //now s_nlines[ni] contains the number of rows we will need to sent to processor ni. 
    s_disp[0]=0; //we compute displacements
    for (int in=0; in<pc.mysize; ++in) s_disp[in+1]=s_disp[in]+s_nlines[in];
    std::valarray<unsigned long> s_lines(s_disp[pc.mysize]);
    
    
    //a little summary: now in r_lines we hold all the lines we need from other nodes,
    //and the section of node in starts at r_lines[r_disp[in]]. we have allocated
    //the space necessary to get the requests of other nodes. let's make the folks know each other!
    mpirv=MPI_Alltoallv(&r_lines[0],&r_nlines[0],&r_disp[0],MPI_UNSIGNED_LONG, 
                        &s_lines[0],&s_nlines[0],&s_disp[0],MPI_UNSIGNED_LONG, 
                        pc.mycomm);
#ifdef PCFS_MDBG
        if (mpirv!=0) efstream<<"ERROR "<<mpirv<<" on node "<<pc.myrank<<" on MPI_Alltoallv (1)\n";
        efstream<<"A) first handshaking performed"<<std::endl;
#endif
    
    /*!! WRITING OUT THE RESULTS
    if (pc.myrank==0) std::cerr<<"DATA RECEIVED\n";
    MPI_Barrier(pc.mycomm);
    for (int in=0; in<pc.mysize; ++in)
    {
        if (pc.myrank==in)
        {
            std::cerr<<"node : "<<pc.myrank<<" requires\n";
            for (unsigned long k=0; k<pc.mysize; ++k) 
            {
                std::cerr<<"from "<<k<<" : "<<s_nlines[k]<<" rows, ";
                for (unsigned long i=s_disp[k]; i<s_disp[k+1]; ++i) 
                    std::cerr<<s_lines[i]<<" ";
                std::cerr<<"\n";
            }
            std::cerr<<"node : "<<pc.myrank<<" will send\n";
            for (unsigned long k=0; k<pc.mysize; ++k) 
            {
                std::cerr<<"to "<<k<<" : "<<r_nlines[k]<<" rows, ";
                for (unsigned long i=r_disp[k]; i<r_disp[k+1]; ++i) 
                    std::cerr<<r_lines[i]<<" ";
                std::cerr<<"\n";
            }
            std::cerr<<"+++++++++++++++++++++++++++++++++++++\n";
        }
        MPI_Barrier(pc.mycomm);
    }
    //*/
    
    //well. now we know what the others want. we must tell them how many elements
    //we will send them, so that they can prepare the storage. at the same time, we 
    //must query and allocate ourselves!
    
#ifdef PCFS_MDBG
//       efstream<< pc.mysize<< " ::  "<< s_disp[pc.mysize]<<std::endl;
                       std::valarray<int> s_nels(s_disp[pc.mysize]);
                       i=0;
                       for (int in=0; in<pc.mysize; ++in)
                       {
//                           efstream<< s_disp[in]<< "->"<< s_disp[in+1]<<"("<<i<<") : ";
//                       efstream.flush();
    //we must convert from global to local row indices when reading in c
                           for (int j=s_disp[in]; j<s_disp[in+1]; ++j)
                               s_nels[i++]=c.rpoints[s_lines[j]-myrowbase+1]-c.rpoints[s_lines[j]-myrowbase];
                       }
                       efstream<<std::endl;
#else
    std::valarray<int> s_nels(s_disp[pc.mysize]);
    i=0;
    for (int in=0; in<pc.mysize; ++in)
    //we must convert from global to local row indices when reading in c
    for (int j=s_disp[in]; j<s_disp[in+1]; ++j)
        s_nels[i++]=c.rpoints[s_lines[j]-myrowbase+1]-c.rpoints[s_lines[j]-myrowbase];
#endif

    //allocate and do the communication
    std::valarray<int>  r_nels(0,r_disp[pc.mysize]);
    mpirv=MPI_Alltoallv(&s_nels[0],&s_nlines[0],&s_disp[0],MPI_INTEGER, 
                   &r_nels[0],&r_nlines[0],&r_disp[0],MPI_INTEGER, 
                pc.mycomm);
    
#ifdef PCFS_MDBG
        if (mpirv!=0) efstream<<"ERROR "<<mpirv<<" on node "<<pc.myrank<<" on MPI_Alltoallv (2)\n";
#endif
    
#ifdef PCFS_MDBG
    //!CHECK BACK!!! for some reason, every now and then, the elements mismatch!
    std::valarray<int>  c_nels(0,s_disp[pc.mysize]);
    mpirv=MPI_Alltoallv(&r_nels[0],&r_nlines[0],&r_disp[0],MPI_INTEGER, 
                         &c_nels[0],&s_nlines[0],&s_disp[0],MPI_INTEGER, 
                         pc.mycomm);
    for(int k=0; k<s_nels.size(); ++k)
        if (s_nels[k]!=c_nels[k]) std::cerr<<"FW-BW MISMATCH @"<<k<<"!!!!"<<std::endl;
#endif

    //ok. now in r_lines[] we have the row indices, and in r_nels[] the corresponding n.of elements.
    //we have already enough info to set up row pointers and allocate the index and values arrays
    wc.rpoints[0]=0;  int ik=0;
    for(i=0; i<wc.wr; ++i)
    {
        wc.rpoints[i+1]=wc.rpoints[i];
        if (ik<r_disp[pc.mysize] && i==r_lines[ik]) wc.rpoints[i+1]+=r_nels[ik++];
    }

#ifdef PCFS_MDBG
    efstream<<"B) second handshaking performed (size "<<wc.rpoints[wc.wr]<<")"<<std::endl;
#endif
    wc.presize(wc.rpoints[wc.wr]);
    
#ifdef BENCHMARK
    hrt.stop();
    TBBenchmarks["pcrs_mult"].tot_avg["time_handshaking"]+=hrt;
    hrt.start(); 
#endif
    
    /*!! WRITING OUT THE RESULTS
    MPI_Barrier(pc.mycomm);
    for (int in=0; in<pc.mysize; ++in)
    {
        if (pc.myrank==in)
        {
            std::cerr<<"node : "<<pc.myrank<<" holds index list\n";
            for (unsigned long k=0; k<=wc.wr; ++k) 
            {
                std::cerr<<wc.rpoints[k]<<" ";
            }
            std::cerr<<"\n";
            std::cerr<<"+++++++++++++++++++++++++++++++++++++\n";
        }
        MPI_Barrier(pc.mycomm);
    }
    //*/
    
    //wunderbar. now we have allocated room for the off-node elements,
    //we know where they are and where they should go. we can send requests
    //and issue receives.
    //we use sends of non-contiguous data, so that 
    //no buffering or copying is needed!
    
#ifndef PCRS_BLOCKING
    //first we collect the indices, then the values.
    //here again, r_ stands for incoming buffers, and s_ for outgoing ones
    std::valarray<MPI_Request> s_req_index(pc.mysize), s_req_data(pc.mysize);
    std::valarray<MPI_Request> r_req_index(pc.mysize), r_req_data(pc.mysize);
    //std::valarray<MPI_Datatype> rb_type(pc.mysize), sb_type(pc.mysize);
    //std::valarray<MPI_Datatype> rv_type(pc.mysize), sv_type(pc.mysize);
    MPI_Datatype r_type, s_type;
    std::valarray<int> r_edisp;  //these contain the pointers to the beginning of row slices
    std::valarray<int> s_edisp;
            
#ifdef PCFS_MDBG
    efstream<<"RECEIVING BUFFER DATA: ";
    efstream.flush();
#endif
    
    unsigned long nsend=0, nrcv=0, ic;
    ic=0;
    for (int in=0; in<pc.mysize; ++in)
    {
        if (r_nlines[in]==0) continue;
        nrcv++;
        r_edisp.resize(r_nlines[in]); 
        int tdisp=0;
        for (long k=r_disp[in]; k<r_disp[in+1]; ++k)
        {
            r_edisp[k-r_disp[in]]=wc.rpoints[r_lines[k]];
            tdisp+=r_nels[k];
        }
#ifdef PCFS_MDBG
        efstream<<in<<":"<<r_nlines[in]<< "  ";//<<std::endl;
#endif
        int rv=MPI_Type_indexed(r_nlines[in],&r_nels[r_disp[in]],&r_edisp[0],pc.mpi_index,&r_type);
        if(rv!=0) 
        { std::cerr<<"ERROR ON NODE: "<<pc.myrank <<" WHILE MAKING INDEX RECV TYPE IN COMM. WITH NODE "<<in<<" RETURN CODE "<<rv<<">>"<<r_nlines[in]<<","<<r_nels[r_disp[in]]<<","<<r_edisp[0] << std::endl;
#ifdef PCFS_MDBG
            efstream<<"n. of chunks: "<<r_nlines[in]<<", buffer size "<<wc.indices.size()<<"\n";
            for (long k=0; k<r_nlines[in]; ++k)
                efstream<<r_nels[r_disp[in]+k]<<" elements starting @ "<<r_edisp[k]<<"\n";
            efstream<<"********************************************************************\n";
            efstream.close();
#endif
            exit(1);
        }
        
        int rc=MPI_Type_commit(&r_type);
#ifdef PCFS_MDBG
        if(rc!=0) 
        {  std::cerr<<"ERROR ON NODE: "<<pc.myrank <<" WHILE COMMITTING INDEX RECV TYPE IN COMM. WITH NODE "<<in<<" RETURN CODE "<<rc<<std::endl; }
#endif
        MPI_Irecv(&wc.indices[0],1,r_type,in,101,pc.mycomm,&r_req_index[ic]);
        MPI_Type_free(&r_type);
        ic++;
    }
    
#ifdef PCFS_MDBG
    efstream<<std::endl;
    efstream<<"SENDING BUFFER DATA: ";
    efstream.flush();
#endif
    ic=0;
    for (int in=0; in<pc.mysize; ++in)
    {
        if (s_nlines[in]==0) continue;
        nsend++;
        int tdisp=0;
        s_edisp.resize(s_nlines[in]); 
        for (long k=s_disp[in]; k<s_disp[in+1]; ++k)
        {
            s_edisp[k-s_disp[in]]=c.rpoints[s_lines[k]-pc.nroots[pc.myrank]];
            tdisp+=s_nels[k];
        }
        
        int rv=MPI_Type_indexed(s_nlines[in],&s_nels[s_disp[in]],&s_edisp[0],pc.mpi_index,&s_type);
        if(rv!=0) 
        { std::cerr<<"ERROR ON NODE: "<<pc.myrank <<" WHILE MAKING INDEX SEND TYPE IN COMM. WITH NODE "<<in<<" RETURN CODE "<<rv<<std::endl;      }

        MPI_Type_commit(&s_type);
        MPI_Isend(&(const_cast<std::valarray<unsigned long> &>(c.indices)[0]),1,s_type,in,101,pc.mycomm,&s_req_index[ic]);
        MPI_Type_free(&s_type);
        ic++;
    }
    
    //go with data
#ifdef PCFS_MDBG
    efstream<<"\nRECEIVING VALUES"<<std::endl;
#endif
    ic=0;
    for (int in=0; in<pc.mysize; ++in)
    {
        if (r_nlines[in]==0) continue;
        r_edisp.resize(r_nlines[in]); 
        for (long k=r_disp[in]; k<r_disp[in+1]; ++k)
            r_edisp[k-r_disp[in]]=wc.rpoints[r_lines[k]];
        int rv=MPI_Type_indexed(r_nlines[in],&r_nels[r_disp[in]],&r_edisp[0],pc.mpi_data,&r_type);
        if(rv!=0) 
        { std::cerr<<"ERROR ON NODE: "<<pc.myrank <<" WHILE MAKING DATA RECV TYPE IN COMM. WITH NODE "<<in<<" RETURN CODE "<<rv<<std::endl;      }
        
        MPI_Type_commit(&r_type);
        MPI_Irecv(&wc.values[0],1,r_type,in,202,pc.mycomm,&r_req_data[ic]);
        MPI_Type_free(&r_type);
        ic++;
    
    }
#ifdef PCFS_MDBG
    efstream<<"ic: "<<ic<<" nrecv " <<nrcv<<"\n";
    efstream<<"SENDING VALUES"<<std::endl;
#endif
    ic=0;
    for (int in=0; in<pc.mysize; ++in)
    {
        if (s_nlines[in]==0) continue;
        s_edisp.resize(s_nlines[in]); 
        for (long k=s_disp[in]; k<s_disp[in+1]; ++k)
            s_edisp[k-s_disp[in]]=c.rpoints[s_lines[k]-pc.nroots[pc.myrank]];
        int rv=MPI_Type_indexed(s_nlines[in],&s_nels[s_disp[in]],&s_edisp[0],pc.mpi_data,&s_type);
        if(rv!=0) 
        { std::cerr<<"ERROR ON NODE: "<<pc.myrank <<" WHILE MAKING DATA SEND TYPE IN COMM. WITH NODE "<<in<<" RETURN CODE "<<rv<<std::endl;      }
        
        MPI_Type_commit(&s_type);
        MPI_Isend(&(const_cast<std::valarray<U>& >(c.values)[0]),1,s_type,in,202,pc.mycomm,&s_req_data[ic]);
        MPI_Type_free(&s_type);
        ic++;
    }
    
#ifdef PCFS_MDBG
    efstream<<"ic: "<<ic<<" nsend " <<nsend<<"\n";
    efstream<<"C) send/receive issued"<<std::endl;
#endif
    
#ifdef BENCHMARK
    hrt.stop();
    TBBenchmarks["pcrs_mult"].tot_avg["time_issuecomm"]+=hrt;
    hrt.start(); 
#endif

    //now we should be having both the on-site and off-node matrices. 
    //the hard part is done, and we can do the product!
    
    //NOW we have to wait for all the communication to be over.
    //!it should be sufficient to wait for successful receive, as this implies that the send is over as well...
    std::valarray<MPI_Status> status(pc.mysize);
    
    mpirv=MPI_Waitall(nsend,&s_req_index[0],&status[0]);
#ifdef PCFS_MDBG
    if (mpirv!=0) efstream<<"ERROR "<<mpirv<<" on node "<<pc.myrank<<" on MPI_Waitall (1)\n";
#endif
    mpirv=MPI_Waitall(nsend,&s_req_data[0],&status[0]);
#ifdef PCFS_MDBG
    if (mpirv!=0) efstream<<"ERROR "<<mpirv<<" on node "<<pc.myrank<<" on MPI_Waitall (2)\n";
#endif
    
    mpirv=MPI_Waitall(nrcv,&r_req_index[0],&status[0]);
#ifdef PCFS_MDBG
    if (mpirv!=0) 
    {
        efstream<<"ERROR "<<mpirv<<" on node "<<pc.myrank<<" on MPI_Waitall (3)\n";
        for (int in=0; in<pc.mysize;in++)
            efstream << " + status on node "<<in<<": "<< status[in].MPI_ERROR<<std::endl;
    }
#endif
    
    mpirv=MPI_Waitall(nrcv,&r_req_data[0],&status[0]);
#ifdef PCFS_MDBG
    if (mpirv!=0) efstream<<"ERROR "<<mpirv<<" on node "<<pc.myrank<<" on MPI_Waitall (4)\n";
    efstream<<"D) waitall is over"<<std::endl;
#endif
    
    
#ifdef BENCHMARK
    hrt.stop();
    TBBenchmarks["pcrs_mult"].tot_avg["time_waitcomm"]+=hrt;
    hrt.start(); 
#endif
    
#else // #ifdef PCRS_BLOCKING
    std::valarray<MPI_Request> req(pc.rows());
    std::valarray<MPI_Status> stat(pc.rows());
    //blocking send implementation, where the nonblocking version is slower or buggy
    unsigned long nreq=0;
    for (int in=0; in<pc.mysize; ++in)
    {
        if (r_nlines[in]==0) continue;
        unsigned long k=r_disp[in];
                
        for (long i=0; i< r_nlines[in]; ++i)
        {
            MPI_Irecv(&(wc.indices[wc.rpoints[r_lines[k]]]),r_nels[k],pc.mpi_index,in,100+i,pc.mycomm,&(req[nreq]));
            nreq++;
            ++k; 
        }
    }
        
#ifdef PCFS_MDBG
    efstream<<"NBRECV INDICES"<<std::endl;
#endif

    MPI_Barrier(pc.mycomm);
    
    for (int in=0; in<pc.mysize; ++in)
    {
        if (s_nlines[in]==0) continue;
        unsigned long k=s_disp[in];
        
        for (long i=0; i<s_nlines[in]; ++i)
        {
            MPI_Send(&(const_cast<std::valarray<unsigned long> &>(c.indices)[c.rpoints[s_lines[k]-pc.nroots[pc.myrank]]]),
                       s_nels[k],pc.mpi_index,in,100+i,pc.mycomm);
            
            k++;
        }
    }
    
    MPI_Waitall(nreq,&req[0],&stat[0]);
    
#ifdef PCFS_MDBG
    efstream<<"RSEND INDICES"<<std::endl;
#endif
    nreq=0;
    for (int in=0; in<pc.mysize; ++in)
    {
        if (r_nlines[in]==0) continue;
        unsigned long k=r_disp[in];
        for (long i=0; i< r_nlines[in]; ++i)
        {
            MPI_Irecv(&wc.values[wc.rpoints[r_lines[k]]],r_nels[k],pc.mpi_data,in,100+r_nlines[in]+i,pc.mycomm,&(req[nreq]));
            nreq++;
            k++;
        }
    }

#ifdef PCFS_MDBG
    efstream<<"NBRECV DATA"<<std::endl;
#endif

    MPI_Barrier(pc.mycomm);
    
    for (int in=0; in<pc.mysize; ++in)
    {
        if (s_nlines[in]==0) continue;
        unsigned long k=s_disp[in];
        for (long i=0; i<s_nlines[in]; ++i)
        {
            MPI_Send(&(const_cast<std::valarray<U> &>(c.values)[c.rpoints[s_lines[k]-pc.nroots[pc.myrank]]]),
                       s_nels[k],pc.mpi_data,in,100+s_nlines[in]+i,pc.mycomm);
            k++;
        }
    }
    
    
#ifdef PCFS_MDBG
    efstream<<"RSEND DATA"<<std::endl;
#endif

    MPI_Waitall(nreq,&req[0],&stat[0]);
    //send indices (blocking)

    
#endif // #ifdef PCRS_BLOCKING
    //std::cerr<<"testing sanity on node " <<pb.myrank<<"\n";
    //wc.sanity();
    //std::cerr<<"sanity check passed on node " <<pb.myrank<<"\n";
//    std::ofstream wcs((std::string("wc.")+int2str(pb.myrank)).c_str());
//    wcs<<wc;
    mult(b,wc,a);
           
#ifdef PCFS_MDBG
    efstream<<"E) MULT DONE"<<std::endl;
    efstream<<"last norm "<<a.LASTNORM<<"\n";
#endif
    
    MPI_Barrier(pa.mycomm);
    
#ifdef BENCHMARK
    hrt.stop();
    TBBenchmarks["pcrs_mult"].tot_avg["time_multiply"]+=hrt;
    hrt.start(); 
#endif
    
    //!we need the TRUE norm, which is the max from what's scattered across processors
    double gthr=0; 
    
    mpirv=MPI_Allreduce(&a.LASTNORM, &gthr, 1, MPI_DOUBLE, MPI_MAX, pa.mycomm); 
#ifdef PCFS_MDBG
    if (mpirv!=0) efstream<<"ERROR "<<mpirv<<" on node "<<pc.myrank<<" on MPI_Allreduce (3)"<<std::endl;
#endif
    if (a.pops.atthresh>0.) {a.trunc(a.pops.atrel?(gthr):a.pops.atthresh);}
    
#ifdef BENCHMARK
    hrt.stop();
    TBBenchmarks["pcrs_mult"].tot_avg["time_trunc"]+=hrt;
    hrt.start(); 
#endif

    MPI_Barrier(pc.mycomm);
#ifdef BENCHMARK
    hrt.stop(); TBBenchmarks["pcrs_mult"].timer.stop(); 
    TBBenchmarks["pcrs_mult"].tot_time+=TBBenchmarks["pcrs_mult"].timer;
    TBBenchmarks["pcrs_mult"].tot_avg["res_sparsity"]+=(pa.size()*1./pa.wr);  //mean number of el. per row of result
    if  ((pa.size()*1./pa.wr)>TBBenchmarks["pcrs_mult"].tot_prop["min_res_sparsity"]) TBBenchmarks["pcrs_mult"].tot_prop["min_res_sparsity"]=(pa.size()*1./pa.wr);
#endif
#ifdef PCFS_MDBG
    efstream<<"F) CLOSING STREAM"<<std::endl;
    efstream.close();
#endif
}

//by now these do nothing but allowing linking :-)
template <class U, class V> void transpose(const PCrsMatrix<U>& pa, PCrsMatrix<V>& pb) 
{
    //THIS IS BROKEN! WE BYPASS IT TO AVOID MESSING UP
    pb=pa;
    return;
    //wow, this is a simple function but requires a lot of communication. first we transpose 
    //the chunks on each processor, then we send them around
    if (pa.mycomm!=pb.mycomm) ERROR("Unable to compute transposition between different communicators.");
    pb.resize(pa.wc,pa.wr);
    
    const CrsMatrix<U>& a=pa.pmat;
    CrsMatrix<V>& b=pb.pmat;
    
    //first, transpose on-site
    CrsMatrix<U> ta; 
    transpose(a,ta);
    //we correct the column indices in ta!
    ta.indices+=pb.nroots[pb.myrank];
    
    //then we tell everybody how they should allocate space for the new elements.
    //we get the elements per row
    std::valarray<int> l_rp(ta.wr);
    for (unsigned long i=0; i<ta.wr; ++i) l_rp[i]=ta.rpoints[i+1]-ta.rpoints[i];

    //ok, now on every processor we have parts of the rows which should go to the others. 
    //every processor should know how many it should be expecting from the others for every
    //of its own rows.
    unsigned long mynrows=(pb.nroots[pb.myrank+1]-pb.nroots[pb.myrank]);
    std::valarray<int> nl_rp(mynrows*pb.mysize);
    
    MPI_Barrier(pb.mycomm);  //!
    for(int in=0; in<pb.mysize; ++in)
    {
        MPI_Gather(&l_rp[pb.nroots[in]],pb.nroots[in+1]-pb.nroots[in],MPI_INTEGER,
                    &nl_rp[0],pb.nroots[in+1]-pb.nroots[in],MPI_INTEGER,in,pb.mycomm);
    }
    
    /*!
    MPI_Barrier(pb.mycomm);
    for (unsigned long in=0; in<pb.mysize; ++in)
    {
        if (in==pb.myrank)
        {
            if (in==0) std::cerr<<"\n DATA: \n"; 
            std::cerr<<"node: "<<pb.myrank<<"\n";
            for(unsigned long in2=0; in2<pb.mysize; ++in2)
            {
                for (unsigned long i=0; i< mynrows; ++i)
                  std::cerr<<nl_rp[in2*mynrows+i]<<" ";
                std::<<"\n";
            }
        }
        MPI_Barrier(pb.mycomm);
    }
    //*/
            
    unsigned long nel=0;
    std::valarray<int> r_pos(mynrows);
    for (unsigned long i=0; i<mynrows; ++i) 
    { 
        r_pos[i]=b.rpoints[i]=nel; 
        for(int in=0; in<pb.mysize; ++in) nel+=nl_rp[in*mynrows+i]; 
    }
    b.rpoints[mynrows]=nel;
    b.presize(nel);
    
    //now we know how the transpose looks like, and we have allocated it as well.
    //we also know where chunks from differents nodes should go. we can therefore
    //make the communication to collect the data
    MPI_Datatype r_type, s_type;
    
    std::valarray<MPI_Request> send_req_index(pb.mysize), send_req_data(pb.mysize);
    std::valarray<MPI_Request> recv_req_index(pb.mysize), recv_req_data(pb.mysize);
    
    //FIRST THE INDICES
    for (int in=0; in<pb.mysize; ++in)
    {
        if (in==pb.myrank)
        {
            for (unsigned long bi=0; bi<mynrows; ++bi) 
                for(long bj=0; bj<nl_rp[in*mynrows+bi]; ++bj)
                    b.indices[r_pos[bi]+bj]=ta.indices[ta.rpoints[pb.nroots[pb.myrank]+bi]+bj];
        } 
        else 
        {
            MPI_Type_indexed(mynrows,&nl_rp[in*mynrows],&r_pos[0],pb.mpi_index,&s_type);
            MPI_Type_commit(&s_type);
            MPI_Irecv(&b.indices[0],1,s_type,in,101,pb.mycomm,&recv_req_index[in]);
            MPI_Type_free(&s_type);
        }
        for (unsigned long bi=0; bi<mynrows; ++bi) 
            r_pos[bi]+=nl_rp[in*mynrows+bi];
    }
    
    std::valarray<int> ta_disp(ta.wr);
    for (unsigned long i=0;  i<ta.wr; ++i) ta_disp[i]=ta.rpoints[i];
    for (int in=0; in<pb.mysize; ++in)
    {
        if (in==pb.myrank)  continue;
        MPI_Type_indexed(pb.nroots[in+1]-pb.nroots[in],&l_rp[pb.nroots[in]],&ta_disp[pb.nroots[in]],pb.mpi_index,&r_type);
        MPI_Type_commit(&r_type);
        MPI_Isend(const_cast<void *>((const void *)&ta.indices[0]),1,r_type,in,101,pb.mycomm,&send_req_index[in]);
        MPI_Type_free(&r_type);
    }
    
    //SECOND THE VALUES
    for (unsigned long bi=0; bi<mynrows; ++bi) r_pos[bi]=b.rpoints[bi];
    
    for (int in=0; in<pb.mysize; ++in)
    {
        if (in==pb.myrank)
        {
            for (unsigned long bi=0; bi<mynrows; ++bi) 
                for(long bj=0; bj<nl_rp[in*mynrows+bi]; ++bj)
                    b.values[r_pos[bi]+bj]=ta.values[ta.rpoints[pb.nroots[pb.myrank]+bi]+bj];
        } 
        else 
        {
            MPI_Type_indexed(mynrows,&nl_rp[in*mynrows],&r_pos[0],pb.mpi_data,&s_type);
            MPI_Type_commit(&s_type);
            MPI_Irecv(&b.values[0],1,s_type,in,101,pb.mycomm,&recv_req_index[in]);
            MPI_Type_free(&s_type);
        }
        for (unsigned long bi=0; bi<mynrows; ++bi) 
            r_pos[bi]+=nl_rp[in*mynrows+bi];
    }
    
    for (int in=0; in<pb.mysize; ++in)
    {
        if (in==pb.myrank)  continue;
        MPI_Type_indexed(pb.nroots[in+1]-pb.nroots[in],&l_rp[pb.nroots[in]],&ta_disp[pb.nroots[in]],pb.mpi_data,&r_type);
        MPI_Type_commit(&r_type);
        MPI_Isend(const_cast<void *>((const void *)&ta.values[0]),1,r_type,in,101,pb.mycomm,&send_req_index[in]);
        MPI_Type_free(&r_type);
    }
    
    //NOW we have to wait for all the communication to be over
    MPI_Status status;
    for (int in=0; in<pb.mysize; ++in)
    {
        if (in==pb.myrank) continue;
        MPI_Wait(&send_req_index[in], &status);
        MPI_Wait(&recv_req_index[in], &status);
    }
    for (int in=0; in<pb.mysize; ++in)
    {
        if (in==pb.myrank) continue;
        MPI_Wait(&send_req_data[in], &status);
        MPI_Wait(&recv_req_data[in], &status);
    }
    
}

template <class U> double norm1(const PCrsMatrix<U>& a, const bool nodiag) 
{
    ERROR("Reached a stub!");
}

template <class U> double norminf(const PCrsMatrix<U>& a, const bool nodiag)
{
    double lval, gval;
    //, itr;
    /*if (nodiag) //we must do it here as the local matrix doesn't know where is the diagonal
    {
        unsigned long bi, k=0;
        for (unsigned long i=0; i<a.pmat.wr; ++i)
        {
            itr=0.; bi=i+a.nroots[a.myrank];
            for(; k<a.pmat.rpoints[i+1]; ++k)
                if (bi!=a.pmat.indices[k]) itr+=abs(a.pmat.values[k]); 
            if (itr>lval) lval=itr;
        }
    }
    else lval=norminf(a.pmat, false);
    */
    lval=norminf(a.pmat,nodiag);
    //MPI_Barrier(a.mycomm); //we need to make sure that everybody has collected its term
    MPI_Allreduce(&lval, &gval, 1, MPI_DOUBLE, MPI_MAX, a.mycomm);  //this is also blocking
    return gval;
}


template <class U> double normfrob(const PCrsMatrix<U>& a, const bool nodiag) 
{
    double lval=0., gval, av;
    
    unsigned long k=0, bi;
    if (nodiag)
    {
        for (unsigned long i=0; i<a.pmat.wr; ++i)
        {
            bi=i+a.nroots[a.myrank];
            for(; k<a.pmat.rpoints[i+1]; ++k)
                if (bi!=a.pmat.indices[k]) { av=abs(a.pmat.values[k]); lval+=av*av; }
        }
    }
    else for(k=0; k<a.pmat.rpoints[a.pmat.wr]; ++k) { av=abs(a.pmat.values[k]); lval+=av*av; }
    
    MPI_Barrier(a.mycomm);
    MPI_Allreduce(&lval, &gval, 1, MPI_DOUBLE, MPI_SUM, a.mycomm);  //this is also blocking
    return sqrt(gval);
}


}; //ends namespace toolbox
#endif //ends #ifndef __MATRIX_PCRS_H
#endif //ends #ifdef MPI
