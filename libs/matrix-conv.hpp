#ifndef __MATRIX_CONV_H
#define __MATRIX_CONV_H
#include "tbdefs.hpp"
#include "matrices.hpp"
#include "matrix-crs.hpp"
#include "matrix-full.hpp"
#include "matrix-coord.hpp"

namespace toolbox {
    
#ifdef __MATRIX_PCRS_H
    /* Gets PCrsMatrix data from the node where serial CrsMatrix has more than zero rows */
template <class U>
template <class T> PCrsMatrix<U>::PCrsMatrix(const CrsMatrix<T>& s, const  MPI_Comm ncomm)
{
    mpi_init(ncomm);
    
    int iamsender=(s.rows()>0?1:0);    //! we rely on the fact that only one node has this cond. fulfilled.

    // now we let everyone know who is the sender
    std::valarray<int> slist(mysize);
    MPI_Allgather(&iamsender, 1, MPI_INTEGER, &slist[0], 1, MPI_INTEGER, mycomm);
    int sender=-1;
    
    for (int i=0; i< mysize; ++i) 
        if (slist[i]==1) if(sender==-1) sender=i; else ERROR("More than one process qualified as sender!");
    if (sender==-1) ERROR("No process qualified as sender!");
    
    // now we get the matrix size and resize it.
    typename PCrsMatrix<U>::index_type dim[2];
    if (iamsender) { dim[0]=s.rows(); dim[1]=s.cols(); }
    MPI_Bcast(dim,2,mpi_index,sender,mycomm);
    resize(dim[0],dim[1]);
    
    // now we copy the options, as if it were an array of char...
    MPI_Bcast(&(const_cast<CrsMatrix<T> &>(s).pops),sizeof(pops),MPI_CHAR,sender,mycomm);
    setops(s.pops);
    
    // now we can send out the row indices to the nodes.
    unsigned long nmyrows=nroots[myrank+1]-nroots[myrank];
    MPI_Request rreq;
    MPI_Irecv(&pmat.rpoints[0],nmyrows+1,mpi_index,sender,101,mycomm,&rreq);
    
    if (iamsender) {
        for (int i=0; i<mysize; ++i)
            MPI_Send(&(const_cast<CrsMatrix<T> &>(s).rpoints[nroots[i]]),nroots[i+1]-nroots[i]+1,mpi_index,i,101,mycomm);
    };
    
    //wait for receive
    MPI_Status rstatus;
    MPI_Wait(&rreq, &rstatus);
    //then shift the indices as necessary, since we are getting chunks of data
    for (typename PCrsMatrix<U>::index_type i=1;i<=nmyrows; ++i)
        pmat.rpoints[i]-=pmat.rpoints[0];
    pmat.rpoints[0]=0;
    pmat.presize(pmat.rpoints[nmyrows]);
    
    //very well. now we can share the column indices and the data!
    MPI_Request rreq_i, rreq_d;
    MPI_Irecv(&pmat.indices[0],pmat.rpoints[nmyrows],mpi_index,sender,202,mycomm,&rreq_i);
    
    if (iamsender) {
        for (int i=0; i<mysize; ++i)
            MPI_Send(&(const_cast<CrsMatrix<T> &>(s).indices[s.rpoints[nroots[i]]]),
                     s.rpoints[nroots[i+1]]-s.rpoints[nroots[i]],mpi_index,i,202,mycomm);
    };
    
    MPI_Irecv(&pmat.values[0],pmat.rpoints[nmyrows],mpi_data,sender,303,mycomm,&rreq_d);
    
    if (iamsender) {
        for (int i=0; i<mysize; ++i)
            MPI_Send(&(const_cast<CrsMatrix<T> &>(s).values[s.rpoints[nroots[i]]]),
                     s.rpoints[nroots[i+1]]-s.rpoints[nroots[i]],mpi_data,i,303,mycomm);
    };
    MPI_Wait(&rreq_i, &rstatus);
    MPI_Wait(&rreq_d, &rstatus);
}

template <class U>
template <class T> CrsMatrix<U>::CrsMatrix(const PCrsMatrix<T>& s) 
{
    RBASE=0;
}

#endif //ends #ifdef __MATRIX_PCRS_H

template <class U>
template <class T> CoordMatrix<U>::CoordMatrix(const CrsMatrix<T>& s) 
{
    resize(s.wr, s.wc);
    data.resize(s.rpoints[s.wr]);
    index_type i=0;
    for (index_type k=0; k<s.rpoints[s.wr]; ++k)
    { 
        while (k==s.rpoints[i+1]) ++i; 
        data[k].i=i;
        data[k].j=s.indices[k];
        data[k].val=(U) s.values[k];
    }
}

template <class U>
template <class T> CoordMatrix<U>::CoordMatrix(const FMatrix<T>& s) 
{
    /*OPTIMIZATION NEEDED*/
    resize(s.wr, s.wc);
    data.resize(0);
    
    unsigned long i=0, j=0;
    CMPoint<U> nv;
    for (index_type ks=0; ks<s.wr*s.wc; ++ks)
    {
        if(j>=s.wc) { ++i; j=0; }
        if (abs(s.data[ks])>0) { nv.val=(U) s.data[ks]; nv.i=i; nv.j=j; data.push_back(nv); }
        ++j;
    }
    sz=data.size();
}

template <class U>
template <class T> CrsMatrix<U>::CrsMatrix(const CoordMatrix<T>& ns) 
{
    CoordMatrix<T> s(ns);
    s.sort();
    RBASE=0;
    resize(s.wr, s.wc);
    values.resize(s.size());
    indices.resize(s.size());
    
    index_type i=0,k; rpoints[0]=0;
    for (k=0; k<s.size(); ++k)
    { 
        while (s.data[k].i>i) { rpoints[i+1]=k; ++i; }
        indices[k]=s.data[k].j;
        values[k]=(U) s.data[k].val; 
    }
    
    while (i<s.wr) { rpoints[i+1]=k; ++i; }
}


template <class U>
template <class T> CrsMatrix<U>::CrsMatrix(const FMatrix<T>& s) 
{
    /*OPTIMIZATION NEEDED*/
    //pass through a conversion to CoordMatrix, we are lazy!
    RBASE=0;
    CoordMatrix<T> cm(s);
    *this=cm;
}

template <class U>
template <class T> FMatrix<U>::FMatrix(const CrsMatrix<T>& s) 
{
    resize(s.wr,s.wc);
        
    index_type i=0, k=s.rpoints[i], h=0; 
        //gets to first non-empty rows
    while (i<wr && s.rpoints[i+1]==s.rpoints[i]) {++i; h+=wc;}
        
    //writes first element (if the matrix is not empty!)
    if (i==wr) return;
    h+=s.indices[k]; data[h]=s.values[k]; ++k; 
        
    while (true)
    {
        if (k==s.rpoints[i+1]) 
        { 
            ++i; if (i==wr) break;
            h+=wc;
                
            continue; 
        }
            
        h+=s.indices[k]-s.indices[k-1];
        data[h]=(U) s.values[k];
        ++k;
    }
}

template <class U>
template <class T> FMatrix<U>::FMatrix(const std::valarray<std::valarray<T> >& s)
{
#ifdef DEBUG
    if (s.size()<1 || s[0].size()<1) ERROR("Trying to initialize FMatrix from empty valarray")
#endif 
    resize(s.size(), s[0].size());
    index_type k=0;
    for (index_type i=0; i<wr; ++i)
    {
#ifdef DEBUG
        if (s[i].size()!=wc) ERROR("Element "<<i<<" in valarray has wrong dimension")
#endif 
        for (index_type j=0; j<wc; ++j) data[k++]=(U)s[i][j];
    }
}

template <class U>
template <class T> FMatrix<U>::FMatrix(const CoordMatrix<T>& s) 
{
    resize(s.wr,s.wc);
    data=0;
    for (index_type k=0; k<s.size(); ++k)
        data[s.data[k].i*s.wc+s.data[k].j]=(U) s.data[k].val;
}

};//ends namespace toolbox
#endif //ends ifndef __MATRIX_CONV_H
