/* Coordinate matrix storage.
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/

#ifndef __MATRIX_COORD_H
#define __MATRIX_COORD_H
#include "tbdefs.hpp"
#include "matrices.hpp"
#include <vector>
#include <algorithm>

/****************************************************************
 *  This is just a stub to implement conversion between matrix  *
 *  formats: no important matrix-matrix operation is implemented*
 ****************************************************************/
namespace toolbox{
    
template <typename U> 
class CMPoint
{
public:
    unsigned long i, j;
    U val;
    bool operator<(const CMPoint<U>& other) const
    {
        return (i<other.i || (i==other.i && j<other.j));
    }
};

template <typename U> std::ostream& operator<<(std::ostream& ostr, const CMPoint<U>& cm)
{ ostr<<cm.i<<" "<<cm.j<<" "<<cm.val<<"\n"; return ostr;} //we insert line breaks, so that output of vector parser is cleaner
template <typename U> std::istream& operator>>(std::istream& istr, CMPoint<U>& cm)
{ istr>>cm.i>>cm.j>>cm.val; return istr; }
template <typename U>
inline bool compare(const CMPoint<U>&a, const CMPoint<U>& b)
{  return (a.i<b.i?true:a.i==b.i && a.j<b.j); }

template <typename U>
class CoordMatrix
{
    template <class T> friend class CoordMatrix;
    template <class T> friend class CrsMatrix;
    template <class T> friend class FMatrix;
    
    friend class IField<CoordMatrix<U> >;
    friend class IField<CoordMatrix<const U> >;
    friend std::ostream& operator<< <> (std::ostream& ostr, const CoordMatrix& cm);

public:
    typedef unsigned long index_type;
    typedef U data_type;
    
private:
    index_type wr, wc, sz; 
    MatrixOptions<CoordMatrix<U> > pops;
    std::vector<CMPoint<U> > data;   //we use a vector, because we are lazy!
     
     
public:
    //matrix properties
    index_type rows() const { return wr; }
    index_type cols() const { return wc; }    
    index_type size() const { return sz; }
    void resize(const index_type& r, const index_type& c)
    {   wr=r; wc=c; data.resize(0); sz=0; }
    
    //non-destructive resize (automatic, as we are using a vector!);
    void presize(const index_type& nel)
    { data.resize(nel); sz=nel; }
    
    void sort()
    { std::sort(data.begin(), data.end()); }
    
    CoordMatrix(const index_type& rows=0, const index_type& cols=0) : wr(rows), wc(cols) {}
    CoordMatrix(const index_type& rows, const index_type& cols, const index_type& nel, 
                index_type* li, index_type* lj, data_type* ld)
    {
        resize(rows,cols);
        data.resize(nel); sz=nel;
        for (index_type k=0; k<nel; ++k)
        { data[k].i=li[k]; data[k].j=lj[k]; data[k].val=ld[k]; }
    }
    
    //conversion operators. these are just declarations, definitions are in the file matrix-conv.hpp
    //which must be included separately
    template <class T> CoordMatrix(const FMatrix<T>&);
    template <class T> CoordMatrix(const CrsMatrix<T>&);
};



}; //ends namespace toolbox
#endif //ends __MATRIX_COORD_H

