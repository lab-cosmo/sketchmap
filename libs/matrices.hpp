/* Function definitions for the matrix classes.
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/

#ifndef __MATRICES_H
#define __MATRICES_H
#include "tbdefs.hpp"
#include "ioparser.hpp"

namespace toolbox {
    
//pre-definitions of classes
template<class U> class FMatrix;
template<class U> class CrsMatrix;
template<class U> class CoordMatrix;
template<class U> class PCrsMatrix;

template<class MType> class MatrixOptions
{
    //the default class is just empty: we fill it up with specialization
};
    
//pre-definitions of matrix operations which needs to be friends: CrsMatrix class
template <class U, class V> void scale(CrsMatrix<U>& a, const V& b);
template <class T> void neg(CrsMatrix<T>& a);
template <class T> void map(CrsMatrix<T>& a, T (*f)(const T));
template <class T> void map(CrsMatrix<T>& a, T (*f)(const T&));
template <class T, class U> void copymap(const CrsMatrix<T>& b, CrsMatrix<U>&a, U (*f)(const T));
template <class T, class U> void copymap(const CrsMatrix<T>& b, CrsMatrix<U>&a, U (*f)(const T&));
template <class U, class V> void transpose(const CrsMatrix<U>& a, CrsMatrix<V>& b);
template <class U, class V> void incr(CrsMatrix<U>& a, const CrsMatrix<V>& b);
template <class U, class V> void incr(CrsMatrix<U>& a, const V& b);
template <class U, class V> void incr(CrsMatrix<U>& a, const std::valarray<V>& b);
template <class U, class V, class T> void incr(CrsMatrix<U>& a, const CrsMatrix<V>& b, const T& s);
template <class U, class V> void add(const CrsMatrix<U>& b, const CrsMatrix<V>& c, CrsMatrix<U>& a);
template <class U, class V> void mult(const CrsMatrix<U>& b, const CrsMatrix<V>& c, CrsMatrix<U>& a);
template <class U> double norm1(const CrsMatrix<U>& a, const bool nodiag=false);
template <class U> double norminf(const CrsMatrix<U>& a, const bool nodiag=false);
template <class U> double normfrob(const CrsMatrix<U>& a, const bool nodiag=false);
template <class U> U trace(const CrsMatrix<U>& a);
//template <class U> double srad_ub(const CrsMatrix<U>& a);

//pre-definitions of matrix operations which needs to be friends: FMatrix class
template <class U, class V> void scale(FMatrix<U>& a, const V& b);
template <class T> void neg(FMatrix<T>& a);
template <class T> void map(FMatrix<T>& a, T (*f)(const T));
template <class T> void map(FMatrix<T>& a, T (*f)(const T&));
template <class T, class U> void copymap(const FMatrix<T>& b, FMatrix<U>&a, U (*f)(const T));
template <class T, class U> void copymap(const FMatrix<T>& b, FMatrix<U>&a, U (*f)(const T&));
template <class U, class V> void transpose(const FMatrix<U>& a, FMatrix<V>& b);
template <class U, class V> void incr(FMatrix<U>& a, const FMatrix<V>& b);
template <class U, class V> void incr(FMatrix<U>& a, const V& b);
template <class U, class V> void incr(FMatrix<U>& a, const std::valarray<V>& b);
template <class U, class V, class T> void incr(FMatrix<U>& a, const FMatrix<V>& b, const T& s);
template <class U, class V> void add(const FMatrix<U>& b, const FMatrix<V>& c, FMatrix<U>& a);
template <class U, class V> void mult(const FMatrix<U>& b, const FMatrix<V>& c, FMatrix<U>& a);
template <class U> double norm1(const FMatrix<U>& a);
template <class U> double norminf(const FMatrix<U>& a);
template <class U> double normfrob(const FMatrix<U>& a);
template <class U>  U trace(const FMatrix<U>& a);

//pre-definitions of matrix operations which needs to be friends: PCrsMatrix class
template <class U, class V> void scale(PCrsMatrix<U>& a, const V& b);
template <class T> void neg(PCrsMatrix<T>& a);
template <class T> void map(PCrsMatrix<T>& a, T (*f)(const T));
template <class T> void map(PCrsMatrix<T>& a, T (*f)(const T&));
template <class T, class U> void copymap(const PCrsMatrix<T>& b, PCrsMatrix<U>&a, U (*f)(const T));
template <class T, class U> void copymap(const PCrsMatrix<T>& b, PCrsMatrix<U>&a, U (*f)(const T&));
template <class U, class V> void transpose(const PCrsMatrix<U>& a, PCrsMatrix<V>& b);
template <class U, class V> void incr(PCrsMatrix<U>& a, const PCrsMatrix<V>& b);
template <class U, class V> void incr(PCrsMatrix<U>& a, const V& b);
template <class U, class V> void incr(PCrsMatrix<U>& a, const std::valarray<V>& b);
template <class U, class V, class T> void incr(PCrsMatrix<U>& a, const PCrsMatrix<V>& b, const T& s);
template <class U, class V> void add(const PCrsMatrix<U>& b, const PCrsMatrix<V>& c, PCrsMatrix<U>& a);
template <class U, class V> void mult(const PCrsMatrix<U>& b, const PCrsMatrix<V>& c, PCrsMatrix<U>& a);
template <class U> double norm1(const PCrsMatrix<U>& a, const bool nodiag=false);
template <class U> double norminf(const PCrsMatrix<U>& a, const bool nodiag=false);
template <class U> double normfrob(const PCrsMatrix<U>& a, const bool nodiag=false);
template <class U> U trace(const PCrsMatrix<U>& a);
template <class T, class V> bool samenodes(const PCrsMatrix<T>& a, const PCrsMatrix<V>& b);

//io-related function declarations
template<typename U> class CMPoint;
template<typename U> class IField<CMPoint<U> >;
template<typename U> class IField<CoordMatrix<U> >;
template<typename U> std::ostream& operator<< (std::ostream& ostr, const CoordMatrix<U>& cm);

template<typename U> class IField<CrsMatrix<U> >;
template<typename U> class IField<const CrsMatrix<U> >;
template<typename U> std::ostream& operator<< (std::ostream& ostr, const CrsMatrix<U>& cm);

template<typename U> class IField<PCrsMatrix<U> >;
template<typename U> class IField<const PCrsMatrix<U> >;
template<typename U> std::ostream& operator<< (std::ostream& ostr, const PCrsMatrix<U>& cm);

template<typename U> class IField<FMatrix<U> >;
template<typename U> class IField<const FMatrix<U> >;
template<typename U> std::ostream& operator<< (std::ostream& ostr, const FMatrix<U>& cm);

//auto truncation mode
enum ATNormMode { at_norminf=0, at_normone=1, at_normfrob=2 };
} //ends namespace toolbox

#endif //ends ifdef__MATRICES_H
