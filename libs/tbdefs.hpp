/* Main include file for libtoolbox.
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License
*/

#ifndef __TBDEFS_H
#define __TBDEFS_H

/**************************************
 *            INCLUDES                *
 **************************************/
#ifdef TB_MPI
#include "mpi.h"
#endif

#include <valarray>
#include <stdlib.h>
#include <iostream>
#include <cstdarg>
#include <complex>
#include <time.h>

/**************************************
 *          MACROS & DEFINES          *
 **************************************/
#define __TB_STD_EPS 1e-5
#define __TB_STD_MAXSTEP 100

#if __STDC_VERSION__ < 199901L
#   if __GNUC__ >= 2
#       define __func__ __FUNCTION__
#   else
#       define __func__ "<unknown>"
#   endif
#endif

#ifdef __DEBUG
#   define ERROR( msg ) { std::cerr << "Error in " << __func__ << ":\n" << msg <<"\n"; }
#else
#   define ERROR( msg ) { std::cerr << "Error in " << __func__ << ":\n" << msg <<"\n"; exit(-1); }
#endif
#define WARNING( msg ) { std::cerr << "Warning in " << __func__ << ":\n" << msg <<"\n"; }

/**************************************
 *        GENERAL FUNCTIONS           *
 **************************************/
inline double conj(const double& g) {return g;}
inline std::complex<double> conj(const std::complex<double>& g) {return std::conj(g);}
#ifndef __PGI
//inline double abs(const double& g) {return fabs(g);}
#endif
namespace toolbox{
#if defined(__i386__)
    extern "C" {
        inline unsigned long long rdtsc() {
            unsigned long lo, hi;
            /* We cannot use "=A", since this would use %rax on x86_64 */
            __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
            return (unsigned long long)hi << 32 | lo;
        }
    }
#endif
}

namespace toolbox {
    //sort template. keep here until we find a better place for this
    template<class U> void heapsort(std::valarray<U>& a)
    {
        unsigned long n=a.size(); long s, r, e, c;
        U swp;
        if (n<=1) return; //does not sort empty or single-element array
    //heapify the array
        s=(n-2)/2; e=n-1;
        while (s>=0)
        {
            r=s;
            while ((c=r*2+1)<=e)
            {
                if (c<e && a[c]<a[c+1]) ++c;
                if (a[r]<a[c]) { swp=a[r]; a[r]=a[c]; a[c]=swp; r=c; } else break;
            }
            --s;
        }

        while (e>0)
        {
            { swp=a[e]; a[e]=a[0]; a[0]=swp; }
            r=0; --e;
            while ((c=r*2+1)<=e)
            {
                if (c<e && a[c]<a[c+1]) ++c;
                if (a[r]<a[c]) { swp=a[r]; a[r]=a[c]; a[c]=swp; r=c; } else break;
            }
        }
    }
}

namespace toolbox{
    double str2float(const std::string& str);
    int str2int(const std::string& str);
    std::string int2str(const long& ival);
    std::string float2str(const double& ival);
    void csv2floats(const std::string&  istr, std::valarray<double>& vl);
};

/**************************************
 *      TYPEDEFS & ENUMS              *
 **************************************/
namespace toolbox {
template <class T> class StaticError {};

class HRTimer {
private:
    typedef unsigned long long hrttype;
    hrttype base, curr;
#if defined(__i386__)
    inline hrttype gettime() { return rdtsc(); }
#else
#ifdef TB_MPI
    inline hrttype gettime()
    {
        double mptime=MPI_Wtime();
        return (hrttype) (mptime*1e6);
    }
#else
    inline hrttype gettime()
    {
        return (hrttype) (clock())*((hrttype) CLOCKS_PER_SEC);
    }
#endif
#endif
    bool fpaused;
public:
    HRTimer() { curr=base=(hrttype) 0; fpaused=false;}
    void start() { if (fpaused) fpaused=false; else base=gettime(); }
    void stop() { curr=gettime()-base; base=0; fpaused=false;}
    void hold() { curr=(base==(hrttype) 0?base:gettime()-base); fpaused=true; }

    inline operator double()
    {  return (double) curr; }
};

template <class U, unsigned long N=1u>
class fixarray
{
private:
    U elems[N];

public:
    fixarray() {
        for (unsigned long i=0; i<N; ++i)
        {
            elems[i]=U();
        }
    }

    template<class V> fixarray(const fixarray<V,N>& iarr) { for (unsigned long i=0; i<N; ++i) elems[i]=(U) iarr[i]; }

    template<class V> fixarray(const V iarr[N]) { for (unsigned long i=0; i<N; ++i) elems[i]=(U) iarr[i]; }

    template<class V> fixarray(const std::valarray<V>& iarr) {
#ifdef DEBUG
        if (iarr.size()!=N) ERROR("Size mismatch in copy constructor")
#endif
        for (unsigned long i=0; i<N; ++i) elems[i]=(U) iarr[i]; }

    //variadic function. beware, NO ERROR CHECK ON N. of args CAN BE DONE
    fixarray(const U& only) { for (unsigned long i=0; i<N; ++i) elems[i]=only; }
    fixarray(const U first, const U second, ...) {
        elems[0]=first;
        if (N<=1) return;
        elems[1]=second;
        if (N<=2) return;
        va_list ap; va_start(ap, second);
        for (unsigned long i=2; i<N; ++i) elems[i]=va_arg(ap,U);
        va_end(ap);
    }

    template<class V> fixarray<U,N>& operator= (const fixarray<V,N>& iarr)
    {
        if ((void*) &iarr== (void *) this) return *this;
        for (unsigned long i=0; i<N; ++i) elems[i]=iarr[i];
        return *this;
    }

    template<class V> fixarray<U,N>& operator= (const V& iel)
    { for (unsigned long i=0; i<N; ++i) elems[i]=iel;  return *this; }

    inline U& operator[] (unsigned long i)
    {
#ifdef DEBUG
        if (i>=N) ERROR("Out of bounds\n")
#endif
        return elems[i];
    }

    inline const U& operator[] (unsigned long i) const
    {
#ifdef DEBUG
        if (i>=N) ERROR("Out of bounds\n")
#endif
        return elems[i];
    }

    operator double() {
        return elems[0];
    }

    inline unsigned long size() { return N; }
};

template <class U>
class fixarray<U ,0u>
{
    //do something which raise an error!
    typename StaticError<U>::Cannot_define_fixarrays_of_size_zero a;
};

//default is NECESSARY, ABSOLUTE, VALUE. flags serve to change these defaults
enum IChkFlags{ ichk_default=0, ichk_sufficient=1, ichk_relative=2, ichk_change=4};

template<class U, unsigned long N=1u> class IterOptions;
template<class U, unsigned long N>
class IterOptions
{
private:
    fixarray<U,N> val;
    fixarray<U,N> oval;
    fixarray<unsigned long,N> nstep;
    bool fmaxstep_reached;

public:
    unsigned long maxstep;
    fixarray<double,N> thresh;
    fixarray<U,N> target;
    fixarray<unsigned long,N> flags;

    void setval(const U& v, unsigned long i=0)
    {
        oval[i]=val[i];
        val[i]=v;
        nstep[i]++;
    }

    unsigned long iter() { unsigned long mx=0; for (unsigned long i=0; i<N; ++i) if (nstep[i]>mx) mx=nstep[i]; return mx; }

    bool maxstep_reached() {return fmaxstep_reached;}

    operator bool() {
        bool rt=true;

        for (unsigned long i=0; i<N; ++i)
        {
            double d;
            if (nstep[i]>=maxstep) { fmaxstep_reached=true; if (flags[i] & ichk_sufficient) return true;  else continue; }
            if (nstep[i]==0) { rt=false; continue; }
            if (nstep[i]==1 && (flags[i] & ichk_change)) { rt=false; continue; }
            if (flags[i] & ichk_change) d=abs(oval[i]-val[i]);
            else d=abs(val[i]-target[i]);

            if (flags[i] & ichk_relative)
            {
                if ((flags[i] &ichk_change) && abs(val[i])>0. ) d*=1/abs(val[i]);
                else if (abs(target[i])>0. ) d*=1./abs(target[i]);
            }
            if (d>thresh[i] && !(flags[i] & ichk_sufficient)) rt=false;
            if (d<thresh[i] && (flags[i] & ichk_sufficient)) { return true;}
        }
        return rt;
    }

    IterOptions(unsigned long nmaxstep=1,  const fixarray<double,N>& nthresh=fixarray<double,N>(0.),
                   const fixarray<U,N>& ntarget=fixarray<U,N>(),
                   const fixarray<IChkFlags,N>& nflags=fixarray<IChkFlags,N>(ichk_default));
                   
    IterOptions(unsigned long nmaxstep, const double& nthresh, const U& ntarget=0., const IChkFlags& nflags=ichk_default)
    {
        thresh=nthresh;
        val=0.;
        oval=0.;
        target=ntarget;
        flags=nflags;
        maxstep=nmaxstep;
        nstep=0;
        fmaxstep_reached=false; }
};

template <class U, unsigned long N>
IterOptions<U,N>::IterOptions(unsigned long nmaxstep,  const fixarray<double,N>& nthresh,
                   const fixarray<U,N>& ntarget,
                   const fixarray<IChkFlags,N>& nflags)
    : val(), oval(), nstep(), maxstep(nmaxstep), thresh(nthresh), target(ntarget), flags(nflags)
{ nstep=0; val=0.; oval=0.;  fmaxstep_reached=false; }

}; //ends namespace toolbox


/**************************************
 *       TIMING & BENCHMARK           *
 **************************************/
#ifdef BENCHMARK
#include <map>
#include <string>
#include <iomanip>
namespace toolbox {

typedef struct _BenchData {
    HRTimer timer;
    unsigned long n_calls;
    double tot_time;
    std::map<std::string, double>  tot_avg;
    std::map<std::string, double>  tot_prop;
    _BenchData() : n_calls(0), tot_time(0.) {}
} BenchData;

class CTBBD : public std::map<std::string,BenchData> {
public:
    void print(std::ostream& ostr=std::cout) {
        CTBBD::iterator it; std::map<std::string, double>::iterator it2;
        ostr<<"********************************************************************************\n";
        ostr<<"  BENCHMARK RESULTS FOR CALLED FUNCTIONS:                                       \n";
        ostr<<"  FUNCTION                 N. CALLS            AV. TIME                         \n";
        for (it=this->begin(); it!=this->end(); ++it)
        {
            ostr<<"* "<<std::setw(19)<<std::left<<(*it).first
                    <<" "<<std::setw(19)<<std::right<<(*it).second.n_calls
                    <<" "<<std::setw(19)<<(*it).second.tot_time/(*it).second.n_calls<<"\n";
            if ((*it).second.tot_avg.begin()!=(*it).second.tot_avg.end())
            {
                for (it2=(*it).second.tot_avg.begin(); it2!=(*it).second.tot_avg.end(); ++it2)
                {
                    ostr<<"  > "<<std::setw(19)<<std::left<<(*it2).first
                            <<" "<<std::setw(19)<<std::left<<(*it2).second/(*it).second.n_calls<<"\n";
                }
            }
            if ((*it).second.tot_prop.begin()!=(*it).second.tot_prop.end())
            {
                for (it2=(*it).second.tot_prop.begin(); it2!=(*it).second.tot_prop.end(); ++it2)
                {
                    ostr<<"  > "<<std::setw(19)<<std::left<<(*it2).first
                            <<" "<<std::setw(19)<<std::left<<(*it2).second<<"\n";
                }
            }
        }
        ostr<<"********************************************************************************\n";
    }

    CTBBD() : std::map<std::string,BenchData>() {}
};
};
#endif

/************************
  NULL OUTPUT STREAM
*************************/
namespace toolbox{
struct nullstream:
        std::ostream {
    struct nullbuf: std::streambuf {
        int overflow(int c) { return traits_type::not_eof(c); }
    } m_sbuf;
    nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf) {}
        };
};

#ifndef __EXTERNALS
namespace std{
    extern toolbox::nullstream cnull;
};

namespace toolbox{
    namespace constant{
        const double pi=3.1415926535897932385;
        const double sqrt2=1.41421356237309504880;
    };
#ifdef BENCHMARK
    extern CTBBD TBBenchmarks;
#endif
};

#else
namespace std{
    toolbox::nullstream cnull;
};
#endif


/************************
  MPI SAFE STREAMS
*************************/
namespace toolbox{
//MPI-SAFE OUTPUT (execute on every node but output only on node 0)
class mpiostream {
private:
    std::ostream& os;
#ifdef TB_MPI
    MPI_Comm mycomm;
#endif

public:
#ifdef TB_MPI
    mpiostream(std::ostream& ros, const MPI_Comm& rcomm=MPI_COMM_WORLD) : os(ros), mycomm(rcomm) {}
#else
    mpiostream(std::ostream& ros) : os(ros) {}
#endif

    template<class T> std::ostream& operator << (T data)
    {
#ifdef TB_MPI
        static int myrank=-1;
        if (myrank==-1) MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        if (myrank==0) return os<<data; else return std::cnull;
#else
        return os<<data;
#endif
    }

    operator std::ostream&()
    {
#ifdef TB_MPI
        static int myrank=-1;
        if (myrank==-1) MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        if (myrank==0) return os; else return std::cnull;
#else
        return os;
#endif
    }
};
};

#ifndef __EXTERNALS
namespace std{
    extern toolbox::mpiostream pout, perr;
};
#endif




#endif //matches #define __TBDEFS_H
