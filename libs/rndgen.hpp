#ifndef __RNDGEN_DEC_H
#define __RNDGEN_DEC_H

#include "tbdefs.hpp"
/****************************************************
 * Declarations for templates for random number gen *
 * More complex generators (gaussian or correlated) *
 * can be connected to several kinds of uniform gen *
 * Each must specitfy the type of the random number *
 * returned by operator(), and the type of the      *
 * state vector.                                    *
 ****************************************************/

namespace toolbox {
/****************************************************
 * Any uniform rnd gen class must be a template depend
 * Abstract random numbers generator class, holding *
 * options of two different data types ass, holding *
 ****************************************************
template <class U, typename STATUS>
class _RndGen {
protected:
    STATUS pstate;
    virtual U extract()=0;
    
public:
    typedef STATUS status_type;
    virtual void setstatus(const STATUS& state)=0;
    
    
    U operator() () { return extract(); }
};
*/

/***************************************************
 * Standard generator for uniformly distributed    *
 * doubles within [0:1] using stdlib drand48.      *
 * Any unif. rnd. gen. must expose at least these  *
 * member functions and properties.                *
 ***************************************************/
class StdRndUniform {
public:
    typedef double rnd_type;
    typedef long status_type;
private:    
    inline rnd_type extract() { return drand48(); }
    status_type pstate;
    
public:
    void setstatus(const status_type& state);
    void getstatus(status_type& state) const { state=pstate; }
    
    //constructors
    StdRndUniform(const status_type& init=0) { pstate=init; }
    StdRndUniform(const StdRndUniform& ru) { setstatus(ru.pstate);}
    
    StdRndUniform& operator=(const StdRndUniform& ru) { setstatus(ru.pstate); return *this; }    
    inline rnd_type operator() () { return extract(); }
};


/****************************************************
 * Mersenne twister random longs, code copied from  *
 * Gnu Scientific Library.                          *
 ****************************************************/
#define MT_N 624 
#define MT_M 397
#define MT_UPPER_MASK 0x80000000UL
#define MT_LOWER_MASK 0x7fffffffUL   
#define MT_MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

struct MTStatus {
    std::valarray<unsigned long> mtvec;
    unsigned short pos;
    
    MTStatus() : mtvec(MT_N), pos(0) {}
    MTStatus& operator =(const MTStatus& mts2) {if (&mts2==this) return *this; mtvec.resize(MT_N); mtvec=mts2.mtvec; pos=mts2.pos; return *this; }
};

class MTRndUniform;
class MTRndUniformL
{
    friend class MTRndUniform;
public:
    typedef unsigned long rnd_type;
    typedef MTStatus status_type;
private:
    rnd_type extract();
    status_type pstate;
    
public:
    void setstatus(const status_type& nstat) { pstate=nstat; } //we don't need to perform any initialization
    void getstatus(status_type& state) const { state=pstate; }
    void seed(const unsigned long& seed=4537);
    
    inline rnd_type operator() () { return extract(); }
    
    // constructors
    MTRndUniformL(const status_type& nstat) { pstate=nstat;}
    MTRndUniformL(const unsigned long& nseed=4357) {  seed(nseed); }
    MTRndUniformL(const MTRndUniformL& cp) { pstate=cp.pstate;}
    MTRndUniformL& operator=(const MTRndUniformL& ru) { setstatus(ru.pstate); return *this; }    
};

/****************************************************
 * Mersenne twister random double, minimal dupl.    * 
 * of the above long version. WITHIN [0,1)          *
 ****************************************************/
class MTRndUniform 
{
public:
    typedef double rnd_type;
    typedef MTStatus status_type;
    
private:
    MTRndUniformL LMT;
    rnd_type extract() { return ((rnd_type) (LMT.extract()))*2.3283064365386962890625e-10; }

public:
    //the status is just a duplicate of the long MT 
    void setstatus(const status_type& nstat) { LMT.setstatus(nstat); }
    void getstatus(status_type& state) const { state=LMT.pstate; }
    void seed(const unsigned long& seed=4537) { LMT.seed(seed); } 
    MTRndUniform(const status_type& nstat) { setstatus(nstat);}
    MTRndUniform(const unsigned long& nseed=4357) { seed(nseed); }
    MTRndUniform(const MTRndUniform& cp) { setstatus(cp.LMT.pstate);}
    MTRndUniform& operator=(const MTRndUniform& ru) { setstatus(ru.LMT.pstate); return *this; } 
    
    inline rnd_type operator() () { return extract(); }
};

/***************************************************
 * Template for gaussian random numbers generation *
 * it requires the returned type, which must def   *
 * some conversion from floating point (unless an  *
 * overload of some of the functions is defined)   *
 * and a uniform number generator returning the    *
 * same type.                                      *
 ***************************************************/
template <class U> 
struct RGPars {
    U mean, sigma;
    RGPars(): mean((U) 0.), sigma((U) 1.) {}
    RGPars(const U& nmean, const U& nsigma) : mean(nmean), sigma(nsigma) {}
};
    
template <class U, class UGType> 
struct RGStatus {
    UGType rustate;
    RGPars<U> pars; 
    U gauss2;
    bool fg2ready;
    
    RGStatus(): rustate(), pars(), gauss2(), fg2ready(false) {};
};

template <class U, class UGEN=StdRndUniform> 
class RndGaussian {
public:
    typedef U rnd_type;
    typedef RGStatus<U, typename UGEN::status_type> status_type; 
private:
    status_type pstate;
    UGEN rugen;
    rnd_type extract();
    
public:    
    //constructors
    RndGaussian(const UGEN& ru=UGEN()): rugen(ru) {pstate=status_type(); rugen.setstatus(pstate.rustate); }
    RndGaussian(const status_type& nstat): rugen() { setstatus(nstat); }    
    RndGaussian(const RGPars<U>& pars): rugen()
        { pstate=status_type(); rugen.getstatus(pstate.rustate); pstate.pars=pars; }
    
    void getstatus(status_type& state) const { state=pstate; }
    void setstatus(const status_type& nstat) { pstate=nstat; rugen.setstatus(nstat.rustate); }
    void setpars(const RGPars<U>& npars) { pstate.pars=npars; }
    
    U mean() {return pstate.pars.mean; }
    U sigma() {return pstate.pars.sigma; }
    UGEN& RUGenerator() { return rugen; } 
    inline rnd_type operator() () { return extract(); }
};


/***************************************************
 * Template for gaussian random numbers generation *
 * with given correlation over a finite timespan   *
 ***************************************************/
template <class U>
struct RCGPars {
    U mean, sigma;
    std::valarray<U> corr;
    RCGPars() : mean(), sigma(), corr() {}
    RCGPars(const U& nmean, const U& nsigma, const std::valarray<U>& ncorr) :
        mean(nmean), sigma(nsigma), corr(ncorr) {}
    RCGPars<U> operator=(const RCGPars<U>& np) 
    { if (this!=&np) { mean=np.mean; sigma=np.sigma; corr.resize(np.corr.size()); corr=np.corr; } return *this; }
};

template <class U, class GGType> 
struct RCGStatus {
    GGType rgstate;
    std::valarray<U> mem;
    std::valarray<U> alpha;
    RCGPars<U> pars;
    RCGStatus(): rgstate(), alpha(), pars() {};
};

#define RCG_ALPHA_ITER 1000
#define RCG_ALPHA_TOL ((U)1.e-8)

template <class U, class GGEN=RndGaussian<U,StdRndUniform> >
class RndCorrGaussian {
public:
    typedef U rnd_type;
    typedef RCGStatus<U, typename GGEN::status_type> status_type;
    
private:
    status_type pstate;
    void check_corr(const std::valarray<U>& corr);
    void get_alpha(const std::valarray<U>& corr, std::valarray<U>& alpha, U tol, unsigned long max_step);
    GGEN rggen;
    rnd_type extract();

public:
    RndCorrGaussian(const GGEN& rg=GGEN()): rggen(rg) {this->pstate=status_type(); rggen.getstatus(this->pstate.rustate); }
    RndCorrGaussian(const status_type& nstat=status_type()): rggen() { setstatus(nstat); setpars(nstat.pars); }    
    RndCorrGaussian(const RCGPars<U>& pars, U tol=RCG_ALPHA_TOL, unsigned long max_step=RCG_ALPHA_ITER): rggen()
        { this->pstate=status_type(); rggen.getstatus(this->pstate.rgstate); setpars(pars,tol, max_step); init();}
    
    
    void setstatus(const status_type& nstat) 
    {   this->pstate=nstat; 
        this->pstate.rgstate.pars.mean=nstat.pars.mean; 
        this->pstate.rgstate.pars.sigma=nstat.pars.sigma; 
        rggen.setstatus(this->pstate.rgstate); }
    void init(); //initialize the memory
    void setpars(const RCGPars<U>& npars, U tol=RCG_ALPHA_TOL, unsigned long max_step=RCG_ALPHA_ITER); 
    
    U mean() {return pstate.pars.mean; }
    U sigma() {return pstate.pars.sigma; }
    
    inline rnd_type operator() () { return extract(); }
    GGEN& RGGenerator() { return rggen; } 
};

/* RndGaussian */
template<class U, class UGEN>
U RndGaussian<U,UGEN>::extract()
{
    if (this->pstate.fg2ready) 
    {
        this->pstate.fg2ready=false;
        return this->pstate.gauss2*this->pstate.pars.sigma+this->pstate.pars.mean;
    }
    
    U x1,x2,w;
	
    do
    {
        x1=2.0*rugen()-1;
        x2=2.0*rugen()-1;
        w=x1*x1+x2*x2;
    } while (w>=1.0);
	
    w=sqrt(-2.0*log(w)/w);
    this->pstate.fg2ready=true;
    this->pstate.gauss2=x2*w;
    
    return x1*w*this->pstate.pars.sigma+this->pstate.pars.mean;
}


/* Mersenne Twister (long version)*/
inline unsigned long MTRndUniformL::extract()
{
    unsigned long k ,y;
    std::valarray<unsigned long>& mt=pstate.mtvec;

    if (pstate.pos >= MT_N)
    {   
        /* generate N words at one time */
        unsigned short int kk;
        for (kk=0; kk < MT_N - MT_M; ++kk)
        {
            y = (mt[kk] & MT_UPPER_MASK) | (mt[kk + 1] & MT_LOWER_MASK);
            mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ MT_MAGIC(y);
        }
        for (; kk < MT_N - 1; ++kk)
        {
            y = (mt[kk] & MT_UPPER_MASK) | (mt[kk + 1] & MT_LOWER_MASK);
            mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ MT_MAGIC(y);
        }

        {
            y = (mt[MT_N - 1] & MT_UPPER_MASK) | (mt[0] & MT_LOWER_MASK);
            mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ MT_MAGIC(y);
        }

        pstate.pos = 0;
    }
    /* Tempering */
  
    k = mt[pstate.pos];
    k ^= (k >> 11);
    k ^= (k << 7) & 0x9d2c5680UL;
    k ^= (k << 15) & 0xefc60000UL;
    k ^= (k >> 18);

    ++pstate.pos;

    return k;
}


/* Correlated gaussian numbers */
template <class U, class GGEN>
void RndCorrGaussian<U,GGEN>::init() 
{ 
    unsigned long sz=this->pstate.pars.corr.size();
    this->pstate.mem.resize(sz);
    for (unsigned long k=0; k<sz; ++k)
        this->pstate.mem[k]=rggen();
} 

//indexing macros
#define CK_M2V(a, b, n) (n*a+b-a*(a+1)/2)
#define CK_L2V(i, j, n) ((n*(n+1)-i*(i+1))/2-j-1)


template <class U, class GGEN>
void RndCorrGaussian<U,GGEN>::get_alpha(const std::valarray<U>& corr, std::valarray<U>& alpha, U tol, unsigned long max_step)
{
    unsigned long sz=corr.size();
    std::valarray<U> V(sz*(sz-1)/2+sz);
    U tv1, tv2;
    /* first we compute M by Cholesky decomposition 
    of the minimal correlation block */  
    
    //init first column
    tv1=V[CK_L2V(0,0,sz)]=sqrt(corr[0]);
    for (unsigned long i=1; i<sz;++i)
        V[CK_L2V(i,0,sz)]=corr[i]/tv1;
    
    for (unsigned long i=1; i<sz;++i)
    {
        for (unsigned long j=1; j<i; ++j)
        {
            V[CK_L2V(i,j,sz)]=corr[i-j];
            for (unsigned long k=0; k<j; ++k)
                V[CK_L2V(i,j,sz)]-=V[CK_L2V(i,k,sz)]*V[CK_L2V(j,k,sz)];
            V[CK_L2V(i,j,sz)]/=V[CK_L2V(j,j,sz)];
        }
        V[CK_L2V(i,i,sz)]=corr[0];
        for (unsigned long k=0; k<i; ++k)
            V[CK_L2V(i,i,sz)]-=V[CK_L2V(i,k,sz)]*V[CK_L2V(i,k,sz)];
        V[CK_L2V(i,i,sz)]=sqrt(V[CK_L2V(i,i,sz)]);
    }
    
    /* then we go on with a recursive formula, 
    requiring some shifting that could be
    probably replaced with clever circular formulation */
    
    for (unsigned long istep=0; istep< max_step; ++istep)
    {
        //first we shitft the matrix. this can be done directly
        //on the collapsed vector, as it should be faster
        unsigned long ks, kc, kt;
        ks=sz*(sz-1)/2+sz-1; kc=kt=1;
        for (unsigned long k=ks; k>=sz; --k)
        {
            --ks; --kc;
            if (kc==0) {kc=kt; ++kt; --ks;}
            V[k]=V[ks];
        }
        
        //then we compute the new row with recursion formula 
        //we use M indexing
        tv1=0;
        for (unsigned long b=sz-1; b>0; --b)
        {
            V[CK_M2V(0,b,sz)]=corr[b];
            for (unsigned long k=b+1; k<sz; ++k)
                V[CK_M2V(0,b,sz)]-=V[CK_M2V(0,k,sz)]*V[CK_M2V(b,k,sz)];
            V[CK_M2V(0,b,sz)]/=V[CK_M2V(b,b,sz)];
            tv1+=V[CK_M2V(0,b,sz)]*V[CK_M2V(0,b,sz)];
        }
        V[CK_M2V(0,0,sz)]=sqrt(corr[0]-tv1);
        
        //we check how far we are from solving the limit equation
        tv1=0;
        for (unsigned long b=0; b<sz;++b) 
        {
            tv2=0;
            for (unsigned long k=b; k<sz; ++k) 
                tv2+=V[k]*V[k-b];
            tv1+=(tv2-corr[b])*(tv2-corr[b]);
        }
        
        if (tv1<=tol*tol) break;
    }   
    
    if (tv1>tol*tol) ERROR("Unable to converge alpha values within required accuracy")
                alpha.resize(sz);
    //we copy the relevant part of converged M to alpha array
    for (unsigned long i=0; i<sz;++i) alpha[i]=V[i];
}

template <class U, class GGEN>
U RndCorrGaussian<U,GGEN>::extract()
{
    U rv=(U) 0.;
    unsigned long sz=this->pstate.mem.size();

    //instead of shifting, we could make up a cyclic array through a modulus operation
    //but this will do, by now....
    for (unsigned long k=0; k<sz-1; ++k)
        this->pstate.mem[k]=this->pstate.mem[k+1];
    
    this->pstate.mem[sz-1]=rggen();
    
    for (unsigned long k=0; k<sz; ++k)
        rv+=this->pstate.mem[k]*this->pstate.alpha[k];
    
    return rv*this->pstate.pars.sigma+this->pstate.pars.mean;
}

template <class U, class GGEN>
void RndCorrGaussian<U,GGEN>::setpars(const RCGPars<U>& npars, U tol, unsigned long max_step)
{
    //make sure that the gaussian random numbers generator has mean zero and variance one: the actual values will be computed a-posteriori
    rggen.setpars(RGPars<U>((U) 0,(U) 1));
    this->pstate.pars=npars;
    //this is tough, we need to check that the correlation is feasible, and set the alphas to the correct values
    //check_corr(npars.corr);
    get_alpha(this->pstate.pars.corr, this->pstate.alpha,tol,max_step);
}

} //ends namespace toolbox

#endif // ends ifdef __RNDGEN_DEC_H
