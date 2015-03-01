/* Minimization library. Mostly implemented as templates.
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/
   
#ifndef __MINSEARCH_H
#define __MINSEARCH_H 0
#include "tbdefs.hpp"
#include "rndgen.hpp"
#include <iomanip>
#include "ioparser.hpp"
#include <fstream>

namespace toolbox {
/***********************************************************************
  COLLECTION OF FUNCTIONS FOR MINIMIZATION. the functions are templates
  acting on a "function object" which must define a series of interface
  member functions: 
   size() which returns the size of the search space
   set/get_vars(valarray) which set and retrieve the values of the 
                          parameters which must be optimized
   get_value() which returns the value of the function for the parameters
   get_gradiens() which returns the gradient
   get_hessian() which return the hessian
   bailout() which the minimizer might call at the end of every loop, 
             breaking out if the return value is true
i.e.
    class minf_prototype {
    public:
        unsigned long size();
        void set_vars(cost std::valarray<double>& vars);
        void get_vars(std::valarray<double>& vars) const;
        void get_value(double &rv);
        void get_gradient(std::valarray<double> &rg);
        void get_hessian(std::valarray<std::valarray<double<> &rh);
        bool bailout();
    };
    get_gradient and get_hessian are required only by the methods
    using gradients and hessian respectively.
************************************************************************/

template <double (*f) (const std::valarray<double>& ), 
            void (*gf) (const std::valarray<double>&, std::valarray<double>&),
            void (*hf) (const std::valarray<double>&, std::valarray<std::valarray<double> >&)>
class StdMinFunction {
private:
    std::valarray<double> coords;
    
public:
    inline unsigned long size() { return coords.size(); }
    
    inline void set_vars(const std::valarray<double>& rv) { if (coords.size()!=rv.size()) coords.resize(rv.size()); coords=rv; }
    inline void get_vars(std::valarray<double>& rv) const { if (coords.size()!=rv.size()) rv.resize(coords.size()); rv=coords; }
    inline void get_value(double& rv) { rv=f(coords); }
    inline void get_gradient(std::valarray<double>& rv) const 
    { if (coords.size()!=rv.size()) rv.resize(coords.size()); gf(coords,rv); }
    inline void get_hessian(std::valarray<std::valarray<double> >& rv) const 
    {
#ifdef DEBUG
        if(gf==NULL) ERROR("Trying to get the hessian, but no function have been provided to the template");
#endif
        if (coords.size()!=rv.size()) rv.resize(coords.size()); hf(coords,rv); 
    } 
    bool bailout() { return false; }
};

template <double (*f) (const std::valarray<double>& )>
class FOnlyMinFunction
{
protected:
    std::valarray<double> coords;

public:
    FOnlyMinFunction<f>& operator=(const FOnlyMinFunction<f>& rf)
    { if (*rf==this) return *this; coords.resize(rf.coords.size()); coords=rf.coords; }
    
    inline unsigned long size() { return coords.size(); }

    inline void set_vars(const std::valarray<double>& rv) { if (coords.size()!=rv.size()) coords.resize(rv.size()); coords=rv; }
    inline void get_vars(std::valarray<double>& rv) const { if (coords.size()!=rv.size()) rv.resize(coords.size()); rv=coords; }
    inline void get_value(double& rv) { rv=f(coords); }
    inline void get_gradient(std::valarray<double>& rv) const 
    { ERROR("Minimization function class does not define a gradient function"); }
    inline void get_hessian(std::valarray<std::valarray<double> >& rv) const 
    { ERROR("Minimization function class does not define a hessian function"); } 
    bool bailout() { return false; }
};

inline double FOFMDummy(const std::valarray<double>& x) { return 0; }
typedef FOnlyMinFunction<FOFMDummy> FOnlyMinFunctionBase; 
    
#define __SMPLX_RHO 1.
#define __SMPLX_CHI 2.
#define __SMPLX_GAMMA 0.5.
#define __SMPLX_SIGMA 0.5.
typedef struct _SimplexOpts {
    double rho, chi, gamma, sigma;
    _SimplexOpts(): rho(1.), chi(2.), gamma(0.5), sigma(0.5) {}
} SimplexOptions;


/***********************************************************************
 DOWNHILL SIMPLEX MINIMIZER
 two convergency parameters, default is both conditions must be
 satisfied:
 * parameter one is the diff. of the value of the function, at the highest 
   and lowest point of the simplex. default is bailout when 
   the change is smaller than __TB_STD_EPS
 * parameter two is the simplex spread, default is bailout when its 
   absolute value goes below __TB_STD_EPS
************************************************************************/

std::valarray<std::valarray<double> > make_simplex(std::valarray<double> ip, double dr, double fr=1.);

template<class FCLASS>
void min_simplex (
                  FCLASS& f,
                  const std::valarray<std::valarray<double> > & initial_simplex,
                  std::valarray<double>& rpos, double& rvalue,
                  IterOptions<double,2> iops=IterOptions<double,2>(
                          __TB_STD_MAXSTEP,
                          fixarray<double,2>(__TB_STD_EPS, __TB_STD_EPS),
                         fixarray<double,2>(0.,0.),
                        fixarray<double,2>(ichk_change, ichk_default)),
                  const SimplexOptions so=SimplexOptions()
                 )
{
#ifdef DEBUG
    if (initial_simplex.size()==0 || initial_simplex.size()-1!=initial_simplex[0].size())
        ERROR("Simplex size should be one more than the search space dimension")
#endif
    unsigned long sdim=initial_simplex.size()-1;
    std::valarray<std::valarray<double> > s_pos=initial_simplex;
    std::valarray<double> s_val(sdim+1);
    std::valarray<double> x_center(sdim), x_new(sdim), x_new2(sdim); 

    double f_new, f_new2,s_spread;
    std::valarray<double> va_tmp(sdim); double d_tmp;
    std::cerr<<"Initialize simplex: \n";
    for (unsigned long i=0; i<sdim+1; i++)
    {
        f.set_vars(s_pos[i]);
        f.get_value(s_val[i]);
        std::cerr<<i<<":"<<s_val[i]<<"  ";
    }
    std::cerr<<"\n";

    //stupid initial sort for the s_val and s_pos arrays
    for (unsigned long i=0; i<sdim+1; i++)
        for (unsigned long j=i+1; j<sdim+1; j++)		
            if (s_val[i]>s_val[j])
            {
                d_tmp=s_val[i]; va_tmp=s_pos[i];
                s_val[i]=s_val[j]; s_pos[i]=s_pos[j];
                s_val[j]=d_tmp; s_pos[j]=va_tmp;
            }
    
    //compute the "corrected" center of mass of the simplex
    double os=s_val[0];
    while (!iops)
    {
        std::cerr<<s_val[0]<<" simplex value ("<<s_val[0]-os<<")\n";
        os=s_val[0];
        x_center=0.;
        for (unsigned long j=0; j<sdim; ++j)
        {
            for (unsigned long i=0; i<sdim; ++i) x_center[j]+=s_pos[i][j];
            x_center[j]/=sdim;
        }
        for (unsigned long j=0; j<sdim; ++j)
            x_new[j]=x_center[j]+so.rho*(x_center[j]-s_pos[sdim][j]);
        f.set_vars(x_new); f.get_value(f_new);
        
        if (f_new<s_val[0]) 
        {
            //expansion
            for (unsigned long j=0; j<sdim; ++j)
                x_new2[j]=x_center[j]+so.chi*(x_new[j]-x_center[j]);
            f.set_vars(x_new2); f.get_value(f_new2);
            
            if (f_new2>=f_new) {x_new2=x_new; f_new2=f_new;}
        }
        else if (f_new>=s_val[sdim-1])
        {
            //contraction
            if (f_new<s_val[sdim])
            {
                //outside
                for (unsigned long j=0; j<sdim; ++j)
                    x_new2[j]=x_center[j]+so.gamma*(x_new[j]-x_center[j]);
                f.set_vars(x_new2); f.get_value(f_new2);
                s_val[sdim]=f_new;
            }
            else 
            {
                //inside
                for (unsigned long j=0; j<sdim; ++j)
                    x_new2[j]=x_center[j]+so.gamma*(s_pos[sdim][j]-x_center[j]);
                f.set_vars(x_new2); f.get_value(f_new2);
            }
        }
        else 
        {
            //reflection
            x_new2=x_new; f_new2=f_new;
        }
        if (f_new2>s_val[sdim])
        {
            //shrink
            for (unsigned long i=1; i<sdim+1; i++)
            {
                s_pos[i]-=s_pos[0]; s_pos[i]*=so.sigma; s_pos[i]+=s_pos[0];
                f.set_vars(s_pos[i]); f.get_value(s_val[i]);
            }
            
            //resort the simplex
            for (unsigned long i=0; i<sdim+1; i++)
                for (unsigned long j=i+1; j<sdim+1; j++)		
                    if (s_val[i]>s_val[j])
                    {
                        d_tmp=s_val[i]; va_tmp=s_pos[i];
                        s_val[i]=s_val[j]; s_pos[i]=s_pos[j];
                        s_val[j]=d_tmp; s_pos[j]=va_tmp;
                    }
        }
        else
        {
            //insert the new point
            unsigned long i;
            for (i=0; s_val[i]<f_new2; i++);
            for (unsigned long j=sdim;j>i;--j)
            { s_val[j]=s_val[j-1]; s_pos[j]=s_pos[j-1]; }
            s_val[i]=f_new2; s_pos[i]=x_new2;			
        }	
        s_spread=0;
        for (unsigned long i=1; i<sdim+1; ++i)
            for (unsigned long j=0; j<sdim; ++j)
                s_spread+=(s_pos[1][j]-s_pos[i][j])*(s_pos[1][j]-s_pos[i][j]);
        s_spread=sqrt(s_spread);
        iops.setval(fabs(s_val[sdim]-s_val[0]),0); 
        iops.setval(s_spread,1);
        if (f.bailout()) break;
    }
#ifdef DEBUG
    if (iops.maxstep_reached() ) ERROR("Minimization finished before reaching required accuracy")
#endif
    rpos.resize(s_pos[0].size());
    rpos=s_pos[0]; rvalue=s_val[0];
}



/***********************************************************************
 NESTED SAMPLING MINIMIZER
 set number of walkers K, and generate K uniformly random samples {x_i, 1 <=i <=K}, evaluate objective function for all of them:  y_i = f(x_i)

loop 
  sort vector y_i, pick i with highest y_i
  set limit<- y_i
  discard sample x_i, replace it with a clone of randomly selected j: x_i <- x_j
  perform random walk on x_i, only subject to constraint f(x_i) < limit, adjust step size to maintain acceptance ratio
end

************************************************************************/
typedef struct _NestSampOpts {
    double mc_step, mc_wall;
    unsigned long steps;
    double adapt, adapt_target;
    _NestSampOpts(): mc_step(1.0), mc_wall(0), steps(10), adapt(1.5), adapt_target(0.3) {}
} NestSampOptions;

std::valarray<std::valarray<double> > make_walkers(unsigned long ndim, unsigned long nwalker, double dr, MTRndUniform rngu=MTRndUniform());

template<class FCLASS>
void min_nestsamp (
                  FCLASS& f,
                  const std::valarray<std::valarray<double> > & initial_walkers,
                  std::valarray<double>& rpos, double& rvalue,
                  IterOptions<double,2> iops=IterOptions<double,2>(
                          __TB_STD_MAXSTEP,
                          fixarray<double,2>(__TB_STD_EPS, __TB_STD_EPS),
                         fixarray<double,2>(0.,0.),
                        fixarray<double,2>(ichk_change, ichk_default)),
                  const NestSampOptions nso=NestSampOptions(),
                  MTRndUniform rngu=MTRndUniform()
                 )
{
#ifdef DEBUG
    if (initial_walkers.size()<2 )
        ERROR("Need at least two walkers to perform nested sampling.")
#endif
    unsigned long nwalk=initial_walkers.size(), wdim=initial_walkers[0].size(), imax, ispawn;
    std::valarray<std::valarray<double> > w_pos(initial_walkers);
    std::valarray<double> w_val(nwalk), npos(wdim);
    std::valarray<unsigned long> w_nok(nwalk);
    double maxw, nval, mcstep;
    
    
    std::cerr<<"Initialize walker values: \n";
    maxw=-1e-100;
    for (unsigned long i=0; i<nwalk; i++)
    {
        f.set_vars(w_pos[i]);
        f.get_value(w_val[i]);        
        if (maxw<w_val[i]) maxw=w_val[i];
    }
    w_nok=0; rpos=w_pos[0]; rvalue=w_val[0];
    std::cerr<<"\n";
    rpos.resize(w_pos[0].size());
    rpos=w_pos[0]; rvalue=w_val[0]; 
    w_nok=0; mcstep=nso.mc_step;
    while (!iops)
    {
        // does some random walking
        double accept=0;
        for (unsigned long k=0; k<nso.steps; ++k)
        {
           std::cerr<<"NS step "<<k<<"\n";
            for (unsigned long i=0; i<nwalk; ++i)
            {               
               for (unsigned long j=0; j<wdim; ++j) npos[j] = rngu()-0.5;
               npos *=  mcstep;
               npos+=w_pos[i];
                              
               if (nso.mc_wall>0) 
               { for (unsigned long j=0; j<wdim; ++j) if (fabs(npos[j])*2>nso.mc_wall) {npos[j]=w_pos[i][j]; accept--;} }
               
               
               f.set_vars(npos); f.get_value(nval);
               if (nval<=maxw) { w_pos[i]=npos; w_val[i]=nval; accept+=wdim;} // also keeps track of the accepted moves
               if (nval<rvalue) { 
                  rvalue=nval; rpos=w_pos[i]; 
                  std::cerr<<"New low energy structure "<<rvalue<<"\n";
                  std::cout<<"# New round\n";
                  for (unsigned long j=0; j<wdim; j+=2) std::cout<<j<<" "<<rpos[j]<<" "<<rpos[j+1]<<"\n";
               } // stores the lowest-energy walker
            }            
        }
        
        // selects the hottest replica
        maxw=w_val[0]; imax=0;        
        for (unsigned long i=1; i<nwalk; i++)
        { if (w_val[i]>maxw) { maxw=w_val[i]; imax=i; } }
        
        // kills the hot replica and spawns it on a random walker
        ispawn = imax;
        while (ispawn==imax) ispawn=rngu()*nwalk;
        w_pos[imax] = w_pos[ispawn];
        w_val[imax] = w_val[ispawn];
        
        // adaptive mc step
        accept/=(nwalk*nso.steps*wdim);        
        if (accept>nso.adapt_target) mcstep*=nso.adapt;
        else mcstep/=nso.adapt;
        std::cerr<<"Setting Nested Sampling threshold "<<maxw<<"\n";
        std::cerr<<"Mean acceptance "<<accept<<", new step "<<mcstep<<"\n";
        
        // sets the checks for bailing out of the loop
        iops.setval(fabs(maxw-rvalue),0); 
        iops.setval(mcstep,1);
        
        
     }
#ifdef DEBUG
    if (iops.maxstep_reached() ) ERROR("Minimization finished before reaching required accuracy")
#endif
}




/***********************************************************************
             LINE SEARCH (no gradients)
 ***********************************************************************/
template<class FCLASS> class Vec2Lin {
    private:
    std::valarray<double> p, d;
    FCLASS pf;
    public:
        Vec2Lin(const std::valarray<double>& np, const std::valarray<double>& nd) :
            p(np), d(nd) 
            {
                if (p.size()<1 || d.size()!=p.size() ) ERROR("Error initializing Vec2Lin with given point & dir");
            };
        Vec2Lin(const FCLASS& npf, const std::valarray<double>& nd) :
                pf(npf), d(nd) 
                {
                    pf.get_vars(p);
                    if (p.size()<1 || d.size()!=p.size() ) ERROR("Error initializing Vec2Lin with given point & dir");
                };
        Vec2Lin& operator=(const Vec2Lin& ov)
        {
            if (this==&ov) return *this;
            p.resize(ov.p.size()); p=ov.p;
            d.resize(ov.d.size()); d=ov.d;
            pf=ov.pf;
            return *this;
        }
        double operator() (double x)
        {
            std::valarray<double> lp(d);
            lp*=x; lp+=p; pf.set_vars(lp);
            double rv; pf.get_value(rv);
            return rv;
        }
};


class LineSearchOpts {
public:
    double istep, maxstep;
    unsigned long maxiter;
    double lsgold, lsigold, lstiny, lstol;
    std::valarray<double> dir;
    LineSearchOpts(): istep(1.), maxstep(100.), maxiter(100), 
                    lsgold(1.618034), lsigold(1./lsgold), lstiny(1e-30), lstol(1e-10) {}
};
 

template<class FCLASS> void ls_bracket(
          Vec2Lin<FCLASS>& v2l, double& a, double& b, double& c,
          double& fa, double& fb, double& fc, const LineSearchOpts& op=LineSearchOpts())
{
    double tmp,r,q,u,ul,fu;
    fa=v2l(a); fb=v2l(b);  //makes sure fb<fa
    if (fa<fb) {tmp=a; a=b; b=tmp; tmp=fa; fa=fb; fb=tmp; }
    
    c=b+op.lsgold*(b-a);
    std::cerr<<"im brack "<<a<<","<<b<<","<<c<<std::endl;
    fc=v2l(c);
    std::cerr<<"fm brack "<<fa<<","<<fb<<","<<fc<<std::endl;
    while (fb>=fc)
    {
        std::cerr<<"bracketing ("<<a<<","<<fa-fb<<") ("<<b<<","<<fb<<") ("<<c<<","<<fc-fb<<")\n";
        r=(b-a)*(fb-fc); q=(b-c)*(fb-fa);
        u=b-((b-c)*q-(b-a)*r)/(2.*(fabs(q-r)>op.lstiny?fabs(q-r):op.lstiny)*(q-r>=0.?1.:-1.));
        ul=b+op.maxstep*(c-b);
        
        if((b-u)*(u-c)>0) 
        {
            fu=v2l(u);
            if (fu<fc) { a=b; fa=fb; b=u; fb=fu; return; }
            else if (fu>fb) {c=u; fc=fu; return; }
            u=c+op.lsgold*(c-b);
            fu=v2l(u);
        }
        else if ((c-u)*(u-ul)>0)
        {
            fu=v2l(u);
            if (fu>fc) {b=c; c=u; u=c+op.lsgold*(c-b); fb=fc; fc=fu; fu=v2l(u);}
        }
        else if ((u-ul)*(ul-c)>=0.) {u=ul; fu=v2l(u); }
        else { u=c+op.lsgold*(c-b); fu=v2l(u); }
        
        a=b; b=c; c=u; fa=fb; fb=fc; fc=fu;
        if (fa==fb &&fb==fc) { std::cerr<<"bracketing failed at a flat point!\n"; return; }
    }
    return;
}

template<class FCLASS> void ls_brent(
   Vec2Lin<FCLASS>& v2l, const double& ax, const double& bx, const double& cx,
   double& x, double& fx, const LineSearchOpts& op=LineSearchOpts())
{
    double a,b,d,e,et,fu,fv,fw,p,q,r,t1,t2,u,v,w,xm;
    a=(ax<cx?ax:cx);
    b=(ax>=cx?ax:cx);
    x=w=v=bx; 
    e=0;
    fv=fw=fx=v2l(x);

    for (unsigned long liter=0; liter<op.maxiter; ++liter)
    {
        xm=0.5*(a+b);
        t1=op.lstol*fabs(x)+op.lstiny; 
        t2=t1+t1;

        if (fabs(x-xm)<=t2-0.5*(b-a)) return;
        if (fabs(e)>t1)
        {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2*(q-r);
            if (q>0) p=-p;
            q=fabs(q);
            et=e;
            e=d;
                                                        
            if ((abs(p)>=fabs(0.5*q*et) || p<=q*(a-x)
                 || p>=q*(b-x))) goto one;
            d=p/q; u=x+d;
            if ((u-a)<t2 || (b-u) <t2) d=(xm-x>=0?fabs(t1):-fabs(t1));
            goto two;
        }
one: 
        if (x>=xm) e=a-x; else e=b-x;
        d=op.lsigold*e;
two:
        if (fabs(d)>=t1) u=x+d; else u=x+(d>=0?fabs(t1):-fabs(t1));
        fu=v2l(u);
        if(fu<=fx)
        {
            if (u>=x) a=x; else b=x;
            v=w; fv=fw;
            w=x; fw=fx;
            x=u; fx=fu;
        }
        else
        {
            if (u<x) a=u; else b=u;
            if (fu<=fw || w==x)
            { v=w; fv=fw; w=u; fw=fu; }
            else if (fu<-fv || v==x || v==w) {v=u; fv=fu; }
        }
    }
    return;
}

template<class FCLASS>
void min_linesearch (
    FCLASS& f,
    const std::valarray<double> & initial_pos,
    std::valarray<double>& rpos, double& rvalue,
    double &rstep,
    const LineSearchOpts& op
)
{
    double f0;
    std::cerr<<"getting initial value\n";
    f.set_vars(initial_pos); f.get_value(f0);
    Vec2Lin<FCLASS> v2l(f,op.dir);
    double t,ft, a,b,c,fa,fb,fc;
    a=0; b=op.istep; 
    std::cerr<<"bracketing"<<std::endl;
    ls_bracket(v2l,a,b,c,fa,fb,fc,op);
    std::cerr<<"bracketed ("<<a<<","<<fa-f0<<") ("<<b<<","<<fb-f0<<") ("<<c<<","<<fc-f0<<")"<<std::endl;
    ls_brent(v2l,a,b,c,t,ft,op);
    std::cerr<<"converged to "<<t<<","<<ft-f0<<std::endl;
    rpos.resize(initial_pos.size()); rpos=initial_pos;
    if (ft<f0)   //something may go wrong. we DON'T go uphill in any case!
    {
        rpos+=t*op.dir; rvalue=ft; rstep=fabs(t);
    }
    else { rvalue=f0; rstep=0.; }
        
}

/***********************************************************************
             POWELL MINIMIZATION (no gradients) NR 10.5
 ***********************************************************************/
 class PowellOpts{
     public:
         LineSearchOpts linesearch;
         double tol, drnd;  bool fmxdiscard, fadjstep;
         unsigned long maxiter;
         PowellOpts(): tol(1e-5), drnd(-1.), fmxdiscard(false), fadjstep(true), maxiter(100) {}
 };

template<class FCLASS>
void min_powell (
      FCLASS& f,
      const std::valarray<double> & initial_pos,
      std::valarray<double>& rpos, double& rvalue,
      const PowellOpts& op,
      std::valarray<std::valarray<double> >& u=std::valarray<double>(0)
     )
{
    RndGaussian<double,StdRndUniform> rng;
    std::valarray<double> pos(initial_pos), npos(pos), step(pos);
    step=op.linesearch.istep;
    unsigned long sz=pos.size();
    LineSearchOpts ls(op.linesearch);
    ls.dir.resize(sz);
    //initialize set of directions
    if (u.size()!=sz) 
    {
        u.resize(sz);
        for (int i=0;i<sz;++i) { u[i].resize(sz); u[i]=0.;  u[i][i]=1.; }
    }
    double fp, fn, fo, mde=0., lsstep;
    unsigned long mid=0;
    f.set_vars(pos); f.get_value(fn);
    //std::cerr<<"INITIAL VALUE IN MINPOWELL "<<fn<<"\n";
    //std::cerr<<"INITIAL PARS IN MINPOWELL "<<pos<<"\n";
    for (unsigned long iter=0; iter<op.maxiter; ++iter)
    {
        std::cerr<<" * STARTING POWELL ITERATION "<<iter<<"   *\n";
        for (unsigned long id=0; id<sz; ++id)
        {
            std::cerr<<id<<"("<<step[id]<<")::";
            for (int i=0; i<sz;++i) std::cerr<<u[id][i]<<" ";
            std::cerr<<"\n";
        }
        
        //ONE POWELL METHOD ITERATION
        std::valarray<double>po(pos);  fp=fn; mde=0.; mid=0;
        std::cerr<<" Initial value: "<<fp<<"\n";
        for (unsigned long id=0; id<sz; ++id)
        {
            std::cerr<<" * minsearch along u["<<id<<"]\n";
            ls.dir=u[id]; ls.istep=step[id]; fo=fn; 
            min_linesearch(f,pos,npos,fn,lsstep,ls);
            pos=npos; 
            if (op.fadjstep) step[id]=(step[id]*0.5+lsstep);
            std::cerr<<" step size: "<<lsstep<<" function val: "<<fn-fp<<"\n";
            if(fo-fn>mde)
            { mde=fo-fn;  mid=id; }
        }
        std::cerr<<" increment "<<fn-fp<<"\n";
        if(2.*(fp-fn)<=op.tol*(fabs(fp)+fabs(fn))) break; //bailout
        //reasons not to update direction set
        if (op.fmxdiscard)
        {
            npos=pos; npos*=2.; npos-=po;
            f.set_vars(npos); f.get_value(fo);
            std::cerr<<" checup "<<fo-fp<<"\n";
            if (fo>=fp) continue;
            if (2.*pow((fp-2*fn+fo)*(fp-fn-mde),2.)-mde*pow(fp-fo,2.) >=0) continue;
            std::cerr<<"SETTING NEW SEARCH DIRECTION!\n";
            ls.dir=pos; ls.dir-=po;
            double nn=0.; 
            for (int i=0; i<sz; ++i) nn+=ls.dir[i]*ls.dir[i]; 
            nn=1./sqrt(nn); ls.dir*=nn;
            min_linesearch(f,pos,npos,fn,lsstep,ls);
            pos=npos; u[mid]=u[sz-1]; step[mid]=step[sz-1]; u[sz-1]=ls.dir; step[sz-1]=op.linesearch.istep;
        }
        else
        {
            for (int i=0; i<sz-1; ++i) { u[i]=u[i+1]; step[i]=step[i+1]; }
            ls.dir=pos; ls.dir-=po;
            double nn=0.; 
            for (int i=0; i<sz; ++i) nn+=ls.dir[i]*ls.dir[i]; 
            nn=1./sqrt(nn); ls.dir*=nn;
            u[sz-1]=ls.dir; step[sz-1]=op.linesearch.istep;
            min_linesearch(f,pos,npos,fn,lsstep,ls);
            pos=npos;
        }
        if (op.drnd>0.) 
        {
            std::cerr<<"RANDOMIZING DIRECTIONS\n";
            for (unsigned long id=0; id<sz; ++id)
            {
                double nn=0.;
                for (int i=0; i<sz; ++i) 
                { u[id][i]+=rng()*op.drnd/sz; nn+=u[id][i]*u[id][i]; }
                nn=1./sqrt(nn);  u[id]*=nn;
            }
        }
    }
}
        
/***********************************************************************
SIMULATED ANNEALING MINIMIZER
************************************************************************/
typedef struct _AnnealingOpts {
    double temp_init, temp_final, mc_step;
    unsigned long steps;
    double adapt, drnd; bool fpowell;
    _AnnealingOpts(): temp_init(1.), temp_final(1e-3), mc_step(0.1) , steps(1000), adapt(1.), drnd(0.), fpowell(false) {}
} AnnealingOptions;

template<class FCLASS, class RNG>
void sim_annealing (
        FCLASS& f,
        const std::valarray<double> & init_vars,
        std::valarray<double>& rpos, double& rvalue,
        const AnnealingOptions ao=AnnealingOptions(),
        std::valarray<std::valarray<double> > u=std::valarray<std::valarray<double> >(0),
        RNG rngen=RNG()
        )
{
    unsigned long nv=init_vars.size(), nu=u.size();
    std::valarray<double> pos(init_vars), upos(nv), npos(nv);
    
    if (nu==0)
    {
        u.resize(nv); nu=nv;
        for (int i=0; i<nv; ++i) { u[i].resize(nv); u[i]=0.; u[i][i]=1.; }
    }
    else for (int i=0; i<nu; ++i) if (u[i].size()!=nv) ERROR("Direction vector size mismatch with state vector");
    
    std::valarray<double> step(ao.mc_step,nu);
    std::valarray<long> accept(0.,nu), tstep(0.,nu);
    double nrg, nnrg, temp=ao.temp_init, ts=std::exp(log(ao.temp_final/ao.temp_init)/ao.steps);
    
    f.set_vars(pos);
    f.get_value(nrg);
    //std::cerr<<"Starting off VALUE: "<<nrg<<"\n";
    //std::cerr<<"Starting off POS: "<<pos<<"\n";
    bool fmoved; 
    unsigned long bu; double be;
    for (unsigned long is=0; is< ao.steps; ++is)
    {
        upos=pos; tstep+=1; fmoved=false; bu=0; be=0.;
        for (unsigned long iu=0; iu<nu; ++iu)
        {
            npos=pos; npos+=u[iu]*step[iu]*(rngen()-0.5);
            f.set_vars(npos); f.get_value(nnrg);
            if (nnrg-nrg<be) { be=(nnrg-nrg); bu=iu; }
            //metropolis step
            if (rngen()<=std::exp((nrg-nnrg)/temp))
            {
                accept[iu]++; fmoved=true;
                pos=npos; nrg=nnrg;
            }
            if (accept[iu]>tstep[iu]/2) 
            {
                step[iu]*=ao.adapt;
            }
            else
            {
                step[iu]/=ao.adapt;
            }
            //std::cerr<<":: "<<nrg<<","<<nnrg<<" "<< accept[iv]<<" out of "<<in<<" for "<<iv<<"\n";
        }
        std::perr<<"********* SIM ANNEALING STATS ("<<std::setw(7)<< is <<") *************\n";
        std::perr<<" energy: "<<nrg<<"\n";
        std::perr<<" temperature: "<<temp<<"\n";
        std::perr<<" mean acceptance: "<<accept.sum()*1.0/tstep.sum()<<"\n";
        std::perr<<" mean step: "<<step.sum()*1.0/nu<<"\n";
        std::perr<<step;
        
        if (ao.drnd>0.)
        {
            //randomize/normalize
            double wu;
            for (unsigned long iu=0; iu<nu; ++iu) 
            {
                wu=0.; for (int i=0; i<nv; ++i) { u[iu][i]+=ao.drnd*(rngen()-0.5)/nv; wu+=u[iu][i]*u[iu][i]; }
                u[iu]*=1./sqrt(wu);
            }
        }
        if (ao.fpowell && fmoved)
        {
            //saves best energy direction;
            std::valarray<double> bd=u[bu];
            for (unsigned long iu=0; iu<nu-1; ++iu) 
            { u[iu]=u[iu+1]; step[iu]=step[iu+1]; accept[iu]=accept[iu+1]; tstep[iu]=tstep[iu+1]; }
            u[nu-1]=pos-upos; double wu=0.;
            for (int i=0; i<nv; ++i) { wu+=u[nu-1][i]*u[nu-1][i]; } wu=sqrt(wu);
            std::perr<<" total displacement: "<<wu<<"\n";
            if (wu>0.) {
                tstep[nu-1]=accept[nu-1]=0; step[nu-1]=2.*wu; u[nu-1]*=1./wu;
                //mix with best energy direction
                u[nu-1]+=bd; wu=0.; for (int i=0; i<nv; ++i) { wu+=u[nu-1][i]*u[nu-1][i]; } 
                wu=sqrt(wu); u[nu-1]*=1./wu; 
            }
        }
        temp*=ts;
    }
    
    f.set_vars(pos);   //reset status to last accepted position
    rvalue=nrg;
/*    std::perr<<"SIM ANN STATS:\n";
    for (iv=0; iv<nv; ++iv)
    {
        std::cerr<<(accept[iv]*nv*1./ao.steps)<<" ";
    }
    std::cerr<<"\n*******************\n";
    */
}

template<class FCLASS>
void sim_annealing (
     FCLASS& f,
     const std::valarray<double> & init_vars,
     std::valarray<double>& rpos, double& rvalue,
     const AnnealingOptions ao=AnnealingOptions(),
     std::valarray<std::valarray<double> > u=std::valarray<std::valarray<double> >(0)
 )
{
    sim_annealing<FCLASS,MTRndUniform>(f,init_vars,rpos,rvalue,ao,u,MTRndUniform());
}

/***********************************************************************
             STEEPEST DESCENT
 ***********************************************************************/
 class SteepestOpts{
     public:
         LineSearchOpts linesearch;
         double tol;          unsigned long maxiter;
         SteepestOpts(): tol(1e-5), maxiter(100) {}
 };

 template<class FCLASS>
 void min_steepest (
      FCLASS& f, 
      const std::valarray<double> & initial_pos,
      std::valarray<double>& rpos, double& rvalue,
      const SteepestOpts& op
      )
{
    std::valarray<double> pos(initial_pos), npos(pos);
    
    unsigned long sz=pos.size(); double fo, fn, lsstep;
    LineSearchOpts ls(op.linesearch);
    ls.dir.resize(sz);
    
    for (unsigned long istep=0; istep<op.maxiter; ++istep)
    {
        f.set_vars(pos); 
        f.get_gradient(ls.dir);
        
        min_linesearch(f,pos,npos,fn,lsstep,ls);
        ls.istep=lsstep*0.5;
        pos=npos;
        std::cerr<<"Step size: "<<ls.istep<<"\n";
        std::cerr<<"Step dir: "<<ls.dir[0]<<"\n";
        std::cerr<<"Current value is "<<fn<<"\n";
    }
    rpos=pos; rvalue=fn;
}

/***********************************************************************
             CONJUGATE GRADIENT
 ***********************************************************************/
 class ConjGradOpts{
     public:
         LineSearchOpts linesearch;
         double tol;  unsigned long maxiter;
         ConjGradOpts(): tol(1e-5), maxiter(100) {}
 };

#define _CG_DEF_STEP 1e-3
template<class FCLASS>
void min_conjgrad (
                    FCLASS& f, 
                    const std::valarray<double> & initial_pos,
                    std::valarray<double>& rpos, double& rvalue,
                    ConjGradOpts& op
                   )
{
    std::valarray<double> pos(initial_pos), npos(pos);
    
    unsigned long sz=pos.size(); double fo, fn, lsstep, gamma, gg;
    LineSearchOpts ls(op.linesearch);
    std::valarray<double> g(sz), og(sz);
    if (ls.dir.size()!=sz) 
    {
        ls.dir.resize(sz); 
    
        f.set_vars(pos); 
        f.get_gradient(ls.dir); //ls.dir*=-1.0;
    }

    og=ls.dir;
    for (unsigned long istep=0; istep<op.maxiter; ++istep)
    {
        
        min_linesearch(f,pos,npos,fn,lsstep,ls);
        f.set_vars(npos); 
        f.get_gradient(g); //g*=-1.0;
        
        gamma=gg=0.0; for (unsigned long i=0; i<sz; ++i) { gg+=og[i]*og[i]; gamma+=(g[i]-og[i])*g[i]; } gamma/=gg;
        
        og=g;
        ls.dir*=gamma; ls.dir+=g;
        
        if (lsstep<=0.0) lsstep=_CG_DEF_STEP;
        ls.istep=lsstep*0.5;
        pos=npos;
        std::cerr<<"Step size: "<<ls.istep<<"\n";
        //std::cerr<<"Step dir: "<<ls.dir<<"\n";
        std::cerr<<"Current value is "<<fn<<"\n";
    }
    op.linesearch.dir.resize(sz); op.linesearch.dir=ls.dir;  op.linesearch.istep=ls.istep;
    rpos=pos; rvalue=fn;
}


/***********************************************************************
PARALLEL TEMPERING ANNEALING MINIMIZER
************************************************************************/
typedef struct _ParaOpts {
    double temp_init, temp_final, temp_factor;
    double dt, tau, tau_avg;
    unsigned long steps, replica;
    _ParaOpts(): temp_init(1.), temp_final(1e-3), temp_factor(1.5), steps(10000), replica(5), dt(0.1), tau(100.0),
                 tau_avg(100.0) {}
} ParaOptions;


template<class FCLASS, class RNG>
void para_temp (
        FCLASS& f,
        const std::valarray<double> & init_vars,
        std::valarray<double>& rpos, double& rvalue,
        const ParaOptions po=ParaOptions(),
        RNG rngen=RNG()
        )
{
    unsigned long nv=init_vars.size(), nr=po.replica;
    double swp, nrgmin;
    std::valarray<double> temp(nr), nrg(nr), kin(nr), cns(nr), pmin(nv), pswp(nv), vdir(nv), wte(nr), wte_dw(nr), wte_up(nr);
    std::valarray<std::valarray<double> > pos(init_vars,nr), vel(init_vars,nr), grad(std::valarray<double>(nv),nr);
    std::valarray<std::vector<double> >  
      wte_history(std::vector<double>(0),nr), wte_heights(std::vector<double>(0),nr), wte_widths(std::vector<double>(0),nr);
    std::valarray<double> wte_gamma(po.temp_factor,nr);
    std::valarray<long> ireplica(nr); std::valarray<double> vreplica(nr);
    double wtev, wtedv;
    
    RndGaussian<double,RNG> grngen(rngen());
    
    //initialize replica temperatures
    double c1, c2, c3, c4;

    /* we compute "smart" running averages with exponentially decaying memory: let f=exp(-dt/tau_m)
      then at each time step tot<--f*tot+y   n<--f*n+1 so after T steps, we will have
      tot=sum_i<T y_i exp(-dt*i/tau_m) ; n=sum_i<T exp(-dt*i/tau_m) 
      and tot/n will yield a sort of running average of the quantity y
    */
    std::valarray<double> totv(0.0,nr), totv2(0.0,nr);  double totn=0.0, totf;
    totf=std::exp(-1.0*po.dt/po.tau_avg);
 
//    temp[0]=po.temp_init; for (unsigned long i=1; i<nr; ++i) temp[i]=temp[i-1]*po.temp_factor;
    temp[0]=po.temp_init; for (unsigned long i=1; i<nr; ++i) temp[i]=temp[i-1]*(1+po.temp_factor/sqrt(nv));
    double ts=std::exp(log(po.temp_final/po.temp_init)/po.steps);

    std::cerr<<"para temp optimizer. steps: "<<po.steps<<" replica: "<<nr<<"\n";
    std::cerr<<"temperatures "<<temp;
    std::cerr<<"scaling "<<ts<<std::endl;

    std::ofstream lowr, highr, ptmd, ptwte, ptreplica; 
    cns=0.0; 
    
    ptmd.open("ptmd.dat");  ptwte.open("ptwte.dat"); ptreplica.open("ptreplica.dat"); 
    double minv, maxv;

//    std::cerr<<"TESTING GRADIENT\n";
//       f.set_vars(pos[0]);
//       f.get_value(nrg[0]);
//       f.get_gradient(grad[0]);
//   for (unsigned long i=0; i<nv; ++i)  
//    {
//      
//      pos[1][i]=pos[0][i]+1e-4; 
//       f.set_vars(pos[1]);
//       f.get_value(nrg[1]);
//      std::cerr<<"coord: "<<i<<" analytical: "<<grad[0][i]<<" numeric: "<<(nrg[1]-nrg[0])/1e-4<<"\n";
//      
//      pos[1][i]=pos[0][i];
//    }

    double vx=0.0, vy=0.0; 
    for (unsigned long ir=0; ir<nr; ++ir) 
    { 
       c2=sqrt(temp[ir]);  
       for (unsigned long i=0; i<nv; ++i) vel[ir][i]=c2*rngen(); 
       vx=vy=0.0; for (unsigned long i=0; i<nv; i+=2) {vx+=vel[ir][i]; vy+=vel[ir][i+1]; }
       vx*=2.0/nv;  vy*=2.0/nv; for (unsigned long i=0; i<nv; i+=2) {vel[ir][i]-=vx; vel[ir][i+1]-=vy; }

       f.set_vars(pos[ir]);
       f.get_value(nrg[ir]);
       f.get_gradient(grad[ir]);
       totv[ir]=nrg[ir]; totv2[ir]=nrg[ir]*nrg[ir]; totn=1.0;
       ireplica[ir]=ir;

    }
    minv=maxv=nrg[0];
    nrgmin=nrg[0]; pmin=pos[0]; wtedv=0.0;
    for (unsigned long is=0; is< po.steps; ++is)
    {    
        std::cerr<<"step: "<<std::setw(6)<<is<<" t0: "<<std::setw(8)<<temp[0]<< "  ";
        ptmd<<is<<" ";     ptreplica<<is<<" ";     
        totn=totf*totn+1.0;
        for (unsigned long ir=0; ir<nr; ++ir)
        {
            // does one MD step ( no trotter splitting, since we don't really care that much about the MD )
            c1=std::exp(-0.5*po.dt/po.tau), c2=sqrt(temp[ir]*(1.0-c1*c1)); 
            c3=ts, c4=sqrt(temp[ir]*(1.0-c3*c3)); 
            
            kin[ir]=0.0;  for (unsigned long i=0; i<nv; ++i) kin[ir]+=vel[ir][i]*vel[ir][i]; kin[ir]*=0.5;
            cns[ir]+=kin[ir];            
            
            //BDP thermostat...
            vdir=vel[ir]*sqrt(0.5/kin[ir]); 
            vel[ir]*=c1;  for (unsigned long i=0; i<nv; ++i) vel[ir][i]+=c2*grngen();  //white noise step
            kin[ir]=0.0;  for (unsigned long i=0; i<nv; ++i) kin[ir]+=vel[ir][i]*vel[ir][i]; // project on the old velocity direction
            vel[ir]=vdir*sqrt(kin[ir]); kin[ir]*=0.5;
            cns[ir]-=kin[ir];            
            //AND ALSO A BIT OF WN LANGEVIN...
            cns[ir]+=kin[ir];
            vel[ir]*=c3;  for (unsigned long i=0; i<nv; ++i) vel[ir][i]+=c4*grngen();  //white noise step            
            vx=vy=0.0; for (unsigned long i=0; i<nv; i+=2) {vx+=vel[ir][i]; vy+=vel[ir][i+1]; }
            vx*=2.0/nv;  vy*=2.0/nv; for (unsigned long i=0; i<nv; i+=2) {vel[ir][i]-=vx; vel[ir][i+1]-=vy; } // followed by com removal
            cns[ir]-=kin[ir];

            vel[ir]+=grad[ir]*(1.0+wtedv)*(-po.dt*0.5);
            pos[ir]+=vel[ir]*po.dt;
            f.set_vars(pos[ir]);
            f.get_value(nrg[ir]);
            f.get_gradient(grad[ir]);

            //wte evaluation
            wtev=wtedv=0.0; double dwte;
            //also evaluate wte[ir](nrg[ir-1])
            if (ir>1) for (unsigned long i=0; i<wte_heights[ir].size(); ++i) 
            {  dwte=(nrg[ir-1]-wte_history[ir][i])/wte_widths[ir][i]; dwte=std::exp(-dwte*dwte)*wte_heights[ir][i];  wtev+=dwte;  }            
            if (ir>1) wte_up[ir-1]=wtev;
            //also evaluate wte[ir-1](nrg[ir])
            wtev=0.0;
            if (ir>1) for (unsigned long i=0; i<wte_heights[ir-1].size(); ++i) 
            {  dwte=(nrg[ir]-wte_history[ir-1][i])/wte_widths[ir-1][i]; dwte=std::exp(-dwte*dwte)*wte_heights[ir-1][i];  wtev+=dwte;  }            
            wte_dw[ir]=wtev;

            wtev=0.0;
            for (unsigned long i=0; i<wte_heights[ir].size(); ++i) 
            {   
               dwte=(nrg[ir]-wte_history[ir][i])/wte_widths[ir][i]; dwte=std::exp(-dwte*dwte)*wte_heights[ir][i]; 
               wtev+=dwte; wtedv+=dwte*(-2.0)*(nrg[ir]-wte_history[ir][i])/(wte_widths[ir][i]*wte_widths[ir][i]);
               wte_heights[ir][i]*=ts;  //must scale down hills, because we are also annealing!
            }            
            wte[ir]=wtev;

            cns[ir]-=wtev*(1.0-ts);
            

            if (nrg[ir]<nrgmin) 
            { 
               nrgmin=nrg[ir]; pmin=pos[ir]; 
               std::cerr<<"saving new min\n";
               lowr.open("min.dat"); lowr.precision(12);
               for (unsigned long i=0; i<nv; ++i) {lowr<<std::setw(10)<<pmin[i]<<"  "; if (i%2==1) lowr<<std::endl;}
               lowr.close();
            }
            // accumulates running averages of the energy
            totv[ir]=totf*totv[ir]+nrg[ir]; totv2[ir]=totf*totv2[ir]+nrg[ir]*nrg[ir];

            vel[ir]+=grad[ir]*(1.0+wtedv)*(-po.dt*0.5);

            kin[ir]=0.0;  for (unsigned long i=0; i<nv; ++i) kin[ir]+=vel[ir][i]*vel[ir][i]; kin[ir]*=0.5;
            cns[ir]+=kin[ir];            

            //BDP thermostat...
            vdir=vel[ir]*sqrt(0.5/kin[ir]); 
            vel[ir]*=c1;  for (unsigned long i=0; i<nv; ++i) vel[ir][i]+=c2*grngen();  //white noise step
            kin[ir]=0.0;  for (unsigned long i=0; i<nv; ++i) kin[ir]+=vel[ir][i]*vel[ir][i]; // project on the old velocity direction
            vel[ir]=vdir*sqrt(kin[ir]); kin[ir]*=0.5;
            cns[ir]-=kin[ir];            
            
            std::cerr<<"| "<<std::setw(10)<<nrg[ir]/nv*2<<" "<<std::setw(10)<<(kin[ir]+nrg[ir]+cns[ir])/nv*2<< " |";

            vreplica[ireplica[ir]]=nrg[ir];
            ptmd<<nrg[ir]/nv*2<<" "<<kin[ir]/nv*2<<" "<<wtev<<" "<<(nrg[ir]+kin[ir]+cns[ir]+wtev)/nv*2<<" ";            
        }
        std::cerr<<"\n";

        for (unsigned long ir=0; ir<nr; ++ir) 
            ptreplica<<vreplica[ir]/nv*2<<" ";
                
        // does PT swaps
//        if (rngen()<po.dt/(2.0*po.tau*(po.temp_factor-1.0)) ) 
        {
           std::cerr<<"attempting swap\n";
           for (unsigned long ir=0; ir<nr-1; ++ir) 
           { 

               // this is not "thermodynamic" PT, since we make the energies size-intensive.
//               if (rngen()<exp(-(nrg[ir]-nrg[ir-1])/(nv*0.5*(temp[ir-1]+temp[ir]) ))   )
               double deltah=(1.0/temp[ir+1]-1.0/temp[ir])*(nrg[ir+1]-nrg[ir]) +
                              (wte[ir]-wte_dw[ir+1])/nrg[ir] +
                              (wte[ir+1]-wte_up[ir])/nrg[ir+1];
               std::cerr<<ir<<" "<<deltah<<std::endl;
               if (rngen()<std::exp(deltah)   )  
               { 
                  std::cerr<<"Swapping replicas "<<ir<<" and " <<ir+1<<std::endl;
                  pswp=pos[ir];  pos[ir]=pos[ir+1];   pos[ir+1]=pswp;
                  pswp=grad[ir]; grad[ir]=grad[ir+1]; grad[ir+1]=pswp;        
                  swp=nrg[ir];   nrg[ir]=nrg[ir+1];   nrg[ir+1]=swp;                         
                  cns[ir]+=nrg[ir+1]-nrg[ir];   cns[ir+1]+=nrg[ir]-nrg[ir+1];

                  pswp=vel[ir]; vel[ir]=vel[ir+1]*sqrt(temp[ir]/temp[ir+1]); vel[ir+1]=pswp*sqrt(temp[ir+1]/temp[ir]);
                  swp=kin[ir]; kin[ir]=kin[ir+1]*(temp[ir]/temp[ir+1]); kin[ir+1]=swp*(temp[ir+1]/temp[ir]); 
                  cns[ir]+=kin[ir+1]*(temp[ir]/temp[ir+1])-kin[ir];   
                  cns[ir+1]+=kin[ir]*(temp[ir+1]/temp[ir])-kin[ir+1];
//                  pswp=vel[ir];   vel[ir+1]*=sqrt(temp[ir]/temp[ir+1]); vel[ir]=vel[ir+1];   pswp*=sqrt(temp[ir+1]/temp[ir]); vel[ir+1]=pswp;
//                  cns[ir]-=kin[ir+1]*temp[ir]/temp[ir+1]-kin[ir];  cns[ir+1]-=kin[ir]*temp[ir+1]/temp[ir]-kin[ir+1];

                  swp=ireplica[ir]; ireplica[ir]=ireplica[ir+1]; ireplica[ir+1]=swp;
               }
            }
        }        
      

        // ADD NEW HILLS 
        if (false && is>po.tau_avg/po.dt && rngen()<0.5*po.dt/po.tau)  // uses the same time constant as the thermostat, for a start
        {
            for (unsigned long ir=0; ir<nr; ++ir) 
            {
               if (minv>nrg[ir]) minv=nrg[ir];               if (maxv<nrg[ir]) maxv=nrg[ir];
               std::cerr<<"ADDING NEW HILL, nrg: "<<nrg[ir]<<" temp "<<temp[ir]<<" gamma "<<wte_gamma[ir]<<"\n"; 
               wte_history[ir].push_back(nrg[ir]); 
               wte_heights[ir].push_back(temp[ir]*0.5*std::exp(-wte[ir]/(temp[ir]*(wte_gamma[ir]-1))));
               wte_widths[ir].push_back(0.5*sqrt(totv2[ir]/totn-(totv[ir]/totn)*(totv[ir]/totn)) );
               std::cerr<<" new height: "<<wte_heights[ir][wte_heights[ir].size()-1]
                        << "  new width: "<<wte_widths[ir][wte_heights[ir].size()-1]<<"\n";
               wtev=wte[ir];
               wtev+=wte_heights[ir][wte_heights[ir].size()-1];
               cns[ir]-=wte_heights[ir][wte_heights[ir].size()-1];  // updates conserved quantity to keep it in order
               //TODO also updates gamma to try to attain the optimal variance 
               unsigned long ih=wte_heights[ir].size()-1;
               ptwte<<ih<<"  "<<wte_history[ir][ih]<<"  "<<wte_heights[ir][ih]<<"  "<<wte_widths[ir][ih]<<"  ";
            }
            ptwte<<std::endl;

            //prints out the bias
            lowr.open("bias.dat"); lowr.precision(7); double dwte;
            for (double x=minv-0.5*(maxv-minv); x<maxv+0.5*(maxv-minv); x+=1e-3*(maxv-minv) )
            {
                lowr<<std::setw(10)<<x<<" ";
                for (unsigned long ir=0; ir<nr; ++ir) 
                { wtev=0.0;
                  for (unsigned long i=0; i<wte_heights[ir].size(); ++i) 
                 {  dwte=(x-wte_history[ir][i])/wte_widths[ir][i]; dwte=std::exp(-dwte*dwte)*wte_heights[ir][i];  wtev+=dwte;  }  
                 lowr<<std::setw(10)<<wtev<<" ";
                 }
                  lowr<<std::endl;
              }
            lowr.close();
         }

        ptmd<<std::endl;        ptreplica<<std::endl;
        if (is%100==0) {
         std::cerr<<"printing snapshots\n";
         for (unsigned long ir=0; ir<nr; ir++)
         {
            lowr.open((std::string("replica")+int2str(ir)+std::string(".dat")).c_str()); lowr.precision(12);
            for (unsigned long i=0; i<nv; ++i) {lowr<<std::setw(10)<<pos[ir][i]<<"  "; if (i%2==1) lowr<<std::endl;}
            lowr.close();            
         }
//         highr.open("highr.dat"); highr.precision(12);
//         for (unsigned long i=0; i<nv; ++i) {highr<<std::setw(10)<<pos[nr-1][i]<<"  "; if (i%2==1) highr<<std::endl;}
//         highr.close();        
        }
        
        temp*=ts;
    }
    rpos=pmin;
    std::cerr<<"done & returning\n";    
    f.set_vars(pmin);   //reset status to minimum position
    f.get_value(rvalue);
    std::cerr<<"done & returning (again)\n";        
}

template<class FCLASS>
void para_temp (
     FCLASS& f,
     const std::valarray<double> & init_vars,
     std::valarray<double>& rpos, double& rvalue,
     const ParaOptions po=ParaOptions()
 )
{
    para_temp<FCLASS,StdRndUniform >(f,init_vars,rpos,rvalue,po,StdRndUniform());
}



/***********************************************************************
             HESSIAN-AWARE MINIMIZER
 ***********************************************************************/
/*class HessOpts{
public:
    double tol;          unsigned long maxiter;
    HessOpts(): tol(1e-5), maxiter(100) {}
};

template<class FCLASS>
void min_hessian (
       FCLASS& f, 
       const std::valarray<double> & initial_pos,
       std::valarray<double>& rpos, double& rvalue,
       const HessOpts& op
      )
{
    std::valarray<double> pos(initial_pos), npos(pos);
   
    unsigned long sz=pos.size(); double fo, fn, lsstep, gamma, gg;
    LineSearchOpts ls(op.linesearch);
    ls.dir.resize(sz); std::valarray<double>g(ls.dir), og(g);
   
    f.set_vars(pos); 
    f.get_gradient(ls.dir); ls.dir*=-1.0;
    og=ls.dir;
    for (unsigned long istep=0; istep<op.maxiter; ++istep)
    {
        min_linesearch(f,pos,npos,fn,lsstep,ls);
        f.set_vars(npos); 
        f.get_gradient(g); g*=-1.0;
       
        gamma=gg=0.0; for (unsigned long i=0; i<sz; ++i) { gg+=og[i]*og[i]; gamma+=(g[i]-og[i])*g[i]; } gamma/=gg;
       
        og=g;
        ls.dir*=gamma; ls.dir+=g;
       
        ls.istep=lsstep*0.5;
        pos=npos;
        std::cerr<<"Step size: "<<ls.istep<<"\n";
        std::cerr<<"Step dir: "<<ls.dir<<"\n";
        std::cerr<<"Current value is "<<fn<<"\n";
    }
    rpos=pos; rvalue=fn;
}*/

}; //ends namespace toolbox
#endif //ends ifdef __MINSEARCH_H
