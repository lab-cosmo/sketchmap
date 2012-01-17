#ifndef __TOOLS_AUTOCORR_H
#define __TOOLS_AUTOCORR_H

#include "tbdefs.hpp"
#include <vector>

#define TB_FFTAC_BUF 1024
#define TB_FFTAC_SCL 1.25
#ifdef TB_FFTAC
#include "fftw3.h"
#endif

/*************************************************
 Class to compute autocorrelation function on the 
 fly. Can read in single data points or valarrays.
 Compute on the fly several quantities of interest.
 *************************************************/
namespace toolbox {
template <class U> class AutoCorrelation;
template <class U> class ACOptions;

template <class U> class ACOptions<AutoCorrelation<U> > {
public:
    U timestep; 
    U xmean, xsigma, xtau, xtau2;
    bool f_exact_mean, f_exact_sigma, f_exact_tau, f_exact_tau2;
    ACOptions(double ntimestep=1., double nxmean=0., double nxsigma=1., 
            bool nf_exact_mean=false, bool nf_exact_sigma=false, 
            double nxtau=1., double nxtau2=1., 
            bool nf_exact_tau=false, bool nf_exact_tau2=false
             ) :
            timestep(ntimestep), xmean(nxmean), xsigma(nxsigma),
            xtau(nxtau), xtau2(nxtau2),
            f_exact_mean(nf_exact_mean), f_exact_sigma(nf_exact_sigma),
            f_exact_tau(nf_exact_tau), f_exact_tau2(nf_exact_tau2){}
};
        
template <class U> class AutoCorrelation {
private:
    std::vector<std::valarray<U> > history;
#ifdef TB_FFTAC
    std::vector<std::vector<U> > autodata;
#else
    std::vector<std::valarray<U> > autodata;
#endif
    std::vector<std::valarray<U> > headdata;
    
    std::vector<U> seriestot;
//    std::vector<U> seriestot2;
    std::vector<double> serieslengths;
    
    bool f_fresh, f_acfresh;
    unsigned long n_collapsed;
    unsigned long iseries;
    unsigned long p_ndata;
    unsigned long p_nseries;
    unsigned long p_ncorr;
    std::valarray<U> p_acf;
    double p_mean, p_mean2;
    U tot, tot2;
    
    void mkmeans();
    void mkacf();
    ACOptions<AutoCorrelation<U> > opts;

public:
    void get_options(ACOptions<AutoCorrelation<U> >& rop)
    { rop=opts; }
    
    void set_options(const ACOptions<AutoCorrelation<U> >& rop)
    { opts=rop; f_fresh=f_acfresh=false; }
    
    
    AutoCorrelation(const unsigned long& ncorr=0, 
                ACOptions<AutoCorrelation<U> > nopts=ACOptions<AutoCorrelation<U> >()) : 
            f_fresh(false), n_collapsed(0), p_ndata(0), p_nseries(1), 
            p_ncorr(ncorr), tot(0.), tot2(0.),
            opts(nopts)
    { reset(ncorr); }
    
    //copy & assign...
    AutoCorrelation(const AutoCorrelation<U>& ac);
    AutoCorrelation<U>& operator=(const AutoCorrelation<U>& ac);

    void reset(const unsigned long& ncorr=0, const unsigned long& nsets=1);
    
    double actime() 
    { 
    //we compute AC time as the integral of the AC function
    //we do so by calling the predefined function, even if this is
    //not optimally efficient, we keep the code cleaner
        double actime=0;
        for (unsigned long h=0; h<p_ncorr; ++h) actime+=(*this)[h];
        actime-=(*this)[0]*0.5;
        return actime*opts.timestep; 
    }
    
    double actime2() 
    { 
    //we compute the square time as the integral of the AC function
    //we do so by calling the predefined function, even if this is
    //not optimally efficient, we keep the code cleaner
        double actime=0, ach;
        for (unsigned long h=0; h<p_ncorr; ++h) { ach=(*this)[h]; actime+=ach*ach; }
        actime-=(*this)[0]*0.5;
        return actime*opts.timestep; 
    }

    double mean() { mkmeans(); return p_mean; }
    double mean2() { mkmeans(); return p_mean2; }
    //ok, beware. this does not return sqrt(<a2>-<a>2), 
    //instead it returns acf(0), i.e. taking into account the 'exact mean'
    //and 'exact sigma2' options 
    double sigma() 
    { 
        mkmeans(); 
        if (opts.f_exact_sigma) return opts.xsigma;
        if (opts.f_exact_mean) return sqrt(p_mean2+opts.xmean*(opts.xmean-2.*p_mean));
        return sqrt(p_mean2-p_mean*p_mean); 
    }
    unsigned long samples() const { return p_ndata; }
    unsigned long maxlag() const { return p_ncorr; }
    
    void add(const U& nel);
    void series_next();
    void series_prev();
    void series_collapse(); //collapse all the data onto dataset 0, without affecting the overall acf

//syntactic sugar to insert a new element into the series, to move to next or previous series
    inline void operator << (const U& nel) { add(nel); }
    inline void operator << (const std::valarray<U>& nseries) { for (unsigned long i=0; i<nseries.size(); ++i) add(nseries[i]); }
    void operator ++() { series_next(); }
    void operator --() { series_prev(); }    

    void getpoint(const unsigned long& h, double& t, U& v);
    void fullanalysis(std::valarray<double>& t, std::valarray<U>& v, std::valarray<U>& b, std::valarray<U>& ev, std::valarray<U>& eb);

    U operator [] (const long& i) {double t,v; getpoint(i,t,v); return v;}
};


/***************************************************************
 *          AutoCorrelation member functions                  *
 ***************************************************************/
template <class U>
void AutoCorrelation<U>::reset(const unsigned long& ncorr, const unsigned long& nsets)
{
    f_fresh=f_acfresh=false;
    n_collapsed=0;
    p_ndata=0;
    p_ncorr=ncorr;
    p_nseries=nsets;
    iseries=0;

    autodata.resize(nsets); 
    headdata.resize(nsets); 
    history.resize(nsets); 
//    seriestot2.resize(nsets);
    seriestot.resize(nsets);
    serieslengths.resize(nsets);
    
    for (unsigned long k=0; k<nsets;++k)
    {        
        autodata[k].resize(ncorr); 
        for (unsigned long i=0; i<autodata[k].size(); ++i)
            autodata[k][i]=0.;
        
        headdata[k].resize(ncorr); headdata[k]=0.;
        history[k].resize(ncorr); history[k]=0.;
        seriestot[k]=serieslengths[k]=0;
        //seriestot2[k]=0
    }
    
    tot=tot2=0.;
    p_acf.resize(ncorr);
}



//copy and assign
template <class U>
AutoCorrelation<U>::AutoCorrelation(const AutoCorrelation<U>& ac)
{
    reset(ac.p_ncorr, ac.p_nseries);
    
    history=ac.history;
    autodata=ac.autodata;
    headdata=ac.headdata;
//    seriestot2=ac.seriestot2;
    seriestot=ac.seriestot;
    serieslengths=ac.serieslengths;
    
    tot=ac.tot; tot2=ac.tot2;
    
    f_fresh=ac.f_fresh; f_acfresh=ac.f_acfresh;
    n_collapsed=0;
    p_ndata=ac.p_ndata; p_nseries=ac.p_nseries; iseries=ac.iseries;
    p_mean=ac.p_mean; p_mean2=ac.p_mean2; p_acf=ac.p_acf;
    opts=ac.opts;
}

template <class U>
AutoCorrelation<U>& AutoCorrelation<U>::operator=(const AutoCorrelation<U>& ac)
{
    if (this!=&ac)
    {
        reset(ac.p_ncorr,ac.p_nseries);
        history=ac.history;
        autodata=ac.autodata;
        headdata=ac.headdata;
//        seriestot2=ac.seriestot2;
        seriestot=ac.seriestot;
        serieslengths=ac.serieslengths;
    
        tot=ac.tot; tot2=ac.tot2;
    
        f_fresh=ac.f_fresh; f_acfresh=ac.f_acfresh;
        n_collapsed=0;
        p_ndata=ac.p_ndata; p_nseries=ac.p_nseries; iseries=ac.iseries;
        p_mean=ac.p_mean; p_mean2=ac.p_mean2; p_acf=ac.p_acf;
        opts=ac.opts;
    }
    return *this;
}

/*************************************************************************
 we compute the AC function defined as 
 C_n(k)=1/(n-k)*1/s^2 \sum_i=k^n-1 (A_i -mu)*(A_i-k -mu)
 mu and sigma can be given as known parameters, or can be computed from
 the data set. in order to compute everything on the fly, without storing
 the complete time series, we rearrange the expression: let <a> be the 
 mean computed from the time series
 then 
 C_n(k)=1/s^2 (1/(n-k)\sum_i=k^n-1 A_i*A_i-k + mu^2 -2<a>*mu*n/(n-k) 
        + mu/(n-k) sum_i=0^k-1[A_i+A(n-i-1)] 
*************************************************************************/
template <class U>
void AutoCorrelation<U>::mkmeans()
{
    // if computed data are "fresh", we do nothing!
    if (!f_fresh)
    {
        if (p_ndata>0)
        {
            p_mean=tot/p_ndata;
            p_mean2=tot2/p_ndata;
        }
        else
        {
            p_mean=0; p_mean2=0;
        }
        f_fresh=true;
    }
}

#ifdef TB_FFTAC
void do_acfft(std::valarray<double>& sacf)
{
    //ok. this is the time to do FFT computation of the AC. 
    unsigned long ftsize=sacf.size()/2+1;
    fftw_complex * ftacf=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*ftsize);
    fftw_plan fwplan=fftw_plan_dft_r2c_1d(sacf.size(),&sacf[0],ftacf,FFTW_ESTIMATE);
    fftw_plan bwplan=fftw_plan_dft_c2r_1d(sacf.size(), ftacf, &sacf[0],FFTW_ESTIMATE);

    fftw_execute(fwplan);
    for(unsigned long i=0; i<ftsize; ++i)
        {  ftacf[i][0]=ftacf[i][0]*ftacf[i][0]+ftacf[i][1]*ftacf[i][1];  ftacf[i][1]=0;  }
    fftw_execute(bwplan);
    sacf*=(1./sacf.size());
    fftw_destroy_plan(fwplan);
    fftw_destroy_plan(bwplan);
    fftw_free(ftacf);
}
#endif

template <class U>
void AutoCorrelation<U>::mkacf()
{
    HRTimer timing;
    if (!f_acfresh)
    {
        mkmeans();
        double mu=(opts.f_exact_mean ? opts.xmean : p_mean);
        double sigma2;
        p_acf=0;
        
        double v;
        
        //here we compute acf of different series, and average them according to the 
            //length of the series themselves
        double nsmp=0;
        for (unsigned long s=0; s<p_nseries; ++s)
        if (serieslengths[s]>p_ncorr) //only series contributing meaningful data are considered
        {
            std::cerr<<"Working through series "<<s<<"\n";
            //in order to get proper renormalization of accumulated series, the 
            //series zero (the collapsed one) must be treated slightly differently,
            //as all the data sets contain the accumulated acf from the other series.
            
            std::valarray<double> sacf;
#ifdef TB_FFTAC
            if (s>0 || n_collapsed==0) 
            {
                sacf.resize(serieslengths[s]+p_ncorr);  sacf=0.;
                for(unsigned long i=0; i<serieslengths[s]; ++i) sacf[i]=autodata[s][i];
                do_acfft(sacf);
            }
            else
            {
                sacf.resize(p_ncorr);  
                for(unsigned long i=0; i<p_ncorr; ++i) sacf[i]=autodata[s][i];
            }
#else
            sacf.resize(p_ncorr);
            for(unsigned long i=0; i<p_ncorr; ++i) sacf[i]=autodata[s][i];
#endif

#ifdef TB_FFTAC
            //if collapsed series, autodata[0] already contains the autocorrelation
            
            
#endif
            std::cerr<<"ACF computed\n";
            timing.start();
            //treats head and tail effects...
            v=0;
            if (s==0 && n_collapsed>0) for (unsigned long h=0; h<p_ncorr; ++h)
            {
                p_acf[h]+=(sacf[h]+mu*(v-2.*seriestot[s]));
                v+=headdata[s][h]+history[s][h];
            }
            else for (unsigned long h=0; h<p_ncorr; ++h)
            {
                p_acf[h]+=(sacf[h]+mu*(v-2.*seriestot[s]));
                v+=headdata[s][h]+history[s][(((unsigned long) serieslengths[s])-1-h)%p_ncorr];
            }
            timing.stop();
            std::cerr<<"Tails computed, time "<<timing<<" \n";
            nsmp+=serieslengths[s];
        }
        for (unsigned long h=0; h<p_ncorr; ++h)
            p_acf[h]*=1./(nsmp-h*(p_nseries+n_collapsed));
        p_acf+=mu*mu;
            /*!
            double nsmp=0;
            
            //here we compute acf of different series, and average them according to the 
            //length of the series themselves
            for (unsigned long s=0; s<p_nseries; ++s)
                if (serieslengths[s]>h) //only series contributing meaningful data are considered
            {
                //in order to get proper renormalization of accumulated series, the 
                //series zero (the collapsed one) must be treated slightly differently,
                //as all the data sets contain the accumulated acf from the other series.
                
                v=0; //compute the tail correction
                for (unsigned long i=0; i<h; ++i) 
                    v+=headdata[s][i]+
                        history[s][
                            s==0 && n_collapsed>0?
                                i:
                                ((unsigned long)serieslengths[s]-i-1)%p_ncorr
                        ];
                std::<<s<<","<<h<<","<<v<<","<<serieslengths[s]<<","<<((autodata[s][h]+mu*(v-2.*seriestot[s]))/
                        (serieslengths[s]-h)+mu*mu)/(sigma()*sigma())<<"\n";
                p_acf[h]+=(autodata[s][h]+mu*(v-2.*seriestot[s]));///(1.-h*1./serieslengths[s]);
                
                nsmp+=serieslengths[s];
            }
            p_acf[h]*=1./(nsmp-h*(p_nseries+n_collapsed));
        }
        p_acf+=mu*mu; 
            */
        //normalization of the ACF. if we gave exact mean but no sigma, we have to compute 
        //acf(0) in order to get the normalization
        
        if (opts.f_exact_sigma) sigma2=opts.xsigma*opts.xsigma;
        else if (! opts.f_exact_mean) sigma2=p_mean2-p_mean*p_mean;
        else sigma2=p_acf[0];
        
        p_acf *= 1./sigma2;
        f_acfresh=true;
    }
}

//see notes of the compute() function!
//averages are computed over the different series
template <class U>void AutoCorrelation<U>::getpoint(const unsigned long& h, double& t, U& v)
{ 
    mkacf();
    t=h*opts.timestep; v=p_acf[h];
}

template <class U>
void AutoCorrelation<U>::fullanalysis(std::valarray<double>& t, std::valarray<U>& v, std::valarray<U>& b, std::valarray<U>& ev, std::valarray<U>& eb)
{
    
    t.resize(p_ncorr); v.resize(p_ncorr); b.resize(p_ncorr); ev.resize(p_ncorr); eb.resize(p_ncorr);
        
    /**********************************************************************
    everything is very easy, here! all the nasty collapsed series things
    have been dealt with in the autocorrelation computation, so we just need 
    to integrate the ac function.
    we compute the "block error" b[m]=sum_k=0^m(2-d_k0) (1-k/m) acf(k)
    ***********************************************************************/
    for (unsigned long h=0; h<p_ncorr; ++h) getpoint(h,t[h],v[h]);
    b[0]=eb[0]=0.;
    b[1]=v[0]*0.5; //acf[0]
    double vtimesm=v[1];
    for (unsigned long m=2; m<p_ncorr; ++m) 
    {   
        b[m]=b[m-1]+vtimesm*(1./(m-1)-1./(m));
        vtimesm+=v[m]*m;
    }
    b*=(2*opts.timestep);
    /*for (unsigned long m=1; m<p_ncorr; ++m) 
    {   
        b[m]=0;
        for (unsigned long h=1; h<m; ++h) b[m]+=(1-h*1./m)*v[h];
        b[m]=(1.+b[m]+b[m])*opts.timestep;
    }*/
    
    /**********************************************************************
    here we compute the errors on ACF, assuming that the AC is a gaussian
    distributed random variable. see Zwanzig and Ailawadi¸ PR 182 (1969)
    ***********************************************************************/
    //to keep code short, we compute actime with the routine, even if in 
    //practice in this way we compute another time the ac function
    double act=(opts.f_exact_tau2?opts.xtau2:actime2());
    act*=2./(p_ndata*opts.timestep); act=(act>0?sqrt(act):0);

    for (unsigned long h=0; h<p_ncorr; ++h) ev[h]=act*(1-v[h]);
    
    //the errors on blocking averages can be computed in a similar way 
    //!BEWARE THIS FORMULA NEED CHECK!!!!!!!!!!!!!!!!
    for (unsigned long h=0; h<p_ncorr; ++h) eb[h]=act*fabs(b[h]-h*opts.timestep/2);
}


/*************************************************************************
 a couple of words to explain how this is done: we store the data in a 
 cyclic buffer, and we add the contributions to the AC at different time
 lags on the fly.  we also need to store for corrections the first p_ncorr 
*************************************************************************/
template <class U> 
void AutoCorrelation<U>::add(const U& nel)
{
    if (n_collapsed>0 && iseries==0) ERROR("Cannot add data to series zero after it has been collapsed"); 
    
    history[iseries][(unsigned long)serieslengths[iseries]%p_ncorr]=nel;
    seriestot[iseries]+=nel;
//    seriestot2[iseries]+=nel*nel;
    tot+=nel;
    tot2+=nel*nel;
    
    unsigned long id=(unsigned long) serieslengths[iseries];
    
#ifdef TB_FFTAC
    //if we compute by FFT, autodata holds the full data series
    if (id>=autodata[iseries].size()) 
        autodata[iseries].resize(autodata[iseries].size()*TB_FFTAC_SCL+TB_FFTAC_BUF);
    autodata[iseries][id]=nel;
#else
    for (unsigned long k=0; k<(id<p_ncorr?id+1:p_ncorr); ++k)
        autodata[iseries][k]+=nel*history[iseries][(id-k)%p_ncorr];
#endif
    if (id<p_ncorr)
    {
        //for (unsigned long k=0; k<=id; ++k)
        //    autodata[iseries][k]+=nel*history[iseries][(id-k)];
        headdata[iseries][id]=nel;
    }
        
    ++p_ndata; ++serieslengths[iseries];
    f_fresh=false;
}

/***********************************************************************************
 this routine consolidate all the series available so far into a single series, 
 retaining the same average autocorrelation which would have been obtained with
 all the series.  since we also consider the multiple tail and head contributions,
 series 0 cannot be touched anymore.
************************************************************************************/
template <class U> 
void AutoCorrelation<U>::series_collapse()
{
    if (n_collapsed==0 && p_nseries>1)
    { 
        //basically now we want history[0] to hold an accumulation of tail elements, 
        //so that history[0][0] holds the last elements and so on. the first time, we need
        //to rearrange the elements in history[0] to comply with this bookkeeping
        //this could be done with just one 'working slot' extra storage, but
        //I do it with a full working array cuz I'm lazy
        std::valarray<double> swp(p_ncorr);
        for(unsigned long h=0; h<p_ncorr; ++h)
            swp[h]=history[0][((unsigned long)serieslengths[0]-h-1)%p_ncorr];
        history[0]=swp;
#ifdef TB_FFTAC
        std::valarray<double> sac(serieslengths[0]+p_ncorr);
        sac=0.;
        for(unsigned long i=0; i<serieslengths[0]; ++i)
            sac[i]=autodata[0][i];
        do_acfft(sac);
        autodata[0].resize(p_ncorr);
        for(unsigned long h=0; h<p_ncorr; ++h) 
            autodata[0][h]=sac[h];
#endif
    }
    
    for (unsigned long s=1;s<p_nseries;++s)
    {
#ifdef TB_FFTAC
        //since data of series will be deleted anyway, we put acf inside it and
        //collapse the data onto series 0
        std::valarray<double> sac(serieslengths[s]+p_ncorr);
        sac=0.;
        for(unsigned long i=0; i<serieslengths[s]; ++i)
            sac[i]=autodata[s][i];
        do_acfft(sac);
        autodata[s].resize(p_ncorr);
        for(unsigned long h=0; h<p_ncorr; ++h) 
            autodata[s][h]=sac[h];
#endif
        //first, we make collapse autodata
        for(unsigned long h=0; h<p_ncorr; ++h)
        {
            autodata[0][h]+=autodata[s][h];
            headdata[0][h]+=headdata[s][h];
            //slightly trickier is collapsing the tails
            history[0][h]+=history[s][((unsigned long)serieslengths[s]-h-1)%p_ncorr];
        }
        serieslengths[0]+=serieslengths[s];
        seriestot[0]+=seriestot[s];
        n_collapsed++;
    }
    
    
    //mkacf(); //updates headdata
    
    //shrink all the arrays
    p_nseries=1;
    history.resize(1); 
    autodata.resize(1); headdata.resize(1);
    serieslengths.resize(1); seriestot.resize(1); 
    iseries=0;
    f_fresh=f_acfresh=false; //things need to be recomputed!
}

template <class U> 
void AutoCorrelation<U>::series_next()
{ 
    ++iseries; 
    if (iseries>=p_nseries)  //automatically insert a new series if needed
    {
        p_nseries=iseries+1;
        history.resize(p_nseries); history[iseries].resize(p_ncorr); history[iseries]=0;
        autodata.resize(p_nseries); autodata[iseries].resize(p_ncorr); 
        for (unsigned long i=0; i<autodata[iseries].size(); ++i)
            autodata[iseries][i]=0.;
        headdata.resize(p_nseries); headdata[iseries].resize(p_ncorr); headdata[iseries]=0;
        serieslengths.resize(p_nseries); serieslengths[iseries]=0;
        seriestot.resize(p_nseries); seriestot[iseries]=0;
        //seriestot2.resize(p_nseries); seriestot2[iseries]=0;
    }
}

template <class U> 
void AutoCorrelation<U>::series_prev()
{
    iseries--; 
    if (iseries<0) iseries=p_nseries-1; //deals gracefully with overflow
}

/*************************************************************************************
  VERY SIMPLE (non-fft, non cumulative) class for cross-correlation. to be improved....
**************************************************************************************/
template <class U> class CrossCorrelation;
template <class U> class CCOptions;

template <class U> class CCOptions<CrossCorrelation<U> > {
    public:
        U timestep; 
        U xmean_a, xmean_b;
        bool f_exact_mean;
        CCOptions(double ntimestep=1., double nxmean_a=0., double nxmean_b=0.,
                  bool nf_exact_mean=false
                 ) :
                timestep(ntimestep), xmean_a(nxmean_a), xmean_b(nxmean_b),
                          f_exact_mean(nf_exact_mean) {}
};
        
template <class U> class CrossCorrelation {
    private:
        std::valarray<U> history;
        std::valarray<U> headdata;
        std::valarray<U> autodata;
        
        bool f_fresh, f_acfresh;
        unsigned long p_ndata;
        unsigned long p_ncorr;
        std::valarray<U> p_acf;
        double p_mean_a, p_mean_a2, p_mean_b, p_mean_b2, p_mean_ab;
        U tot_a, tot_b, tot_a2, tot_b2, tot_ab;
    
        void mkmeans();
        void mkacf();
        CCOptions<CrossCorrelation<U> > opts;

    public:
        void get_options(CCOptions<CrossCorrelation<U> >& rop)
        { rop=opts; }
    
        void set_options(const CCOptions<CrossCorrelation<U> >& rop)
        { opts=rop; f_fresh=f_acfresh=false; }
    
    
        CrossCorrelation(const unsigned long& ncorr=0, 
                CCOptions<CrossCorrelation<U> > nopts=CCOptions<CrossCorrelation<U> >()) : 
                f_fresh(false), p_ndata(0), p_ncorr(ncorr), tot_a(0.), tot_b(0.), tot_a2(0.), tot_b2(0.),
                        tot_ab(0.), opts(nopts)
        { reset(ncorr); }
    
    //copy & assign...
        CrossCorrelation(const CrossCorrelation<U>& ac);
        CrossCorrelation<U>& operator=(const CrossCorrelation<U>& ac);

        void reset(const unsigned long& ncorr=0);
    
        double actime() 
        { 
    //we compute AC time as the integral of the AC function
    //we do so by calling the predefined function, even if this is
    //not optimally efficient, we keep the code cleaner
            double actime=0;
            for (unsigned long h=1; h<p_ncorr; ++h) actime+=(*this)[h];
            actime+=(*this)[0]*0.5;
            return actime*opts.timestep; 
        }
    
        double actime2() 
        { 
    //we compute the square time as the integral of the AC function
    //we do so by calling the predefined function, even if this is
    //not optimally efficient, we keep the code cleaner
            double actime=0, ach;
            for (unsigned long h=0; h<p_ncorr; ++h) { ach=(*this)[h]; actime+=ach*ach; }
            actime-=(*this)[0]*0.5;
            return actime*opts.timestep; 
        }

        double mean_a() { mkmeans(); return p_mean_a; }
        double mean_a2() { mkmeans(); return p_mean_a2; }
        double mean_b() { mkmeans(); return p_mean_b; }
        double mean_b2() { mkmeans(); return p_mean_b2; }
    //ok, beware. this does not return sqrt(<a2>-<a>2), 
    //instead it returns acf(0), i.e. taking into account the 'exact mean'
    //and 'exact sigma2' options 
        double sigma_a() 
        { 
            mkmeans(); 
            if (opts.f_exact_mean) return sqrt(p_mean_a2+opts.xmean_a*(opts.xmean_a-2.*p_mean_a));
            return sqrt(p_mean_a2-p_mean_a*p_mean_a); 
        }
        double sigma_b() 
        { 
            mkmeans(); 
            if (opts.f_exact_mean) return sqrt(p_mean_b2+opts.xmean_b*(opts.xmean_b-2.*p_mean_b));
            return sqrt(p_mean_b2-p_mean_b*p_mean_b); 
        }
    
        double cross_ab() 
        { 
            mkmeans(); 
            if (opts.f_exact_mean) return sqrt(p_mean_ab+opts.xmean_a*opts.xmean_b-
                        (opts.xmean_a*p_mean_b+opts.xmean_b*p_mean_a));
            return sqrt(p_mean_ab-p_mean_a*p_mean_b); 
        }
        
        unsigned long samples() const { return p_ndata; }
        unsigned long maxlag() const { return p_ncorr; }
    
        void add(const U& na, const U& nb);
                                
//syntactic sugar to insert a new element into the series, to move to next or previous series
        void getpoint(const unsigned long& h, double& t, U& v);
        void fullanalysis(std::valarray<double>& t, std::valarray<U>& v);

        U operator [] (const long& i) {double t,v; getpoint(i,t,v); return v;}
};


template <class U>
void CrossCorrelation<U>::reset(const unsigned long& ncorr)
{
    f_fresh=f_acfresh=false;
    p_ndata=0;
    p_ncorr=ncorr;

    autodata.resize(ncorr); history.resize(ncorr); headdata.resize(ncorr); 
    for (unsigned long i=0; i<ncorr; ++i)
        autodata[i]=history[i]=headdata[i]=0.;
    
    tot_a=tot_a2=tot_b=tot_b2=tot_ab=0.;
    p_acf.resize(ncorr);
}



//copy and assign
template <class U>
CrossCorrelation<U>::CrossCorrelation(const CrossCorrelation<U>& ac)
{
    reset(ac.p_ncorr);
    
    history=ac.history;
    autodata=ac.autodata;
    headdata=ac.headdata;
    
    tot_a=ac.tot_a; tot_a2=ac.tot_a2;
    tot_b=ac.tot_b; tot_b2=ac.tot_b2;
    tot_ab=ac.tot_ab;
    
    f_fresh=ac.f_fresh; f_acfresh=ac.f_acfresh;
    p_ndata=ac.p_ndata; 
    p_mean_a=ac.p_mean_a; p_mean_a2=ac.p_mean_a2; 
    p_mean_b=ac.p_mean_b; p_mean_b2=ac.p_mean_b2;
    p_mean_ab=ac.p_mean_ab;
     
    p_acf=ac.p_acf;
    opts=ac.opts;
}

template <class U>
CrossCorrelation<U>& CrossCorrelation<U>::operator=(const CrossCorrelation<U>& ac)
{
    if (this!=&ac)
    {
        reset(ac.p_ncorr);
        history=ac.history;
        autodata=ac.autodata;
        headdata=ac.headdata;
    
        tot_a=ac.tot_a; tot_a2=ac.tot_a2;
        tot_b=ac.tot_b; tot_b2=ac.tot_b2;
        tot_ab=ac.tot_ab;
    
        f_fresh=ac.f_fresh; f_acfresh=ac.f_acfresh;
        p_ndata=ac.p_ndata; 
        p_mean_a=ac.p_mean_a; p_mean_a2=ac.p_mean_a2; 
        p_mean_b=ac.p_mean_b; p_mean_b2=ac.p_mean_b2;
        p_mean_ab=ac.p_mean_ab;
     
        p_acf=ac.p_acf;
        opts=ac.opts;
    }
    return *this;
}

template <class U>
void CrossCorrelation<U>::mkmeans()
{
    // if computed data are "fresh", we do nothing!
    if (!f_fresh)
    {
        if (p_ndata>0)
        {
            p_mean_a=tot_a/p_ndata;
            p_mean_a2=tot_a2/p_ndata;
            p_mean_b=tot_b/p_ndata;
            p_mean_b2=tot_b2/p_ndata;
            p_mean_ab=tot_ab/p_ndata;
        }
        else
        {
            p_mean_a=0; p_mean_a2=0;
            p_mean_b=0; p_mean_b2=0;
            p_mean_ab=0; 
        }
        f_fresh=true;
    }
}

/*************************************************************************
 we compute the AC function defined as 
 C_n(k)=1/(n-k)*1/s^2 \sum_i=k^n-1 (A_i -mu)*(A_i-k -mu)
 mu and sigma can be given as known parameters, or can be computed from
 the data set. in order to compute everything on the fly, without storing
 the complete time series, we rearrange the expression: let <a> be the 
 mean computed from the time series
 then 
 C_n(k)=1/s^2 (1/(n-k)\sum_i=k^n-1 A_i*A_i-k + mu^2 -2<a>*mu*n/(n-k) 
        + mu/(n-k) sum_i=0^k-1[A_i+A(n-i-1)] 
*************************************************************************/
template <class U>
void CrossCorrelation<U>::mkacf()
{
    if (!f_acfresh)
    {
        mkmeans();
        double mu_a=(opts.f_exact_mean ? opts.xmean_a : p_mean_a);
        double mu_b=(opts.f_exact_mean ? opts.xmean_b : p_mean_b);
        double ccab0;
        p_acf=0;
        
        double va,vb;
        
        std::valarray<double> sacf;
        sacf.resize(p_ncorr);
        for(unsigned long i=0; i<p_ncorr; ++i) sacf[i]=autodata[i];

        for (unsigned long h=0; h<p_ncorr; ++h)
        {
            va=vb=0.; //compute the tail correction
            for (unsigned long i=0; i<h; ++i) 
            {    
                va+=headdata[i];
                vb+=history[((unsigned long)p_ndata-i-1)%p_ncorr];
            }
            p_acf[h]=(autodata[h]+mu_b*(va-tot_a)+mu_a*(vb-tot_b))/(p_ndata-h);
        }
        p_acf+=mu_a*mu_b;

        
        //how should cross-correlation be normalized? as for now, we don't normalize at all...
        p_acf *= 1./(sigma_a()*sigma_b());
       f_acfresh=true;
    }
}

//see notes of the compute() function!
//averages are computed over the different series
template <class U>void CrossCorrelation<U>::getpoint(const unsigned long& h, double& t, U& v)
{ 
    mkacf();
    t=h*opts.timestep; v=p_acf[h];
}

template <class U>
        void CrossCorrelation<U>::fullanalysis(std::valarray<double>& t, std::valarray<U>& v)
{
    
    t.resize(p_ncorr); v.resize(p_ncorr);
        
    for (unsigned long h=0; h<p_ncorr; ++h) getpoint(h,t[h],v[h]);
}


/*************************************************************************
 a couple of words to explain how this is done: we store the data in a 
 cyclic buffer, and we add the contributions to the AC at different time
 lags on the fly.  we also need to store for corrections the first p_ncorr 
*************************************************************************/
template <class U> 
void CrossCorrelation<U>::add(const U& na, const U& nb)
{
    history[(unsigned long)p_ndata%p_ncorr]=nb;
    tot_a+=na; tot_b+=nb; tot_a2+=na*na; tot_b2+=nb*nb; tot_ab+=na*nb;
    
    unsigned long id=(unsigned long) p_ndata;
    
    for (unsigned long k=0; k<(id<p_ncorr?id+1:p_ncorr); ++k)
        autodata[k]+=na*history[(id-k)%p_ncorr];
    if (id<p_ncorr)
    {
        headdata[id]=na;
    }
    
    ++p_ndata; 
    f_fresh=false;
}


} //ends namespace toolbox (autocorr class)
#endif  //ends #ifndef __TOOLS_AUTOCORR_H


