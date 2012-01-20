/* A library to compute histograms in one and many dimensions
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/

#ifndef __TOOLS_HISTOGRAM_H
#define __TOOLS_HISTOGRAM_H

#include "tbdefs.hpp"
#include <vector>
#include <limits>
namespace toolbox {
    template <class U> class Histogram;
    template <class U> class HGOptions;
    
    enum HGWindowMode {HGWDelta,  HGWBox, HGWTriangle};
    template <class U> class HGOptions<Histogram<U> > {
    public:
        HGWindowMode window;
        U window_width;
        double (*window_function)(const U& a, const U& b);
        std::valarray<U> boundaries;
        
        HGOptions(const HGWindowMode& nwindow=HGWDelta, const U nwindow_width=(U)0., const std::valarray<U>& nbnd=std::valarray<U>(0)) :
            window(nwindow), window_width(nwindow_width) { boundaries.resize(nbnd.size()); boundaries=nbnd; }
        HGOptions(const HGOptions& hgo) : window(hgo.window), window_width(hgo.window_width),
                  window_function(hgo.window_function), boundaries(hgo.boundaries) {}
        HGOptions& operator=(const HGOptions& hgo) 
        {
            if (&hgo==this) return *this;
            window=hgo.window;
            window_width=hgo.window_width;
            window_function=hgo.window_function;
            boundaries.resize(hgo.boundaries.size()); boundaries=hgo.boundaries;
            return *this;
        }
    };
    
    template <class U> class Histogram {
    private:
        std::valarray<double> bins;
        double below, above;
        double ndata;
        HGOptions<Histogram<U> > opts;

    public:
        void reset() 
        {
            bins.resize(opts.boundaries.size()-1);
            bins=0.; ndata=0; above=below=0.;
        }
        
        double samples() { return ndata; }
        
        void get_options(HGOptions<Histogram<U> >& rop)
        { rop=opts; }
    
        void set_options(const HGOptions<Histogram<U> >& rop)
        { opts=rop; reset(); }
    
        Histogram(const HGOptions<Histogram<U> >& rop=HGOptions<Histogram<U> >())
        { set_options(rop); }
        
        Histogram(const U& a, const U&b, unsigned long n, HGWindowMode hgw=HGWDelta, double ww=0.)
        {
            opts.window=hgw; opts.window_width=ww;
            opts.boundaries.resize(n+1);
            for (unsigned long i=0; i<=n; ++i) opts.boundaries[i]=a+(b-a)*(1.*i)/n;
            reset();
        }
        
        Histogram(const Histogram<U>& ac) { set_options(ac.opts); }
        
        Histogram<U>& operator=(const Histogram<U>& ac)
        {
            if (&ac==this) return *this;
            set_options(ac.opts);
        }
        
    //copy & assign...
        template<class T> Histogram(const Histogram<T>& ac);
        template<class T> Histogram<U>& operator=(const Histogram<T>& ac);

        void add(const U& nel, double weight=1.0);

    //syntactic sugar to insert a new element into the series, to move to next or previous series
        inline void operator << (const U& nel) { add(nel); }
        inline void operator << (const std::valarray<U>& nseries) { for (unsigned long i=0; i<nseries.size(); ++i) add(nseries[i]); }
    
        void get_bins(std::valarray<double>& rbins) const
        {
            rbins.resize(bins.size());
            rbins=bins*(1./ndata);
        }
        
        void get_bins(std::valarray<U>& cent, std::valarray<U>& ws, std::valarray<double>& rbins) const
        {
            rbins.resize(bins.size()); cent.resize(bins.size()); ws.resize(bins.size());
            rbins=bins*(1./ndata);
            for (unsigned long i=0; i<bins.size();++i)
            {
                cent[i]=(opts.boundaries[i]+opts.boundaries[i+1])/2.;
                ws[i]=(opts.boundaries[i+1]-opts.boundaries[i]);
            }
        }
        
        void get_outliers(double& rabove, double& rbelow) const
        {
            rabove=above/ndata; rbelow=below/ndata;
        }
    };
                
    template<class U>
    std::ostream& operator<<(std::ostream& os, const Histogram<U>& his)
    {
        std::valarray<U> wx, ww, wf;
        his.get_bins(wx,ww,wf);
        os.precision(12);
        os.setf(std::ios::scientific);
        os.width(14);
        for (unsigned int i=0; i<wx.size(); ++i)
            os<<wx[i]<<" "<<wf[i]<<" "<<ww[i]<<"\n";
        return os;
    }
    
    template<class U>
    double __hgwfbox(const U& a, const U& b)
    {
        //std::cerr<<a<<","<<b;
        if (b<=-0.5 || a >=0.5) return 0.;
        
        U ia, ib;
        ia=(a>-0.5?a:-0.5);
        ib=(b< 0.5?b:0.5);
        return ib-ia;
    }
    
    template<class U>
    double __hgwftri(const U& a, const U& b)
    {
        //std::cerr<<a<<","<<b;
        if (b<=-1. || a >=1.) return 0.;
        
        U ia, ib;
        ia=(a>-1.?a:-1.);
        ib=(b< 1.?b:1.);
        
        return (ib*(2.-fabs(ib))-ia*(2.-fabs(ia)))*0.5;
    }

    template<class U>
    void Histogram<U>::add(const U& nel, double weight)
    {
        long bs=bins.size(),ia=0, ib=bs/2, ic=bs;
        //First, finds in which bin lies the center
        while (ib>ia && ib<ic)
        {
            if (nel>opts.boundaries[ib+1])
            {   ia=ib;  ib=(ib+ic+1)/2; }
            else
            {
                if (nel>opts.boundaries[ib]) break;
                else { ic=ib; ib=(ia+ib)/2;  }
            }
        }
        //now, nel is between boundaries[ib] and boundaries[ib+1] OR below;
        double (*wf)(const U& a, const U& b);
        switch(opts.window)
        {
        case HGWDelta:
            if (ib==0) 
            {
                if (nel>=opts.boundaries[0]) bins[0]+=weight/(opts.boundaries[1]-opts.boundaries[0]);
                else below+=weight;
            }
            else if(ib==bs) above+=weight;
            else bins[ib]+=weight/(opts.boundaries[ib+1]-opts.boundaries[ib]);
            break;
        case HGWBox:
            wf=__hgwfbox;
            break;
        case HGWTriangle:
            wf=__hgwftri;
            break;
        default:
            ERROR("Windowing mode not implemented yet!\n");
        }
        if (opts.window!=HGWDelta) 
        {
            double nb=0.;
            for (ia=ib-1; ia>=0; --ia)
            {
                nb=wf((opts.boundaries[ia]-nel)/opts.window_width,(opts.boundaries[ia+1]-nel)/opts.window_width)
                        /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                //std::cerr<<":"<<nb<<" ? "<<bins[ia]<<" \n";
                if (nb==0) break; else bins[ia]+=nb*weight;
            }
            if (ia<0) {
                below+=weight*wf(-std::numeric_limits<U>::max(),(opts.boundaries[0]-nel)/opts.window_width);
            }
            for (ia=ib; ia<bs; ++ia)
            {
                nb=wf((opts.boundaries[ia]-nel)/opts.window_width,(opts.boundaries[ia+1]-nel)/opts.window_width)
                        /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                //std::cerr<<":"<<nb<<" ? "<<bins[ia]<<" \n";
                if (nb==0) break; else bins[ia]+=weight*nb;
            }
            if (ia==bs) {
                above+=weight*wf((opts.boundaries[bs]-nel)/opts.window_width,std::numeric_limits<U>::max());
            }
        }
        ndata+=weight;
    }
    
/***********************************************************
 *             N-DIMENSIONAL HISTOGRAM                     *
 ***********************************************************/

template <class U> class NDHistogram {
private:
    unsigned long dim;
    std::valarray<double> bins, vols;
    std::valarray<long> nbin;
    double outliers;
    double ndata;
    std::valarray<HGOptions<Histogram<U> > > opts;
    long c2b(const std::valarray<long>& vl) 
    {
        long ts=1,rp=0;
        for (int i=0; i<dim;++i)
        { rp+=vl[i]*ts; ts*=nbin[i]; }
        return rp;
    }

public:
    void reset() 
    {
        long tsz=1; nbin.resize(dim);
        for (int i=0; i<dim; ++i) tsz*=(nbin[i]=(opts[i].boundaries.size()-1));
        if (tsz<0) tsz=0;
        
        bins.resize(tsz); 
        bins=0.; ndata=0; outliers=0.;
        
        std::valarray<long> cp(dim); cp=0;
        long k=0; vols.resize(tsz); 
        while (cp[dim-1]<nbin[dim-1])
        {
            vols[k]=1.;
            for (int i=0; i<dim; ++i) vols[k]*=(opts[i].boundaries[cp[i]+1]-opts[i].boundaries[cp[i]]);
            k++; cp[0]++;
            for (int i=0; i<dim-1; ++i) if (cp[i]>=nbin[i]) {cp[i]=0; ++cp[i+1];}
        }
    }

    double samples() { return ndata; }

    void get_options(std::valarray<HGOptions<Histogram<U> > >& rop)
    { rop.resize(opts.size); rop=opts; }

    void set_options(const std::valarray<HGOptions<Histogram<U> > >& rop)
    { opts.resize(dim=rop.size()); opts=rop; reset(); }

    NDHistogram(const std::valarray<HGOptions<Histogram<U> > >& rop)
    { set_options(rop); }

//copy & assign...
    template<class T> NDHistogram(const NDHistogram<T>& ac);
    template<class T> NDHistogram<U>& operator=(const NDHistogram<T>& ac);

    void add(const std::valarray<U>& nel, double weight=1.0);

//syntactic sugar to insert a new element into the series, to move to next or previous series
    inline void operator << (const std::valarray<U>& nel) { add(nel); }
    
    double max() const { return bins.max()/ndata; }
    double min() const { return bins.min()/ndata; }
    
    
    void get_bins(std::valarray<double>& rbins) const
    {
        rbins.resize(bins.size());
        rbins=bins*(1./ndata);
    }

    void get_bins(std::valarray<std::valarray<U> >& cent, std::valarray<std::valarray<U> >& ws, std::valarray<double>& rbins) const
    {
        rbins.resize(bins.size()); cent.resize(dim); ws.resize(dim);
        rbins=bins*(1./ndata);
        for (unsigned long i=0; i<dim;++i)
        {
            cent[i].resize(nbin[i]); ws[i].resize(nbin[i]);
            for (unsigned long k=0; k<nbin[i];++k)
            {
                cent[i][k]=(opts[i].boundaries[k]+opts[i].boundaries[k+1])/2.;
                ws[i][k]=(opts[i].boundaries[k+1]-opts[i].boundaries[k]);
            }
        }
    }
    
    void get_bin(const std::valarray<long>& index, std::valarray<double>& center, double& val)
    {
        long ibin=c2b(index);
        center.resize(dim); val=bins[ibin]*(1./ndata);
        for (unsigned long i=0; i<dim;++i) 
            center[i]=(opts[i].boundaries[index[i]]+opts[i].boundaries[index[i]+1])/2.;
    }
    
    void get_outliers(double& routliers) const
    {
        routliers=outliers/ndata; 
    }
};

template<class U>
std::ostream& operator<<(std::ostream& os, const NDHistogram<U>& his)
{
    std::valarray<U> wf;
    std::valarray<std::valarray<U> > wx, ww;
    his.get_bins(wx,ww,wf);
    os.precision(12);
    os.setf(std::ios::scientific);
    os.width(14);
    
    U outliers; his.get_outliers(outliers);
    os<<"# Fraction of outliers: "<<outliers<<"\n";
    for (unsigned int i=0; i<wx.size(); ++i)
    {
        os<<"# x("<<i<<"): ";
        for (unsigned int j=0; j<wx[i].size(); ++j) os<<wx[i][j]<<" ";
        os<<"\n# w("<<i<<"): ";
        for (unsigned int j=0; j<ww[i].size(); ++j) os<<ww[i][j]<<" ";
        os<<"\n";
    }
    
    for (unsigned int i=0; i<wf.size(); ++i)
        os<<wf[i]<<"\n";
    return os;
}

template<class U>
void NDHistogram<U>::add(const std::valarray<U>& nel, double weight)
{
    std::valarray<long> p(dim);
    long bs, ia, ib, ic;
    //First, finds the "coordinates" of the center
    for (int i=0; i<dim; ++i)
    {
        bs=nbin[i]; ia=0; ic=bs; ib=(ic+1)/2;
        while (ib>ia && ib<ic)
        {
            if (nel[i]>opts[i].boundaries[ib+1])
            {   ia=ib;  ib=(ib+ic+1)/2; }
            else
            {
                if (nel[i]>opts[i].boundaries[ib]) break;
                else { ic=ib; ib=(ia+ib)/2;  }
            }
        }
        p[i]=ib;
    }
    ndata+=weight;
    //now, nel is between boundaries[p[i]] and boundaries[p[i]+1] OR below in each dim
    double (*wf)(const U& a, const U& b);
    std::valarray<std::valarray<double> > tbins(dim);
    int i;
    for (i=0; i<dim; ++i) 
    {
        tbins[i].resize(nbin[i]); tbins[i]=0.;
        //std::cerr<<"running through dimension "<<i<<":"<<tbins[i].size()<<"\n";
        std::valarray<long> ap(p);
        /*only delta-function, by now...*/
        switch(opts[i].window)
        {
        case HGWDelta:
            if ((p[i]==0 && nel[i]< opts[i].boundaries[0]) || p[i]>=nbin[i])  break;
            tbins[i][p[i]]=weight;///(opts[i].boundaries[p[i]+1]-opts[i].boundaries[p[i]]);
            break;
        case HGWBox:
            wf=__hgwfbox;
            break;
        case HGWTriangle:
            wf=__hgwftri;
            break;
        default:
            ERROR("Windowing mode not implemented yet!\n");
        }
        if (opts[i].window!=HGWDelta) 
        {
            double nb=0.;
            for (ia=p[i]-1; ia>=0; --ia)
            {
                nb=wf((opts[i].boundaries[ia]-nel[i])/opts[i].window_width,(opts[i].boundaries[ia+1]-nel[i])/opts[i].window_width);///(opts[i].boundaries[ia+1]-opts[i].boundaries[ia]);
                //std::cerr<<nb<<" ** \n";
                if (nb==0) break; else tbins[i][ia]+=nb;
            }
            for (ia=p[i]; ia<nbin[i]; ++ia)
            {
                nb=wf((opts[i].boundaries[ia]-nel[i])/opts[i].window_width,(opts[i].boundaries[ia+1]-nel[i])/opts[i].window_width);///(opts[i].boundaries[ia+1]-opts[i].boundaries[ia]);
                //std::cerr<<nb<<" ** \n";
                if (nb==0) break; else tbins[i][ia]+=nb;
            }
        }
    }
    //now we must make all the products and increment
    if (i<dim) { outliers+=weight; return; }
    std::valarray<long> minp(dim), maxp(dim); minp=0; maxp=0;
    double outs=1.;
    //lower boundary of nonzero region
    for (i=0; i<dim; ++i)
    {
        int j;
        for (j=0; j<nbin[i] && tbins[i][j]==0.;) ++j;
        if (j==nbin[i]) { outliers+=weight; return; }
        minp[i]=j;
        for (;j<nbin[i]&&tbins[i][j]!=0.;) ++j;
        maxp[i]=j;
    }
    std::valarray<long> cp(minp);
    while (cp[dim-1]<maxp[dim-1])
    {
        int j;
        
        long k=c2b(cp);
        //for (i=0; i<dim; ++i) std::cerr<<cp[i]<<" ";
        //std::cerr<<std::endl;
        double tv=1.; for (i=0; i<dim; ++i) tv*=tbins[i][cp[i]];
        bins[k]+=tv/vols[k]*weight;
        outs-=tv;//*vols[k]; 
        cp[0]++;
        for (i=0; i<dim-1; ++i) if (cp[i]>=maxp[i]) {cp[i]=minp[i]; ++cp[i+1];}
    }
    //std::cerr<<"outliers: "<<outs<<"\n";
    outliers+=outs*weight;
}

}//ends namespace toolbox
#endif  //ends ifdef __TOOLS_HISTOGRAM_H
