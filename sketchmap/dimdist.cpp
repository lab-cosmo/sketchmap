#include "dimreduce.hpp"
#include "clparser.hpp"
#include "matrix-io.hpp"
#include "matrix-conv.hpp"
#include "tools-histogram.hpp"

using namespace toolbox;

void banner() 
{
    std::cerr
            << " USAGE: dimdist -D hi-dim -d low-dim -P hd-file -p ld-file -pi period [-w]      \n"
            << "               [-maxd maxd[ maxd2]] [-nbin nbin[ nbin2]] [-wbin wbin]           \n"
            << "               [-gnuplot] [-h]                                                  \n"
            << "                                                                                \n"
            << " computes a histogram of D-d fits based on the lists of points given in input.  \n"
            << " dimension is set by -D option, and the projection is performed down to the     \n"
            << " dimensionality specified by -d. Optionally, high-dimensional data may be       \n"
            << " assumed to lie in a hypertoroidal space with period -pi.                       \n"
            << " Optionally, function fun [identity, sigmoid, compress] can be applied to the   \n"
            << " distances, using the scale parameter sigma.                                    \n"
            << " Data must be provided in the files set by -P and -p in the form:               \n"
            << " X1_1, X1_2, ... X1_D                                                           \n"
            << " X2_1, X2_2, ... X2_D                                                           \n"
            << " Output is given as Di di n(Di,di) dD dd, optionally [-gnuplot] with breaks     \n"
            << " to make it compatible with gnuplot.                                            \n";
}
typedef struct { double D, d; } dpair;

int main(int argc, char**argv)
{
    CLParser clp(argc,argv);
    unsigned long D, d; double peri, speri, wbin; bool fgnuplot, fhelp, fweight, flowmem;
    std::vector<double> maxd; std::vector<unsigned long> nbin;
    std::string fmds, fhd, fld;
    bool fok=clp.getoption(D,"D",(unsigned long) 3) && 
            clp.getoption(d,"d",(unsigned long) 0) &&
            clp.getoption(peri,"pi",0.0) &&
            clp.getoption(speri,"spi",0.0) &&
            clp.getoption(fhelp,"h",false) &&
            clp.getoption(fweight,"w",false) &&
            clp.getoption(fgnuplot,"gnuplot",false) && 
            clp.getoption(flowmem,"lowmem",false) && 
            clp.getoption(fld,"p",std::string("")) &&
            clp.getoption(fhd,"P",std::string("")) &&
            clp.getoption(nbin,"nbin",std::vector<unsigned long>(1,100u)) &&
            clp.getoption(maxd,"maxd",std::vector<double>(0)) &&
            clp.getoption(wbin,"wbin",0.0);
    
    if (fhelp || !fok) { banner(); exit(1); }
    if (fhd=="" || (d>0 && fld=="")) ERROR("Hi-dim and low-dim points must be provided by the -P and -p options"); 
    
    std::ifstream sP(fhd.c_str()), sp(fld.c_str()); 
    
    // reads points from standard input
    FMatrix<double> HP, lp;
    std::vector<std::vector<double> > plist; std::vector<double> point(D), pweight;
    while (sP.good())
    {
        for (int i=0; i<D; i++) sP>>point[i];
        double nw;  if (fweight) sP>>nw; else nw=1;
        
        if (sP.good()) { plist.push_back(point); pweight.push_back(nw); }
    }
    
    HP.resize(plist.size(), D);
    for (unsigned long i=0; i<plist.size(); ++i) for (unsigned long j=0; j<D; ++j) HP(i,j)=plist[i][j];
    
    if (d>0) {
        point.resize(d); plist.clear();
        while (sp.good())
        {
            for (int i=0; i<d; i++) sp>>point[i];
            if (sp.good()) plist.push_back(point);
        }
        
        lp.resize(plist.size(), d);
        for (unsigned long i=0; i<plist.size(); ++i) for (unsigned long j=0; j<d; ++j) lp(i,j)=plist[i][j];
        if (lp.rows()!=HP.rows()) ERROR("Mismatch between hi-d and low-d point sets");
    }
    
    NLDRMetricPBC nperi; NLDRMetricEuclid neuclid; NLDRMetricSphere nsphere;
    nperi.periods.resize(D); nperi.periods=peri;
    nsphere.periods.resize(D); nsphere.periods=speri;
    NLDRMetric *metric;
    
    if (peri==0.0 && speri==0.0) metric=&neuclid;
    else if (speri==0) { metric=&nperi; }
    else { metric=&nsphere; }
    
    //computes distances. this could be made on the fly, but would make impossible
    //set automatically the size of the histogram
    std::valarray<dpair> distances; 
    std::valarray<double> dweight; 
    if (!flowmem) {distances.resize((HP.rows()*(HP.rows()-1))/2); dweight.resize((HP.rows()*(HP.rows()-1))/2);}
    
    
    unsigned long k=0; double mxD, mxd; mxD=mxd=0.0;
    if (!flowmem)
    for(unsigned long i=0; i<HP.rows(); ++i) for(unsigned long j=0; j<i; ++j)
    { 
        distances[k].D=metric->dist(&HP(i,0),&HP(j,0),D);
        dweight[k]=pweight[i]*pweight[j];
        if (distances[k].D>mxD) mxD=distances[k].D;
        if (d>0)
        {
            distances[k].d=neuclid.dist(&lp(i,0),&lp(j,0),d);
            if (distances[k].d>mxd) mxd=distances[k].d;
        }
        ++k;
    }
    else
    for(unsigned long i=0; i<HP.rows(); ++i) for(unsigned long j=0; j<i; ++j)
    {   
        double kD, kd;
        kD=metric->dist(&HP(i,0),&HP(j,0),D);
        if (kD>mxD) mxD=kD;
        if (d>0)
        {
            kd=neuclid.dist(&lp(i,0),&lp(j,0),d);
            if (kd>mxd) mxd=kd;
        }
        ++k;
    }
    
    if (d>0)
    {
        //creates n-dimensional histogram and bins it
        std::valarray<HGOptions<Histogram<double> > > hgo(2);
        hgo[0].window=hgo[1].window=(wbin==0.0?HGWDelta:HGWTriangle);
        hgo[0].window_width=hgo[1].window_width=wbin;
        if (nbin.size()==1) nbin.push_back(nbin[0]);
        if (nbin.size()==0 || nbin.size() >2) ERROR("-nbin options requires one or two arguments.");
        nbin[0]+=1; nbin[1]+=1;
        hgo[0].boundaries.resize(nbin[0]); hgo[1].boundaries.resize(nbin[1]);
        if (maxd.size() >2) ERROR("-maxd options requires one or two arguments, if given.");
        if (maxd.size()==0) { maxd.push_back(mxD); maxd.push_back(mxd); }
        if (maxd.size()==1) maxd.push_back(maxd[0]);
        for(unsigned long i=0; i<2; ++i)  for (k=0; k<hgo[i].boundaries.size();k++)
            hgo[i].boundaries[k]=k*maxd[i]/(hgo[i].boundaries.size()-1);
        
        NDHistogram<double> ndh(hgo); std::valarray<double> dval(2);
    #define verysmall 1e-3
        if (!flowmem)
        for (k=0; k<distances.size(); ++k) 
        {
            dval[0]=distances[k].D; dval[1]=distances[k].d;  ndh.add(dval,dweight[k]); 
        }
        else 
        for(unsigned long i=0; i<HP.rows(); ++i) for(unsigned long j=0; j<i; ++j)
        {
            dval[0]=metric->dist(&HP(i,0),&HP(j,0),D);
            dval[1]=neuclid.dist(&lp(i,0),&lp(j,0),d);
            ndh.add(dval,pweight[i]*pweight[j]); 
        }
        //outputs
        double outliers;  ndh.get_outliers(outliers);
        std::cout<<"# Fraction outside: "<<outliers<<std::endl;
        if (!fgnuplot) std::cout <<ndh; 
        else
        {
            
            std::valarray<long> ind(2); std::valarray<double> cen(2); double val;
            for (int i=0; i<hgo[0].boundaries.size()-1; ++i) 
            {
                for (int j=0; j<hgo[1].boundaries.size()-1; ++j)
                {
                    ind[0]=i; ind[1]=j;
                    ndh.get_bin(ind,cen,val);
                    std::cout<<cen[0]<<"\t"<<cen[1]<<"\t"<<val<<"\n";
                }
                std::cout<<std::endl;
            }
        }
    }
    else
    {
        // just bin distances in 1D
        HGOptions<Histogram<double> > hgo;
        hgo.window=(wbin==0.0?HGWDelta:HGWTriangle);
        hgo.window_width=wbin;
        if (nbin.size()==0 || nbin.size() >2) ERROR("-nbin options requires one or two arguments.");
        nbin[0]+=1; 
        hgo.boundaries.resize(nbin[0]);
        if (maxd.size() >2) ERROR("-maxd options requires one or two arguments, if given.");
        if (maxd.size()==0) { maxd.push_back(mxD); }

        for (k=0; k<hgo.boundaries.size();k++) hgo.boundaries[k]=k*maxd[0]/(hgo.boundaries.size()-1);
        
        Histogram<double> hh(hgo); double dval; double weight;
#define verysmall 1e-3
        if (!flowmem)
            for (k=0; k<distances.size(); ++k) hh.add(distances[k].D,dweight[k]);
        else
            for(unsigned long i=0; i<HP.rows(); ++i) for(unsigned long j=0; j<i; ++j) hh.add(metric->dist(&HP(i,0),&HP(j,0),D),pweight[i]*pweight[j]);
        //outputs
        double outd,outu;  hh.get_outliers(outd,outu);
        std::cout<<"# Fraction outside: "<<outd<<" - "<<outu<<std::endl;
        std::cout <<hh;
    }
}