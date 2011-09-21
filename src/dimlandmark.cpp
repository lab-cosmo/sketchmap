#include "dimreduce.hpp"
#include "clparser.hpp"
#include "matrix-io.hpp"

using namespace toolbox;
void banner() 
{
    std::cerr
            << " USAGE: landmark -D dim -n nlandmark [-pi period] [-spi period] [-wi]           \n"
            << "                [-i] [-w] [-mode stride|rnd|minmax] [-lowmem]                   \n"
            << "                                                                                \n"
            << " selects n landmark points from N given in input, with dimensionality D, using  \n"
            << " a greedy MinMax approach (minmax, default), a constant-stride or random        \n"
            << " criterion.                                                                     \n"
            << " Points can be read with weights [-wi], and distances may be computed with      \n"
            << " thoroidal [-pi] or spherical [-spi] periodic boundary conditions.              \n"
            << " Optionally [-w] assigns a weight to each landmark by performing a voronoi      \n"
            << " tassellation and assigning the original N points. [-i] print indices of points.\n"
            << " Use -lowmem when -n and N are both large and the nxN distance matrix won't     \n"
            << " fit memory.                                                                    \n";
}

int main(int argc, char**argv)
{
    CLParser clp(argc,argv);
    unsigned long D,n,N; 
    bool fhelp, fweight, findex, flowm, finw;
    double peri, speri; std::string smode;
    bool fok=clp.getoption(D,"D",(unsigned long) 3) && 
            clp.getoption(n,"n",(unsigned long) 100) &&
            clp.getoption(fhelp,"h",false) &&
            clp.getoption(fweight,"w",false) &&
            clp.getoption(findex,"i",false) &&
            clp.getoption(finw,"wi",false) &&
            clp.getoption(flowm,"lowmem",false) &&
            clp.getoption(smode,"mode",std::string("minmax")) &&
            clp.getoption(speri,"spi",0.0) &&
            clp.getoption(peri,"pi",0.0);
    
    if (fhelp || !fok) { banner(); exit(1); }
    // reads points from standard input
    std::vector<std::vector<double> > plist; std::vector<double> point(D), wlist;
    double neww;
    while (std::cin.good())
    {
        for (int i=0; i<D; i++) std::cin>>point[i];
        if (finw) std::cin>>neww; else neww=1.0;
        if (std::cin.good()) { plist.push_back(point); wlist.push_back(neww); }
    }
    
    FMatrix<double> HP, LP, LD;
    HP.resize(N=plist.size(), D);
    
    for (unsigned long i=0; i<N; ++i) for (unsigned long j=0; j<D; ++j) HP(i,j)=plist[i][j];
    NLDRMetricPBC nperi; NLDRMetricEuclid neuclid; NLDRMetricSphere nsphere;
    nperi.periods.resize(D); nperi.periods=peri;
    nsphere.periods.resize(D); nsphere.periods=speri;
    
    NLDRMetric *metric;
    if (peri==0.0 && speri==0.0) metric=&neuclid;
    else if (speri==0) { metric=&nperi; }
    else { metric=&nsphere; }

    LP.resize(n, D); 
    if (flowm) LD.resize(0,0); else LD.resize(n, N);
    
    LP.row(0)=HP.row(0);
    double maxd; unsigned long maxj; std::valarray<double> mdlist(N); 
    std::valarray<unsigned long> isel(n); 
    isel[0]=0;
    std::cerr<<"picking "<<n <<" points out of "<<N<<"\n";
    if (flowm)
    {
        double dij;
        if (smode=="minmax") 
        {
            for (unsigned long j=0; j<N; ++j) { mdlist[j]=metric->dist(&LP(0,0),&HP(j,0),D); } 
            for (unsigned long i=1; i<n; ++i)
            {
                //std::cerr<<mdlist<<"\n";
                maxd=0.;  for (unsigned long j=0; j<N; ++j) if (mdlist[j]>maxd) {maxd=mdlist[j]; maxj=j;}
                std::cerr<<"selecting point "<<i<<" : "<<maxj<<"("<<maxd<<")\n";
                isel[i]=maxj;
                LP.row(i)=HP.row(maxj);
                for (unsigned long j=0; j<N; ++j) 
                { dij=metric->dist(&LP(i,0),&HP(j,0),D); if (mdlist[j]>dij) mdlist[j]=dij; }
            }
        }
        else if (smode=="stride")
        {
            unsigned long stride=(N/n);
            for (unsigned long j=0; j<N; ++j) { mdlist[j]=metric->dist(&LP(0,0),&HP(j,0),D); } 
            for (unsigned long i=1; i<n; i++)
            {
                LP.row(i)=HP.row(i*stride);
                isel[i]=i*stride;
                maxd=0.0;
                for (unsigned long j=0; j<N; ++j) 
                { dij=metric->dist(&LP(i,0),&HP(j,0),D); if (mdlist[j]> dij) mdlist[j]=dij; }

                maxd=0.;  for (unsigned long j=0; j<N; ++j) if (mdlist[j]>maxd) maxd=mdlist[j];
                std::cerr<<"selecting point "<<i<<" : "<<i*stride<<"("<<maxd<<")\n";
            }
        }
        else ERROR("Selection mode "<<smode<<" not implemented yet\n");
    }
    else
    {
        if (smode=="minmax") 
        {
            LD*=0.0;
            for (unsigned long j=0; j<N; ++j) { LD(0,j)=metric->dist(&LP(0,0),&HP(j,0),D); } 
            mdlist=LD.row(0);
            for (unsigned long i=1; i<n; ++i)
            {
                
                maxd=0.;  for (unsigned long j=0; j<N; ++j) if (mdlist[j]>maxd) {maxd=mdlist[j]; maxj=j;}
                std::cerr<<"selecting point "<<i<<" : "<<maxj<<"("<<maxd<<")\n";
                LP.row(i)=HP.row(maxj);
                isel[i]=maxj;
                for (unsigned long j=0; j<N; ++j) 
                { LD(i,j)=metric->dist(&LP(i,0),&HP(j,0),D); if (mdlist[j]> LD(i,j)) mdlist[j]=LD(i,j); }
            }
        }
        else if (smode=="stride")
        {
            unsigned long stride=(N/n); mdlist=-1.0;
            for (unsigned long i=0; i<n; i++)
            {
                LP.row(i)=HP.row(i*stride);
                isel[i]=i*stride;
                maxd=0.0;
                for (unsigned long j=0; j<N; ++j) { LD(i,j)=metric->dist(&LP(i,0),&HP(j,0),D); if (mdlist[j]<0.0 || mdlist[j]> LD(i,j)) mdlist[j]=LD(i,j); }
                
                maxd=0.;  for (unsigned long j=0; j<N; ++j) if (mdlist[j]>maxd) maxd=mdlist[j];
                std::cerr<<"selecting point "<<i<<" : "<<i*stride<<"("<<maxd<<")\n";
            }
        }
        else ERROR("Selection mode "<<smode<<" not implemented yet\n");
    }
    std::valarray<unsigned long> weights(n); weights=0;
    if (fweight) 
    {  
        double mind; unsigned long mini;
        if (flowm)
        {
            double dij;
            for (unsigned long j=0; j<N; ++j)
            {
                mind=metric->dist(&LP(0,0),&HP(j,0),D); mini=0;
                for (unsigned long i=1; i<n; ++i)
                {
                    if ((dij=metric->dist(&LP(i,0),&HP(j,0),D))<mind) { mind=dij; mini=i; }
                }
                weights[mini]+=wlist[j];
            }
        }
        else
        {
            for (unsigned long j=0; j<N; ++j)
            {
                mind=LD(0,j); mini=0;
                for (unsigned long i=1; i<n; ++i)
                {
                    if (LD(i,j)<mind) { mind=LD(i,j); mini=i; }
                }
                weights[mini]+=wlist[j];
            }
        }
    }
    std::cout<<std::scientific; std::cout.precision(6); 
    std::cout<<"# "<<n<<" landmark points selected out of "<<N<<" and chosen by "<<smode<<"\n";
    std::cout<<"# Max distance detected: "<<maxd<<"\n";
    for (unsigned long i=0; i<n; ++i)
    {
        if (findex) std::cout<<isel[i]<<" ";
        for (unsigned long j=0; j<D; ++j) std::cout<<LP(i,j)<<" ";
        if (fweight) std::cout<<weights[i];
        std::cout<<"\n";
    }
}