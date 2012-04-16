/* Strided or min-max selection of landmarks
   --------------------------------------------------
   Author: Michele Ceriotti, 2011
   Distributed under the GNU General Public License  
*/

#include "dimreduce.hpp"
#include "clparser.hpp"
#include "matrix-io.hpp"
#include "rndgen.hpp"

using namespace toolbox;
void banner() 
{
    std::cerr
            << " USAGE: landmark -D dim -n nlandmark [-pi period] [-spi period] [-wi]           \n"
            << "                [-i] [-w] [-mode stride|minmax|resample|staged] [-lowmem] [-rnd nrnd]  \n"
            << "                                                                                \n"
            << " selects n landmark points from N given in input, with dimensionality D, using  \n"
            << " a greedy MinMax approach (minmax, default), a constant-stride or random        \n"
            << " criterion.                                                                     \n"
            << " Points can be read with weights [-wi], and distances may be computed with      \n"
            << " thoroidal [-pi] or spherical [-spi] periodic boundary conditions.              \n"
            << " Optionally [-w] assigns a weight to each landmark by performing a voronoi      \n"
            << " tassellation and assigning the original N points. [-i] print indices of points.\n"
            << " Use -lowmem when -n and N are both large and the nxN distance matrix won't     \n"
            << " fit memory.                                                                    \n"
            << " -rnd selecs nrnd random points before starting the usual selection scheme.    \n";
}

int main(int argc, char**argv)
{
    CLParser clp(argc,argv);
    unsigned long D,n,N,nrand,seed; 
    bool fhelp, fweight, findex, flowm, finw;
    double peri, speri, gamma, alpha; std::string smode;
    bool fok=clp.getoption(D,"D",(unsigned long) 3) && 
            clp.getoption(n,"n",(unsigned long) 100) &&
            clp.getoption(nrand,"rnd",(unsigned long) 0) &&
            clp.getoption(seed,"seed",(unsigned long) 12345) &&            
            clp.getoption(fhelp,"h",false) &&
            clp.getoption(fweight,"w",false) &&
            clp.getoption(findex,"i",false) &&
            clp.getoption(finw,"wi",false) &&
            clp.getoption(flowm,"lowmem",false) &&
            clp.getoption(smode,"mode",std::string("minmax")) &&
            clp.getoption(gamma,"gamma",1.0) &&
            clp.getoption(alpha,"alpha",1.0) &&            
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
    
    double maxd, dij; unsigned long maxj; std::valarray<double> mdlist(N); 
    std::valarray<unsigned long> isel(n); 
    isel[0]=0;
    if (nrand>0) 
    {
        StdRndUniform rndgen(seed);
        std::cerr<<"picking "<<nrand<<" random points (duplicates are possible)\n";
        isel[0]=maxj=rndgen()*N;
        LP.row(0)=HP.row(maxj);
        for (unsigned long j=0; j<N; ++j) { mdlist[j]=metric->dist(&LP(0,0),&HP(j,0),D); }
        for (unsigned long i=1; i<nrand; ++i)
        {
            isel[i]=maxj=rndgen()*N;
            LP.row(i)=HP.row(maxj);
            for (unsigned long j=0; j<N; ++j) 
            { dij=metric->dist(&LP(i,0),&HP(j,0),D); if (mdlist[j]>dij) mdlist[j]=dij; }
        }
    }
    else
    {
        isel[0]=0;
        LP.row(0)=HP.row(0);
        for (unsigned long j=0; j<N; ++j) { mdlist[j]=metric->dist(&LP(0,0),&HP(j,0),D); }         
    }
    
    std::cerr<<"picking "<<n <<" points out of "<<N<<"\n";
    if (flowm)
    {
        if (smode=="minmax") 
        {
            for (unsigned long i=(nrand==0?1:nrand); i<n; ++i)
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
            for (unsigned long i=(nrand==0?1:nrand); i<n; i++)
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
        else if (smode=="resample")
        {
            StdRndUniform rndgen(seed);
            std::cerr<<"initialize with seed"<<seed<<"first number is "<<rndgen()<<"\n";
            mdlist=0.0;
            isel[0]=maxj=rndgen()*N;   //first picks a point randomly
            LP.row(0)=HP.row(maxj);
            double totpi=0.0, selpi, lij, lmax;
            double a2=alpha*alpha*2;
            for (unsigned long i=1; i<n; i++)
            {
               //compute a new set of distances
               totpi=0.0;
               for (unsigned long j=0; j<N; ++j) 
//               { dij=metric->dist(&LP(i-1,0),&HP(j,0),D); mdlist[j]+=(dij==0.0?1.0e200:pow(dij,-gamma)); totpi+=pow(mdlist[j],-gamma); }
//               { dij=metric->dist(&LP(i-1,0),&HP(j,0),D); mdlist[j]+=(dij==0.0?1.0e200:pow(dij,-2.0)); totpi+=wlist[j]*pow(mdlist[j],-gamma); }
//               { dij=metric->dist(&LP(i-1,0),&HP(j,0),D); mdlist[j]+=exp(-dij*dij/(2*0.25)); totpi+=wlist[j]*pow(mdlist[j],-gamma); }
//             // mdlist[j] stores the logarithm of the KMC weight for point j
               // since the individual points contribution can be tiny, it is important to accumulate the logs
               { 
                  dij=metric->dist(&LP(i-1,0),&HP(j,0),D); 
                  lij=dij*dij/(a2);
                  if (i==1) mdlist[j]=gamma*lij; // in the first step there is little to think. just weight in the 
                  else
                  {
                     if (mdlist[j]<gamma*lij)
                        mdlist[j]=mdlist[j]-gamma*log(1.0+exp(-lij+mdlist[j]/gamma));
                     else
                        mdlist[j]=gamma*lij-gamma*log(1.0+exp(lij-mdlist[j]/gamma));         
                 }
                  if (j==0 || mdlist[j]>lmax) lmax=mdlist[j];
//                  std::cerr<<"mdlist  "<<dij<<" :: "<<j<<" , "<<mdlist[j]<<"\n";
               }
               for (unsigned long j=0; j<N; ++j) totpi+=wlist[j]*exp(mdlist[j]-lmax);
               std::cerr<<i<<"  : lmax= "<<lmax<< "  totpi "<< totpi<< "\n";
//               exit(1);
               selpi=rndgen()*totpi;
               

               //a better search could be used here but hey....               
               unsigned long j=0;
               for (j=0; j<N; ++j) 
//               {  selpi-=wlist[j]*pow(mdlist[j],-gamma); if(selpi<0.0) { j++; break; } }
               {  selpi-=wlist[j]*exp(mdlist[j]-lmax); if(selpi<0.0) { j++; break; } }
               
               isel[i]=j-1; LP.row(i)=HP.row(isel[i]);      
               std::cerr<<"selecting point "<<i<<" : "<<isel[i]<<"("<<wlist[j]*exp(mdlist[j]-lmax)/totpi<<")\n";
            }
        }
        else if (smode=="staged")        
        {
            StdRndUniform rndgen(seed);
            FMatrix<double> MP;
            unsigned long m=sqrt(n*N);
            std::valarray<unsigned long> msel(m);
            MP.resize(m, D); 
            std::cerr<<"initialize with seed"<<seed<<"first number is "<<rndgen()<<"\n";
            mdlist=0.0;
            msel[0]=maxj=rndgen()*m;   //first picks a point randomly
            MP.row(0)=HP.row(maxj);
            for (unsigned long j=0; j<N; ++j) { mdlist[j]=metric->dist(&MP(0,0),&HP(j,0),D); }       
        
            std::cerr<<"now selecting "<<m<< " minmax points\n";
            for (unsigned long i=1; i<m; ++i)
            {
                //std::cerr<<mdlist<<"\n";
                maxd=0.;  for (unsigned long j=0; j<N; ++j) if (mdlist[j]>maxd) {maxd=mdlist[j]; maxj=j;}
                std::cerr<<"selecting point "<<i<<" : "<<maxj<<"("<<maxd<<")\n";
                msel[i]=maxj;
                MP.row(i)=HP.row(maxj);
                for (unsigned long j=0; j<N; ++j) 
                { dij=metric->dist(&MP(i,0),&HP(j,0),D); if (mdlist[j]>dij) mdlist[j]=dij; }
            }
            
            //now gets Voronoi weights for the m points
            std::valarray<double> mweights(m); mweights=0.0;
            std::cerr<<"now running Voronoi weight assignment\n";
            double mind; unsigned long mini;
            std::vector<std::vector<unsigned long> > vplist(m);
            for (unsigned long j=0; j<N; ++j)
            {
                mind=metric->dist(&LP(0,0),&HP(j,0),D); mini=0;
                for (unsigned long i=1; i<m; ++i)
                {
                    if ((dij=metric->dist(&MP(i,0),&HP(j,0),D))<mind) { mind=dij; mini=i; }
                }
                mweights[mini]+=wlist[j];
                vplist[mini].push_back(j);
                
            }

            double tw=0.0;
            for (unsigned long i=0; i<m; ++i) {  mweights[i]=pow(mweights[i],gamma); tw+=mweights[i]; }
             
            //now picks n points
            double selpi; unsigned long subsel;
            std::cerr<<"now running picking "<<n<<" final points total weight is "<<tw <<" \n";
            for (unsigned long i=0; i<n; ++i)
            {
               selpi=rndgen()*tw;
               unsigned long j=0;
               for (j=0; j<m; ++j) 
               {  selpi-=mweights[j]; if(selpi<0.0) { j++; break; } }
               j--;
               
               std::cerr<<"picked "<<j<<" weight: "<<mweights[j]<<"\n";
 //             isel[i]=msel[j]; LP.row(i)=MP.row(j);      

               //actually, picks a random point from the voronoi polihedron
               subsel=vplist[j][rndgen()*vplist[j].size()];
               std::cerr<<"  subpoint "<<subsel<<" selected\n";
               isel[i]=subsel; LP.row(i)=HP.row(subsel);
            }
            
            //TODO the "voronoi" weight of these landmarks is now NOT related to the original probability distribution
        }
        else ERROR("Selection mode "<<smode<<" not implemented yet\n");
    }
    else
    {  //!TODO DOUBLE CHECK THIS CODE PATH
        if (smode=="minmax") 
        {
            LD*=0.0;
            for (unsigned long i=0; i<nrand; ++i)
            for (unsigned long j=0; j<N; ++j) { LD(i,j)=metric->dist(&LP(i,0),&HP(j,0),D); } 
            mdlist=LD.row(0);
            for (unsigned long i=(nrand==0?1:nrand); i<n; ++i)
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
            for (unsigned long i=(nrand==0?1:nrand); i<n; i++)
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
    std::cout<<"#"; if (nrand>0) std::cout<<" "<<nrand<<" points selected randomly.  ";
    std::cout<<" Max distance detected: "<<maxd<<"\n";
    for (unsigned long i=0; i<n; ++i)
    {
        if (findex) std::cout<<isel[i]<<" ";
        for (unsigned long j=0; j<D; ++j) std::cout<<LP(i,j)<<" ";
        if (findex && finw) std::cout<<wlist[isel[i]]<<"  ";
        if (fweight) std::cout<<weights[i];
        
        std::cout<<"\n";
    }
}
