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
            << "                [-i] [-w] [-ifirst i] [-seed s] [-unique]                       \n"
            << "                [-mode stride|random|minmax|resample|staged]                    \n"
            << "                [-gamma g] [-wgamma wg]                                         \n"
            << " tassellation and assigning the original N points. [-i] print indices of points.\n"
            << " selects n landmark points from N given in input, with dimensionality D, using  \n"
            << " a greedy MinMax approach (minmax, default), a constant-stride or random        \n"
            << " criterion.                                                                     \n"
            << " Points can be read with weights [-wi], and distances may be computed with      \n"
            << " thoroidal [-pi] or spherical [-spi] periodic boundary conditions.              \n"
            << " Optionally [-w] assigns a weight to each landmark by performing a voronoi      \n"
            << " tassellation and assigning the original N points. [-wgamma] specifies that     \n"
            << " weights must be manipolated by making them match the probability density       \n"
            << " elevated at power wg. [-i] print indices of points.                            \n"
            << " Generally the first point is chosen randomly unless otherwise specified with   \n"
            << " the option [-ifirst]. The PRNG seed can be set by [-seed] and one can make sure\n"
            << " that the selected landmarks are unique by [-unique].                           \n"
            << " [-mode mode] specifies how the landmarks should be selected:                   \n"
            << "   stride: n points are simply taken equally spaced along the trajectories.     \n"            
            << "   random: n points are picked at random.                                       \n"            
            << "   minmax: farthest point sampling.                                             \n"            
            << "   resample: mixed method which interpolates between random and minmax.         \n"            
            << "             it is controlled by the parameter gamma, where gamma~0 --> random  \n"            
            << "             and gamma -> infty --> minmax                                      \n"            
            << "   staged: picks the landmarks by a two-stage procedure, where first one takes  \n"            
            << "           several points to estimate the probability distribution of the       \n"            
            << "           input points P(x), then picks n with a scaled probability which is   \n"            
            << "           approximately P^gamma. In practice, gamma = 1 means that points are  \n" 
            << "           sampled randomly, gamma<1 means that the points are sampled at a     \n"            
            << "           'higher temperature' and gamma>1 that they are sampled at lower T.   \n";

            
}

int main(int argc, char**argv)
{
    CLParser clp(argc,argv);
    unsigned long D,n,N,seed; 
    long ifirst;
    bool fhelp, fweight, findex, flowm, finw, funique;
    double peri, speri, gamma, wgamma; std::string smode;
    bool fok=clp.getoption(D,"D",(unsigned long) 3) && 
            clp.getoption(n,"n",(unsigned long) 100) &&
            clp.getoption(ifirst,"first",(long) -1) &&
            clp.getoption(seed,"seed",(unsigned long) 12345) &&            
            clp.getoption(fhelp,"h",false) &&
            clp.getoption(fweight,"w",false) &&
            clp.getoption(findex,"i",false) &&
            clp.getoption(funique,"unique",false) &&            
            clp.getoption(finw,"wi",false) &&
            clp.getoption(smode,"mode",std::string("minmax")) &&
            clp.getoption(gamma,"gamma",1.0) &&     
            clp.getoption(wgamma,"wgamma",1.0) &&            
            clp.getoption(speri,"spi",0.0) &&
            clp.getoption(peri,"pi",0.0);
    
    if (fhelp || !fok) { banner(); exit(1); }

    // reads points from standard input
    std::vector<std::vector<double> > plist; std::vector<double> point(D), wlist;
    FMatrix<double> HP, LP;
    std::valarray<double> weights(n); weights=0;


    double neww;
    while (std::cin.good())
    {
        for (int i=0; i<D; i++) std::cin>>point[i];
        if (finw) std::cin>>neww; else neww=1.0;  // if requested, reads in weights for the points
        if (std::cin.good()) { plist.push_back(point); wlist.push_back(neww); }
    }
    
    // points are read in a vector of vectors (push_back is just too simple!). actually, we want them to be in a NxD matrix
    HP.resize(N=plist.size(), D);    
    for (unsigned long i=0; i<N; ++i) for (unsigned long j=0; j<D; ++j) HP(i,j)=plist[i][j];
    
    // metric objects initialized as per input    
    NLDRMetricPBC nperi; NLDRMetricEuclid neuclid; NLDRMetricSphere nsphere;
    nperi.periods.resize(D); nperi.periods=peri;
    nsphere.periods.resize(D); nsphere.periods=speri;
    
    NLDRMetric *metric;
    if (peri==0.0 && speri==0.0) metric=&neuclid;
    else if (speri==0) { metric=&nperi; }
    else { metric=&nsphere; }

    //allocates space for the landmarks
    LP.resize(n, D); 
           
    double maxd, dij; unsigned long maxj; std::valarray<double> mdlist(N); 
    std::valarray<unsigned long> isel(n); 
    StdRndUniform rndgen(seed);
    
//    isel[0]=0;
//    if (nrand>0) 
//    {
//        
//        std::cerr<<"picking "<<nrand<<" random points (duplicates are possible)\n";
//        isel[0]=maxj=rndgen()*N;
//        LP.row(0)=HP.row(maxj);
//        for (unsigned long j=0; j<N; ++j) { mdlist[j]=metric->dist(&LP(0,0),&HP(j,0),D); }
//        for (unsigned long i=1; i<nrand; ++i)
//        {
//            isel[i]=maxj=rndgen()*N;
//            LP.row(i)=HP.row(maxj);
//            for (unsigned long j=0; j<N; ++j) 
//            { dij=metric->dist(&LP(i,0),&HP(j,0),D); if (mdlist[j]>dij) mdlist[j]=dij; }
//        }
//    }
//    else
//    {
//        isel[0]=0;
//        LP.row(0)=HP.row(0);
//        for (unsigned long j=0; j<N; ++j) { mdlist[j]=metric->dist(&LP(0,0),&HP(j,0),D); }         
//    }

    // First point is picked at random -- except if otherwise specified     
    if (ifirst<0) isel[0]=maxj=rndgen()*N; else isel[0]=maxj=ifirst; 
    LP.row(0)=HP.row(maxj);
    for (unsigned long j=0; j<N; ++j) { mdlist[j]=metric->dist(&LP(0,0),&HP(j,0),D); }
    
    std::cerr<<"picking "<<n <<" points out of "<<N<<"\n";
    bool isunique;

    if (smode=="stride")
    {
        // select points with a fixed stride. In this case, we obviously start from zero so -ifirst 
        // request is not honoured. 
        
        unsigned long stride=(N/n);
        for (unsigned long i=0; i<n; i++)
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
    else if (smode=="random")
    {
        // select points randomly from the provided data. again, no point in honouring -ifirst
        
        for (unsigned long i=0; i<n; i++)
        {
            std::cerr<<"picking point "<<i<<std::endl;
            
            isunique=false;
            while (!isunique)
            {     
                isel[i]=maxj=rndgen()*N; 
                isunique=true; // if requested, makes sure that unique points are selected
                if (funique) for (unsigned long j=0; j<i; ++j) if (isel[j]==isel[i]) { isunique=false; break; }
            }
                    
            LP.row(i)=HP.row(isel[i]);
            maxd=0.0;
            for (unsigned long j=0; j<N; ++j) 
            { dij=metric->dist(&LP(i,0),&HP(j,0),D); if (mdlist[j]> dij) mdlist[j]=dij; }
            maxd=0.;  for (unsigned long j=0; j<N; ++j) if (mdlist[j]>maxd) maxd=mdlist[j];
        }
    }
    else if (smode=="minmax") 
    {
        // farthest point sampling selection of the points. no risk on getting duplicates here
        
        for (unsigned long i=1; i<n; ++i)
        {

            maxd=0.;  for (unsigned long j=0; j<N; ++j) if (mdlist[j]>maxd) {maxd=mdlist[j]; maxj=j;}
            std::cerr<<"selecting point "<<i<<" : "<<maxj<<"("<<maxd<<")\n";

            isel[i]=maxj;
            LP.row(i)=HP.row(maxj);
            for (unsigned long j=0; j<N; ++j) 
            { dij=metric->dist(&LP(i,0),&HP(j,0),D); if (mdlist[j]>dij) mdlist[j]=dij; }
        }
    }
    else if (smode=="resample")
    {

        // picks points with a random strategy which tries to interpolate between farthest sampling 
        // and random sampling. 

        double totpi=0.0, selpi, lij, lmax;
        for (unsigned long i=1; i<n; i++)
        {
           //compute a new set of distances
           totpi=0.0;         
           // computes w_i = (sum_j<i d_ij^-gamma)^-gamma
           for (unsigned long j=0; j<N; ++j) 
           { dij=metric->dist(&LP(i-1,0),&HP(j,0),D); mdlist[j]+=(dij==0.0?1.0e200:pow(dij,-gamma)); totpi+=pow(mdlist[j],-gamma); }

           //these are alternative ways to compute the distance-dependent weight
           
//               { dij=metric->dist(&LP(i-1,0),&HP(j,0),D); mdlist[j]+=(dij==0.0?1.0e200:pow(dij,-2.0)); totpi+=wlist[j]*pow(mdlist[j],-gamma); }
//               { dij=metric->dist(&LP(i-1,0),&HP(j,0),D); mdlist[j]+=exp(-dij*dij/(2*0.25)); totpi+=wlist[j]*pow(mdlist[j],-gamma); }
//             // mdlist[j] stores the logarithm of the KMC weight for point j
           // since the individual points contribution can be tiny, it is important to accumulate the logs
//           { 
//              dij=metric->dist(&LP(i-1,0),&HP(j,0),D); 
//              lij=dij*dij/(a2);
//              if (i==1) mdlist[j]=gamma*lij; // in the first step there is little to think. just weight in the 
//              else
//              {
//                 if (mdlist[j]<gamma*lij)
//                    mdlist[j]=mdlist[j]-gamma*log(1.0+exp(-lij+mdlist[j]/gamma));
//                 else
//                    mdlist[j]=gamma*lij-gamma*log(1.0+exp(lij-mdlist[j]/gamma));         
//             }
//              if (j==0 || mdlist[j]>lmax) lmax=mdlist[j];
////                  std::cerr<<"mdlist  "<<dij<<" :: "<<j<<" , "<<mdlist[j]<<"\n";
//           }
//           for (unsigned long j=0; j<N; ++j) totpi+=wlist[j]*exp(mdlist[j]-lmax);
//           std::cerr<<i<<"  : lmax= "<<lmax<< "  totpi "<< totpi<< "\n";
// this is for the -- deprecated -- exp(dist) method, which requires using log weights (commented section above)
//           {  selpi-=wlist[j]*exp(mdlist[j]-lmax); if(selpi<0.0) { j++; break; } }


           //picks a point randomly
           isunique = false;
           unsigned long j=0;
           while (!isunique) { 
               selpi=rndgen()*totpi;

               //a better search could be used here but hey....               

               for (j=0; j<N; ++j) 
               {  selpi-=wlist[j]*pow(mdlist[j],-gamma); if(selpi<0.0) { j++; break; } }           
               isel[i]=j-1;  
               isunique=true;
               if (funique) for (j=0; j<i; ++j) if (isel[i]==isel[j]) { isunique=false; break; }                
           }    
               
           LP.row(i)=HP.row(isel[i]);      

           std::cerr<<"selecting point "<<i<<" : "<<isel[i]<<"("<<wlist[j]*exp(mdlist[j]-lmax)/totpi<<")\n";
        }
    }
    else if (smode=="staged")        
    {
        // "staged" landmark selection. two rounds of selections are run. 
        // first, sqrt(nN) points are picked by farthest point sampling to get 
        // uniform coverage of the accessible areas. then, their weight is determined
        // by Voronoi polyhedra. then, a subset is picked at a -- possibly -- rescaled temperature
        
        FMatrix<double> MP;
        unsigned long m=sqrt(n*N);
        std::valarray<unsigned long> msel(m);

        MP.resize(m, D);         
        MP.row(0)=HP.row(maxj);
        for (unsigned long j=0; j<N; ++j) { mdlist[j]=metric->dist(&MP(0,0),&HP(j,0),D); }       
    
        std::cerr<<"now selecting "<<m<< " minmax points\n";
        for (unsigned long i=1; i<m; ++i)
        {
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
        unsigned long j=0;
        for (unsigned long i=0; i<n; ++i)
        {
            isunique=false;
            while (!isunique) {
               selpi=rndgen()*tw;
               for (j=0; j<m; ++j) 
               {  selpi-=mweights[j]; if(selpi<0.0) { j++; break; } }
               j--;

               std::cerr<<"picked "<<j<<" weight: "<<mweights[j]<<"\n";
               subsel=vplist[j][rndgen()*vplist[j].size()];
               
               isunique=true; 
               if (funique) for (unsigned long k=0; k<i; ++k) if (subsel==isel[k]) { isunique=false; break; }            
            }
           
           //actually, picks a random point from the voronoi polihedron
           std::cerr<<"  subpoint "<<subsel<<" selected "<< pow(mweights[j],wgamma/gamma)<<" \n";
           
           isel[i]=subsel; LP.row(i)=HP.row(subsel); 
           
           //so, mweights contain an estimate of the original probability density raised to power gamma.
           //now, we want to make it the probability to power wgamma, so must also correct for the original distortion
           weights[i]=pow(mweights[j],wgamma/gamma);
        }
        weights*=1.0/weights.sum();
    }
    else ERROR("Selection mode "<<smode<<" not implemented yet\n");



    if (fweight && smode!="staged") 
    {  
        double mind; unsigned long mini;
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
        double tw;
        for (unsigned long i=0; i<n; ++i) { weights[i]=pow(weights[i],wgamma); tw+=weights[i]; }
        weights*=1.0/tw;
        
//        }
//        else
//        {
//            for (unsigned long j=0; j<N; ++j)
//            {
//                mind=LD(0,j); mini=0;
//                for (unsigned long i=1; i<n; ++i)
//                {
//                    if (LD(i,j)<mind) { mind=LD(i,j); mini=i; }
//                }
//                weights[mini]+=wlist[j];
//            }
//        }
    }
    
    
    std::cout<<std::scientific; std::cout.precision(10); 
    std::cout<<"# "<<n<<" landmark points selected out of "<<N<<" and chosen by "<<smode<<"\n";
    std::cout<<"# Max distance detected: "<<maxd<<"\n";
    for (unsigned long i=0; i<n; ++i)
    {
        if (findex) std::cout<<isel[i]<<" ";
        for (unsigned long j=0; j<D; ++j) std::cout<<LP(i,j)<<" ";
        if (findex && finw) std::cout<<wlist[isel[i]]<<"  ";
        if (fweight) std::cout<<weights[i];
        
        std::cout<<"\n";
    }
}
