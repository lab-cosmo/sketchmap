#include "dimreduce.hpp"
#include "clparser.hpp"
#include "matrix-io.hpp"

using namespace toolbox;
void banner() 
{
    std::cerr
            << " USAGE: dimproj -D hi-dim -d low-dim -P hd-file -p ld-file [-pi period] [-w]    \n"
            << "               [-g1 grid1] [-g2 grid2] [-gw width] [-gt exp.factor] [-cgmin st] \n"
            << "               [-fun function] [-sigma sigma] [-h] [-print] < input             \n"
            << "                                                                                \n"
            << " computes the projection of the points given in input, given landmark points.   \n"
            << " dimension is set by -D option, and the projection is performed down to the     \n"
            << " dimensionality specified by -d. Optionally, high-dimensional data may be       \n"
            << " assumed to lie in a hypertoroidal space with period -pi.                       \n"
            << " If -cgmin is specified, a global minimization is attempted followed by st      \n"
            << " steps of conjugate gradient minimization.                                      \n"
            << " Optionally, function fun [identity, sigmoid, compress] can be applied to the   \n"
            << " distances, using the scale parameter sigma.                                    \n"
            << " Data must be provided in the files set by -P and -p and in input in the form:  \n"
            << " X1_1, X1_2, ... X1_D                                                           \n"
            << " X2_1, X2_2, ... X2_D                                                           \n"
            << " -print specifies that the cost functions for each point must be printed out    \n";
}

int main(int argc, char**argv)
{
    CLParser clp(argc,argv);
    unsigned long D,d,n,grid1,grid2,cgsteps; double gtemp, gwidth;
    double hdsigma, hdexpa, hdexpb; double ldsigma, ldexpa, ldexpb; 
    std::string fP, fp, fdist; bool fhelp, fprint, fweight;
    double peri, speri;
    bool fok=clp.getoption(D,"D",(unsigned long) 3) && 
            clp.getoption(d,"d",(unsigned long) 2) &&
            clp.getoption(fP,"P",std::string("")) &&
            clp.getoption(fp,"p",std::string("")) &&
            clp.getoption(grid1,"g1",(unsigned long) 21) &&
            clp.getoption(grid2,"g2",(unsigned long) 201) &&
            clp.getoption(cgsteps,"cgmin", (unsigned long) 0) &&
            clp.getoption(gwidth,"gw", 4.0) &&
            clp.getoption(gtemp,"gt", 1e-3) &&
            clp.getoption(hdsigma,"sigma",2.0) &&
            clp.getoption(hdexpa,"expa",4.0) &&
            clp.getoption(hdexpb,"expb",4.0) &&
            clp.getoption(ldsigma,"lsigma",-1.0) &&
            clp.getoption(ldexpa,"lexpa",-1.0) &&
            clp.getoption(ldexpb,"lexpb",-1.0) &&
            clp.getoption(fdist,"fun",std::string("identity")) &&
            clp.getoption(fhelp,"h",false) &&
            clp.getoption(fweight,"w",false) &&
            clp.getoption(fprint,"print",false) &&
            clp.getoption(peri,"pi",0.0) &&
            clp.getoption(speri,"spi",0.0);
    
    if (ldsigma<0) ldsigma=hdsigma;    if (ldexpa<0) ldexpa=hdexpa;    if (ldexpb<0) ldexpb=hdexpb;
    
    if (fhelp || !fok) { banner(); exit(1); }
    if (fP=="" || fp=="") ERROR("Hi-dim and low-dim points must be provided by the -P and -p options"); 
   
    std::ifstream sP(fP.c_str()), sp(fp.c_str()); 
    // reads points from standard input
    FMatrix<double> HP, lp;
    std::vector<std::vector<double> > plist; std::vector<double> point(D), pweight;
    while (sP.good())
    {
        double nw;
        for (int i=0; i<D; i++) sP>>point[i];
        if (fweight) sP>>nw; else nw=1;
        if (sP.good()) { plist.push_back(point); pweight.push_back(nw); }
    }
    
    HP.resize(plist.size(), D);
    for (unsigned long i=0; i<plist.size(); ++i) for (unsigned long j=0; j<D; ++j) HP(i,j)=plist[i][j];
    
    point.resize(d); plist.clear();
    while (sp.good())
    {
        for (int i=0; i<d; i++) sp>>point[i];
        if (sp.good()) plist.push_back(point);
    }
    
    lp.resize(plist.size(), d);
    for (unsigned long i=0; i<plist.size(); ++i) for (unsigned long j=0; j<d; ++j) lp(i,j)=plist[i][j];
    
    
    if ((n=lp.rows())!=HP.rows()) ERROR("HD and LD point list mismatch");
    
    NLDRProjection nlproj; NLDROptions opts;
    NLDRMetricPBC nperi; NLDRMetricEuclid neuclid;
    NLDRMetricSphere nsphere;
    nperi.periods.resize(D); nperi.periods=peri;
    nsphere.periods.resize(D); nsphere.periods=speri;
    
    if (peri==0.0 && speri==0.0) opts.nopts.ometric=&neuclid;
    else if (speri==0) { opts.nopts.ometric=&nperi; }
    else { opts.nopts.ometric=&nsphere; std::cerr<<"Spherical geodesic distances\n"; }
    
    std::valarray<double> tfpars;
    if (fdist=="identity")
    { tfpars.resize(0); opts.tfunH.set_mode(NLDRIdentity,tfpars); opts.tfunL=opts.tfunH; }
    else if (fdist=="compress")
    { 
        tfpars.resize(1); 
        tfpars[0]=hdsigma; opts.tfunH.set_mode(NLDRCompress,tfpars); 
        tfpars[0]=ldsigma; opts.tfunL.set_mode(NLDRCompress,tfpars); 
    }
    else if (fdist=="sigmoid")
    { 
        tfpars.resize(1); 
        tfpars[0]=hdsigma; opts.tfunH.set_mode(NLDRSigmoid,tfpars); 
        tfpars[0]=ldsigma; opts.tfunL.set_mode(NLDRSigmoid,tfpars); 
    }
    else if (fdist=="xsigmoid")
    {
        tfpars.resize(3); 
        tfpars[0]=hdsigma; tfpars[1]=hdexpa; tfpars[2]=hdexpb; opts.tfunH.set_mode(NLDRXSigmoid,tfpars); 
        tfpars[0]=ldsigma; tfpars[1]=ldexpa; tfpars[2]=ldexpb; opts.tfunL.set_mode(NLDRXSigmoid,tfpars); 
        std::cerr<<"XSIGMOID: "<<hdsigma<<" - "<<hdexpa<<" - "<<hdexpb<<
                ldsigma<<" - "<<ldexpa<<" - "<<ldexpb<<"\n";
    }
    else ERROR("Undefined map function for iterative distance matching");
    
    opts.grid1=grid1; opts.grid2=grid2; opts.gtemp=gtemp; opts.gwidth=gwidth; opts.cgsteps=cgsteps;
    nlproj.set_options(opts);
    std::valarray<double> nw(n); for (int i=0; i<n; i++)nw[i]=pweight[i];
    nlproj.set_points(HP,lp,nw);
    
    std::valarray<double> NP(D), PP(D), pp(d);
    unsigned long ip=0;
    while (std::cin.good())
    {
        std::cerr<<"Projecting "<<ip++<<"\n";
        for (int i=0; i<D; i++) std::cin>>NP[i];
        if (! std::cin.good()) break;
        double mind;
        if (fprint) nlproj.interp_out=std::string("interpolant.")+int2str(ip); 
        double perr=nlproj.project(NP, PP, pp, mind);
        
        for (int i=0; i<d; i++) std::cout<<pp[i]<<" "; std::cout<<perr<<" "<<mind<<"\n";
    }
}