/* Performs out-of-sample sketch-map embedding
   --------------------------------------------------
   Author: Michele Ceriotti, 2011
   Distributed under the GNU General Public License  
*/

#include "dimreduce.hpp"
#include "clparser.hpp"
#include "matrix-io.hpp"

using namespace toolbox;
void banner() 
{
    std::cerr
            << " USAGE: dimproj -D hi-dim -d low-dim -P hd-file -p ld-file [-pi period] [-w]    \n"
            << "               [-grid gw,g1,g2 ] [-cgmin st] [-gt temp]                         \n"
            << "               [-fun-hd s,a,b] [-fun-ld s,a,b] [-h] [-print]       < input      \n"
            << "                                                                                \n"
            << " computes the projection of the points given in input, given landmark points.   \n"
            << " dimension is set by -D option, and the projection is performed down to the     \n"
            << " dimensionality specified by -d. Optionally, high-dimensional data may be       \n"
            << " assumed to lie in a hypertoroidal space with period -pi.                       \n"
            << " A global minimization is performed on a grid ranging between -gw and +gw in d  \n"            
            << " dimensions. First, the stress function is computed on g1 points per dim, then  \n"
            << " an interpolated grid with g2 points is evaluated, and the min selected.        \n"                        
            << " If -gt is used, an exponential averaging is used instead of the min.           \n"                        
            << " If -cgmin is specified, the global minimization is followed by st steps of     \n"
            << " conjugate gradient minimization.                                               \n"
            << " Optionally, a sigmoid function can be applied in the D-dim space (-fun-hd)     \n"
            << " and/or in the d-dim space (-fun-ld). Both parameters must be followed by three \n"
            << " comma separated reals corresponding to sigma, a, b.                            \n"
            << " Data must be provided in the files set by -P and -p and in input in the form:  \n"
            << " X1_1, X1_2, ... X1_D                                                           \n"
            << " X2_1, X2_2, ... X2_D                                                           \n"
            << " If [-w] is present, an extra element per row is expected in the high-dimension \n"
            << " file, which is meant to be the weighting factor for the landmark point.        \n"
            << " -print specifies that the stress functions for each point must be printed out  \n";
}

int main(int argc, char**argv)
{
    CLParser clp(argc,argv);
    unsigned long D,d,n,cgsteps; double gtemp;
    std::string fP, fp, fdhd, fdld, gpars; bool fhelp, fprint, fweight;
    double peri, speri;
    bool fok=clp.getoption(D,"D",(unsigned long) 3) && 
            clp.getoption(d,"d",(unsigned long) 2) &&
            clp.getoption(fP,"P",std::string("")) &&
            clp.getoption(fp,"p",std::string("")) &&
            clp.getoption(gpars,"grid",std::string("1.0,21,201")) &&
            clp.getoption(cgsteps,"cgmin", (unsigned long) 0) &&
            clp.getoption(gtemp,"gt",0.0) &&            
            clp.getoption(fdhd,"fun-hd",std::string("identity")) &&
            clp.getoption(fdld,"fun-ld",std::string("identity")) &&
            clp.getoption(fhelp,"h",false) &&
            clp.getoption(fweight,"w",false) &&
            clp.getoption(fprint,"print",false) &&
            clp.getoption(peri,"pi",0.0) &&
            clp.getoption(speri,"spi",0.0);
        
    if (fhelp || !fok) { banner(); exit(1); }
    if (fP=="" || fp=="") ERROR("Hi-dim and low-dim points must be provided by the -P and -p options"); 
   
    std::ifstream sP(fP.c_str()), sp(fp.c_str()); 
    if (sP.fail()) ERROR("Unable to open high-dim file.");
    if (sp.fail()) ERROR("Unable to open low-dim file.");

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
    std::valarray<double> fhdpars(0.0,3), fldpars(0.0,3), fgrid(0.0,3);    
    if (fdhd=="identity")
    { tfpars.resize(0); opts.tfunH.set_mode(NLDRIdentity,tfpars); }
    else 
    {
      csv2floats(fdhd,tfpars);  fhdpars=tfpars; 
      std::cerr<<"high-dim pars"<<tfpars<<"\n";     
      if (tfpars.size()==2) 
      {
        opts.tfunH.set_mode(NLDRGamma,tfpars);       
      }
      else if (tfpars.size()==3) 
      {
        opts.tfunH.set_mode(NLDRXSigmoid,tfpars);     
      }
      else
      {  ERROR("-fun-hd argument must be of the form sigma,a,b or sigma,n");  }
    }
    
    if (fdld=="identity")
    { tfpars.resize(0); opts.tfunL.set_mode(NLDRIdentity,tfpars); }
    else 
    {
      csv2floats(fdld,tfpars); fldpars=tfpars;      
      std::cerr<<"lo-dim pars"<<tfpars<<"\n";     
      if (tfpars.size()==2) 
      {
        opts.tfunL.set_mode(NLDRGamma,tfpars);  
      }
      else if (tfpars.size()==3) 
      {
        opts.tfunL.set_mode(NLDRXSigmoid,tfpars);
      }
      else
      {  ERROR("-fun-ld argument must be of the form sigma,a,b or sigma,n");  }
    }
    
    csv2floats(gpars,tfpars); if (tfpars.size()<3) ERROR("-grid argument requires gw,g1,g2")    
    opts.grid1=tfpars[1]; opts.grid2=tfpars[2]; opts.gwidth=tfpars[0]; opts.gtemp=gtemp; opts.cgsteps=cgsteps;
    
    nlproj.set_options(opts);
    std::valarray<double> nw(n); for (int i=0; i<n; i++)nw[i]=pweight[i];
    nlproj.set_points(HP,lp,nw);
    
    std::valarray<double> NP(D), PP(D), pp(d);
    std::cout.precision(12); std::cout.setf(std::ios::scientific);
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
