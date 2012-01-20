/* Program to perform sketch-map nonlinear dimensionality reduction
   --------------------------------------------------
   Author: Michele Ceriotti, 2011
   Distributed under the GNU General Public License  
*/

#include "dimreduce.hpp"
#include "clparser.hpp"
#include "matrix-io.hpp"
#include "matrix-conv.hpp"

using namespace toolbox;

void banner() 
{
    std::cerr
            << " USAGE: dimred -D hi-dim -d low-dim -pi period [-v|-vv] [-h] [-w] [-init file]  \n"
            << "               [-center] [-plumed] [-fun-hd s,a,b] [-fun-ld s,a,b] [-imix mix]  \n"
            << "               [-preopt steps] [-grid gw,g1,g2] [-gopt steps]                   \n"            
            << "                                                                                \n"
            << " compute the dimensionality reduction of data points given in input. The high   \n"
            << " dimension is set by -D option, and the projection is performed down to the     \n"
            << " dimensionality specified by -d. Optionally, high-dimensional data may be       \n"
            << " assumed to lie in a hypertoroidal space with period -pi.                       \n"
            << " Data must be provided in input in the format                                   \n"
            << " X1_1, X1_2, ... X1_D [w1]                                                      \n"
            << " X2_1, X2_2, ... X2_D [w2]                                                      \n"
            << " where wi's are optional weights to be given if -w is chosen.                   \n"
            << " Verbosity of output is controlled by -v and -vv options, and optionally        \n"
            << " output can be made compatible with the PLUMED implementation of bespoke CVs    \n"
            << " by the -plumed option. -center weight-centers points around the origin.        \n"
            << " The mode of operation is that initial low-dim points are loaded from -init     \n"
            << " file. If absent, multi-dimensional scaling is performed to get starting pos.   \n"
            << " Then, iterative optimization starts, for -preopt steps of conjugate gradient.  \n"
            << " The stress function is given by chi=mix*chi_id+(1-mix) chi_fun where chi_id    \n"
            << " is the quadratic discrepancy of distances, and chi_fun is computed applying    \n"
            << " -fun-hd and -fun-ld (both default to identity, otherwise sigma,a,b must be     \n"
            << " given, which control the shape of the sigmoid function.                        \n"
            << " If -grid is given, after this pre-optimization a pointwise-global optimizer    \n"
            << " is run, which minimized one point at a time on a grid over [-gw:gw], with g1   \n"
            << " points on the coarse grid and g2 points on the fine grid. -gopt steps of CG    \n"
            << " optimizer are then performed.                                                  \n";
}

int main(int argc, char**argv)
{
    CLParser clp(argc,argv);    
    double sat1, sat2, irnd, imix;
    unsigned long D,d,dts,nn, presteps, gsteps, pluneigh; 
    double sm, neps,peri,speri; bool fverb, fveryverb, fplumed, fhelp,  fweight, fcenter;
    std::string fmds, finit, fdhd, fdld, gpars, itermode;
    bool fok=clp.getoption(D,"D",(unsigned long) 3) && 
            clp.getoption(d,"d",(unsigned long) 2) &&
            clp.getoption(peri,"pi",0.0) &&
            clp.getoption(speri,"spi",0.0) &&
            clp.getoption(fverb,"v",false) &&  
            clp.getoption(fweight,"w",false) &&
            clp.getoption(fveryverb,"vv",false) &&
            clp.getoption(fhelp,"h",false) &&
            clp.getoption(fplumed,"plumed",false) && 
            clp.getoption(fcenter,"center",false) && 
            clp.getoption(gpars,"grid",std::string("")) &&
            clp.getoption(gsteps,"gopt", (unsigned long) 0) &&
            clp.getoption(presteps,"preopt",(unsigned long) 0) &&
            clp.getoption(fdhd,"fun-hd",std::string("identity")) &&
            clp.getoption(fdld,"fun-ld",std::string("identity")) &&
            clp.getoption(itermode,"imode",std::string("conjgrad")) &&
            clp.getoption(finit,"init",std::string("")) &&
            clp.getoption(irnd,"randomize",0.0) &&
            clp.getoption(imix,"imix",0.0) &&
            clp.getoption(sat1,"sat1",0.1) &&
            clp.getoption(sat2,"sat2",1e-8) &&

            clp.getoption(sm,"smooth",-1e-3) &&  
            clp.getoption(nn,"neigh",(unsigned long) 4) &&
            clp.getoption(neps,"ncut",0.0) && 
            clp.getoption(pluneigh,"nplumed",(unsigned long) 10);

    if (fhelp || !fok) { banner(); exit(1); }
    std::vector<std::vector<double> > plist; std::vector<double> point(D), weights;
    
    // reads points from standard input
    while (std::cin.good())
    {
        double nw;
        for (int i=0; i<D; i++) std::cin>>point[i];
        if (fweight) std::cin>>nw; else nw=1.0; 
        if (std::cin.good()) { plist.push_back(point); weights.push_back(nw); }
    }

    std::valarray<std::valarray<double> > hplist, lplist; 
    FMatrix<double> mpoints(plist.size(),D);
    for (int i=0; i<plist.size(); i++) for (int j=0; j<D; j++) mpoints(i,j)=plist[i][j];

    NLDRProjection nlproj;
    NLDRMetricPBC nperi; NLDRMetricEuclid neuclid; NLDRMetricSphere nsphere;
    nperi.periods.resize(D); nperi.periods=peri;
    nsphere.periods.resize(D); nsphere.periods=speri;
            
    NLDRMDSReport mdsreport;
    NLDRMDSOptions mdsopts; mdsopts.lowdim=d; mdsopts.verbose=fveryverb;
    if (peri==0.0 && speri==0.0) mdsopts.metric=&neuclid;
    else if (speri==0) { mdsopts.metric=&nperi; }
    else { mdsopts.metric=&nsphere; }
    
    NLDRITEROptions iteropts;
    NLDRITERReport iterreport;
    std::valarray<double> tfpars, fhdpars(0.0,3), fldpars(0.0,3), fgrid(0.0,3);
    iteropts.lowdim=d; iteropts.verbose=fveryverb; iteropts.metric=mdsopts.metric;

    if (fdhd=="identity")
    { tfpars.resize(0); iteropts.tfunH.set_mode(NLDRIdentity,tfpars); }
    else 
    {
      csv2floats(fdhd,tfpars); if (tfpars.size()<3) ERROR("-fun-hd argument must be of the form sigma,a,b")
      std::cerr<<"high-dim pars"<<tfpars<<"\n";
      iteropts.tfunH.set_mode(NLDRXSigmoid,tfpars);  fhdpars=tfpars;
    }
    
    if (fdld=="identity")
    { tfpars.resize(0); iteropts.tfunL.set_mode(NLDRIdentity,tfpars); }
    else 
    {
      csv2floats(fdld,tfpars); if (tfpars.size()<3) ERROR("-fun-ld argument must be of the form sigma,a,b")
      std::cerr<<"lo-dim pars"<<tfpars<<"\n";      
      iteropts.tfunL.set_mode(NLDRXSigmoid,tfpars);  fldpars=tfpars; 
    }
    
    bool doglobal;
    if (gpars!="")
    {
      csv2floats(gpars,tfpars); if (tfpars.size()<3) ERROR("-grid argument requires gw,g1,g2");    
      fgrid=tfpars;
      iteropts.grid1=tfpars[1]; iteropts.grid2=tfpars[2]; iteropts.gridw=tfpars[0];  doglobal=true;
    }  else doglobal=false;
    
    std::cerr<<"Initialization done, running dim. reduction\n";
    
    if (itermode=="conjgrad") iteropts.minmode=NLDRCGradient;
    else if (itermode=="simplex") iteropts.minmode=NLDRSimplex;
    else if (itermode=="anneal") iteropts.minmode=NLDRAnnealing;
    iteropts.saopts.temp_init=sat1; iteropts.saopts.temp_final=sat2;
    iteropts.weights.resize(weights.size()); for (unsigned long i=0; i<weights.size();++i) iteropts.weights[i]=weights[i]; iteropts.imix=imix;
    
    iteropts.ipoints.resize(mpoints.rows(),d);
    if (finit!="")
    {
        //reads initial values of LD points(might be just useless, unless iterative method is requested)
        std::ifstream fip(finit.c_str());
        for (unsigned long i=0; i<mpoints.rows(); i++)
            for (unsigned long j=0; j<d; j++) fip>>iteropts.ipoints(i,j);
        RndGaussian<double> prng;
        if (irnd>0) for (unsigned long i=0; i<mpoints.rows(); i++)
            for (unsigned long j=0; j<d; j++) iteropts.ipoints(i,j)+=prng()*irnd;
    }
    else
    {
      //initialize from classical MDS
      NLDRMDS(mpoints,nlproj,mdsopts,mdsreport);
      nlproj.get_points(hplist,lplist);
      
      for (unsigned long i=0; i<mpoints.rows(); i++)
         for (unsigned long j=0; j<d; j++) iteropts.ipoints(i,j)=lplist[i][j];
    }
    
    iteropts.global=false; iteropts.steps=presteps;
    if (presteps>0) 
    {  
      NLDRITER(mpoints,nlproj,iteropts,iterreport);
      nlproj.get_points(hplist,lplist);  
    }
        
    if (doglobal)
    {
      iteropts.global=true; iteropts.steps=gsteps; 
      for (unsigned long i=0; i<mpoints.rows(); i++)
         for (unsigned long j=0; j<d; j++) iteropts.ipoints(i,j)=lplist[i][j];      
      NLDRITER(mpoints,nlproj,iteropts,iterreport);
    }
    
    if (fplumed)
    {
        std::cout << "NLANDMARKS " <<lplist.size()<<"\n\n";
        //std::cout << "SIGMOID "<<hdsigma<<"\n\n"; //!TODO output for sigma in low-dim as well
        std::cout<<"LOW_D_FUNCTION TYPE "<<
                (fdld=="identity"?"distance":"xsigmoid"
                )<<" SIGMA "<<fldpars[0]<<" POWERS "<<fldpars[1]<<" "<<fldpars[2]<<"\n";
        std::cout<<"HIGH_D_FUNCTION TYPE "<<
                (fdhd=="identity"?"distance":"xsigmoid"
                )<<" SIGMA "<<fhdpars[0]<<" POWERS "<<fhdpars[1]<<" "<<fhdpars[2]<<"\n";

        std::cout << "LIMITS> \n"<<-fgrid[0]<<" "<<fgrid[0]<<"\n"<<-fgrid[0]<<" "<<fgrid[0]<<"\nLIMITS<\n";
        
        std::cout << "HIGH_D>\n";
        for (int i=0; i<hplist.size(); i++)
        { for (int h=0; h<D; h++) std::cout <<hplist[i][h]<<" "; std::cout<<"\n"; }
        std::cout << "HIGH_D<\n\n";
        std::cout << "LOW_D>\n";
        for (int i=0; i<lplist.size(); i++)
        { for (int h=0; h<d; h++) std::cout <<lplist[i][h]<<" "; std::cout<<"\n"; }
        std::cout << "LOW_D<\n\n";
        if (fweight)
        {
            std::cout << "WEIGHTS>\n";
            for (int i=0; i<lplist.size(); i++)
            { std::cout<<weights[i]<<"\n"; }
            std::cout << "WEIGHTS<\n\n";
        }
        /*
        std::cerr<<"Computing " <<pluneigh<<" Plumed neighbours\n";
        FMatrix<double> MP; 
        std::cerr<<"matrix conversion\n";
        MP=hplist; 
        std::cerr<<"compute neighbors\n";
        NLDRNeighborOptions pnopts; pnopts.maxneigh=pluneigh;  pnopts.greediness=NLDRAsym;
        pnopts.ometric=nopts.ometric;
        NLDRNeighborList nlist(MP,pnopts);
        
        std::cerr<<"printing neighbors\n";
        std::cout<<"MAXNEIGHBOURS "<<pnopts.maxneigh<<"\n";
        
        std::cout << "NEIGHBOURS>\n";
        for (int i=0; i<lplist.size(); i++) 
        { for (int h=0; h<nlist.nneigh(i); h++) std::cout <<nlist.index(i,h)<<" "; std::cout<<"\n"; }
        std::cout << "NEIGHBOURS<\n";
        */
    }
    else 
    {
    if (fverb || fveryverb)
    {
        if (presteps==0 && ! doglobal)
        {
            std::cout << " ######################## MDS REPORT ###################\n";
            std::cout << " # Large Eigenvalues of M: \n # ";
            for (int i=0; i<mdsreport.deval.size(); ++i) std::cout<<mdsreport.deval[i]<<" "; 
            if (fveryverb) std::cout<<"("<<mdsreport.dp1eval<<")"; std::cout <<"\n";
            std::cout << " # Error in fitting LD points: "<<mdsreport.ld_error<<"\n";
            std::cout << " # y1 .. yd "<<(fveryverb?" ld_error ":"")<<"\n";
        }
        else 
        {
            std::cout << " ################### ITERATIVE"<<(doglobal?" gMDS REPORT #############\n":" MDS REPORT ##############\n");
            std::cout << " # Computed with function: "<<fmds<<" with pars ?????\n";
            std::cout << " # Conjugate gradient steps: "<<presteps+gsteps<<"\n";
            std::cout << " # Error in fitting LD points: "<<iterreport.ld_error<<"\n";
            std::cout << " # y1 .. yd "<<(fveryverb?" ld_error ":"")<<"\n";
        }
    }
    
    std::valarray<double> com(d); com=0.0;
    if (fcenter)
    { for (int i=0; i<lplist.size(); i++) com+=lplist[i]; com*=1.0/lplist.size(); }
    for (int i=0; i<lplist.size(); i++)
    {
        for (int h=0; h<d; h++)  std::cout<<lplist[i][h]-com[h]<<" ";
        if (fveryverb) 
        {
            if (presteps==0 && ! doglobal) std::cout<<mdsreport.ld_errors[i]<<" ";
            else std::cout<<iterreport.ld_errors[i]<<" ";
        }
        std::cout<<std::endl;
    }
    }
    return 0;
}
