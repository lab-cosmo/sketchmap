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
    double hdsigma, hdexpa, hdexpb; double ldsigma, ldexpa, ldexpb;
    double sat1, sat2, irnd, imix;
    unsigned long D,d,dts,nn, nsteps, pluneigh; double sm, neps,peri,speri; bool fverb, fveryverb, fplumed, fhelp,  fweight, fcenter;
    unsigned long g1,g2; double gw;
    std::string nlrmode, fmds, finit, itermode;
    bool fok=clp.getoption(D,"D",(unsigned long) 3) && 
            clp.getoption(d,"d",(unsigned long) 2) &&
            clp.getoption(nlrmode,"mode",std::string("LLE")) &&  
            clp.getoption(peri,"pi",0.0) &&
            clp.getoption(speri,"spi",0.0) &&
            clp.getoption(fverb,"v",false) &&  
            clp.getoption(fweight,"w",false) &&
            clp.getoption(fveryverb,"vv",false) &&
            clp.getoption(fhelp,"h",false) &&
            clp.getoption(fplumed,"plumed",false) && 
            clp.getoption(fcenter,"center",false) && 
            clp.getoption(hdsigma,"sigma",2.0) &&
            clp.getoption(hdexpa,"expa",4.0) &&
            clp.getoption(hdexpb,"expb",4.0) &&
            clp.getoption(ldsigma,"lsigma",-1.0) &&
            clp.getoption(ldexpa,"lexpa",-1.0) &&
            clp.getoption(ldexpb,"lexpb",-1.0) &&
            clp.getoption(fmds,"fun",std::string("identity")) &&
            clp.getoption(itermode,"imode",std::string("conjgrad")) &&
            clp.getoption(sat1,"sat1",0.1) &&
            clp.getoption(sat2,"sat2",1e-8) &&
            clp.getoption(nsteps,"steps",(unsigned long) 100) &&
            clp.getoption(finit,"init",std::string("")) &&
            clp.getoption(irnd,"randomize",0.0) &&
            clp.getoption(imix,"imix",0.0) &&
            clp.getoption(g1,"g1",(unsigned long) 11) &&
            clp.getoption(g2,"g2",(unsigned long) 101) &&
            clp.getoption(gw,"gw",20.0) &&
            clp.getoption(dts,"dts",(unsigned long) 0) &&
            clp.getoption(sm,"smooth",-1e-3) &&  
            clp.getoption(nn,"neigh",(unsigned long) 4) &&
            clp.getoption(neps,"ncut",0.0) && 
            clp.getoption(pluneigh,"nplumed",(unsigned long) 10);

    if (ldsigma<0) ldsigma=hdsigma;    if (ldexpa<0) ldexpa=hdexpa;    if (ldexpb<0) ldexpb=hdexpb;
    
    if (fhelp || !fok) { banner(); exit(1); }
    if (dts==0) dts=d; 
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
    NLDRNeighborOptions nopts; nopts.greediness=NLDRAsym; nopts.maxneigh=nn; nopts.cutoff=neps; 
    NLDRMetricPBC nperi; NLDRMetricEuclid neuclid; NLDRMetricSphere nsphere;
    nperi.periods.resize(D); nperi.periods=peri;
    nsphere.periods.resize(D); nsphere.periods=speri;
    
    if (peri==0.0 && speri==0.0) nopts.ometric=&neuclid;
    else if (speri==0) { nopts.ometric=&nperi; }
    else { nopts.ometric=&nsphere; std::cerr<<"spherical distances\n"; }
        
    NLDRLLEReport llereport;
    NLDRLLEOptions lleopts; lleopts.nlopts=nopts; lleopts.smooth=sm;  lleopts.lowdim=d; lleopts.dimts=dts;
    lleopts.rmlonesome=true; lleopts.verbose=fveryverb; 
    
    NLDRMDSReport mdsreport;
    NLDRMDSOptions mdsopts; mdsopts.lowdim=d; mdsopts.verbose=fveryverb;
    if (peri==0.0 && speri==0.0) mdsopts.metric=&neuclid;
    else if (speri==0) { mdsopts.metric=&nperi; }
    else { mdsopts.metric=&nsphere; }
    
    NLDRITEROptions iteropts;
    NLDRITERReport iterreport;
    std::valarray<double> tfpars;
    iteropts.lowdim=d; iteropts.verbose=fveryverb; iteropts.steps=nsteps;
    if (fmds=="identity")
    { tfpars.resize(0); iteropts.tfunH.set_mode(NLDRIdentity,tfpars); iteropts.tfunL=iteropts.tfunH; }
    else if (fmds=="compress")
    { 
        tfpars.resize(1); 
        tfpars[0]=hdsigma; iteropts.tfunH.set_mode(NLDRCompress,tfpars); 
        tfpars[0]=ldsigma; iteropts.tfunL.set_mode(NLDRCompress,tfpars); 
    }
    else if (fmds=="sigmoid")
    { 
        tfpars.resize(1); 
        tfpars[0]=hdsigma; iteropts.tfunH.set_mode(NLDRSigmoid,tfpars); 
        tfpars[0]=ldsigma; iteropts.tfunL.set_mode(NLDRSigmoid,tfpars); 
    }
    else if (fmds=="xsigmoid")
    {
        tfpars.resize(3); 
        tfpars[0]=hdsigma; tfpars[1]=hdexpa; tfpars[2]=hdexpb; iteropts.tfunH.set_mode(NLDRXSigmoid,tfpars); 
        tfpars[0]=ldsigma; tfpars[1]=ldexpa; tfpars[2]=ldexpb; iteropts.tfunL.set_mode(NLDRXSigmoid,tfpars); 
        std::cerr<<"XSIGMOID: "<<hdsigma<<" - "<<hdexpa<<" - "<<hdexpb<<
                ldsigma<<" - "<<ldexpa<<" - "<<ldexpb<<"\n";
    }
    else ERROR("Undefined map function for iterative distance matching");
    
    std::cerr<<iteropts.tfunH.f(5.0)<<","<<iteropts.tfunH.df(5.0)<<"  test function\n";
    iteropts.metric=mdsopts.metric;
    iteropts.grid1=g1; iteropts.grid2=g2; iteropts.gridw=gw; 
    
    /*
    std::valarray<double> nr(3), nz(3); nr=1.0; nz=0.0;
    std::cerr<<lleopts.nlopts.ometric->dist(nr,nz)<<" euclid\n";
    return 0;*/
    std::cerr<<"Initialization done, running dim. reduction\n";
    
    enum { mMDS, mLLE, mITER } modes;
    if (nlrmode=="LLE") { modes=mLLE; lleopts.mode=LLE;}
    else if (nlrmode=="LLTE")  { modes=mLLE; lleopts.mode=LLTE; }
    else if (nlrmode=="HLLE")  { modes=mLLE; lleopts.mode=HLLE; }
    else if (nlrmode=="WLLE")  { modes=mLLE; lleopts.mode=LLE; lleopts.nlopts.kw=nn; }
    else if (nlrmode=="WHLLE") { modes=mLLE; lleopts.mode=HLLE; lleopts.nlopts.kw=nn; }
    else if (nlrmode=="WLLTE") { modes=mLLE; lleopts.mode=LLTE; lleopts.nlopts.kw=nn; }
    else if (nlrmode=="MDS")   { modes=mMDS; mdsopts.mode=MDS; }
    else if (nlrmode=="SMDS")  { modes=mMDS; mdsopts.mode=SMDS; }
    else if (nlrmode=="TMDS")  { modes=mMDS; mdsopts.mode=TMDS; }
    else if (nlrmode=="IMDS")  { modes=mITER;  iteropts.global=false; }
    else if (nlrmode=="GMDS")  { modes=mITER;  iteropts.global=true; }
    else ERROR("Unsupported NLDR mode. Use LLE or LLTE.");
    
    if (itermode=="conjgrad") iteropts.minmode=NLDRCGradient;
    else if (itermode=="simplex") iteropts.minmode=NLDRSimplex;
    else if (itermode=="anneal") iteropts.minmode=NLDRAnnealing;
    iteropts.saopts.temp_init=sat1; iteropts.saopts.temp_final=sat2;
    iteropts.weights.resize(weights.size()); for (unsigned long i=0; i<weights.size();++i) iteropts.weights[i]=weights[i]; iteropts.imix=imix;
    
    if (finit!="")
    {
        //reads initial values of LD points(might be just useless, unless iterative method is requested)
        iteropts.ipoints.resize(mpoints.rows(),d);
        std::ifstream fip(finit.c_str());
        for (unsigned long i=0; i<mpoints.rows(); i++)
            for (unsigned long j=0; j<d; j++) fip>>iteropts.ipoints(i,j);
        RndGaussian<double> prng;
        if (irnd>0) for (unsigned long i=0; i<mpoints.rows(); i++)
            for (unsigned long j=0; j<d; j++) iteropts.ipoints(i,j)+=prng()*irnd;
    }
    
    if (modes==mLLE) NLDRLLE(mpoints,nlproj,lleopts,llereport);
    else if (modes==mMDS) NLDRMDS(mpoints,nlproj,mdsopts,mdsreport);
    else if (modes==mITER) NLDRITER(mpoints,nlproj,iteropts,iterreport);
    
    nlproj.get_points(hplist,lplist);
    if (fplumed)
    {
        std::cout << "NLANDMARKS " <<lplist.size()<<"\n\n";
        //std::cout << "SIGMOID "<<hdsigma<<"\n\n"; //!TODO output for sigma in low-dim as well
        std::cout<<"LOW_D_FUNCTION TYPE "<<
                (fmds=="identity"?"distance":
                 fmds=="compress"?"compress":
                 fmds=="sigmoid"?"sigmoid":
                 fmds=="xsigmoid"?"general":"unknown"
                )<<" SIGMA "<<ldsigma<<" POWERS "<<ldexpa<<" "<<ldexpb<<"\n";
        std::cout<<"HIGH_D_FUNCTION TYPE "<<
                (fmds=="identity"?"distance":
                 fmds=="compress"?"compress":
                 fmds=="sigmoid"?"sigmoid":
                fmds=="xsigmoid"?"general":"unknown"
                )<<" SIGMA "<<hdsigma<<" POWERS "<<hdexpa<<" "<<hdexpb<<"\n";

        std::cout << "LIMITS> \n"<<-gw<<" "<<gw<<"\n"<<-gw<<" "<<gw<<"\nLIMITS<\n";
        
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
        if (modes==mLLE) 
        {
            std::cout << " ######################## LLE REPORT ###################\n";
            std::cout << " # Error in fitting HD points: "<<llereport.hd_error<<"\n";
            std::cout << " # Small Eigenvalues of M: \n # ";
            for (int i=0; i<llereport.deval.size(); ++i) std::cout<<llereport.deval[i]<<" "; 
            if (fveryverb) std::cout<<"("<<llereport.dp1eval<<")"; std::cout <<"\n";
            std::cout << " # Error in fitting LD points: "<<llereport.ld_error<<"\n";
            std::cout << " # y1 .. yd "<<(fveryverb?" hd_error ld_error ":"")<<"\n";
        }
        else if (modes==mMDS)
        {
            std::cout << " ######################## MDS REPORT ###################\n";
            std::cout << " # Large Eigenvalues of M: \n # ";
            for (int i=0; i<mdsreport.deval.size(); ++i) std::cout<<mdsreport.deval[i]<<" "; 
            if (fveryverb) std::cout<<"("<<mdsreport.dp1eval<<")"; std::cout <<"\n";
            std::cout << " # Error in fitting LD points: "<<mdsreport.ld_error<<"\n";
            std::cout << " # y1 .. yd "<<(fveryverb?" ld_error ":"")<<"\n";
        }
        else if (modes==mITER)
        {
            std::cout << " ################### ITERATIVE"<<(iteropts.global?" gMDS REPORT #############\n":" MDS REPORT ##############\n");
            std::cout << " # Computed with function: "<<fmds<<" with pars ?????\n";
            std::cout << " # Conjugate gradient steps: "<<nsteps<<"\n";
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
            if (modes==mLLE) std::cout<<llereport.hd_errors[i]<<" "<<llereport.ld_errors[i]<<" ";
            else if (modes==mMDS) std::cout<<mdsreport.ld_errors[i]<<" ";
            else if (modes==mITER) std::cout<<iterreport.ld_errors[i]<<" ";
        }
        std::cout<<std::endl;
    }
    }
    return 0;
}
