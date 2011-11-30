#include "dimreduce.hpp"
#include "matrix-io.hpp"
#include "matrix-conv.hpp"
#include "tensor-full.hpp"
#include "interpol.hpp"
#include "linalg.hpp"

#ifdef ARPACK
namespace tblapack {
extern "C" {
#define dsaupd dsaupd_
#define dseupd dseupd_
    int dsaupd( int* IDO, char* BMAT, int*N, char* WHICH, int* NEV, 
             double* TOL, double *RESID, int* NCV, double* V, int *LDV, 
             int* IPARAM, int* IPNTR, double *WORKD, 
             double *WORKL, int*LWORKL, int*INFO );
    int dseupd( int* RVEC, char* HOWMNY, int *SELECT, 
                double *D, double *Z, int *LDZ, double *SIGMA, 
                char* BMAT, int*N, char* WHICH, int* NEV, 
                double* TOL, double *RESID, int* NCV, double* V, int *LDV, 
                int* IPARAM, int* IPNTR, double *WORKD, 
                double *WORKL, int*LWORKL, int*INFO);
    
}};
#endif

namespace toolbox {
#define __NLDR_SIGMA

double NLDRFunction::nldr_identity(double x) const { return x; }
double NLDRFunction::nldr_didentity(double x) const { return 1.0; }
void NLDRFunction::nldr_identity(double x, double &rf, double& rdf) const { rf=x; rdf=1.0; }

inline double NLDRFunction::nldr_sigmoid(double x) const 
{ double sx=x*pars[0];  return 1.0-1.0/(1.0+sx*sx); }
inline double NLDRFunction::nldr_dsigmoid(double x) const 
{  double sx=x*pars[0]; sx=1.0/(1.0+sx*sx); return x*(sx*sx)*pars[1]; }
inline void NLDRFunction::nldr_sigmoid(double x, double &rf, double& rdf) const  
{  double sx=x*pars[0]; sx=1.0/(1.0+sx*sx); rf=1.0-sx; rdf=x*(sx*sx)*pars[1]; }

inline double NLDRFunction::nldr_compress(double x) const
{  double sx=x*pars[0];  return 1.0-1.0/(1.0+sx);  }
inline double NLDRFunction::nldr_dcompress(double x) const
{  double sx=x*pars[0];   sx=1.0/(1.0+sx);  return (sx*sx)*pars[0]; }
inline void NLDRFunction::nldr_compress(double x, double &rf, double& rdf)  const
{  double sx=x*pars[0];   sx=1.0/(1.0+sx); rf=1.0-sx;  rdf=(sx*sx)*pars[0]; }

//this must return 1-(1+(2^(a/b)-1)(x/s)^a)^-b/a, needs cleaning up (precompute constants)
inline double NLDRFunction::nldr_xsigmoid(double x) const
{ double sx=x*pars[0];  return 1.0-pow(1.0+pars[1]*pow(sx,pars[2]),pars[4]); }
inline double NLDRFunction::nldr_dxsigmoid(double x) const
{ double sx=x*pars[0];  sx=pars[1]*pow(sx,pars[2]);  return pars[3]*sx/x*pow(1.+sx,pars[4]-1.0); }
inline void NLDRFunction::nldr_xsigmoid(double x, double &rf, double& rdf) const
{ double sx=x*pars[0];  sx=pars[1]*pow(sx,pars[2]); rf=pow(1.0+sx,pars[4]); rdf=pars[3]*sx/x*rf/(1.0+sx); rf=1.0-rf; }


NLDRFunction::NLDRFunction(const NLDRFunction& nf)
{
    pmode=nf.pmode;
    pars.resize(nf.pars.size()); pars=nf.pars;
    switch(pmode)
    { 
        case NLDRIdentity: pf=&NLDRFunction::nldr_identity; pdf=&NLDRFunction::nldr_didentity; pfdf=&NLDRFunction::nldr_identity;  break;
        case NLDRCompress: pf=&NLDRFunction::nldr_compress; pdf=&NLDRFunction::nldr_dcompress; pfdf=&NLDRFunction::nldr_compress;  break;
        case NLDRSigmoid: pf=&NLDRFunction::nldr_sigmoid; pdf=&NLDRFunction::nldr_dsigmoid; pfdf=&NLDRFunction::nldr_sigmoid;  break;
        case NLDRXSigmoid: pf=&NLDRFunction::nldr_xsigmoid; pdf=&NLDRFunction::nldr_dxsigmoid; pfdf=&NLDRFunction::nldr_xsigmoid;  break;
    }
}

NLDRFunction& NLDRFunction::operator=(const NLDRFunction& nf)
{
    if (&nf==this) return *this;
    pmode=nf.pmode;
    pars.resize(nf.pars.size()); pars=nf.pars;
    switch(pmode)
    { 
        case NLDRIdentity: pf=&NLDRFunction::nldr_identity; pdf=&NLDRFunction::nldr_didentity; pfdf=&NLDRFunction::nldr_identity;  break;
        case NLDRCompress: pf=&NLDRFunction::nldr_compress; pdf=&NLDRFunction::nldr_dcompress; pfdf=&NLDRFunction::nldr_compress;  break;
        case NLDRSigmoid: pf=&NLDRFunction::nldr_sigmoid; pdf=&NLDRFunction::nldr_dsigmoid; pfdf=&NLDRFunction::nldr_sigmoid;  break;
        case NLDRXSigmoid: pf=&NLDRFunction::nldr_xsigmoid; pdf=&NLDRFunction::nldr_dxsigmoid; pfdf=&NLDRFunction::nldr_xsigmoid;  break;
    }
    return *this;
}
   
void NLDRFunction::set_mode(NLDRFunctionMode mode, const std::valarray<double>& npars)
{
    switch (mode) 
    {
    case NLDRIdentity:
        pars.resize(0); 
        pf=&NLDRFunction::nldr_identity; pdf=&NLDRFunction::nldr_didentity; pfdf=&NLDRFunction::nldr_identity;
        break;
    case NLDRCompress: //computes  1-1/(1+x/npars[0])
        if (npars.size()!=1) ERROR("Wrong number of parameters for Compress transfer function");
        pars.resize(1); pars[0]=1.0/npars[0];
        pf=&NLDRFunction::nldr_compress; pdf=&NLDRFunction::nldr_dcompress; pfdf=&NLDRFunction::nldr_compress;
        break;
    case NLDRSigmoid: //computes  1-1/(1+(x/npars[0])^2)
        if (npars.size()!=1) ERROR("Wrong number of parameters for Sigmoid transfer function");
        pars.resize(2); pars[0]=1.0/npars[0]; pars[1]=2.0*pars[0]*pars[0];
        pf=&NLDRFunction::nldr_sigmoid; pdf=&NLDRFunction::nldr_dsigmoid; pfdf=&NLDRFunction::nldr_sigmoid;
        break;
    case NLDRXSigmoid:  //this must return 1-(1+(2^(a/b)-1)(x/s)^a)^-b/a s=np[0] a=np[1] b=np[2]
        if (npars.size()!=3) ERROR("Wrong number of parameters for XSigmoid transfer function");
        pars.resize(5); pars[0]=1.0/npars[0]; pars[1]=pow(2.,npars[1]/npars[2])-1.0; pars[2]=npars[1]; pars[3]=npars[2]; pars[4]=-npars[2]/npars[1];
        pf=&NLDRFunction::nldr_xsigmoid; pdf=&NLDRFunction::nldr_dxsigmoid; pfdf=&NLDRFunction::nldr_xsigmoid;
        break;
    default: ERROR("Unsupported transfer function");
    }
    pmode=mode;
}

/*! A collection of metric functions */
double NLDRMetricEuclid::pdist(const double* a, const double* b, unsigned long n) const
{
    double d=0.0;
    for (unsigned long i=0; i<n; ++i) d+=(b[i]-a[i])*(b[i]-a[i]);
    return sqrt(d);
}

void NLDRMetricPBC::pdiff(const double* a, const double* b, double* c, unsigned long n) const
{
#ifdef DEBUG
    if (n!=periods.size()) ERROR("Periodicity array has wrong dimensions\n");
#endif
    double dx;

    for (unsigned long i=0; i<n; ++i) 
    { 
        dx=(b[i]-a[i]); 
        dx/=periods[i]; 
        dx-=round(dx);  
        dx*=periods[i]; 
        c[i]=dx;
    }
}
double NLDRMetricPBC::pdist(const double* a, const double* b, unsigned long n) const
{
#ifdef DEBUG
    if (n!=periods.size()) ERROR("Periodicity array has wrong dimensions\n");
#endif
    double d=0.0, dx;

    for (unsigned long i=0; i<n; ++i) 
    { 
        dx=(b[i]-a[i]);
        dx/=periods[i];
        dx-=round(dx);
        dx*=periods[i];
        d+=dx*dx;
    }
    return sqrt(d);
}

void NLDRMetricSphere::pdiff(const double* a, const double* b, double* c, unsigned long n) const
{
#ifdef DEBUG
    if (n!=periods.size()) ERROR("Periodicity array has wrong dimensions\n");
#endif
    double dx;
    //!difference between two points in hyperspherical coordinates is not really well defined...
    
    for (unsigned long i=0; i<n-1; ++i) c[i]=(b[i]-a[i]); 
    
    dx=b[n-1]-a[n-1];  dx/=periods[n-1];  dx-=round(dx);  dx*=periods[n-1];  c[n-1]=dx;
}
double NLDRMetricSphere::pdist(const double* a, const double* b, unsigned long n) const
{
#ifdef DEBUG
    if (n!=periods.size()) ERROR("Periodicity array has wrong dimensions\n");
#endif
    double d=0.0;
    //!geodesic distance on the n-sphere
    
    double xs, ys, xi, yi, xy, twopi=toolbox::constant::pi*2.0;
    xy=0.0;
    xs=ys=1.0;
    for (unsigned long i=0; i<n; ++i) 
    {
        xi=xs*cos(a[i]*twopi/periods[i]);
        yi=ys*cos(b[i]*twopi/periods[i]);
        xs=xs*sin(a[i]*twopi/periods[i]);
        ys=ys*sin(b[i]*twopi/periods[i]);
        xy+=xi*yi;
    }
    xy+=xs*ys;
    if (xy>=1.0) d=0.0; else if (xy<=-1.0) d=constant::pi; else d=acos(xy);
    
    /*std::cerr<<"****************************************************\n";
    for (unsigned long i=0; i<n; ++i) std::cerr<<a[i]<<" "; std::cerr<<"\n";
    for (unsigned long i=0; i<n; ++i) std::cerr<<b[i]<<" "; std::cerr<<"\n";
    std::cerr<<"cos(d) "<<xy<<" >> "<<d<<"\n";
    */
    return d;
}

/*Builds neighbor list*/
void NLDRNeighborList::nlbuildup(const FMatrix<double>& points)
{
    unsigned long n=points.rows(), d=points.cols(), tn=0, pn;
    std::valarray<NLDRNeighbor> ni(n); 
    std::valarray<std::vector<NLDRNeighbor> > nn(n);
    
    npoint.resize(n+1); 
    
    //we rather compute distances twice but avoid storing a matrix which is n x n
    for (unsigned long i=0; i<n; i++)
    {
        //builds array with all the distances
        for (unsigned long j=0; j<n; j++) {
            ni[j].j=j; 
            ni[j].d=opts.ometric->pdist(&(const_cast<FMatrix<double>&>(points)(i,0)), &(const_cast<FMatrix<double>&>(points)(j,0)), d); }
        heapsort(ni); //sorts distances
        //picks maxneigh and/or within cutoff neighbors (skips self)
        unsigned long k=0; 
        while (
               ((opts.kw>0 && k<opts.kw) || ( 
               (opts.cutoff==0.0 || ni[k+1].d<opts.cutoff ) &&  
                 (opts.maxneigh==0 || (opts.kw==0 && k<opts.maxneigh)) ))
                && k<n) ++k;
        
        nn[i].resize(k); for (unsigned long j=0; j<k; ++j) nn[i][j]=ni[j+1];
        npoint[i]=k;  tn+=k; 
    }
    if (opts.kw>0)
    {
        std::cerr<<"Rebuilding list with weighted distances.\n";
        //variables for WLLE metric
        double wllec1, wllec2; double wllea, wlleb, wllel; std::valarray<double> wlletau(d);
        wllec2=sqrt(2.0)*exp(lgamma((d+1)*0.5)-lgamma(d*0.5));
        wllec1=wllec2/d; tn=0;
        std::cerr<<"WLLEC1,2  "<<wllec1<<","<<wllec2<<"\n";
        //builds deformed distance à la WLLE - which needs cycling again over the neighbors
        for (unsigned long i=0; i<n; i++)
        {
            //computes for the metric
            wlletau=0.0; wllel=0.0;
            for (unsigned long j=0; j<nn[i].size(); ++j)
            { 
                for (unsigned long k=0; k<d; ++k) wlletau[k]+=points(nn[i][j].j,k); 
                wllel+=nn[i][j].d;
            }
            wlletau*=1.0/nn[i].size(); for (unsigned long k=0; k<d; ++k) wlletau[k]-=points(i,k);
            wllel*=1.0/nn[i].size();
            
            wllea=wllel/wllec2; 
            wlleb=0.0; for (unsigned long k=0; k<d; ++k) wlleb+=wlletau[k]*wlletau[k]; wlleb=sqrt(wlleb);
            //std::cerr<<"L,  |avg| "<< wllel<<" "<<wlleb<<"\n";
            wlletau*=1.0/wlleb; wlleb*=1.0/wllec1;
            //std::cerr<<"WLLE a,b "<<wllea<<","<<wlleb<<"\n";
            if (wlleb>wllea) wlleb=wllea; //avoids getting negative distances just because the data are not CAM-distributed! this method just sucks.
            //builds array with all the distances
            for (unsigned long j=0; j<n; j++) {
                
                ni[j].j=j; 
                if (i==j) { ni[j]=0.0; continue; }
                ni[j].d=opts.ometric->pdist(&(const_cast<FMatrix<double>&>(points)(i,0)), &(const_cast<FMatrix<double>&>(points)(j,0)), d); 
                //std::cerr<<i<<","<<j<<" "<<ni[j].d<<" >> ";
                double wlletx=0.0; for (unsigned long k=0; k<d; ++k) wlletx+=(points(j,k)-points(i,k))*wlletau[k];
                ni[j].d=ni[j].d/(wllea+wlleb*wlletx/ni[j].d);
                //std::cerr<<ni[j].d<<"\n";
            }

            heapsort(ni); //sorts distances
            //picks maxneigh and/or within cutoff neighbors (skips self)
            unsigned long k=0; 
            while ( (opts.cutoff==0.0 || ni[k+1].d<opts.cutoff ) &&  
                     (opts.maxneigh==0 || k<opts.maxneigh ) && k<n) ++k;
    
            nn[i].resize(k); for (unsigned long j=0; j<k; ++j) nn[i][j]=ni[j+1];
            npoint[i]=k;  tn+=k; 
        }
    }
        
    //implements greedy/liberal symmetrization of the neighbor list
    if (opts.greediness!=NLDRAsym)
    for (unsigned long i=0; i<n; ++i)
    {
        for (unsigned long j=0; j<npoint[i]; ++j)
        {
            bool ijsym=false; int ijj=nn[i][j].j;
            for (unsigned long k=0; k<npoint[ijj]; ++k) if (nn[ijj][k].j==i) { ijsym=true; break; } 
            if (!ijsym) {
                switch(opts.greediness) {
                    case NLDRGreedy:
                        nn[ijj].push_back(NLDRNeighbor(i,nn[i][j].d)); 
                        ++npoint[ijj]; ++tn;
                        break;
                    case NLDRLiberal:
                        nn[i].erase(nn[i].begin()+j);
                        --npoint[i]; --j; --tn;
                        break;
                    default:
                        ERROR("Unsupported NL symmetrization scheme");
                }
            }
        }
    }
    
    //collapses the (expensive) storage as vectors to a "compressed" storage
    nlist.resize(tn); unsigned long k=0; pn=tn=0;
    for (unsigned long i=0; i<n; ++i)
    {
        pn=tn; tn+=npoint[i]; std::valarray<NLDRNeighbor> nni(nn[i].size()); 
        for (unsigned long j=0; j<npoint[i]; ++j) nni[j]=nn[i][j]; if (npoint[i]>1) heapsort(nni);
        for (unsigned long j=0; j<npoint[i]; ++j) nlist[k++]=nni[j];
        npoint[i]=pn;
    }
    npoint[n]=tn;
}

NLDRNeighbor& NLDRNeighborList::rneigh(unsigned long i, unsigned long j)
{
#ifdef DEBUG
    if (i>npoint.size()-1 || j>(npoint[i+1]-npoint[i])) ERROR("Trying to access non-existent neighbor.");
#endif 
    return nlist[npoint[i]+j];
}

unsigned long NLDRNeighborList::index(unsigned long i, unsigned long j) const 
{ return const_cast<NLDRNeighborList*>(this)->rneigh(i,j).j; }
double NLDRNeighborList::dist(unsigned long i, unsigned long j) const
{ return const_cast<NLDRNeighborList*>(this)->rneigh(i,j).d; }
NLDRNeighbor NLDRNeighborList::neigh(unsigned long i, unsigned long j) const
{ return const_cast<NLDRNeighborList*>(this)->rneigh(i,j); }
void NLDRNeighborList::RemoveLonesome(std::valarray<unsigned long>& llone)
{
    unsigned long nlone=0;
    for (unsigned long i=0; i<npoint.size()-1; ++i) if ((npoint[i+1]-npoint[i])<=0) ++nlone;
    llone.resize(nlone); if (nlone==0) return;
    std::valarray<unsigned long> newnp(npoint.size()-nlone);
    nlone=0;
    for (unsigned long i=0; i<npoint.size()-1; ++i) 
    {
        newnp[i-nlone]=npoint[i];
        if ((npoint[i+1]-npoint[i])<=0) { llone[nlone]=i; ++nlone; }
    }
    
    //then updates indices
    for (unsigned long i=0; i<nlist.size(); ++i) 
    { 
        unsigned long dni=0;
        for (unsigned long j=0; j<nlone; ++j) if (nlist[i].j>llone[j]) ++dni;  else break;
        nlist[i].j-=dni;
    }
    
    newnp[newnp.size()-1]=npoint[npoint.size()-1];
    npoint.resize(newnp.size()); npoint=newnp;
}

void NLDRProjection::calc_PM()
{
    PM.resize(n); HV.resize(n); LV.resize(n);
    std::cerr<<"Building neighbor list\n";
    if (ftainted) neigh.Build(P, opts.nopts);
    
    std::valarray<double> v1(D), v2(D), v3(D);
    std::valarray<double> w1(d), w2(d), w3(d);
    
    std::cerr<<d<<":"<<D<<" Computing pseudoinverses\n";
    for (unsigned long i=0; i<n; ++i)
    {
        //for (unsigned long j=0; j<neigh.nneigh(i); ++j) std::cout<<neigh.index(i,j)<<" "; std::cout<<"\n";
        std::cerr<<i<<" nneigh: "<<neigh.nneigh(i)<<"\n";
        HV[i].resize(neigh.nneigh(i),D);
        v1=P.row(i);
        for (unsigned long j=0; j<neigh.nneigh(i); ++j)
        { 
            v2=P.row(neigh.index(i,j)); 
            opts.nopts.ometric->diff(v2,v1,v3);
            HV[i].row(j)=v3; 
        }
        
        LV[i].resize(neigh.nneigh(i),d);
        w1=p.row(i);
        
        for (unsigned long j=0; j<neigh.nneigh(i); ++j)
        { 
            w2=p.row(neigh.index(i,j)); 
            w3=w2; w3-=w1; //MAKE IT MORE GENERAL BY USING AN ARBITRARY METRIC 
            LV[i].row(j)=w3; 
        }
        PseudoInverse(HV[i],PM[i]);
    }
    
    ftainted=false;
}

void NLDRProjection::set_vars(const std::valarray<double>& rv)
{
    vx.resize(rv.size()); vx=rv; 
    vv=0.0; vg.resize(rv.size()); vg=0.0;
    std::valarray<double> v1(d), V1(D); double ld, lfd, ldfd;
    
   // vdX=0.0;  //dChi/dXj
   // vdXx=0.0;  //d2Chi/dXj dxk
    double diffdist, distdX1, tw=0.0;
    for (unsigned long i=0; i<n; ++i)
    {
        //!TODO TEST must check and implement curse-of-dimensionality correction
        v1=p.row(i); v1-=vx; ld=0.0; for (unsigned long j=0; j<d; ++j) ld+=v1[j]*v1[j]; ld=sqrt(ld);
        
        if (ld<=0.0) continue;
        opts.tfunL.fdf(ld,lfd,ldfd);
        
        diffdist=(vfd[i]-lfd); 
        vv+=diffdist*diffdist*w[i]; 
        
        /*
        distdX1=2.0*vf1d[i]/nd[i].d*ldfd/ld;
        for (unsigned long j=0; j<D; ++j) for (unsigned long k=0; k<d; ++k) 
                vdXx(j,k)+=vxx(i,j)*v1[k]*distdX1;
        */
        v1*=2.0*diffdist*ldfd/ld*w[i];
        vg+=v1;
        tw+=w[i];
        /*
        distdX1=2.0*diffdist*vf1d[i]/nd[i].d;
        V1=vxx.row(i); V1*=distdX1;
        vdX+=V1;
        */
    }
    vv*=1.0/tw; vg*=1.0/tw; //vdX*=1.0/n;
}
            
double NLDRProjection::project(const std::valarray<double>& np, std::valarray<double>& hp,  std::valarray<double>& lp, double &md)
{
    //Finds closer-by point and optimized cost function starting from there
    std::valarray<double> v1(D), v2(D); 
    nd.resize(n); 
    vxx.resize(n,D); 
    NLDRNeighbor mind; mind.d=-1.0;
    //calc all distances
    vfd.resize(n); vf1d.resize(n);
    std::cerr<<" n. points" <<n<<"\n";
    for (unsigned long i=0; i<n; ++i)
    {
        v1=P.row(i);
        opts.nopts.ometric->diff(np,v1,v2);
        vxx.row(i)=v2;
        nd[i].j=i; 
        //! WRONG DISTANCE (in spherical)
        //nd[i].d=0.0; for (unsigned long k=0; k<D; ++k) nd[i].d+=v2[k]*v2[k]; nd[i].d=sqrt(nd[i].d);
        nd[i].d=opts.nopts.ometric->dist(np,v1);
        opts.tfunH.fdf(nd[i].d,vfd[i],vf1d[i]);
        if (nd[i].d<mind.d || mind.d<0.0)  mind=nd[i];
    }
    vdX.resize(D); vdXx.resize(D,d);
    //IMDS PROJ with global minimum search
    //creates grid of points and gradients
    unsigned long ngrid=opts.grid1; double wgrid=opts.gwidth;
    if (d!=2) ERROR ("Integral projector is implemented only for 2D");
    std::valarray<double> gx(ngrid), gy(ngrid);
    FMatrix<double> gridU(ngrid,ngrid); std::valarray<double> x(d);
    FTensor<double, 3> gridDU(fixarray<unsigned long, 3>(ngrid,ngrid,d));

    for (unsigned long i=0; i<ngrid; i++) { gx[i]=gy[i]=(-1.0+i*2.0/(ngrid-1)); }
    gx*=wgrid; gy*=wgrid; 
    double minf, reff=-1.;
    for (unsigned long i=0; i<ngrid; i++)
        for (unsigned long j=0; j<ngrid; j++)
    { 
        x[0]=gx[i]; x[1]=gy[j];  
        set_vars(x); get_value(gridU(i,j)); get_gradient(x);
        if (reff<0 || reff>gridU(i,j)) minf=reff=gridU(i,j);
        gridDU(i,j,0)=x[0];  gridDU(i,j,1)=x[1];
        //gridU(i,j)=(x[0]/wgrid)*(x[0]/wgrid)+(x[1]/wgrid)*(x[1]/wgrid);
        //gridDU(i,j,0)=2*(x[0]/wgrid)/wgrid; gridDU(i,j,1)=2*(x[1]/wgrid)/wgrid; 
        //gridU(i,j)=sin((constant::pi*x[0]/wgrid))*cos((constant::pi*x[1]/wgrid));
        //gridDU(i,j,0)=constant::pi*cos((constant::pi*x[0]/wgrid))*cos((constant::pi*x[1]/wgrid))/wgrid; gridDU(i,j,1)=-constant::pi*sin((constant::pi*x[0]/wgrid))*sin((constant::pi*x[1]/wgrid))/wgrid; 
        
    }

    InterpolateBicubic bicub; 
    std::cerr<<"create bicubic interpolator\n";
    bicub.set_table(gx,gy,gridU,gridDU); 
    
    //std::cerr<<"printing interpolated function\n";
    //std::ofstream intfile(interp_out.c_str());
    
    ngrid=opts.grid2; gx.resize(ngrid); gy.resize(ngrid); 
    for (unsigned long i=0; i<ngrid; i++) { gx[i]=gy[i]=(-1.0+i*2.0/(ngrid-1)); }
    gx*=wgrid; gy*=wgrid; 
    
    double tr, tx, ty, tt, kt=opts.gtemp, echi; tr=tx=ty=tt=0.0;
    fixarray<double,2> df;
    fixarray<double,2 > mx;

    for (unsigned long i=0; i<ngrid; i++)
    {
        for (unsigned long j=0; j<ngrid; j++)
        {
            double f; fixarray<double,2> x;
            x[0]=gx[i]; x[1]=gy[j];
            bicub.get_ydy(x,f,df);
            //intfile<<gx[i]<<" "<<gy[j]<<" "<<f<<"\n";
            echi=exp(-(f-reff)/kt);
            if (f<=minf) { mx=x; minf=f; }
            tt+=echi; tx+=x[0]*echi; ty+=x[1]*echi; tr+=sqrt(x[0]*x[0]+x[1]*x[1])*echi;
        }
        //intfile<<std::endl;
    }
    //intfile.close();
    
    tr*=1.0/tt; tx*=1.0/tt; ty*=1.0/tt; 
    
    std::ofstream pointfile("minimum");
    pointfile<<tx<<" "<<ty<<"\n";
    pointfile<<tr*cos(atan2(ty,tx))<<" "<<tr*sin(atan2(ty,tx))<<"\n";
    pointfile.close();
    
    
    hp.resize(np.size()); hp=np;
    
    lp.resize(d);
    if (opts.cgsteps>0)
    {
        ConjGradOpts cgop;
        cgop.maxiter=opts.cgsteps; 
        cgop.linesearch.maxiter=4; cgop.linesearch.lstol=1e-9; 
        std::valarray<double> pcoords(2), rpos; pcoords[0]=mx[0]; pcoords[1]=mx[1];  double rvalue;
        set_vars(pcoords);
        std::cerr<<"starting minimization at "<<pcoords<<" ("<<minf<<")\n";
        min_conjgrad(*this,pcoords,rpos, rvalue, cgop );
        get_vars(lp);
    }
    else
    {
        lp[0]=tr*cos(atan2(ty,tx)); lp[1]=tr*sin(atan2(ty,tx));
    }
    double f; fixarray<double,2>  lx; 
    lx[0]=lp[0]; lx[1]=lp[1];
    bicub.get_ydy(lx,f,df);
    md=mind.d;
    return f;
    //exit(1);
    //IMDS PROJECTION : MULTIPLE MINIMA <-> jumpy trajectories
    /*
    static std::valarray<double> xx(0); 
    if (xx.size()==0) { xx.resize(d); xx=0.0; }
    
    set_vars(xx);
    std::valarray<double> rpos; double rvalue;
    SteepestOpts steepop; steepop.maxiter=10;
    steepop.linesearch.maxiter=5; steepop.linesearch.lstol=1e-10; 
    min_conjgrad(*this,xx, rpos, rvalue, steepop );
    
    hp.resize(D); lp.resize(d);
    hp=np; get_vars(lp); xx=lp;
    return mind.d;
    */
    /*
    
    GEOMETRIC PROJECTION : UNSTABLE AND NONDISCRIMINATING
    if (ftainted) calc_PM();
    
    std::valarray<double> v1(D), v2(D); 
    std::valarray<NLDRNeighbor> nd(n); double mind=-1.0;
    //calc all distances
    for (unsigned long i=0; i<n; ++i)
    {
        v1=P.row(i);
        nd[i].j=i; nd[i].d=opts.nopts.ometric->dist(np,v1);
        if (nd[i].d<mind || mind<0.0)  mind=nd[i].d;
    }
    
    std::valarray<double> b, w1(d), pi(d), Pi(D); double tweight=0.0, weight;
    hp.resize(D); lp.resize(d); hp=0.0; lp=0.0;
    for (unsigned long i=0; i<n; ++i)
    {
        if (nd[i].d>mind*opts.acutoff) continue;
        
        v1=P.row(i); 
        opts.nopts.ometric->diff(np,v1,v2);
        b.resize(neigh.nneigh(i)); b=0.0; 
        for (unsigned long j=0; j<neigh.nneigh(i); ++j) 
            for (unsigned long k=0; k<D; ++k) b[j]+=PM[i](k,j)*v2[k];
        
        pi=0.0; Pi=0.0;
        for (unsigned long j=0; j<neigh.nneigh(i); ++j) 
        { v1=HV[i].row(j); v1*=b[j]; Pi+=v1; } 
        
        for (unsigned long j=0; j<neigh.nneigh(i); ++j) 
        { w1=LV[i].row(j); w1*=b[j]; pi+=w1; }
        Pi+=P.row(i); pi+=p.row(i);

        if (nd[i].d==0.0) { hp=Pi; lp=pi; tweight=1.0; break; }
        weight=pow(1.0/nd[i].d,opts.dpower)*taper(nd[i].d,mind,mind*opts.acutoff);
        Pi*=weight; pi*=weight; tweight+=weight;
        hp+=Pi; lp+=pi;
        
    }
    
    hp*=1.0/tweight; lp*=1.0/tweight;
    return mind;
    */
}

double isqrt(double x) { return 1.0/sqrt(x); }
void NLDRLLE(FMatrix<double>& points, NLDRProjection& proj, const NLDRLLEOptions& opts, NLDRLLEReport& report)
{
    unsigned long dts=(opts.dimts==0?opts.lowdim:opts.dimts);
    proj.P=points; proj.D=points.cols(); proj.d=opts.lowdim; proj.n=points.rows();
    proj.p.resize(proj.n,proj.d); 
    
    std::cerr<<"Building neighbor list\n";
    NLDRNeighborList nlist(proj.P,opts.nlopts);
    if (opts.rmlonesome) 
    {
        std::valarray<unsigned long> ilone;
        nlist.RemoveLonesome(ilone);
        //removing points with no connections
        unsigned long nlone=ilone.size();
            
        if (nlone>0) 
        {
            std::cerr<<"Removing isolated points\n";
            proj.P.resize(proj.n-nlone,proj.D);
            unsigned long k=0, j=0;
            for (unsigned long i=0; i<proj.n; ++i) 
            {
                if (j<nlone && ilone[j]==i) { ++j; continue; }
                for (unsigned long h=0; h<proj.D; ++h) proj.P(k,h)=points(i,h);
                ++k;
            }
            proj.n-=nlone;
            proj.p.resize(proj.n,proj.d);
        }
    }
    
    
    std::cerr<<"Building weights\n";
    std::valarray<std::valarray<double> > weights(proj.n);
    std::valarray<std::valarray<unsigned long> > indices(proj.n);
    
    if (opts.verbose) { report.hd_errors.resize(proj.n); report.hd_errors=0.0;  report.ld_errors.resize(proj.n); report.ld_errors=0.0; }
    report.ld_error=report.hd_error=0.0;
    
    for (unsigned long i=0; i<proj.n; ++i)  {  unsigned long m=nlist.nneigh(i); if (m==0) ERROR("Point "<<i<<" is isolated!\n"); weights[i].resize(m); indices[i].resize(m); }
    
    CrsMatrix<double> W, WT; double dp=dts*(dts+1)/2; std::valarray<double> x1(proj.D), x2(proj.D), dx12(proj.D);
    if (opts.mode==HLLE) W.resize(proj.n*dp,proj.n);
    for (unsigned long i=0; i<proj.n; ++i)
    {
        unsigned long m=nlist.nneigh(i);

        FMatrix<double> C(m,m), C1, G(proj.D,m), GT, H, HT;
        x1.resize(proj.D); x2.resize(proj.D); 
        x1=proj.P.row(i);
        //builds matrix with neighbor vector shifts around the central point
        for (unsigned long j=0; j<m; ++j) 
        {
            x2=proj.P.row(nlist.index(i,j));
            opts.nlopts.ometric->diff(x2,x1,dx12);

            G.col(j)=dx12;
        } 
        transpose(G,GT);

        std::cerr<<"POINT "<<i<<": ";
        for (unsigned long j=0; j<m; ++j)  std::cerr<<nlist.index(i,j)<<"("<<nlist.dist(i,j)<<") ";
        std::cerr<<"\n";
            
        if(opts.mode==LLE)
        {
             mult(GT,G,C);
        }
        else if (opts.mode==LLTE)
        {
            //First, finds d principal components
            mult(G,GT,C);
            FMatrix<double> P, PT; std::valarray<double> p;
            EigenSolverSym(C,P,p,LAEMIndex,proj.D-dts,proj.D-1);

            transpose(P,PT); mult(PT,G,H);  transpose(H, HT);
            mult(HT,H,C);  // now C is projected on the most relevant eigenvalues
        }
        else if (opts.mode==HLLE)
        {
            mult(GT,G,C);
            FMatrix<double> U; std::valarray<double> u;
            EigenSolverSym(C,U,u,LAEMIndex,m-dts,m-1);
            
            //builds Hessian estimator
            x1.resize(m); x2.resize(m);
            FMatrix<double> Xi(m,1+dts+dp);
            x1=1.; Xi.col(0)=x1;
            for (unsigned long h=0; h<dts; ++h) Xi.col(1+h)=U.col(h);
            unsigned long k=0;
            
            for (unsigned long h=0; h<dts; ++h) { x1=U.col(h); for (unsigned long h2=h; h2<dts; ++h2) 
            {
                x2=U.col(h2); x2*=x1;
                Xi.col(1+dts+k)=x2; k++;
            } }

            //GS-orthogonalize columns
            FMatrix<double> R(Xi.cols(),Xi.cols(),0.0);
            for (unsigned long h1=0; h1<R.rows(); h1++)
            {
                R(h1,h1)=0.0; for (unsigned long h2=0; h2<Xi.rows(); h2++) R(h1,h1)+=Xi(h2,h1)*Xi(h2,h1); R(h1,h1)=sqrt(R(h1,h1));
                for (unsigned long h2=0; h2<Xi.rows(); h2++) Xi(h2,h1)*=1./R(h1,h1);
                if (h1==R.rows()-1) break;
                for (unsigned long h2=h1+1; h2<R.rows(); h2++) 
                {
                    R(h1,h2)=0.0; for (unsigned long h3=0; h3<Xi.rows(); h3++) R(h1,h2)+=Xi(h3,h1)*Xi(h3,h2);
                    for (unsigned long h3=0; h3<Xi.rows(); h3++) Xi(h3,h2)-=R(h1,h2)*Xi(h3,h1);
                }
            }

            H.resize(dp,m);
            for (unsigned long h=0; h<dp; ++h) H.row(h)=Xi.col(h+1+dts);
        }
        else ERROR("Unsupported LLE mode");

        if (opts.mode==LLE || opts.mode==LLTE) 
        {
            //smoothens the matrix
            double tr=0.0; 
            if (opts.smooth<0) { tr=normfrob(C)*(-opts.smooth); }
            else tr=opts.smooth;
            for (unsigned long j=0; j<m; ++j) C(j,j)+=tr;
            
            MatrixInverse(C,C1);
            
            double beta, lambda;
            beta=0.0;
            for (unsigned long j=0; j<m; ++j) for (unsigned long k=0; k<m; ++k) beta+=C1(j,k);
            lambda=1.0/beta;
            for (unsigned long j=0; j<m; ++j)
            {
                indices[i][j]=nlist.index(i,j);  weights[i][j]=0.0; 
                for (unsigned long k=0; k<m; ++k) weights[i][j]+=C1(j,k);
            }
            weights[i]*=lambda;
            
            //with LLTE we are computing the error in the reconstruction of PROJECTED high-D points
            for (unsigned long h=0; h<(opts.mode==LLE?proj.D:proj.d); ++h) 
            { 
                double ed=0.0;
                if (opts.mode==LLE)
                {
                    for (unsigned long j=0; j<m; ++j) ed+=weights[i][j]*proj.P(nlist.index(i,j),h);
                    ed-=proj.P(i,h); ed*=ed;
                }
                else if (opts.mode==LLTE)
                {
                    for (unsigned long j=0; j<m; ++j) ed+=weights[i][j]*H(h,j);
                    ed*=ed;  
                }
                report.hd_error+=ed; if (opts.verbose) report.hd_errors[i]+=ed;
            }
            if (opts.verbose) report.hd_errors[i]=sqrt(report.hd_errors[i]);
        }
        else if (opts.mode==HLLE)
        {
            if (opts.verbose)
            {
                FMatrix<double> HT, HHT;
                transpose(H,HT); mult(HT,H,HHT);
                report.hd_errors[i]=trace(HHT)/HHT.rows();
                FMatrix<double> U; std::valarray<double> u;
                EigenSolverSym(HHT,U,u);
                std::cerr<<" singular values of HHT for point "<<i<<"\n"<<u;
                std::cerr<<report.hd_errors[i]<<" "<<normfrob(H)<<"\n";
            }
/*            FMatrix<double> HT, HHT, U; std::valarray<double> u;
            transpose(H,HT); mult(HT,H,HHT);
            EigenSolverSym(C,U,u);
            report.hd_errors[i]=u[u.size()-1];
            std::cerr<<" singular values of HHT for point "<<i<<"\n";
            //if (i==1130 || i==833) H*=0.0; 
            std::cerr<<u;
*/
            //builds W straight away
            for (unsigned long h=0; h<dp; ++h)
            {
                for (unsigned long j=0; j<m; ++j)
                    W(i*dp+h,nlist.index(i,j))=H(h,j);
            }
        }
    }
    report.hd_error=sqrt(report.hd_error/proj.n);
    
    
    if (opts.mode==LLE||opts.mode==LLTE)
    {
        //we have the weight matrix as a sparse matrix
        std::cerr<<"Building sparse W matrix\n";
        W=CrsMatrix<double>(proj.n,proj.n,indices,weights);
        W*=-1.0;  for (unsigned long i=0; i<proj.n; ++i) W(i,i)+=1.0;
    }
    //finds eigenvalues by calling ARPACK routines, else uses full matrix algebra...
#ifdef ARPACK
    int IDO, N, NEV, NCV, LDV, IPARAM[11], IPNTR[11], LWORKL, INFO;
    char BMAT, WHICH[2]; 
    BMAT='I'; WHICH[0]='S'; WHICH[1]='A';
    NEV=proj.d+2; N=proj.n;
    NCV=50*(NEV+1); if(NCV>N) NCV=N; 
    LDV=N;
    for (int i=0; i<11; i++) IPARAM[i]=IPNTR[i]=0;
    IPARAM[0]=1; IPARAM[2]=3*N; IPARAM[3]=1; 
    IPARAM[4]=NEV; IPARAM[6]=1; 
    double TOL=1e-5; 
    IDO=0; LWORKL=NCV*(NCV+8); INFO=0; 
    std::valarray<double> V(N*NCV), RESID(N), WORKD(3*N), WORKL(LWORKL); 
    V=0.0; RESID=0.0; WORKD=0.0; WORKL=0.0;
    std::valarray<double> y(N);
    IDO=0;
    
    do {
/*        std::cerr<<"Calling DSAUPD\n";
        std::cerr<<"IDO "<<IDO<<"\n";
        std::cerr<<"N "<<N<<"\n";
        std::cerr<<"NEV "<<NEV<<"\n";
        std::cerr<<"NCV "<<NCV<<"\n";
        std::cerr<<"LDV "<<LDV<<"\n";
        std::cerr<<"LWORKL "<<LWORKL<<"\n";
        std::cerr<<"INFO "<<INFO<<"\n";
        std::cerr<<"IPARAM "; for (int i=0; i<11; i++)std::cerr<<IPARAM[i]<<" "; std::cerr<<"\n";
        std::cerr<<"IPNTR "; for (int i=0; i<11; i++)std::cerr<<IPNTR[i]<<" "; std::cerr<<"\n";
        */
        tblapack::dsaupd(&IDO, &BMAT, &N, &WHICH[0], &NEV, 
               &TOL, &RESID[0], &NCV, &V[0], &LDV, 
               &IPARAM[0], &IPNTR[0], &WORKD[0], 
               &WORKL[0], &LWORKL, &INFO);
        
        std::valarray<double> x(WORKD[std::slice(IPNTR[0]-1,N,1)]); //takes slice
        mult(W,x,y); x=y; Tmult(W,x,y);
        std::cerr<<"Returned IDO: "<<IDO<<" INFO: "<<INFO<<"\n";
        if (INFO!=0) std::cerr<<" IPARAM[4]: "<<IPARAM[4]<<"\n";
        WORKD[std::slice(IPNTR[1]-1,N,1)]=y;
    } while (IDO==1 || IDO==-1);
    
    int RVEC; std::valarray<int> SELECT(NCV); std::valarray<double> D(N), Z(N*NEV);
    int LDZ; char HOWMNY; double SIGMA;
    RVEC=1; HOWMNY='A'; LDZ=N; SIGMA=0.0; SELECT=0; D=0.0; Z=0.0; 
    
    tblapack::dseupd(&RVEC, &HOWMNY, &SELECT[0], &D[0], &Z[0], &LDZ, &SIGMA,
            &BMAT, &N, WHICH, &NEV, 
            &TOL, &RESID[0], &NCV, &V[0], &LDV, 
            IPARAM, IPNTR, &WORKD[0], 
            &WORKL[0], &LWORKL, &INFO);
    for (int i=0; i<NEV; i++) std::cerr<<"EIGV: "<<i<<" is "<<D[i]<<"\n";
    
    for (unsigned long i=0; i<proj.n; ++i) for (unsigned long h=0; h<proj.d; ++h) 
        proj.p(i,h)=Z[i+(h+1)*proj.n];
#else
    std::cerr<<"Finding low-dim points (matrix-matrix mult)\n";
    transpose(W,WT); CrsMatrix<double> MM; mult(WT,W,MM); 
    FMatrix<double> FM(MM), Q; std::valarray<double> q;
    
    report.ld_errors=FM.diag();
    
    std::valarray<double> c(proj.n);
    for (unsigned long i=0; i<proj.n; ++i) if (FM(i,i)<0.5)
    {
        c=FM.col(i); c*=c; std::cerr<<"element "<<i<<" diag "<<FM(i,i)<<" row |sum| "<<sqrt(c.sum())<<"\n";
    }
        
    /*
    double dw=2;
    for (unsigned long i=0; i<proj.n; ++i) if (FM(i,i)<2)
    for (unsigned long j=0; j<proj.n; ++j) 
    { FM(i,j)-=dw/proj.n; if(i==j) FM(i,i)+=dw; else FM(j,i)-=dw/proj.n; }
    */
    
    /*FM.diag()+=wbad;
    wbad*=(1.0/proj.n);
    for (unsigned long i=0; i<proj.n; ++i) FM.row(i)-=wbad;*/
    //std::cerr<<"DIAGONAL: "<<wbad<<"\n";
    //FM(833,833)+=1e-3; FM(1130,1130)+=1e-3;
    /*FM.row(833)*=wbad; FM.col(833)*=wbad; 
    FM.row(1130)*=wbad; FM.col(1130)*=wbad; 
    std::cerr<<"Scaled elements of FM "<<FM(833,0)<<"; "<<FM(1130,1130)<<"\n";
    FM.row(100)*=wbad; FM.col(100)*=wbad; 
    FM.row(101)*=wbad; FM.col(101)*=wbad; 
    FM.row(102)*=wbad; FM.col(102)*=wbad; */
    
    //std::cerr<<FM<<"\n";
    std::cerr<<"Finding low-dim points (diagonalization "<<FM.rows()<<"x"<<FM.cols()<<")\n";
    if (opts.verbose) {
        EigenSolverSym(FM,Q,q,LAEMIndex,proj.d+1,proj.d+1);
        report.dp1eval=q[0];
    }
    //std::cerr<<"Finding zero eigenvec\n";
    //EigenSolverSym(FM,Q,q,LAEMIndex,0,0);
    //std::cerr<<"zero-eigenval: "<<q<<"eigenvec"<<Q<<"\n";
    EigenSolverSym(FM,Q,q,LAEMIndex,1,proj.d);
    //std::cerr<<"Finding zero eigenvec\n";
    report.deval.resize(proj.d); report.deval=q;
#endif
    
    if (opts.mode==LLE||opts.mode==LLTE)
    {
        for (unsigned long i=0; i<proj.n; ++i) for (unsigned long h=0; h<proj.d; ++h) proj.p(i,h)=Q(i,h);

        for (unsigned long i=0; i<proj.n; ++i)
        {
            unsigned long m=nlist.nneigh(i);
            for (unsigned long h=0; h<proj.d; ++h) 
            { 
                double ed=0.0;
                for (unsigned long j=0; j<m; ++j) ed+=weights[i][j]*proj.p(nlist.index(i,j),h);
                ed-=proj.p(i,h); ed*=ed;
                report.ld_error+=ed; if (opts.verbose) report.ld_errors[i]+=ed*ed;
            }
            if (opts.verbose) report.ld_errors[i]=sqrt(report.ld_errors[i]);
        }
    }
    else if (opts.mode==HLLE) 
    {
        std::cerr<<"HLLE Orienting coordinates\n";
        Q*=1.0/sqrt(proj.n);
        FMatrix<double> QT, R, R2; transpose(Q,QT); mult(QT,Q,R); MatrixFunctionSym(R,&isqrt,R2);
        mult(Q,R2,R);
        //R=Q; // no scaling
        for (unsigned long i=0; i<proj.n; ++i) for (unsigned long h=0; h<proj.d; ++h) proj.p(i,h)=R(i,h);
    }
}

void NLDRMDS(FMatrix<double>& points, NLDRProjection& proj, const NLDRMDSOptions& opts, NLDRMDSReport& report, const FMatrix<double>& outd)
{
    
    if (opts.metric==NULL) ERROR("Uninitialized metric pointer\n");
    proj.P=points; proj.D=points.cols(); proj.d=opts.lowdim; proj.n=points.rows();
    proj.p.resize(proj.n,proj.d);
    
    unsigned long n=points.rows(), D=points.cols();
    
    FMatrix<double> dist(n,n);
    FMatrix<double> M(dist);
    FMatrix<double> Q; std::valarray<double> q,sq;
    double sr;
    if (outd.rows()==n && outd.cols()==n) dist=outd;  
    else
    {
        std::cerr<<"Building distance matrix\n";
        for (unsigned long i=0; i<n; i++)
        {
            dist(i,i)=0;
            for (unsigned long j=0; j<i; j++) {
                dist(i,j)=dist(j,i)=opts.metric->dist(&(const_cast<FMatrix<double>&>(points)(i,0)), &(const_cast<FMatrix<double>&>(points)(j,0)), D); 
            }
        }
    }
    
    if (opts.mode==TMDS)
    {
        //prepares a fake call to NLDRMDS
        NLDRProjection tproj; NLDRMDSOptions topts(opts); NLDRMDSReport treport;
        treport=report;
        unsigned long th=0.0;
        topts.mode=SMDS; topts.lowdim=1; topts.metric=new NLDRMetricEuclid();
        report.deval.resize(proj.d);
        
        for (th=0; th<proj.d; th++)
        {
            //computes the great circle size
            sr=0.0;
            for (unsigned long i=0; i<n; i++)
                for (unsigned long j=0; j<i; j++) 
                    if(dist(i,j)>sr) sr=dist(i,j);
            sr/=constant::pi;
            
            //projects onto a circle
            std::cerr<<"projecting onto a circle, radius: "<<sr<<"\n";
            NLDRMDS(points,tproj,topts,treport,dist);
            
            report.deval[th]=treport.deval[0];
            
            for (unsigned long i=0; i<n; i++) proj.p(i,th)=tproj.p(i,0);
            
            //updates distance matrix
            double tdij;
            for (unsigned long i=0; i<n; i++) for (unsigned long j=0; j<i; j++)
            {
                tdij=fabs(tproj.p(i,0)-tproj.p(j,0));
                while (tdij>1) tdij-=2;
                
                tdij=fabs(tdij)*constant::pi*sr;
                tdij=dist(i,j)*dist(i,j)-tdij*tdij;
                if (tdij<0.0) std::cerr<<"negative corrected distance "<<tdij<<" was "<<dist(i,j)*dist(i,j)<<"\n";
                dist(i,j)=dist(j,i)=sqrt(fabs(tdij));
            }
        }
        std::cerr<<"Iterative toroidal projection complete\n";
    }
    else
    {
     
    if (opts.mode==MDS)
    {
        std::cerr<<"Multiplying to find full metric matrix\n";
        for (unsigned long i=0; i<n; i++)
            for (unsigned long j=0; j<n; j++) 
                M(i,j)=dist(i,j)*dist(i,j)*(-0.5);
        
        //we don't really want to do matrix-matrix multiplies, since H is a low-rank matrix (better, it can be decomposed as such)
        //we need to compute HMH with H being 1-n^-1 1 1^T
        FMatrix<double> T1(n,n); std::valarray<double> tv(n); double td;
        
        //first, get HM
        for (unsigned long i=0; i<n; i++) { tv=M.col(i); td=tv.sum(); tv=(-td/n); T1.col(i)=tv; }
        M+=T1;
        //and finish up with (HM)H
        for (unsigned long i=0; i<n; i++) { tv=M.row(i); td=tv.sum(); tv=(-td/n); T1.row(i)=tv; }
        M+=T1;
    }
    else if (opts.mode==SMDS)
    {
        std::cerr<<"Finding spherical metric matrix\n";
        sr=0.0;
        for (unsigned long i=0; i<n; i++)
            for (unsigned long j=0; j<i; j++) 
                if(dist(i,j)>sr) sr=dist(i,j);
        sr/=constant::pi;
        M*=0.0; M+=1.0;
        for (unsigned long i=0; i<n; i++)
            for (unsigned long j=0; j<i; j++) 
                M(i,j)=M(j,i)=cos(dist(i,j)/sr);
        M*=(sr*sr);
    }
    
    std::cerr<<"Finding low-dim points (diagonalization "<<M.rows()<<"x"<<M.cols()<<")\n";
    unsigned long neva=proj.d; if (opts.mode==SMDS) neva++;
    EigenSolverSym(M,Q,q,LAEMIndex,proj.n-1-neva+1,proj.n-1);
    sq.resize(q.size()); sq=sqrt(q);
    for (unsigned long i=0; i<proj.n; ++i) Q.row(i)*=sq;
    
    report.deval.resize(proj.d); 
    for (unsigned long h=0; h<proj.d; ++h) { report.deval[h]=q[neva-1-h]; }
    report.ld_error=report.deval.sum()/trace(M);

    if (opts.mode==MDS)
    {
        for (unsigned long i=0; i<proj.n; ++i) for (unsigned long h=0; h<proj.d; ++h) proj.p(i,h)=Q(i,proj.d-1-h);
    }
    else if (opts.mode==SMDS)
    {
        double tx=0.0;
        for (unsigned long i=0; i<proj.n; ++i) 
        {
            tx=0.0;
            proj.p(i,0)=atan2(Q(i,proj.d),Q(i,proj.d-1))/constant::pi;
            tx+=Q(i,proj.d)*Q(i,proj.d);
            for (unsigned long h=1; h<proj.d; ++h) 
            {
                tx+=Q(i,proj.d-h)*Q(i,proj.d-h);
                proj.p(i,h)=atan2(sqrt(tx),Q(i,proj.d-h-1))/constant::pi;
            }
            sq=1./sqrt(tx);
            Q.row(i)*=sq; //while we are at it, we normalize the vectors so that they really ARE on top of a sphere
        }
    }
    
    if (opts.verbose) {
        std::cerr<< " Computing detailed error report\n";
        EigenSolverSym(M,Q,q,LAEMIndex,proj.n-proj.d-1,proj.n-proj.d-1);
        report.dp1eval=q[0];
        report.ld_errors.resize(n); report.ld_errors=0.0;
        double dij;
        for (unsigned long i=0; i<n; i++)
        {
            for (unsigned long j=0; j<i; j++)
            {
                if (opts.mode==MDS)
                { 
                    dij=opts.metric->dist(&(const_cast<FMatrix<double>&>(proj.p)(i,0)), &(const_cast<FMatrix<double>&>(proj.p)(j,0)), proj.d); 
                }
                else if (opts.mode==SMDS)
                {
                    dij=0.0; for (unsigned long h=0; h<=proj.d; ++h) dij+=Q(i,h)*Q(j,h);
                    dij=acos(dij)*sr;
                }
                
                dij=(dij-dist(i,j));
                dij*=dij;
                report.ld_errors[i]+=dij;
                report.ld_errors[j]+=dij;
            }
        }
        report.ld_errors*=(1.0/n);
        report.ld_error=report.ld_errors.sum()*(1.0/n);
    }
    } //ends if (opts.mode!=TMDS)
}

void NLDRITERChi::set_vars(const std::valarray<double>& rv) 
{ 
    if (coords.size()!=rv.size()) 
    {
        coords.resize(rv.size()); 
        pgrad.resize(rv.size());
    }
    coords=rv;
    
    FMatrix<double> ld(n,n); ld*=0.0;
    for (unsigned long i=0; i<n; i++)
        for (unsigned long j=0; j<i; j++)
    { ld(j,i)=ld(i,j)=metric->dist(&coords[i*d], &coords[j*d],d);  }
    
    //std::cerr<<"LOW-DIM DISTANCES" <<ld<<"\n";
    double fld, dfld, gij, wij, tw=0.0;
    pval=0.0; pgrad=0.0; 
    double dcw=1.0;
    
    for (unsigned long i=0; i<n; i++)
        for (unsigned long j=0; j<i; j++) 
    { 
        tfun.fdf(ld(i,j),fld,dfld);
        
        wij=weights[i]*weights[j]*dcw;
        tw+=wij;
        pval+=((fhd(i,j)-fld)*(fhd(i,j)-fld)*(1.0-imix)+imix*(hd(i,j)-ld(i,j))*(hd(i,j)-ld(i,j)))  *wij;
        gij=((fhd(i,j)-fld)*dfld*(1.0-imix)+imix*(hd(i,j)-ld(i,j))) /ld(i,j)*wij;

        for (unsigned long h=0; h<d; h++)
        {
            pgrad[i*d+h]+=gij*(coords[i*d+h]-coords[j*d+h]);
            pgrad[j*d+h]-=gij*(coords[i*d+h]-coords[j*d+h]);
        }
    }
    pval*=1.0/tw;
    pgrad*=1.0/tw;
}

void NLDRITERChi::get_value(double& rv) const
{
    rv=pval;
}

void NLDRITERChi::get_gradient(std::valarray<double>& rv) const
{
    if (rv.size()!=pgrad.size()) rv.resize(pgrad.size());
    rv=pgrad;
}

void compute_chi1(const std::valarray<double>& x, const std::valarray<double>& w, const double& imix, const NLDRFunction& tfun, FMatrix<double>& md, FMatrix<double>& mfd, FMatrix<double> p, long skipi, double& vv, std::valarray<double>& vg)
{
    unsigned long n=p.rows(), d=x.size();
    std::valarray<double> vx(x);
    vv=0.0; vg.resize(d); vg=0.0;
    std::valarray<double> v1(d), vd(n), vfd(n); double ld, lfd, ldfd, tw=0.0;
    
    double diffdist, imix2=1.0-imix;
    vfd=mfd.row(skipi); vd=md.row(skipi);

    for (unsigned long i=0; i<n; ++i)
    {
        if (i==skipi) continue; 
        //test: ONLY for global min we weight more the closer points. what happens?
        //ndw=1.0/vfd[i]/vd[i];  nothing particularly good!
        
        
        v1=p.row(i); v1-=vx; ld=0.0; for (unsigned long j=0; j<d; ++j) ld+=v1[j]*v1[j]; ld=sqrt(ld);

        if (ld<=0.0) continue;
        
        tfun.fdf(ld, lfd, ldfd); 

        //std::cerr<<vd[i]<<" ("<<vfd[i]<<")"<<" >> "<<ld<<" ("<<lfd<<")"<< "\n";
        diffdist=(vfd[i]-lfd); 
        vv+=(diffdist*diffdist*imix2+imix*(vd[i]-ld)*(vd[i]-ld))*w[i];

        v1*=2.0*((diffdist)*ldfd*imix2+imix*(vd[i]-ld))/ld*w[i];
        tw+=w[i];
        vg+=v1;
    }
    vv*=1.0/tw; vg*=1.0/tw; 
}


void NLDRITER(FMatrix<double>& points, NLDRProjection& proj, const NLDRITEROptions& opts, NLDRITERReport& report)
{
        
    if (opts.metric==NULL) ERROR("Uninitialized metric pointer\n");
    proj.P=points; proj.D=points.cols(); proj.d=opts.lowdim; proj.n=points.rows();
    proj.p.resize(proj.n,proj.d);
    
    NLDRITERChi chiobj; chiobj.hd.resize(proj.n,proj.n);
    chiobj.fhd.resize(proj.n,proj.n); chiobj.n=proj.n; chiobj.d=proj.d;
    chiobj.metric=new NLDRMetricEuclid; chiobj.tfun=opts.tfunL; 
    chiobj.weights.resize(proj.n); if (opts.weights.size()==0) chiobj.weights=1.0; else chiobj.weights=opts.weights; chiobj.imix=opts.imix;
    
    std::valarray<double> pcoords(proj.n*proj.d); 
    FMatrix<double> distHD(chiobj.fhd), distLD(proj.n,proj.n);
    
    std::cerr<<"Building distance matrix\n";
    for (unsigned long i=0; i<proj.n; i++)
    {
        chiobj.fhd(i,i)=0;
        for (unsigned long j=0; j<i; j++) {
            chiobj.fhd(i,j)=chiobj.fhd(j,i)=opts.metric->dist(&(const_cast<FMatrix<double>&>(points)(i,0)), &(const_cast<FMatrix<double>&>(points)(j,0)), proj.D); 
        }
    }
    chiobj.hd=chiobj.fhd; //!TEST
    
    if (opts.ipoints.size()==0)
    {
        std::cerr<<"Initializing low-dim points from MDS\n";
        NLDRMDSOptions mdsopts; NLDRMDSReport mdsreport;
        mdsopts.lowdim=proj.d; mdsopts.verbose=false; 
        mdsopts.mode=MDS; mdsopts.metric=opts.metric;
        NLDRMDS(points, proj, mdsopts, mdsreport, chiobj.fhd);
    }
    else proj.p=opts.ipoints;
    
    for (unsigned long i=0; i<proj.n; i++) 
        for (unsigned long h=0; h<proj.d; h++) pcoords[i*proj.d+h]=proj.p(i,h);

    for (unsigned long i=0; i<proj.n; i++)
        for (unsigned long j=0; j<i; j++)  { chiobj.fhd(j,i)=chiobj.fhd(i,j)=opts.tfunH.f(chiobj.fhd(i,j));
       // std::cerr<<chiobj.hd(i,j)<<", "<<chiobj.fhd(i,j)<<", "<<chiobj.dfd(chiobj.hd(i,j))<<"\n"; 
        }
    chiobj.set_vars(pcoords); 
    IterOptions<double,2> iops=IterOptions<double,2>( 
                               opts.steps,
                               fixarray<double,2>(1e-5, 1e-5),fixarray<double,2>(0.,0.),
                               fixarray<double,2>(ichk_change, ichk_default));
    double ferr; 
    
    AnnealingOptions saop(opts.saopts);
    ConjGradOpts cgop(opts.cgopts);
    if (cgop.maxiter==0) cgop.maxiter=opts.steps; 
    if (saop.steps==0) saop.steps=opts.steps; 
    
    std::valarray<double> rpos; double rvalue;
    
    if (opts.global)
    {
        std::cerr<<"Preliminary opt\n";
        for (unsigned long i=0; i<proj.n; i++) 
            for (unsigned long h=0; h<proj.d; h++) pcoords[i*proj.d+h]=proj.p(i,h);
        min_conjgrad(chiobj,pcoords,rpos, rvalue, cgop );
        chiobj.get_vars(pcoords);
        for (unsigned long i=0; i<proj.n; i++) 
            for (unsigned long h=0; h<proj.d; h++) proj.p(i,h)=pcoords[i*proj.d+h];
        
        std::cerr<<"Pointwise global minimization\n";
        for (unsigned long ip=0; ip<proj.n; ip++)
        {
            std::cerr<<"Minimizing pt. "<<ip<<", weight: "<<chiobj.weights[ip]<<"\n";
            unsigned long ngrid=opts.grid1; double wgrid=opts.gridw;
            if (proj.d!=2) ERROR ("Integral projector is implemented only for 2D");
            std::valarray<double> gx(ngrid), gy(ngrid);
            FMatrix<double> gridU(ngrid,ngrid); std::valarray<double> x(proj.d), rg(proj.d);
            FTensor<double, 3> gridDU(fixarray<unsigned long, 3>(ngrid,ngrid,proj.d));

            for (unsigned long i=0; i<ngrid; i++) { gx[i]=gy[i]=(-1.0+i*2.0/(ngrid-1)); }
            gx*=wgrid; gy*=wgrid; 
            
            /*
            x[0]=15.0; x[1]=15.0;  std::cerr<<" printing at rim\n"; double rf;
            compute_chi1(x,chiobj.weights,opts.tfunL,chiobj.hd,chiobj.fhd,proj.p,ip,rf,rg);
            x[0]=0; x[1]=0;  std::cerr<<" printing at center\n";
            compute_chi1(x,chiobj.weights,opts.tfunL,chiobj.hd,chiobj.fhd,proj.p,ip,rf,rg);
            exit(1);
            */
            
            for (unsigned long i=0; i<ngrid; i++)
                for (unsigned long j=0; j<ngrid; j++)
            {
                x[0]=gx[i]; x[1]=gy[j];
                compute_chi1(x,chiobj.weights,chiobj.imix,opts.tfunL,chiobj.hd,chiobj.fhd,proj.p,ip,gridU(i,j),rg);
                //std::cerr<<x[0]<<","<<x[1]<<">>"<<gridU(i,j)<<":"<<rg[0]<<","<<rg[1]<<"\n";
                gridDU(i,j,0)=rg[0];  gridDU(i,j,1)=rg[1];
            }

            //std::cerr<<gridU<<" GRIDU\n"; 
            InterpolateBicubic bicub; 
            bicub.set_table(gx,gy,gridU,gridDU); 
    
            ngrid=opts.grid2; gx.resize(ngrid); gy.resize(ngrid); 
            for (unsigned long i=0; i<ngrid; i++) { gx[i]=gy[i]=(-1.0+i*2.0/(ngrid-1)); }
            gx*=wgrid; gy*=wgrid;
            
            double f; fixarray<double,2> df,x0,xx,minx; double initf, minchi;
            
            
            x0[0]=x[0]=proj.p(ip,0); x0[1]=x[1]=proj.p(ip,1); 
            compute_chi1(x,chiobj.weights,chiobj.imix,opts.tfunL,chiobj.hd,chiobj.fhd,proj.p,ip,initf,rg);
            minchi=initf;

            
            std::cerr<<"Initial value: "<<proj.p(ip,0)<<","<<proj.p(ip,1)<<"   >>  "<<minchi<<"\n";
            bool fmoved=false;
            for (unsigned long i=0; i<ngrid; i++)  for (unsigned long j=0; j<ngrid; j++)
            {
                xx[0]=gx[i]; xx[1]=gy[j];
                bicub.get_ydy(xx,f,df);
                if (f<minchi) { fmoved=true; minchi=f; minx=xx; }
            }
            
            if (!fmoved) continue; // do not optimize if it was not moved.
            //double checks on uninterpolated functionx[0]=gx[i]; x[1]=gy[j];
            x[0]=minx[0]; x[1]=minx[1];
            compute_chi1(x,chiobj.weights,chiobj.imix,opts.tfunL,chiobj.hd,chiobj.fhd,proj.p,ip,f,rg);
            if (f>=initf) { std::cerr<<"False positive due to interpolation\n"; continue; } // it's possible that a (very small) decrease was due to errors in interpolant
            proj.p(ip,0)=minx[0]; proj.p(ip,1)=minx[1];
            
            std::ofstream restf((std::string("global.")+int2str(ip)).c_str());
            for (unsigned long i=0; i<proj.n; i++) 
            {  restf<<proj.p(i,0)<<" "<<proj.p(i,1)<<"\n"; }
            restf.close();
            std::cerr<<"Globally optimized to: "<<proj.p(ip,0)<<","<<proj.p(ip,1)<<"   >>  "<<minchi<<"\n";
            if (sqrt((x0[0]-minx[0])*(x0[0]-minx[0])+(x0[1]-minx[1])*(x0[1]-minx[1]))>3.0) 
            {
                std::cerr<<" *** LARGE DISPLACEMENT ***\n";
                std::cerr<<"printing interpolated function\n";
/*/                std::ofstream intfile((std::string("interpolation.")+int2str(ip)).c_str()); 
                for (unsigned long i=0; i<ngrid; i++)  { for (unsigned long j=0; j<ngrid; j++)
                    {
                        xx[0]=gx[i]; xx[1]=gy[j];
                        bicub.get_ydy(xx,f,df);
                        intfile<<xx[0]<<" "<<xx[1]<<" "<<f<<"\n";
                    } intfile<<"\n";}  }
                intfile.close(); //*/ 
            }
            //intfile.close(); ERROR("BREALK");
            for (unsigned long i=0; i<proj.n; i++) 
                for (unsigned long h=0; h<proj.d; h++) pcoords[i*proj.d+h]=proj.p(i,h);
            min_conjgrad(chiobj,pcoords,rpos, rvalue, cgop );
            chiobj.get_vars(pcoords);
            for (unsigned long i=0; i<proj.n; i++) 
                for (unsigned long h=0; h<proj.d; h++) proj.p(i,h)=pcoords[i*proj.d+h];
        }
    }
    else
    {
        switch (opts.minmode)
        {
            case NLDRCGradient:
                min_conjgrad(chiobj,pcoords,rpos, rvalue, cgop );
                break;
            case NLDRSimplex:
                min_simplex(chiobj,make_simplex(pcoords,opts.simplex_spread, opts.simplex_mult),pcoords,ferr, iops);
                break;
            case NLDRAnnealing:
                sim_annealing(chiobj,pcoords,rpos,rvalue,saop);
                break;
        }
    }
    
    chiobj.get_vars(pcoords);
    for (unsigned long i=0; i<proj.n; i++) 
        for (unsigned long h=0; h<proj.d; h++) proj.p(i,h)=pcoords[i*proj.d+h];
    
    chiobj.get_value(report.ld_error);
    
    for (unsigned long i=0; i<proj.n; i++)
        for (unsigned long j=0; j<i; j++) distLD(j,i)=distLD(i,j)=chiobj.metric->dist(&proj.p(i,0),&proj.p(j,0),proj.d);
    
    if (opts.verbose) {
        std::cerr<< " Computing detailed error report\n";
        report.ld_errors.resize(proj.n); report.ld_errors=0.0;
        double dij, tw=0.0, tww=0.0; 
        for (unsigned long i=0; i<proj.n; i++)
        {
            for (unsigned long j=0; j<i; j++)
            {
                dij=(chiobj.fhd(i,j)-opts.tfunL.f(distLD(i,j)));
                dij=(dij*dij*(1.0-chiobj.imix)
                        +chiobj.imix*(chiobj.hd(i,j)-distLD(i,j))*(chiobj.hd(i,j)-distLD(i,j)))
                        *chiobj.weights[i]*chiobj.weights[j];
                tww+=chiobj.weights[i]*chiobj.weights[j];
                report.ld_errors[i]+=dij;
                report.ld_errors[j]+=dij;
            }
        }
        tw=chiobj.weights.sum();
        report.ld_error=report.ld_errors.sum()*(0.5/tww);
        for (unsigned long i=0; i<proj.n; i++) report.ld_errors[i]*=(1.0/(tw-chiobj.weights[i])/chiobj.weights[i]);
    }
    delete chiobj.metric;
}



// void NLDRLLTE(FMatrix<double>& points, NLDRProjection& proj, const NLDRLLTEOptions& opts, NLDRLLTEReport& report)
// {
//     proj.P=points; proj.D=points.cols(); proj.d=opts.lowdim; proj.n=points.rows();
//     proj.p.resize(proj.n,proj.d);
//     
//     std::cerr<<"Building neighbor list\n";
//     NLDRNeighborList nlist(proj.P,opts.nlopts);
//     //std::cerr<<proj.P<<"BEFORE \n";
//     if (opts.rmlonesome) 
//     {
//         std::valarray<unsigned long> ilone;
//         nlist.RemoveLonesome(ilone);
//         //removing points with no connections
//         unsigned long nlone=ilone.size();
//             
//         if (nlone>0) 
//         {
//             std::cerr<<"Removing isolated points\n";
//             proj.P.resize(proj.n-nlone,proj.D);
//             unsigned long k=0, j=0;
//             for (unsigned long i=0; i<proj.n; ++i) 
//             {
//                 if (j<nlone && ilone[j]==i) { ++j; continue; }
//                 for (unsigned long h=0; h<proj.D; ++h) proj.P(k,h)=points(i,h);
//                 ++k;
//             }
//             proj.n-=nlone;
//             proj.p.resize(proj.n,proj.d);
//         }
//     }
//     
//     std::cerr<<"Building weights\n";
//     std::valarray<std::valarray<double> > weights(proj.n);
//     std::valarray<std::valarray<unsigned long> > indices(proj.n);
//     
//     if (opts.verbose) { report.hd_errors.resize(proj.n); report.hd_errors=0.0;  report.ld_errors.resize(proj.n); report.ld_errors=0.0; }
//     report.ld_error=report.hd_error=0.0;
//     
//     for (unsigned long i=0; i<proj.n; ++i)
//     {
//         unsigned long m=nlist.nneigh(i);
//         if (m==0) ERROR("Point "<<i<<" is isolated!\n");
//         FMatrix<double> G(proj.D,m), GT, H, HT, C, C1; 
//         
//         for (unsigned long h=0; h<proj.D; ++h) for (unsigned long j=0; j<m; ++j) G(h,j)=points(nlist.index(i,j),h)-points(i,h);
//         
//         //First, finds d principal components
//         transpose(G,GT); mult(G,GT,C);
//         FMatrix<double> P, PT; std::valarray<double> p;
//         EigenSolverSym(C,P,p,LAEMIndex,proj.D-proj.d,proj.D-1);
// 
//         transpose(P,PT); mult(PT,G,H);  transpose(H, HT);
//         mult(HT,H,C);  // now C is projected on the most relevant eigenvalues
//         
//         //smoothens the matrix
//         double tr=0.0; 
//         if (opts.smooth<0) { tr=normfrob(C)*-opts.smooth; }
//         else tr=opts.smooth;
//         for (unsigned long j=0; j<m; ++j) C(j,j)+=tr;
//         
//         MatrixInverse(C,C1); 
//         
//         double beta, lambda;
//         beta=0.0;
//         for (unsigned long j=0; j<m; ++j) for (unsigned long k=0; k<m; ++k) beta+=C1(j,k);
//         lambda=1.0/beta;
// 
//         weights[i].resize(m); indices[i].resize(m);
//         for (unsigned long j=0; j<m; ++j)
//         {
//             indices[i][j]=nlist.index(i,j);  weights[i][j]=0.0; 
//             for (unsigned long k=0; k<m; ++k) weights[i][j]+=C1(j,k);
//         }
//         weights[i]*=lambda;
// 
//         for (unsigned long h=0; h<proj.d; ++h) 
//         { 
//             double ed=0.0;
//             for (unsigned long j=0; j<m; ++j) ed+=weights[i][j]*H(h,j);
//             ed*=ed;  report.hd_error+=ed; if (opts.verbose) report.hd_errors[i]+=ed*ed;
//         }
//         //should also compute TRUE HD errors, not just the projected ones!
//         if (opts.verbose) report.hd_errors[i]=sqrt(report.hd_errors[i]);
//     }
//     report.hd_error=sqrt(report.hd_error/proj.n);
//     
//     //we have the weight matrix as a sparse matrix
//     std::cerr<<"Building sparse W matrix\n";
//     CrsMatrix<double> W(proj.n,proj.n,indices,weights), WT;
//     W*=-1.0;  for (unsigned long i=0; i<proj.n; ++i) W(i,i)+=1.0;
// 
// 
//     std::cerr<<"Finding low-dim points\n";
//     transpose(W,WT); CrsMatrix<double> M; mult(WT,W,M); 
//     FMatrix<double> FM(M), Q; std::valarray<double> q;
//     EigenSolverSym(FM,Q,q,LAEMIndex,1,proj.d);
//     for (unsigned long i=0; i<proj.n; ++i) for (unsigned long h=0; h<proj.d; ++h) 
//             proj.p(i,h)=Q(i,h);
//     
//     report.deval.resize(proj.d); report.deval=q;
//     if (opts.verbose) {
//         EigenSolverSym(FM,Q,q,LAEMIndex,proj.d+1,proj.d+1);
//         report.dp1eval=q[0];
//     }
//     
//     for (unsigned long i=0; i<proj.n; ++i)
//     {
//         unsigned long m=nlist.nneigh(i);
//         for (unsigned long h=0; h<proj.d; ++h) 
//         { 
//             double ed=0.0;
//             for (unsigned long j=0; j<m; ++j) ed+=weights[i][j]*proj.p(nlist.index(i,j),h);
//             ed-=proj.p(i,h); ed*=ed;
//             report.ld_error+=ed; if (opts.verbose) report.ld_errors[i]+=ed*ed;
//         }
//         if (opts.verbose) report.ld_errors[i]=sqrt(report.ld_errors[i]);
//     }
// }

} //ends namespace toolbox
