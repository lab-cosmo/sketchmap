/* Sketch-map library header file
   --------------------------------------------------
   Author: Michele Ceriotti, 2011
   Distributed under the GNU General Public License  
*/

#ifndef __DIMREDUCE_H
#define __DIMREDUCE_H 0
#include "tbdefs.hpp"
#include "minsearch.hpp"
#include "matrix-full.hpp"
#include "matrix-crs.hpp"

namespace toolbox {

enum NLDRFunctionMode { NLDRIdentity, NLDRSigmoid, NLDRCompress, NLDRXSigmoid, NLDRGamma, NLDRWarp };
class NLDRFunction {
    typedef double (NLDRFunction::*NLDRFP)(double) const;
    typedef void (NLDRFunction::*NLDRFPc)(double, double&, double&) const;
private:
    NLDRFP pf, pdf; NLDRFPc pfdf;
    NLDRFunctionMode pmode;
    std::valarray<double> pars;

    double nldr_identity(double x) const;
    double nldr_didentity(double x) const;
    void nldr_identity(double x, double& rf, double& rdf) const;
    double nldr_sigmoid(double x) const;
    double nldr_dsigmoid(double x) const;
    void nldr_sigmoid(double x, double& rf, double& rdf) const;
    double nldr_xsigmoid(double x) const;
    double nldr_dxsigmoid(double x) const;
    void nldr_xsigmoid(double x, double& rf, double& rdf) const;
    double nldr_compress(double x) const;
    double nldr_dcompress(double x) const;
    void nldr_compress(double x, double& rf, double& rdf) const;

    double nldr_gamma(double x) const;
    double nldr_dgamma(double x) const;
    void nldr_gamma(double x, double& rf, double& rdf) const;

    double g(double y) const;
    double dg(double y) const;    
    double nldr_warp(double x) const;
    double nldr_dwarp(double x) const;
    void nldr_warp(double x, double& rf, double& rdf) const;

public:
    
    NLDRFunction(NLDRFunctionMode nmode=NLDRIdentity, const std::valarray<double>& npars=std::valarray<double>()) { set_mode(nmode,npars); }
    NLDRFunction(const NLDRFunction& nf); 
    NLDRFunction& operator=(const NLDRFunction& nf); 
    
    void set_mode(NLDRFunctionMode mode, const std::valarray<double>& pars);
    inline double f(double x) const { return (this->*pf)(x); }
    inline double df(double x) const { return (this->*pdf)(x); };
    inline void fdf(double x, double& rf, double& rdf) const { (this->*pfdf)(x,rf, rdf); }
}; 

double nldr_identity(const double& x);
double nldr_didentity(const double& x);
double nldr_sigmoid(const double& x);
double nldr_dsigmoid(const double& x);
double nldr_xpsigmoid(const double& x);
double nldr_dxpsigmoid(const double& x);
double nldr_compress(const double& x);
double nldr_dcompress(const double& x);

class NLDRLLEOptions; class NLDRLLEReport; 
class NLDRMDSOptions; class NLDRMDSReport; 
class NLDRITEROptions; class NLDRITERReport; 
class NLDRNeighborList; class NLDRNeighborOptions;
class NLDRProjection;

/*! A collection of metric functions */
//null-metric class
class NLDRMetric {
    friend class NLDRNeighborList;
    private: 
        virtual double pdist(const double* a, const double* b, unsigned long n) const { return 0.0; };
        virtual void pdiff(const double* a, const double* b, double* c, unsigned long n) const { for (unsigned long i=0; i<n; ++i) c[i]=a[i]-b[i]; };
    public:
        inline double dist(const double* a, const double* b, unsigned long n)
        {
            return pdist(a,b,n);
        }
        double dist(const std::valarray<double>& a, const std::valarray<double>& b) const
        { 
#ifdef DEBUG
            if (a.size()!=b.size()) ERROR("Vector size mismatch in distance.");
#endif
            return pdist(&a[0],&b[0],a.size()); 
        }
        void diff(const std::valarray<double>& a, const std::valarray<double>& b, std::valarray<double>& c) const
        {
#ifdef DEBUG
            if (a.size()!=b.size()) ERROR("Vector size mismatch in distance.");
#endif
            c.resize(a.size()); pdiff(&a[0], &b[0], &c[0], a.size()); 
        }
};

class NLDRMetricEuclid: public NLDRMetric {
    private: 
        double pdist(const double* a, const double* b, unsigned long d) const;
};

class NLDRMetricPBC: public NLDRMetric {
    private: 
        double pdist(const double* a, const double* b, unsigned long d) const;
        void pdiff(const double* a, const double* b, double* c, unsigned long n) const;
    public:
        std::valarray<double> periods;
        NLDRMetricPBC() : periods() {}
        NLDRMetricPBC(const NLDRMetricPBC& no) : periods(no.periods) {}
        NLDRMetricPBC& operator= (const NLDRMetricPBC& no) 
        { if (&no==this) return *this; periods.resize(no.periods.size()); periods=no.periods; }
};

class NLDRMetricSphere: public NLDRMetric {
    private: 
        double pdist(const double* a, const double* b, unsigned long d) const;
        void pdiff(const double* a, const double* b, double* c, unsigned long n) const;
    public:
        std::valarray<double> periods;
        NLDRMetricSphere() : periods() {}
        NLDRMetricSphere(const NLDRMetricSphere& no) : periods(no.periods) {}
        NLDRMetricSphere& operator= (const NLDRMetricSphere& no) 
        { if (&no==this) return *this; periods.resize(no.periods.size()); periods=no.periods; }
};

enum NLDRNeighborGreediness { NLDRGreedy, NLDRLiberal, NLDRAsym }; 
class NLDRNeighborOptions{
    public:
        NLDRMetric *ometric;
        unsigned long kw;
        unsigned long maxneigh; double cutoff; NLDRNeighborGreediness greediness;
 
        NLDRNeighborOptions() : ometric(NULL), kw(0), maxneigh(4), greediness(NLDRGreedy), cutoff(0.0) {}
};

class NLDRNeighbor { public: unsigned long j; double d; NLDRNeighbor(unsigned long nj=0, double nd=0.0) :j(nj),d(nd) {}};
inline bool operator < (const NLDRNeighbor &lhs, const NLDRNeighbor &rhs) { return lhs.d< rhs.d; }

class NLDRNeighborList {
    friend void NLDRLLE(FMatrix<double>& points, NLDRProjection& proj, const NLDRLLEOptions& opts, NLDRLLEReport& report);
    private:
        NLDRNeighborOptions opts;
        std::valarray<NLDRNeighbor> nlist;
        std::valarray<unsigned long>  npoint;
        void nlbuildup(const FMatrix<double>& points);
        NLDRNeighbor& rneigh(unsigned long i, unsigned long j); 
    
    public:
        NLDRNeighborList& operator=(const NLDRNeighborList& nn)
        {
            if (&nn==this) return *this;
            nlist.resize(nn.nlist.size()); nlist=nn.nlist;
            npoint.resize(nn.npoint.size()); npoint=nn.npoint;
            opts=nn.opts;
        }
        NLDRNeighborList(const NLDRNeighborOptions& nopts=NLDRNeighborOptions()) : opts(nopts) {}
        NLDRNeighborList(const std::valarray<std::valarray<double> >& points, 
                         const NLDRNeighborOptions& nopts=NLDRNeighborOptions()) { Build(points, nopts); }
        NLDRNeighborList(const FMatrix<double>& points, 
                         const NLDRNeighborOptions& nopts=NLDRNeighborOptions()) : opts(nopts) { nlbuildup(points); }
        void Build(const FMatrix<double>& lpoints) { nlbuildup(lpoints); }
        void Build(const FMatrix<double>& lpoints, const NLDRNeighborOptions& nopts) { opts=nopts; nlbuildup(lpoints); }
        void RemoveLonesome(std::valarray<unsigned long>& llone);
    
        unsigned long size() const {return npoint.size()-1;}
        unsigned long nneigh(unsigned long i) const 
        {
#ifdef DEBUG
            if (i>=npoint.size()) ERROR("Index out of bounds in neighbor list");
#endif
            return npoint[i+1]-npoint[i]; 
        } 
        unsigned long index(unsigned long i, unsigned long j) const;
        double dist(unsigned long i, unsigned long j) const;  
        NLDRNeighbor neigh(unsigned long i, unsigned long j) const; 
};

class NLDROptions {
public:
    NLDRNeighborOptions nopts;
    double acutoff, gwidth, gtemp;
    unsigned long grid1, grid2, cgsteps;
    
    NLDRFunction tfunH, tfunL;
    NLDROptions() : nopts(), acutoff(2.0), gwidth(20.0), gtemp(0.001), grid1(21), grid2(201), cgsteps(0), tfunH(), tfunL() {}
}; 

class NLDRProjection {
    friend void NLDRLLE(FMatrix<double>& points, NLDRProjection& proj, const NLDRLLEOptions& opts, NLDRLLEReport& report);
    friend void NLDRMDS(FMatrix<double>& points, NLDRProjection& proj, const NLDRMDSOptions& opts, NLDRMDSReport& report, const FMatrix<double>& outd);
    friend void NLDRITER(FMatrix<double>& points, NLDRProjection& proj, const NLDRITEROptions& opts, NLDRITERReport& report,  const FMatrix<double>& outd);
    friend void NLDRIProj(const NLDRProjection& proj, const std::valarray<double>& X, std::valarray<double>& x); 
    
private:
    
    FMatrix<double> dmu, lmds_mp;
    
    NLDROptions opts;
    NLDRNeighborList neigh;
    unsigned long D, d, n;
    FMatrix<double> P, p; std::valarray<double> w;
    std::valarray<FMatrix<double> > HV, LV, PM;
    void calc_PM();
    bool ftainted;
    
    //new out-of-set
    std::valarray<NLDRNeighbor> nd; 
    FMatrix<double> vxx, vdXx; 
    std::valarray<double> vx, vfd, vf1d, vg, vdX; double vv;
public:
    std::string interp_out;
    void get_vars(std::valarray<double>& rv) const { rv.resize(vx.size()); rv=vx; }
    void set_vars(const std::valarray<double>& rv);  
    void get_value(double& rv) const { rv=vv; }
    void get_gradient(std::valarray<double>& rv) const {rv.resize(vg.size()); rv=vg; }
    

public:
    void set_options(const NLDROptions& nopts)
    { opts=nopts; ftainted=true; }
    void get_options(NLDROptions& nopts)
    { nopts=opts; }
    void set_points(const std::valarray<std::valarray<double> >& nP, const std::valarray<std::valarray<double> >& np, const std::valarray<double>& nw=std::valarray<double>(0))
    {
        if ((n=nP.size())<1 || (D=nP[0].size())<1) ERROR("Hi-dimensional array has inconsistent sizes.");
        if (np.size()!=n || (d=np[0].size())<1 || d>D) ERROR("Low-dimensional array has inconsistent sizes.");
        P=nP; p=np; ftainted=true;
        w.resize(n); if (nw.size()==0) w=1.0; else w=nw;
    }
    
    void set_points(const FMatrix<double>& nP, const FMatrix<double>& np, const std::valarray<double>& nw=std::valarray<double>(0))
    {
        if ((n=nP.rows())<1 || (D=nP.cols())<1) ERROR("Hi-dimensional array has inconsistent sizes.");
        if (np.rows()!=n || (d=np.cols())<1 || d>D) ERROR("Low-dimensional array has inconsistent sizes.");
        P=nP; p=np; ftainted=true;
        w.resize(n); if (nw.size()==0) w=1.0; else w=nw;
    }
    
    void get_points(std::valarray<std::valarray<double> >& nP, std::valarray<std::valarray<double> >& np) 
    { 
        nP.resize(n); np.resize(n); 
        for (unsigned long i=0; i<n; ++i) { nP[i].resize(D); for (unsigned long h=0; h<D; ++h) nP[i][h]=P(i,h); }
        for (unsigned long i=0; i<n; ++i) { np[i].resize(d); for (unsigned long h=0; h<d; ++h) np[i][h]=p(i,h); }
    }
    
    double project(const std::valarray<double>& np, std::valarray<double>& hp, std::valarray<double>& lp, double &md);
};

enum NLDRLLEMode { LLE, LLTE, HLLE }; 
class NLDRLLEOptions
{
public:
    NLDRNeighborOptions nlopts; NLDRLLEMode mode; bool verbose; 
    unsigned long lowdim, dimts; double smooth; bool rmlonesome;
    NLDRLLEOptions() : nlopts(), mode(LLE), verbose(false), 
                   lowdim(2), dimts(2), smooth(-1e-5), rmlonesome(false) {}
};

class NLDRLLEReport
{
public:
    std::valarray<double> hd_errors; double hd_error;
    std::valarray<double> ld_errors; double ld_error;
    std::valarray<double> deval;  double dp1eval;

    NLDRLLEReport& operator=(const NLDRLLEReport& nr)
    {
        if (&nr==this) return *this;
        hd_errors.resize(nr.hd_errors.size()); hd_errors=nr.hd_errors;  hd_error=nr.hd_error;
        ld_errors.resize(nr.ld_errors.size()); ld_errors=nr.ld_errors;  ld_error=nr.ld_error;
        deval.resize(nr.deval.size()); deval=nr.deval;  dp1eval=nr.dp1eval;
        return *this;
    }
};

enum NLDRMDSMode { MDS, SMDS, TMDS }; 
class NLDRMDSOptions
{
public:
    NLDRMetric *metric; NLDRMDSMode mode; bool verbose; 
    unsigned long lowdim; 
    NLDRMDSOptions() : metric(NULL), mode(MDS), verbose(false), lowdim(2) {}
};

class NLDRMDSReport
{
public:
    std::valarray<double> ld_errors; double ld_error;
    std::valarray<double> deval;  double dp1eval;
};


enum NLDRIterMin { NLDRSimplex, NLDRCGradient, NLDRAnnealing, NLDRParatemp };

class NLDRITEROptions
{
public:
    NLDRFunction tfunH, tfunL;
    NLDRMetric *metric; bool verbose, global; 
    unsigned long grid1, grid2; double gridw, imix;
    unsigned long lowdim, steps; 
    std::valarray<double> weights; FMatrix<double> dweights;
    FMatrix<double> ipoints;
    NLDRIterMin minmode;
    AnnealingOptions saopts;
    ConjGradOpts cgopts;
    ParaOptions ptopts;
    double simplex_spread, simplex_mult;
    
    NLDRITEROptions() : tfunH(NLDRIdentity), tfunL(NLDRIdentity), metric(NULL), verbose(false), 
                   lowdim(2), global(false), grid1(11), grid2(101), gridw(20.0), imix(0.0), ipoints(),
                   minmode(NLDRCGradient), weights(0), dweights(0,0) 
                   {
                       saopts.temp_init=1e-4; saopts.temp_final=1e-20;
                       saopts.steps=0; saopts.mc_step=1e-1; saopts.adapt=1.05; saopts.drnd=0.2;

                       ptopts.temp_init=5e-7; ptopts.temp_final=2e-9;  ptopts.temp_factor=2.0; 
                       ptopts.steps=0; ptopts.replica=6; ptopts.dt=1.0; ptopts.tau=10;
                       ptopts.tau_avg=20; 

                       cgopts.maxiter=0;
                       cgopts.linesearch.maxiter=5; cgopts.linesearch.lstol=5e-10; 
                       simplex_spread=1.0; simplex_mult=1.1;
                   }
};

class NLDRITERReport
{
public:
    std::valarray<double> ld_errors; double ld_error;
};

class NLDRITERChi: public toolbox::FOnlyMinFunctionBase {
private:
    double pval; std::valarray<double> pgrad;
public:
    unsigned long n; unsigned long d; 
    double imix;
    NLDRFunction tfun;
    FMatrix<double> hd, fhd;
    NLDRMetric *metric;
    std::valarray<double> weights; FMatrix<double> dweights;
    void set_vars(const std::valarray<double>& rv); 
    void get_value(double& rv) const;
    void get_gradient(std::valarray<double>& rv) const;
    NLDRITERChi() : n(0), d(0), imix(0.0) {}
};

void NLDRLLE(FMatrix<double>& points, NLDRProjection& proj, const NLDRLLEOptions& opts, NLDRLLEReport& report);
void NLDRMDS(FMatrix<double>& points, NLDRProjection& proj, const NLDRMDSOptions& opts, NLDRMDSReport& report, const FMatrix<double>& outd=FMatrix<double>(0,0));
void NLDRITER(FMatrix<double>& points, NLDRProjection& proj, const NLDRITEROptions& opts, NLDRITERReport& report, const FMatrix<double>& outd=FMatrix<double>(0,0));
void NLDRIProj(const NLDRProjection& proj, const std::valarray<double>& X, std::valarray<double>& x);
}; //ends namespace toolbox
#endif //ends #ifndef __DIMREDUCE_H
