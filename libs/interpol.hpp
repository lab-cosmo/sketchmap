#ifndef __INTERPOL_H
#define __INTERPOL_H
#include "tbdefs.hpp"
#include "matrix-full.hpp"
#include "tensor-full.hpp"

namespace toolbox {
    
class InterpolateBicubic;

class OneDGauge {
private:
    std::valarray<double> xlist;
    long search1(const double& x);
    long jold, n;

public:
    double& operator[] (const unsigned long i)      { return xlist[i]; }
    double operator[] (const unsigned long i) const { return xlist[i]; }
    long operator() (const double x) { return search1(x); }
    
    void set_table(const std::valarray<double>& nlist) { xlist.resize(n=nlist.size()); xlist=nlist; jold=0; }
    OneDGauge(const std::valarray<double>& nlist=std::valarray<double>()) : xlist(nlist) { n=xlist.size(); jold=0; }
    
    OneDGauge(const OneDGauge& ng) { xlist.resize(n=ng.xlist.size()); xlist=ng.xlist; jold=ng.jold; }
    OneDGauge& operator =(const OneDGauge& ng) { if (&ng==this) return *this; xlist.resize(n=ng.xlist.size()); xlist=ng.xlist; jold=ng.jold;  return *this; }
};

/* ONE-D SPLINE INTERPOLATION */
class InterpolateSpline {
friend class InterpolateBicubic;
private:
    OneDGauge xgauge;
    std::valarray<double> ylist;
    std::valarray<double> y2list;

    unsigned long n;
    double intrp(long j, const double& x);
    double dintrp(long j, const double& x);
    
public:
    void set_table(std::valarray<double>& nxlist, std::valarray<double>& nylist);
    inline double operator() (const double& x) { return intrp(xgauge(x),x); }
    inline double d(const double& x)           { return dintrp(xgauge(x),x); }
    
    InterpolateSpline() : xgauge(), n(0) {}
    InterpolateSpline(std::valarray<double>& nxlist, std::valarray<double>& nylist) { set_table(nxlist,nylist); }
};

inline double InterpolateSpline::intrp(long j, const double& x)
{
    double xh,xl,h,b,a;
    xh=xgauge[j+1]; xl=xgauge[j];
    h=xh-xl; a=(xh-x)/h; b=(x-xl)/h;
    return a*ylist[j]+b*ylist[j+1]+(a*(a*a-1.0)*y2list[j]+b*(b*b-1.0)*y2list[j+1])*h*h/6.0;
}

inline double InterpolateSpline::dintrp(long j, const double& x)
{
    double xh,xl,h,b,a;
    xh=xgauge[j+1]; xl=xgauge[j];
    h=xh-xl; a=(xh-x)/h; b=(x-xl)/h;
    return (ylist[j+1]-ylist[j])/h+(-(3.0*a*a-1.0)*y2list[j]+(3.0*b*b-1.0)*y2list[j+1])*h/6.0;
}

/* ND BICUBIC SPLINE INTERPOLATION */
class InterpolateBicubic {
    private:
        fixarray<long, 2> n;
        fixarray<OneDGauge, 2u> xlist;
        FTensor<double, 4> clist;  
        InterpolateSpline is1, is2;

        double intrp(const std::valarray<long>& j, const std::valarray<double>& x);
        void dintrp(const fixarray<unsigned long, 2>& j, const fixarray<double, 2>& x, double& y, fixarray<double, 2>& dy);
    
    public:
        
        void set_table(const std::valarray<double>& nx1list, const std::valarray<double>& nx2list, const FMatrix<double>& nylist, const FTensor<double,3>& ny1list=(FTensor<double,3>()));
        inline double operator() (const std::valarray<double>& x) { }//return intrp(search1(x),x); }
        inline void get_ydy(const fixarray<double, 2>& x, double& y, fixarray<double, 2>& dy) 
        { dintrp(fixarray<unsigned long, 2>(xlist[0](x[0]), xlist[1](x[1])), x, y, dy); }

/*        InterpolateSpline() : jold(0), n(0) {}
        InterpolateSpline(std::valarray<double>& nxlist, std::valarray<double>& nylist) { set_table(nxlist,nylist); }*/
};

}; //ends namespace toolbox
#endif //ends __INTERPOL_H

