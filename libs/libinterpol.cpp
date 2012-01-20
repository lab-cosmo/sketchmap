/* A minimal library for function interpolation.
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/
   
#include "interpol.hpp"
#include "matrix-io.hpp"

namespace toolbox {

long OneDGauge::search1(const double& x)
{
    long jl=jold,ju,jm,inc=1;
    if (n<2) ERROR("Invalid size of interpolation arrays");
    if (jl<0 || jl>=n) { jl=0; ju=n-1; }
    else 
    {
        if (x>=xlist[jl]) 
        {
            while(true)
            {
                ju=jl+inc; 
                if (ju>=n-1) {ju=n-1; break;}
                else if(x<xlist[ju]) break;
                jl=ju; inc+=inc;
            }
        }
        else
        {
            ju=jl; 
            while(true)
            {
                jl=jl-inc;
                if (jl<=0) { jl=0; break; }
                else if(x>=xlist[jl]) break;
                ju=jl; inc+=inc;
            }
        }
    }
    while (ju-jl>1) 
    {
        jm=(ju+jl)/2;
        if (x>xlist[jm]) jl=jm; else ju=jm;
    }
    jold=jl;
    return jl;
}

void InterpolateSpline::set_table(std::valarray<double>& nxlist, std::valarray<double>& nylist)
{
    xgauge=OneDGauge(nxlist);  n=nxlist.size(); 
    if (n!=nylist.size()) ERROR("X and Y lists size mismatch in constructor\n");
    ylist.resize(n); ylist=nylist;
    
    std::valarray<double> u(n-1); y2list.resize(n);
    y2list[0]=u[0]=0.0;

    double p, qn, sig, un;
    for (unsigned long i=1; i<n-1; ++i)
    {
        sig=(xgauge[i]-xgauge[i-1])/(xgauge[i+1]-xgauge[i-1]);
        p=sig*y2list[i-1]+2.0;  y2list[i]=(sig-1)/p; 
        u[i]=(ylist[i+1]-ylist[i])/(xgauge[i+1]-xgauge[i])-(ylist[i]-ylist[i-1])/(xgauge[i]-xgauge[i-1]);
        u[i]=(6*u[i]/(xgauge[i+1]-xgauge[i-1])-sig*u[i-1])/p;
    }
    qn=un=0.0;
    y2list[n-1]=(un-qn*u[n-2])/(qn*y2list[n-2]+1.0);
    for (long k=n-2; k>=0; --k)
    {
        y2list[k]=y2list[k]*y2list[k+1]+u[k];
    }
}


void IBicCoeff(const fixarray<double,4>& y, const fixarray<double,4>& dy1, const fixarray<double,4>& dy2, const fixarray<double,4>& d2y12, const double d1, const double d2, FMatrix<double>&c) 
{
    static int wt_d[16*16]=
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
    2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
    0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
    -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
    9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
    -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
    2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
    -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
    4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};
    int l,k,j,i;
    double xx,d1d2=d1*d2;
    std::valarray<double> cl(16),x(16);
    static FMatrix<int> wt(16,16); l=0; for (i=0;i<16;i++)  for (j=0;j<16;j++) wt(i,j)=wt_d[l++];
    for (i=0;i<4;i++)  { x[i]=y[i]; x[i+4]=dy1[i]*d1; x[i+8]=dy2[i]*d2;  x[i+12]=d2y12[i]*d1d2;  }
    for (i=0;i<16;i++) {  xx=0.0;  for (k=0;k<16;k++) xx += wt(i,k)*x[k];  cl[i]=xx;        }
    l=0; for (i=0;i<4;i++)  for (j=0;j<4;j++) c(i,j)=cl[l++];
}


void InterpolateBicubic::set_table(const std::valarray<double>& nx1list, const std::valarray<double>& nx2list, const FMatrix<double>& nylist, const FTensor<double,3>& ny1list)
{
    fixarray<double,4> y,  dy1, dy2, d2y12;
    double d1, d2;
    n[0]=nx1list.size(); n[1]=nx2list.size(); 
    xlist[0]=OneDGauge(nx1list); xlist[1]=OneDGauge(nx2list);
    
    //create lists with derivatives, if not present
    
    FMatrix<double> ylist(nylist), dcross(nylist); 
    FTensor<double,3> dylist(ny1list);
    //FTensor<double,3> dylist;
    if (dylist.size()==0) 
    {
        
        dylist.resize(fixarray<unsigned long,3>(n[0],n[1],2)); 
        dylist(0,0,0)=(ylist(1,0)-ylist(0,0))/(xlist[0][1]-xlist[0][0]);
        dylist(0,0,1)=(ylist(0,1)-ylist(0,0))/(xlist[1][1]-xlist[1][0]);
        for (unsigned long i2=1; i2<n[1]-1; ++i2)
        {
            dylist(0,i2,0)=(ylist(1,i2)-ylist(0,i2))/(xlist[0][1]-xlist[0][0]);
            dylist(0,i2,1)=(ylist(0,i2+1)-ylist(0,i2-1))/(xlist[1][i2+1]-xlist[1][i2-1]);
        }
        dylist(0,n[1]-1,0)=(ylist(1,n[1]-1)-ylist(0,n[1]-1))/(xlist[0][1]-xlist[0][0]);
        dylist(0,n[1]-1,1)=(ylist(0,n[1]-1)-ylist(0,n[1]-2))/(xlist[1][n[1]-1]-xlist[1][n[1]-2]);
        for (unsigned long i1=1; i1<n[0]-1; ++i1)
        {
            dylist(i1,0,0)=(ylist(i1+1,0)-ylist(i1,0))/(xlist[0][i1+1]-xlist[0][i1]);
            dylist(i1,0,1)=(ylist(i1,1)-ylist(i1,0))/(xlist[1][1]-xlist[1][0]);
            for (unsigned long i2=1; i2<n[1]-1; ++i2)
            {
                dylist(i1,i2,0)=(ylist(i1+1,i2)-ylist(i1-1,i2))/(xlist[0][i1+1]-xlist[0][i1-1]);
                dylist(i1,i2,1)=(ylist(i1,i2+1)-ylist(i1,i2-1))/(xlist[1][i2+1]-xlist[1][i2-1]);
            }
            dylist(i1,n[1]-1,0)=(ylist(i1+1,n[1]-1)-ylist(i1,n[1]-1))/(xlist[0][i1+1]-xlist[0][i1]);
            dylist(i1,n[1]-1,1)=(ylist(i1,n[1]-1)-ylist(i1,n[1]-2))/(xlist[1][n[1]-1]-xlist[1][n[1]-2]);
        }
        dylist(n[0]-1,0,0)=(ylist(n[0]-1,0)-ylist(n[0]-2,0))/(xlist[0][n[0]-1]-xlist[0][n[0]-2]);
        dylist(n[0]-1,0,1)=(ylist(n[0]-1,1)-ylist(n[0]-1,0))/(xlist[1][1]-xlist[1][0]);
        for (unsigned long i2=1; i2<n[1]-1; ++i2)
        {
            dylist(n[0]-1,i2,0)=(ylist(n[0]-1,i2)-ylist(n[0]-2,i2))/(xlist[0][n[0]-1]-xlist[0][n[0]-2]);
            dylist(n[0]-1,i2,1)=(ylist(n[0]-1,i2+1)-ylist(n[0]-1,i2-1))/(xlist[1][i2+1]-xlist[1][i2-1]);
        }
        dylist(n[0]-1,n[1]-1,0)=(ylist(n[0]-1,n[1]-1)-ylist(n[0]-2,n[1]-1))/(xlist[0][n[0]-1]-xlist[0][n[0]-2]);
        dylist(n[0]-1,n[1]-1,1)=(ylist(n[0]-1,n[1]-1)-ylist(n[0]-1,n[1]-2))/(xlist[1][n[1]-1]-xlist[1][n[1]-2]);
        
        std::cerr<<"filling up list of second derivatives\n";
        dcross*=0.0; //on the boundaries, assume no curvature
        for (unsigned long i1=1; i1<n[0]-1; ++i1)
        {
            for (unsigned long i2=1; i2<n[1]-1; ++i2)
            {
                dcross(i1,i2)=(ylist(i1+1,i2+1)+ylist(i1-1,i2-1)-ylist(i1+1,i2-1)-ylist(i1-1,i2+1))/
                        ((xlist[0][i1+1]-xlist[0][i1-1])*(xlist[1][i2+1]-xlist[1][i2-1]));
            }
        }
    }
    else
    {
        dcross*=0.0; //on the boundaries, assume no curvature
        for (unsigned long i1=1; i1<n[0]-1; ++i1)
        {
            for (unsigned long i2=1; i2<n[1]-1; ++i2)
            {
                dcross(i1,i2)=(ylist(i1+1,i2+1)+ylist(i1-1,i2-1)-ylist(i1+1,i2-1)-ylist(i1-1,i2+1))/
                        ((xlist[0][i1+1]-xlist[0][i1-1])*(xlist[1][i2+1]-xlist[1][i2-1]));
            }
        }
    }
    
    clist.resize(fixarray<unsigned long,4>(n[0]-1,n[1]-1,4,4));
    
    FMatrix<double> tc(4,4);
    for (unsigned long i1=0; i1<n[0]-1; ++i1)
    {
        d1=xlist[0][i1+1]-xlist[0][i1];
        for (unsigned long i2=0; i2<n[1]-1; ++i2)
        {
            //std::cerr<<"test deriv: "<<dylist(i1,i2,0)<<","<<ny1list(i1,i2,0)<<"\n";
            d2=xlist[1][i2+1]-xlist[1][i2];
            y[0]=ylist(i1,i2); y[1]=ylist(i1+1,i2); y[2]=ylist(i1+1,i2+1); y[3]=ylist(i1,i2+1);
            dy1[0]=dylist(i1,i2,0); dy1[1]=dylist(i1+1,i2,0); dy1[2]=dylist(i1+1,i2+1,0); dy1[3]=dylist(i1,i2+1,0);
            dy2[0]=dylist(i1,i2,1); dy2[1]=dylist(i1+1,i2,1); dy2[2]=dylist(i1+1,i2+1,1); dy2[3]=dylist(i1,i2+1,1);
            d2y12[0]=dcross(i1,i2); d2y12[1]=dcross(i1+1,i2); d2y12[2]=dcross(i1+1,i2+1); d2y12[3]=dcross(i1,i2+1);
            /*std::cerr<<"x   "<<xlist[0][i1]<<","<<xlist[1][i2]<<"\n";
            std::cerr<<"y   "<<y[0]<<","<<y[1]<<","<<y[2]<<","<<y[3]<<"\n";
            std::cerr<<"dy1 "<<dy1[0]<<","<<dy1[1]<<","<<dy1[2]<<","<<dy1[3]<<"\n";
            std::cerr<<"dy2 "<<dy2[0]<<","<<dy2[1]<<","<<dy2[2]<<","<<dy2[3]<<"\n";
            std::cerr<<"d2y "<<d2y12[0]<<","<<d2y12[1]<<","<<d2y12[2]<<","<<d2y12[3]<<"\n";
            std::cerr<<"bicubic coeff"<<tc<<"\n";*/
            
            IBicCoeff(y, dy1, dy2, d2y12, d1, d2, tc);
            
            for (unsigned long i=0; i<4; ++i) for (unsigned long j=0; j<4; ++j) clist(i1,i2,i,j)=tc(i,j);
        }
    }
}

void InterpolateBicubic::dintrp(const fixarray<unsigned long, 2>& j, const fixarray<double, 2>& x, double& y, fixarray<double, 2>& dy)
{
    double x1l=xlist[0][j[0]], x1u=xlist[0][j[0]+1], x2l=xlist[1][j[1]], x2u=xlist[1][j[1]+1];
    /*std::cerr<<x[0]<<" "<<x[1]<<"\n";
    std::cerr<<j[0]<<" "<<j[1]<<"\n";
    std::cerr<<x1l<<":"<<x1u<<"   "<<x2l<<":"<<x2u<<"\n"; */
    if (x1u == x1l || x2u == x2l) ERROR("Wrong boundaries detected in interpolation");
    double d1=x1u-x1l,d2=x2u-x2l,t=(x[0]-x1l)/d1,u=(x[1]-x2l)/d2;
    
    double ansy,ansy1,ansy2; ansy=ansy2=ansy1=0.0;
    for (long i=3;i>=0;i--) {
        ansy=t*ansy+((clist(j[0],j[1],i,3)*u+clist(j[0],j[1],i,2))*u+clist(j[0],j[1],i,1))*u+clist(j[0],j[1],i,0);
        ansy2=t*ansy2+(3.0*clist(j[0],j[1],i,3)*u+2.0*clist(j[0],j[1],i,2))*u+clist(j[0],j[1],i,1);
        ansy1=u*ansy1+(3.0*clist(j[0],j[1],3,i)*t+2.0*clist(j[0],j[1],2,i))*t+clist(j[0],j[1],1,i);
        //if (i==0) break;
    }
    ansy1 /= d1;
    ansy2 /= d2;
    y=ansy;
    dy[0]=ansy1;
    dy[1]=ansy2;
}
}; //ends namespace toolbox
