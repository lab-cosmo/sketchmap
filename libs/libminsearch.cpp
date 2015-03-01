/* Minimization library. Mostly implemented as templates.
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/
   
#include "tbdefs.hpp"
#include "minsearch.hpp"

namespace toolbox {

std::valarray<std::valarray<double> > make_simplex(std::valarray<double> ip, double dr, double fr)
{
    std::valarray<std::valarray<double> > rs(ip,ip.size()+1);
    for (unsigned long i=0; i<ip.size();++i)
    { rs[i+1][i]+=dr; rs[i+1][i]*=fr; }
    return rs;
}

std::valarray<std::valarray<double> > make_walkers(unsigned long ndim, unsigned long nwalker, double dr, MTRndUniform rngu)
{
    std::valarray<std::valarray<double> > rs(std::valarray<double>(ndim), nwalker);
    for (unsigned long i=0; i<nwalker;++i) for (unsigned long j=0; j<ndim;++j)
    { rs[i][j]=(rngu()-0.5)*dr; };
    return rs;
}

}; //ends namespace toolbox
