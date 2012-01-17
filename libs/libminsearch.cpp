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

}; //ends namespace toolbox
