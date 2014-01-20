/* Toolbox library. A few simple functions used all around.
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/

#define __EXTERNALS 1
#include "tbdefs.hpp"
#include <vector>
//externals and very simple general-purpose functions

namespace toolbox {
#ifdef BENCHMARK
    CTBBD TBBenchmarks;
#endif
};

namespace std{
    toolbox::mpiostream pout(std::cout), perr(std::cerr);
};


namespace toolbox{
    double str2float(const std::string& str)
    {
        std::stringstream sstr(str);
        double rval; sstr >> rval;
        return rval;
    }

    int str2int(const std::string& str)
    {
        std::stringstream sstr(str);
        int rval; sstr >> rval;
        return rval;
    }
        
    std::string int2str(const long& ival)
    {
        std::stringstream sstr;
        sstr << ival;
        std::string rval; sstr >> rval;
        return rval;
    }
        
    std::string float2str(const double& ival)
    {
        std::stringstream sstr;
        sstr << ival;
        std::string rval; sstr >> rval;
        return rval;
    }

    void csv2floats(const std::string&  istr, std::valarray<double>& vv)
    {
        std::vector<double> vl(0);
        vl.clear(); std::string ls=istr;
        int pos=0;
        while( (pos = ls.find_first_of(',')) != ls.npos )
        {
            if(pos > 0)
            {
                vl.push_back(str2float(ls.substr(0,pos)));
            }
            ls=ls.substr(pos+1);
        }
        if(ls.length() > 0)
        {
            vl.push_back(str2float(ls));
        }
        vv.resize(vl.size()); //copies onto a valarray
        for (int k=0; k<vl.size(); k++) vv[k]=vl[k];
    }
};
