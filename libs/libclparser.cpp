#include "clparser.hpp"

namespace toolbox{
CLParser::CLParser(int const& argc, char** const& argv)
{
    ppars.resize(argc);
    for (long i=0; i<argc;++i)
        ppars[i]=argv[i];
}

template<>
bool CLParser::getoption(bool& value, const std::string& parname)
{ 
    //for a bool parameter we just check the existence of the key on the commandline
    std::string spar=std::string("-")+parname;
    unsigned long i=0; value=false;
    while (i<ppars.size() && ppars[i]!=spar) ++i;
    if (i==ppars.size()) return false; //key not found!
    value=true; return true;
}

template <>
bool CLParser::getoption(std::string& value, const std::string& parname)
{
    //for a string parameter we don't use the stringstream!
    std::string spar=std::string("-")+parname;
    unsigned long i=0;
    while (i<ppars.size() && ppars[i]!=spar) ++i;
    if (i==ppars.size()) return false; //key not found!

    if (i+1==ppars.size() || ppars[i+1].size()==0) 
    {
#ifdef DEBUG
        ERROR("Key value missing on command line\n");
#endif
        return false;
    }
    
    value=ppars[i+1]; return true;
}

//for double type, we must accept negative numbers!
template <>
bool CLParser::getoption(double& value, const std::string& parname)
{
    std::string spar=std::string("-")+parname;
    unsigned long i=0;
    while (i<ppars.size() && ppars[i]!=spar) ++i;
    if (i==ppars.size()) return false; //key not found!

    if (i+1==ppars.size() || ppars[i+1].size()==0) 
    {
#ifdef DEBUG
        ERROR((std::string("Key value \"")+parname+std::string("\" missing on command line\n")));
#endif
        return false;
    }
    
    std::stringstream ls(ppars[i+1]); 
    ls >> value; 
    if (ls.good() || ls.eof()) return true; 
    else return false;
}

};
