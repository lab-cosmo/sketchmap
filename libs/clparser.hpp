#ifndef __CLPARSER_H
#define __CLPARSER_H
#include "tbdefs.hpp"
#include <vector>
#include <string>
#include <sstream>
namespace toolbox {

class CLParser {
private: 
    std::vector<std::string> ppars;
public: 
    CLParser(int const& argc, char ** const & argv);
    template <typename T> bool getoption(std::vector<T>& value, const std::string& parname);
    template <typename T> bool getoption(std::vector<T>& value, const std::string& parname, const std::vector<T>& defval);
    template <typename T> bool getoption(T& value, const std::string& parname);
    template <typename T> bool getoption(T& value, const std::string& parname, const T& defval);
};

template <typename T>
bool CLParser::getoption(T& value, const std::string& parname)
{
    std::string spar=std::string("-")+parname;
    unsigned long i=0;
    while (i<ppars.size() && ppars[i]!=spar) ++i;
    if (i==ppars.size()) return false; //key not found!

    if (i+1==ppars.size() || ppars[i+1].size()==0 || ppars[i+1][0]=='-') 
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

template <typename T>
bool CLParser::getoption(std::vector<T>& value, const std::string& parname)
{
    std::string spar=std::string("-")+parname;
    unsigned long i=0; T iel;
    while (i<ppars.size() && ppars[i]!=spar) ++i;
    if (i==ppars.size()) return false; //key not found!

    if (i+1==ppars.size() || ppars[i+1].size()==0 || ppars[i+1][0]=='-') 
    {
#ifdef DEBUG
        ERROR("Key value missing on command line\n");
#endif
        return false;
    }
    
    ++i; unsigned long nel=0;
    while (i<ppars.size() && ppars[i][0]!='-') 
    {
        nel++; 
        std::stringstream ls(ppars[i]); 
        ls >> iel;
        value.push_back(iel);

        if (!(ls.good() || ls.eof())) return false;
        ++i;
    }
    
    if (nel>0) return true; else return false;
}


template <typename T>
bool CLParser::getoption(T& value, const std::string& parname, const T& defval)
{
    if (!getoption(value, parname)) 
        value=defval;
    return true;
}

template <typename T>
bool CLParser::getoption(std::vector<T>& value, const std::string& parname, const std::vector<T>& defval)
{
    if (!getoption(value, parname)) 
        value=defval;
    return true;
}

template<>
bool CLParser::getoption(bool& value, const std::string& parname);

template <>
bool CLParser::getoption(std::string& value, const std::string& parname);

}; //ends namespace toolbox
#endif //ends ifndef __CLPARSER_H
