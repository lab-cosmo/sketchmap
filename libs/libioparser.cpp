#include "ioparser.hpp"
/***************************************************************************
 * This contains a bunch of member functions connected with io-lib classes *
 * which can (have to) be compiled in a single translation unit, not being *
 * either templates or inline functions                                    *
 ***************************************************************************/
namespace toolbox{

bool io_label_seek(std::istream& istr, const std::string& label)
{
    std::string rlab;
    while ((istr>>rlab) && rlab!=label);
    if (rlab!=label) return false; else return true;
}

/*  BEGINS IOMap MEMBER FUNCTIONS */
IOMap::~IOMap()
{ 
    std::map<std::string, IFBase*>::iterator it=begin();
    while (it!=end()) { delete it->second; ++it;}
}

IFBase& IOMap::operator[] (const std::string& key)
{
    std::map<std::string, IFBase*>::iterator it;
    
    if ((it=find(key))==end()) ERROR("Trying to access non-existing element!");
    return *(it->second);
}

IOMap& IOMap::operator<<(std::istream& istr)
{
    bool fnested=false, fstarted=false;
    std::string dump, label;
    while (istr>>label)
    {
        //handle comments 
        if (label[0]=='#')  { getline(istr,dump); continue; }
        //get what is supposed to be a label! if the first label we get is a brace, this means this is a sub-input 
        //and that we must look for a terminator as well
        
        if (label=="") continue;
        if (label=="{") 
        {
            if (fnested || fstarted) { ERROR("Dangling open brace in input file."); }
            else { fnested=true; continue; } 
        }

        if (label=="}") 
        { if (fnested) break; else ERROR("Dangling closed brace in input file."); }
        
        if (find(label)!=end()) 
        {
            fstarted=true;
            if (((*this)[label].flags & iff_set) && (flags & if_warn_multiple))
                std::cerr<<"Warning. Option "<<label<<" is set multiple times.\n";
            istr >> (*this)[label]; 
            if ((flags & if_warn_nonvalid) && ((*this)[label].flags & iff_nonvalid))
                std::cerr<<"Warning. Invalid settings for option "<<label<<".\n";
        }
        else if (flags & if_warn_unused)
        std::cerr<<"Warning. I cannot understand what "<<label<<" means.\n";
    }
    
    //checks whether some flags have not been set and are not optional
    if (flags & if_warn_incomplete)
    {
        std::map<std::string, IFBase*>::iterator it=begin();
        while (it!=end())
        {
            if (!((it->second->flags & iff_set) || (it->second->flags & iff_optional)))
            {
                std::cerr<<"Warning. Option "<<it->first<<" has not a default and have not been defined.\n";
            }
            ++it;
        }
    }
    return *this;
} 

const IOMap& IOMap::operator>>(std::ostream& ostr) const
{
    std::map<std::string, IFBase*>::const_iterator it=begin();
    for (; it!=end(); ++it) ostr << (* (it->second));
    return *this;
}
/*  ENDSS IOMap MEMBER FUNCTIONS */

/*  BEGINS IField<IOMap> MEMBER FUNCTIONS */
template <> bool IField<IOMap>::operator>> (std::ostream& ostr) const
{
    ostr<<name<<" {\n";
    value >> ostr;
    ostr<<" }\n";
    return true;
} 

template <> bool IField<IOMap>::operator<< (std::istream& istr)
{
    value << istr; 
    flags|=iff_set; flags&=~iff_nonvalid; 
    std::map<std::string, IFBase*>::const_iterator it=value.begin();
    //if any of the members is nonset or nonvalid sets the flags of the IOMap field accordingly
    
    for (; it!=value.end(); ++it) 
    {  if (it->second->flags & iff_nonvalid==iff_nonvalid) flags&=~iff_nonvalid; }
    
    return true;
}
/*  ENDS IField<IOMap> MEMBER FUNCTIONS */
}; //ends namespace toolbox
