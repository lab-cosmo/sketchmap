/* A simple yet sophisticated input parser
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License
*/

#ifndef __IOPARSER_H
#define __IOPARSER_H

#include "tbdefs.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <valarray>
#include <vector>

/**********************************************************************
 * Definitions of helper classes for structured input/output.         *
 * The core class is IOMap, which contains a map connecting key fields*
 * as strings to IField named objects.                                *
 * IField class is a template, inheriting from a IFBase class in      *
 * order to have polymorphism.                                        *
 * We try to keep polymorphism intricacies as hidden as possible from *
 * the interface, but inside here it is a bloody mess.                *
 **********************************************************************/
namespace toolbox{
//bit maps for flags valid for field types
typedef enum _IFFlags {
    iff_optional=1,  //the field is optional, i.e. a meaningful default is provided
    iff_set=2,       //this bit is set when the field have been read from input
    iff_nonvalid=4   //field have been read but the value was not in valid form
} IFFlags;
//bit maps for IOMap class input check
typedef enum _IOMapFlags {
    if_warn_unused=1,      //warn about fields not recognized in the input
    if_warn_nonvalid=2,    //warn about invalid reads
    if_warn_incomplete=4,  //warn about non-initialized mandatory fields
    if_warn_multiple=8,    //warn if multiple definitions of the same field are encountered
    if_warn_all=~0         //self-explaining :-)
} IOMapFlags;

/*         PRE-DECLARATIONS        */
class IFBase;
template<class T> class IField;
class IOMap;

/*     BEGINS IOMAP CLASS DEFS     */
class IOMap: public std::map<std::string, IFBase *> {
friend class IFBase;   //input fields are friend
template <typename T> friend class IField;
public:
    unsigned long flags;
    //constructors
    IOMap(): std::map<std::string, IFBase *>(), flags(0) {}
    ~IOMap();

    //overloads map operator [] so as to avoid automatic fields creation
    template <class T> void insert(T& var, const std::string& name);
    template <class T> void insert(T& var, const std::string& name, const T& val);

    IFBase& operator[] (const std::string& key);

    IOMap& operator<<(std::istream& istr);
    const IOMap& operator>>(std::ostream& ostr) const;
};

inline std::ostream& operator<< (std::ostream& ostr, const IOMap& iom) { iom >> ostr; return ostr; }
inline std::istream& operator>> (std::istream& istr, IOMap& iom)  { iom << istr; return istr; }

//a couple of template functions for IOMap
template <class T> void IOMap::insert(T& var, const std::string& name)
{
    if (name==std::string("")) ERROR("Cannot insert elements with no name!");
    if (find(name)!=end()) ERROR("Cannot insert two elements with the same name!");
    this->std::map<std::string, IFBase *>::operator[](name)=new IField<T>(var,name);
}

template <class T> void IOMap::insert(T& var, const std::string& name, const T& val)
{
    if (name==std::string("")) ERROR("Cannot insert elements with no name!");
    if (find(name)!=end()) ERROR("Cannot insert two elements with the same name!");
    this->std::map<std::string, IFBase *>::operator[](name)=new IField<T>(var,name,val);
}

/*       ENDS IOMAP CLASS DEFS     */

/*   BEGINS BASE CLASS FOR INPUT FIELD   */
std::ostream& operator<< (std::ostream& ostr, const IFBase& ifield);
std::istream& operator>> (std::istream& istr, IFBase& ifield);
class IFBase {
friend class IOMap;  //IOMap is friend
protected:
    std::string name;    //field is identified by this string
public:
    unsigned long flags;
    //int mikt;
    //return true or false depending from success
    virtual std::ostream& operator>> (std::ostream& ostr) const { return ostr; };
    virtual std::istream& operator<< (std::istream& istr) {return istr;};

    //using stringstream we get easily the ability to read from strings as well as from streams
    virtual bool operator>> (std::string& str) const { std::ostringstream os; std::ostream& rval=(os << (*this)); str=os.str(); return os.good(); }
//    virtual bool operator>> (std::string& str) const { std::ostringstream os; bool rval; os << (*this); rval=bool(os); str=os.str(); return rval; }
    virtual std::istream& operator<< (const std::string& str) { std::istringstream is(str); return ((*this) << is); }

    IFBase(std::string nname) : name(nname), flags(0) {}
    virtual ~IFBase() {}
};

//these are to allow for standard syntax with the streams
inline std::ostream& operator<< (std::ostream& ostr, const IFBase& ifield) { ifield >> ostr; return ostr; }
inline std::istream& operator>> (std::istream& istr, IFBase& ifield)  {ifield << istr; return istr; }

/*     ENDS BASE CLASS FOR INPUT FIELD   */

/*   BEGINS TEMPLATE CLASS FOR INPUT FIELD   */
/**********************************************************************
 * This is a standard input field object. Any object can work with it *
 * provided operators to read & write from streams are defined.       *
 * The syntax is of the form label { ..... }\n. Since we must parse   *
 * just until the terminator "}\n", this means that the input must be   *
 * read and stored in a stringstream for further processing, which    *
 * might pose mermory concers with nested types.  A further limitation*
 * is that no free "}" characters can be present in the nested input. *
 * These are the drawbacks of the generality of the template, and can *
 * be overcome by defining specialized versions of the read-write     *
 * member functions.                                                  *
 **********************************************************************/
template <class T>
class IField : public IFBase {
public:
    //the variable we are reading is stored as a reference, so that the
    //value is updated as it is read. assignment operator is used.
    T& value;
public:
    using IFBase::operator>>;
    using IFBase::operator<<;
    using IFBase::flags;

    IField(T& var, const std::string& nname) : IFBase(nname), value(var) {}
    IField(T& var, const std::string& nname, const T& defval) : IFBase(nname), value(var)
    { value=defval; flags |=iff_optional; }

    std::ostream& operator>>(std::ostream& ostr) const
    {
        ostr<< name << " { ";
        ostr<<value;
        ostr<<" }\n ";
        return ostr;
    }

    std::istream& operator<<(std::istream& istr)
    {
        long nbraces=0; bool fopened=false, fbreak=false;
        std::size_t dpos, ppos;
        std::string s;
        std::stringstream ss;
        while (!fbreak && getline(istr,s))
        {
            //handle comments
            dpos=s.find('#'); ppos=0;
            //parse line for braces
            while ((ppos=s.find_first_of("{}", ppos))!=std::string::npos && (dpos==0 || ppos<dpos))
            {
                if (s[ppos]=='}') --nbraces;
                else {
                    ++nbraces; fopened=true;
                    if (nbraces==1)
                    {
                        //first open brace is the beginning of the enclosing string, so
                        //we must trim the string and reset the search
                        ss.clear(); s=s.substr(ppos+1); ppos=0;
                    }
                }
                if (nbraces==0 && fopened)
                {
                    //the field is over. we trim the string, push it into ss and break out.
                    //since in general istreams don''t allow seeking, the remaining text after } is
                    //trashed, that's why "}\n" is the suggested terminator rather than "}" alone
                    if (ppos>0) s=s.substr(0,ppos);  else s="";
                    fbreak=true; break;
                }
                else { if (nbraces<0) { ERROR("Mismatched braces in input file."); } }
                ++ppos;
            }
            ss<<s<<"\n";
        }

        //we use default reader for type
        ss >> value;
        //we have no mean of checking wether read was successful since we only require type T to define
        //iostream operators << >>. we must be optimistic and assume that everything went all-right
        flags|=iff_set;
        if (ss.good()) { flags &=~iff_nonvalid; } else { flags|=iff_nonvalid; }
        return istr;
    }

    operator T() { return value; }
};

/***********************************************************************************************
 * Every template declaration must have a const counterpart, as obviously one cannot define a  *
 * istream << operator for a field of a const type. This is more than academic, if we want to  *
 * use IOMap for output as well, as we might need to build an IOMap out of a const reference   *
 ***********************************************************************************************/
template <class T>
class IField<const T>: public IFBase {
public:  //!,mikt
    const T& value;
public:
    using IFBase::operator>>;
    using IFBase::operator<<;
    using IFBase::flags;

    IField(const T& var, const std::string& nname) : IFBase(nname), value(var) { flags|=(iff_set|iff_optional);}
    std::ostream& operator>>(std::ostream& ostr) const
    {
        ostr<< name << " { ";
        ostr<< value;
        ostr<<" }\n ";
        return ostr;
    }
    operator T() { return value; }
};

/*   ENDS TEMPLATE CLASS FOR INPUT FIELD   */

/*   BEGINS SPECIALIZATION FOR INTRINSIC TYPES   */
/******************************************************************
  Specializations of template IField<> for intrinsic types,
  which we write as a single-line parameter.
*******************************************************************/
//function for reading intrinsic types (i.e. parameters of the form "label value")
// for which operators << >> for streams are defined. just to make one-line definitions below
template <class U>
inline std::istream& read_udt(std::istream& istr, U& rval, unsigned long& rflags)
{
    istr.clear(); istr >> rval;
    rflags|=iff_set;
    if (istr.good()) { rflags &=(~iff_nonvalid); return istr;} else { rflags|=iff_nonvalid; istr.setstate(std::ios::failbit); return istr; }
}
template<class U>
inline std::ostream& write_udt(std::ostream& ostr, const U& rval, const std::string rname)
{ ostr << rname << "\t"  << rval <<"\n";  return ostr;}
#define __MK_IT_IOFIELD(type) \
  template<> inline std::istream& toolbox::IField<type>::operator<<(std::istream& istr)     \
  { return toolbox::read_udt(istr, value, flags); } \
  template<> inline std::ostream& toolbox::IField<type>::operator>>(std::ostream& ostr) const  \
  { return toolbox::write_udt(ostr, value, name); } \
  template<> inline std::ostream& toolbox::IField<const type>::operator>>(std::ostream& ostr) const \
  { return toolbox::write_udt(ostr, value, name); }

__MK_IT_IOFIELD(std::complex<double>)
__MK_IT_IOFIELD(double)
__MK_IT_IOFIELD(long)
__MK_IT_IOFIELD(float)
__MK_IT_IOFIELD(char)
__MK_IT_IOFIELD(unsigned char)
__MK_IT_IOFIELD(unsigned long)
__MK_IT_IOFIELD(std::string)

template<> inline std::istream& IField<bool>::operator<<(std::istream& istr)
{
    std::string val;
    istr>>val; flags|=iff_set; flags &=(~iff_nonvalid);
    if (val=="true") value=true; else if (val=="false") value=false;
    else { flags|=iff_nonvalid; istr.setstate(std::ios::failbit); return istr; }
    return istr;
}

template<> inline std::ostream& IField<bool>::operator>>(std::ostream& ostr) const       { ostr<<name<<" "<<(value?"true":"false")<<"\n"; return ostr;}
template<> inline std::ostream& IField<const bool>::operator>>(std::ostream& ostr) const { ostr<<name<<" "<<(value?"true":"false")<<"\n"; return ostr;}

/*   ENDS SPECIALIZATION FOR INTRINSIC TYPES   */

/*   BEGINS SPECIALIZATION FOR VALARRAYS   */
/*****************************************************************************
  Specialization for arbitrary valarrays!!! This should work for valarrays
  of any intrinsic type which supports i/o via streams. we must specialize
  the whole class because of some weird thing about c++ standard I didn't
  really understand
*****************************************************************************/
template<typename T>
class IField<std::valarray<T> >: public IFBase
{
private:
    std::valarray<T>& value;
public:
    using IFBase::operator>>;
    using IFBase::operator<<;

    IField(std::valarray<T>& var, const std::string& nname) : IFBase(nname), value(var) {}
    IField(std::valarray<T>& var, const std::string& nname, const std::valarray<T>& defval) : IFBase(nname), value(var) { value.resize(defval.size()); value=defval; flags |=iff_optional; }
    std::ostream& operator>> (std::ostream& ostr) const;
    std::istream& operator<< (std::istream& istr);

    operator std::valarray<T>() { return value; }
};

//valarray specialization member functions
template<typename T>
std::istream& IField<std::valarray<T> >::operator<<(std::istream& istr)
{
    //reads number of elements
    unsigned long nel;
    flags|=iff_set;
    istr.clear(); istr >> nel;
    if(istr.good()) { flags &=(~iff_nonvalid); }
    else { flags|=iff_nonvalid; istr.setstate(std::ios::failbit); }

    value.resize(nel);
    for (unsigned long i=0; i<nel; ++i)
    {
        istr>>value[i];
        if(! istr.good()) { flags|=iff_nonvalid; }
    }
    return istr;
}

template<typename T>
std::ostream& IField<std::valarray<T> >::operator>>(std::ostream& ostr) const
{
    ostr << name << " \t"<< value.size()<<"\n";
    for (unsigned long i=0; i<value.size(); ++i) ostr<<value[i]<<" ";
    ostr<<"\n";
    return ostr;
}

template<typename T>
class IField<const std::valarray<T> >: public IFBase
{
private:
    const std::valarray<T>& value;
public:
    using IFBase::operator>>;
    using IFBase::operator<<;

    IField(const std::valarray<T>& var, std::string nname) : IFBase(nname), value(var) {}
    std::ostream& operator>> (std::ostream& ostr) const;
    operator std::valarray<T>() { return value; }
};

template<typename T>
std::ostream& IField<const std::valarray<T> >::operator>>(std::ostream& ostr) const
{
    ostr << name << " \t"<< value.size()<<"\n";
    for (unsigned long i=0; i<value.size(); ++i) ostr<<value[i]<<" ";
    ostr<<"\n";
    return ostr;
}

template <class T>
std::ostream& operator << (std::ostream& os, const std::valarray<T>& t)
{
    IField<std::valarray<T> > ift(const_cast<std::valarray<T>&>(t),"");
    ift>>os;
    return os;
}

template <class T>
std::istream& operator >> (std::istream& is, std::valarray<T>& t)
{
    IField<std::valarray<T> > ift(t,"",t);
    is>>ift;
    return is;
}
/*   ENDS SPECIALIZATION FOR VALARRAYS   */

/*   BEGINS SPECIALIZATION FOR VECTORS   */
/*****************************************************************************
  Specialization for arbitrary VECTORS!!! just the same syntax as for valarrays
*****************************************************************************/
template<typename T>
class IField<std::vector<T> >: public IFBase
{
private:
    std::vector<T>& value;
public:
    using IFBase::operator>>;
    using IFBase::operator<<;

    IField(std::vector<T>& var, const std::string& nname) : IFBase(nname), value(var) {}
    IField(std::vector<T>& var, const std::string& nname, const std::vector<T>& defval) : IFBase(nname), value(var) { value.resize(var.size()); value=defval; flags |=iff_optional; }
    std::ostream& operator>> (std::ostream& ostr) const;
    std::istream& operator<< (std::istream& istr);

    operator std::vector<T>() { return value; }
};

/*
  vector specialization member functions
  vectors are written out in the same form as valarrays, but for
  being more flexible, they can also be inputted as
  label { \n
   element_1 \n
   element_2 \n
   :
   }
*/

template<typename T>
std::istream& IField<std::vector<T> >::operator<<(std::istream& istr)
{
    //reads number of elements or braces starting the vector
    unsigned long nel;

    std::stringstream ls; std::string word; unsigned char ch; T el;
    flags|=iff_set;
    istr >> word;

    if (word=="{")
    {
        //braces mode
        value.resize(0);
        while (istr.good())
        {
            ch=istr.get();
            switch (ch)
            {
            case '}':
                return istr; break; //finished reading!
            case ' ':
            case '\t':
            case '\n':
                continue; //go on
                break;
            case '#': //comment
                getline(istr,word);
                continue;
                break;
            default:
                istr.putback(ch);
                break;
            }
            istr >> el;
            value.push_back(el);
        }
    }
    else
    {
        //valarray mode
        ls.str(word); ls >> nel;

        //!we should be checkin' we are reading a number!LAZYWARN
        //if(word==dummy)
        { flags &=(~iff_nonvalid); }
        //else { flags|=iff_nonvalid; return false; }

        value.resize(nel);
        for (unsigned long i=0; i<nel; ++i)
        {
            istr>>value[i];
            if(! istr.good()) { flags|=iff_nonvalid; }
        }

    }
    return istr;
}

template<typename T>
std::ostream& IField<std::vector<T> >::operator>>(std::ostream& ostr) const
{
    ostr << name << " \t"<< value.size()<<"\n";
    for (unsigned long i=0; i<value.size(); ++i) ostr<<value[i]<<" ";
    ostr<<"\n";
    return ostr;
}

template<typename T>
class IField<const std::vector<T> >: public IFBase
{
private:
    const std::vector<T>& value;
public:
    using IFBase::operator>>;
    using IFBase::operator<<;

    IField(const std::vector<T>& var, const std::string& nname) : IFBase(nname), value(var) {}
    std::ostream& operator>> (std::ostream& ostr) const;
    operator std::vector<T>() { return value; }
};

template<typename T>
std::ostream& IField<const std::vector<T> >::operator>>(std::ostream& ostr) const
{
    ostr << name << " \t"<< value.size()<<"\n";
    for (unsigned long i=0; i<value.size(); ++i) ostr<<value[i]<<" ";
    ostr<<"\n";
    return ostr;
}

template <class T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& t)
{
    IField<std::vector<T> > ift(const_cast<std::vector<T>&>(t),"");
    ift>>os;
    return os;
}

template <class T>
std::istream& operator >> (std::istream& is, std::vector<T>& t)
{
    IField<std::vector<T> > ift(t,"",t);
    is>>ift;
    return is;
}
/*   ENDS SPECIALIZATION FOR VECTORS   */

/*   BEGINS SPECIALIZATION FOR MAPS   */
/*****************************************************************************
  Specialization for arbitrary MAPS!!!
*****************************************************************************/
template<class K, class T>
class IField< std::map< K, T > >: public IFBase
{
private:
    std::map<K,T>& value;
public:
    using IFBase::operator>>;
    using IFBase::operator<<;

    IField(std::map<K,T>& var, const std::string& nname) : IFBase(nname), value(var) {}
    IField(std::map<K,T>& var, const std::string& nname, const std::map<K,T>& defval) :
      IFBase(nname), value(var) { value.resize(var.size()); value=defval; flags |=iff_optional; }
    std::ostream& operator>> (std::ostream& ostr) const;
    std::istream& operator<< (std::istream& istr);

    operator std::map<K,T>() { return value; }
};

/*
  map specialization member functions
  the input format is
  label { \n
   key1 value1 \n
   key2 value2 \n
   :
   }
*/

template<typename K,typename T>
std::istream& IField<std::map<K,T> >::operator<<(std::istream& istr)
{
    //reads number of elements or braces starting the vector
    unsigned long nel;
    std::stringstream ls;
    std::string word;
    unsigned char ch;
    T el; K key;
    flags|=iff_set;
    istr >> word;

    if (word=="{")
    {
        //braces mode
        value.clear();
        while (istr.good())
        {
            ch=istr.get();
            switch (ch)
            {
            case '}':
                return istr; break; //finished reading!
            case ' ':
            case '\t':
            case '\n':
                continue; //go on
                break;
            case '#': //comment
                getline(istr,word);
                continue;
                break;
            default:
                istr.putback(ch);
                break;
            }
            istr >> key >> el;
            value[key]=el;
        }
    }
    else
    {
        ERROR("Wrong syntax for input of type map.");
        istr.setstate(std::ios::failbit);
    }
    return istr;
}

template<class K, class T>
std::ostream& IField<std::map<K,T> >::operator>>(std::ostream& ostr) const
{
    typename std::map<K,T>::iterator it;
    ostr << name << " \t {\n";
    for ( it=value.begin() ; it != value.end(); it++ )
      ostr<< (*it).first << " \t " << (*it).second << std::endl;
    ostr<<" }\n";
    return ostr;
}

template<typename K,typename T>
class IField<const std::map<K,T> >: public IFBase
{
private:
    const std::map<K,T>& value;
public:
    using IFBase::operator>>;
    using IFBase::operator<<;

    IField(const std::map<K,T>& var, const std::string& nname) : IFBase(nname), value(var) {}
    std::ostream& operator>> (std::ostream& ostr) const;
    operator std::map<K,T>() { return value; }
};

template<typename K,typename T>
std::ostream& IField<const std::map<K,T> >::operator>>(std::ostream& ostr) const
{
    typename std::map<K,T>::const_iterator it;
    ostr << name << " \t {\n";
    for ( it=value.begin() ; it != value.end(); it++ )
        ostr<< (*it).first << " \t " << (*it).second << std::endl;
    ostr<<" }\n";
    return ostr;
}

template<typename K,typename T>
std::ostream& operator << (std::ostream& os, const std::map<K,T>& t)
{
    IField<std::map<K,T> > ift(const_cast<std::map<K,T>&>(t),std::string(""));
    ift>>os;
    return os;
}

template<typename K,typename T>
std::istream& operator >> (std::istream& is, std::map<K,T>& t)
{
    IField<std::map<K,T> > ift(t,std::string(""),t);
    is>>ift;
    return is;
}
/*   ENDS SPECIALIZATION FOR MAPS   */


/*****************************************************************************
  OK. This is tricky. This defines a IOMap field, so it is the key to
  nested input.
*****************************************************************************/
template<> std::ostream& IField<IOMap>::operator>> (std::ostream& ostr) const;
template<> std::istream& IField<IOMap>::operator<< (std::istream& istr);

}; //ends namespace toolbox
#endif //ends #define __IOPARSER_H
