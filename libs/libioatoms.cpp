#include "ioatoms.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace toolbox {

bool ReadXYZFrame(std::istream& istr, AtomFrame& curfr)
{
    std::string line;
    std::stringstream ss;
    AtomData curat;
    unsigned long nat, i; double cprop; 
    if(!getline(istr, line)) return false;
    curfr.index=0; curfr.ats.resize(0);
    ss.clear(); ss<<line;
    ss>>nat;
    if (!getline(istr,curfr.comment)) ERROR("Read failed while getting frame "<<curfr.index<<".");
    for (i=0; i<nat && getline(istr, line); ++i)
    {
        ss.clear(); ss<<line;
        curat.props.resize(0);
        ss>>curat.name>>curat.x>>curat.y>>curat.z;
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        while ((ss>>cprop)) curat.props.push_back(cprop);
        curfr.ats.push_back(curat);
    }
    return true;
}

void ReadXYZ(std::istream& istr, std::vector<AtomFrame>& frames)
{
    frames.resize(0);
    AtomFrame curfr;
    unsigned long nfr=0;
    while(ReadXYZFrame(istr,curfr)) {curfr.index=(++nfr); frames.push_back(curfr);}
}
    
bool ReadDLPFrame(std::istream& istr, AtomFrame& curfr)
{
    std::string line, dummy;
    std::stringstream ss;
    AtomData curat;
    unsigned long nat, ifr, ipbc, ikey, i; double cprop; 
    if(!getline(istr, line)) return false;
    curfr.index=0; curfr.ats.resize(0);
    ss.clear(); ss<<line;
    ss>>dummy>>ifr>>nat>>ikey>>ipbc>>curfr.nprops["timestep"];
    if(ipbc>0)
    {
        istr>>curfr.nprops["axx"]>>curfr.nprops["axy"]>>curfr.nprops["axz"];
        istr>>curfr.nprops["ayx"]>>curfr.nprops["ayy"]>>curfr.nprops["ayz"];
        istr>>curfr.nprops["azx"]>>curfr.nprops["azy"]>>curfr.nprops["azz"];
        getline(istr,dummy);
    }
    for (i=0; i<nat && getline(istr, line); ++i)
    {
        ss.clear(); ss.str(line);
        curat.props.resize(0);
        ss>>curat.name>>dummy>>curat.nprops["mass"]>>curat.nprops["charge"];
        getline(istr, line); ss.clear(); ss.str(line);
        ss>>curat.x>>curat.y>>curat.z;
        
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        if (ikey>0)
        {
            getline(istr, line); ss.clear(); ss.str(line);
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop);
        }
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        if (ikey>1)
        {
            getline(istr, line); ss.clear(); ss.str(line);
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop);
        }
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        curfr.ats.push_back(curat);
    }
    return true;
}
    
bool ReadDLPConf(std::istream& istr, AtomFrame& curfr)
{
    std::string line, dummy;
    std::stringstream ss;
    AtomData curat;
    unsigned long nat, ifr, ipbc, ikey, i; double cprop; 
    if(!getline(istr, line)) return false;  //drops comment line
    if(!getline(istr, line)) return false;  //reads header
    curfr.index=0; curfr.ats.resize(0);
    ss.clear(); ss<<line;
    ss>>ikey>>ipbc>>nat>>dummy;
    if(ipbc>0)
    {
        istr>>curfr.nprops["axx"]>>curfr.nprops["axy"]>>curfr.nprops["axz"];
        istr>>curfr.nprops["ayx"]>>curfr.nprops["ayy"]>>curfr.nprops["ayz"];
        istr>>curfr.nprops["azx"]>>curfr.nprops["azy"]>>curfr.nprops["azz"];
        getline(istr,dummy);
    }
    nat=0;
    for (i=0; getline(istr, line); ++i)
    {
        nat++;
        ss.clear(); ss.str(line);
        curat.props.resize(0);
        ss>>curat.name>>dummy;
        getline(istr, line); ss.clear(); ss.str(line);
        ss>>curat.x>>curat.y>>curat.z;
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        if (ikey>0)
        {
            getline(istr, line); ss.clear(); ss.str(line);
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop);
        }
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        if (ikey>1)
        {
            getline(istr, line); ss.clear(); ss.str(line);
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop);
        }
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        curfr.ats.push_back(curat);
    }
    return true;
}

bool WritePDBFrame(std::ostream& ostr, AtomFrame& frame)
{
    double lx, ly, lz, ca, cb, cc; 
    if (frame.nprops.count("axx")>0) 
    {
        lx=sqrt(frame.nprops["axx"]*frame.nprops["axx"]+frame.nprops["axy"]*frame.nprops["axy"]+frame.nprops["axz"]*frame.nprops["axz"]);
        ly=sqrt(frame.nprops["ayx"]*frame.nprops["ayx"]+frame.nprops["ayy"]*frame.nprops["ayy"]+frame.nprops["ayz"]*frame.nprops["ayz"]);
        lz=sqrt(frame.nprops["azz"]*frame.nprops["azx"]+frame.nprops["azy"]*frame.nprops["azy"]+frame.nprops["azz"]*frame.nprops["azz"]);
        cc=acos((frame.nprops["axx"]*frame.nprops["ayx"]+frame.nprops["axy"]*frame.nprops["ayy"]+frame.nprops["axz"]*frame.nprops["ayz"])/(lx*ly))*180/constant::pi;
        cb=acos((frame.nprops["axx"]*frame.nprops["azx"]+frame.nprops["axy"]*frame.nprops["azy"]+frame.nprops["axz"]*frame.nprops["azz"])/(lx*ly))*180/constant::pi;;
        ca=acos((frame.nprops["azx"]*frame.nprops["ayx"]+frame.nprops["azy"]*frame.nprops["ayy"]+frame.nprops["azz"]*frame.nprops["ayz"])/(lx*ly))*180/constant::pi;;
    }
    else
    {
        lx=ly=lz=0.0; 
        ca=cb=cc=90.0;
    }
    
    //outputs cell parameters (if present)
    ostr<<"CRYST1"
            <<std::setw(9)<<std::fixed<<std::setprecision(3)<<lx
            <<std::setw(9)<<std::fixed<<std::setprecision(3)<<ly
            <<std::setw(9)<<std::fixed<<std::setprecision(3)<<lz
            <<std::setw(7)<<std::fixed<<std::setprecision(2)<<ca
            <<std::setw(7)<<std::fixed<<std::setprecision(2)<<cb
            <<std::setw(7)<<std::fixed<<std::setprecision(2)<<cc
        <<" P 1       "
        <<"   1"
        <<std::endl;
    
    for (unsigned long i=0; i<frame.ats.size(); i++)
    {
        ostr<<"ATOM  "<<std::setw(5)<<i+1;
        ostr<<std::skipws<<std::setw(4)<<std::right<<frame.ats[i].name<<" ";
        ostr<<" "<<std::setw(3)<<std::right<<"X"<<" "<<" "<<"   1"<<" "<<"   ";
        ostr<<std::setw(8)<<std::fixed<<std::setprecision(3)<<frame.ats[i].x
            <<std::setw(8)<<std::fixed<<std::setprecision(3)<<frame.ats[i].y
            <<std::setw(8)<<std::fixed<<std::setprecision(3)<<frame.ats[i].z;
        ostr<<std::setw(6)<<std::fixed<<std::setprecision(2)
            <<(frame.ats[i].nprops.count("occupancy")>0?frame.ats[i].nprops["occupancy"]:0.0)
            <<std::setw(6)<<std::fixed<<std::setprecision(2)
            <<(frame.ats[i].nprops.count("beta")>0?frame.ats[i].nprops["beta"]:0.0);
        ostr<<"          "<<std::setw(2)<<frame.ats[i].name;
        ostr<<std::setw(2)<<(frame.ats[i].nprops.count("charge")>0?(int)frame.ats[i].nprops["charge"]:0);
        ostr<<std::endl;
    }
    ostr<<"END"<<std::endl;
}

}; //ends namespace toolbox
