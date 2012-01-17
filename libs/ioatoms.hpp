#ifndef __IOATOMS_H
#define __IOATOMS_H

#include "tbdefs.hpp"
#include <vector>
#include <map>

namespace toolbox {
struct AtomData {
    std::string name;
    double x, y, z;
    std::vector<double> props;
    std::map<std::string,double> nprops;
    AtomData() : name(""), x(0.), y(0.), z(0.), props(0),  nprops() {}
};

struct AtomFrame {
    std::string comment;
    unsigned long index;
    std::vector<AtomData> ats;
    std::vector<double> props;
    std::map<std::string,double> nprops;
};

void ReadXYZ(std::istream& istr, std::vector<AtomFrame>& frames);
bool ReadXYZFrame(std::istream& istr, AtomFrame& frames);
bool ReadDLPFrame(std::istream& istr, AtomFrame& frames);
bool ReadDLPConf(std::istream& istr, AtomFrame& frame);
bool WritePDBFrame(std::ostream& ostr, AtomFrame& frame);
}; //ends namespace toolbox
#endif //ends ifdef __IOATOMS_H
