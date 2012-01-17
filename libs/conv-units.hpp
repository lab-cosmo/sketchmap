#ifndef __CONV_UNITS_HPP
#define __CONV_UNITS_HPP
#include "tbdefs.hpp"
#include "ioparser.hpp"
namespace toolbox{
class UnitConv {
public:
    enum UCUnit { UEnergy, UMass, UTime, UFrequency, ULength };

    enum UCEnergy    { Joule, Hartree, Electronvolt,  Kelvin };
    enum UCMass      { Kilogram, AMU, AtomicMass};
    enum UCTime      { Second, Picosecond, AtomicTime };
    enum UCFrequency { Hertz, Terahertz, Inversecm, AtomicFrequency};
    enum UCLength    { Meter, AtomicLength, Angstrom };

    UCEnergy Energy, IEnergy;
    UCMass Mass, IMass;
    UCLength Length, ILength;
    UCTime Time, ITime;
    UCFrequency Frequency, IFrequency;
    UnitConv() : Energy(Hartree), IEnergy(Hartree),
                 Mass(AtomicMass), IMass(AtomicMass),
                 Length(AtomicLength), ILength(AtomicLength),
                 Time(AtomicTime), ITime(AtomicTime),
                 Frequency(AtomicFrequency), IFrequency(AtomicFrequency) {}

    double u2i(const UCUnit& unit, const double& val, bool internal=false);
    double i2u(const UCUnit& unit, const double& val);
};

std::ostream& operator<< (std::ostream& ostr, const UnitConv::UCEnergy& oo);
std::istream& operator>> (std::istream& istr, UnitConv::UCEnergy& ii);
std::ostream& operator<< (std::ostream& ostr, const UnitConv::UCMass& oo);
std::istream& operator>> (std::istream& istr, UnitConv::UCMass& ii);
std::ostream& operator<< (std::ostream& ostr, const UnitConv::UCTime& oo);
std::istream& operator>> (std::istream& istr, UnitConv::UCTime& ii);
std::ostream& operator<< (std::ostream& ostr, const UnitConv::UCFrequency& oo);
std::istream& operator>> (std::istream& istr, UnitConv::UCFrequency& ii);
std::ostream& operator<< (std::ostream& ostr, const UnitConv::UCLength& oo);
std::istream& operator>> (std::istream& istr, UnitConv::UCLength& ii);
std::ostream& operator<< (std::ostream& ostr, const UnitConv& oo);
std::istream& operator>> (std::istream& istr, UnitConv& ii);

__MK_IT_IOFIELD(UnitConv::UCTime);
__MK_IT_IOFIELD(UnitConv::UCMass);
__MK_IT_IOFIELD(UnitConv::UCLength);
__MK_IT_IOFIELD(UnitConv::UCEnergy);
__MK_IT_IOFIELD(UnitConv::UCFrequency);
};     //ends namespace toolbox
#endif //ends #define __CONV_UNITS_HPP