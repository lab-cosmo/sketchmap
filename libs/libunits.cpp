#include "conv-units.hpp"

namespace toolbox {
double UnitConv::u2i(const UCUnit& unit, const double& val, bool internal)
{
    double cf;
    switch (unit)
    {
        case UEnergy:
            switch (internal?IEnergy:Energy)
            {
                case Joule:        cf= 2.2936831964769027E17      ; break;
                case Hartree:      cf= 1.                         ; break;
                case Electronvolt: cf= 0.036748872195972296       ; break;
                case Kelvin:       cf= 3.16677469769886E-06       ; break;
            } break;
        case UTime:
            switch (internal?ITime:Time)
            {
                case Second:       cf= 4.134137379e+16        ; break;
                case AtomicTime:   cf= 1.                     ; break;
                case Picosecond:   cf= 4.134137379e+04        ; break;
            } break;
        case UFrequency:
            switch (internal?IFrequency:Frequency)
            {
                case Hertz:       cf= 2.4188843e-17*2.*toolbox::constant::pi ; break;
                case AtomicFrequency:   cf= 1.                                ; break;
                case Terahertz:   cf= 2.4188843e-05*2.*toolbox::constant::pi ; break;
                case Inversecm:   cf= 7.251632747e-7*2.*toolbox::constant::pi; break;
            } break;
        case UMass:
            switch (internal?IMass:Mass)
            {
                case Kilogram:     cf= 9.10938215E-31;            ; break;
                case AtomicMass:   cf= 1.                         ; break;
                case AMU:          cf= 1822.8900409014022         ; break;
            } break;
        case ULength:
            switch (internal?ILength:Length)
            {
                case Meter:        cf= 1.889726132885643E10       ; break;
                case AtomicLength: cf= 1.                         ; break;
                case Angstrom:     cf= 1.889726132885643E00       ; break;
            } break;
    }
    if (!internal) cf=cf/u2i(unit,1.,true);
    return val*cf;
}

double UnitConv::i2u(const UCUnit& unit, const double& val)
{
    return val/u2i(unit,1.);
}

std::ostream& operator<< (std::ostream& ostr, const UnitConv::UCEnergy& oo)
{
    if (oo==UnitConv::Joule) ostr<<" joule \n";
    else if (oo==UnitConv::Hartree) ostr<<" hartree \n";
    else if (oo==UnitConv::Electronvolt) ostr<<" electronvolt \n";
    else if (oo==UnitConv::Kelvin) ostr<<" kelvin\n";
    return ostr;
}
std::istream& operator>> (std::istream& istr, UnitConv::UCEnergy& ii)
{
    std::string is;
    istr>>is;
    if (is=="joule") ii=UnitConv::Joule;
    else if (is=="hartree") ii=UnitConv::Hartree;
    else if (is=="electronvolt") ii=UnitConv::Electronvolt;
    else if (is=="kelvin") ii=UnitConv::Kelvin;
    return istr;
}
std::ostream& operator<< (std::ostream& ostr, const UnitConv::UCMass& oo)
{
    if (oo==UnitConv::AtomicMass) ostr<<" atomic \n";
    else if (oo==UnitConv::AMU) ostr<<" amu \n";
    else if (oo==UnitConv::Kilogram) ostr<<" kilogram \n";
    return ostr;
}
std::istream& operator>> (std::istream& istr, UnitConv::UCMass& ii)
{
    std::string is;
    istr>>is;
    if (is=="atomic") ii=UnitConv::AtomicMass;
    else if (is=="amu") ii=UnitConv::AMU;
    else if (is=="kilogram") ii=UnitConv::Kilogram;
    return istr;
}
std::ostream& operator<< (std::ostream& ostr, const UnitConv::UCTime& oo)
{
    if (oo==UnitConv::AtomicTime) ostr<<" atomic \n";
    else if (oo==UnitConv::Second) ostr<<" second \n";
    else if (oo==UnitConv::Picosecond) ostr<<" picosecond \n";
    return ostr;
}
std::istream& operator>> (std::istream& istr, UnitConv::UCTime& ii)
{
    std::string is;
    istr>>is;
    if (is=="atomic") ii=UnitConv::AtomicTime;
    else if (is=="second") ii=UnitConv::Second;
    else if (is=="picosecond") ii=UnitConv::Picosecond;
    return istr;
}
std::ostream& operator<< (std::ostream& ostr, const UnitConv::UCFrequency& oo)
{
    if (oo==UnitConv::AtomicFrequency) ostr<<" atomic \n";
    else if (oo==UnitConv::Hertz) ostr<<" hertz \n";
    else if (oo==UnitConv::Terahertz) ostr<<" terahertz \n";
    else if (oo==UnitConv::Inversecm) ostr<<" inversecm \n";
    return ostr;
}
std::istream& operator>> (std::istream& istr, UnitConv::UCFrequency& ii)
{
    std::string is;
    istr>>is;
    if (is=="atomic") ii=UnitConv::AtomicFrequency;
    else if (is=="hertz") ii=UnitConv::Hertz;
    else if (is=="terahertz") ii=UnitConv::Terahertz;
    else if (is=="inversecm") ii=UnitConv::Inversecm;
    return istr;
}
std::ostream& operator<< (std::ostream& ostr, const UnitConv::UCLength& oo)
{
    if (oo==UnitConv::AtomicLength) ostr<<" atomic \n";
    else if (oo==UnitConv::Meter) ostr<<" meter \n";
    else if (oo==UnitConv::Angstrom) ostr<<" angstrom \n";
    return ostr;
}
std::istream& operator>> (std::istream& istr, UnitConv::UCLength& ii)
{
    std::string is;
    istr>>is;
    if (is=="atomic") ii=UnitConv::AtomicLength;
    else if (is=="meter") ii=UnitConv::Meter;
    else if (is=="angstrom") ii=UnitConv::Angstrom;
    return istr;
}
std::ostream& operator<< (std::ostream& ostr, const UnitConv& oo)
{
    toolbox::IOMap iom;
    iom.insert(oo.Energy,"energy");
    iom.insert(oo.Time,"time");
    iom.insert(oo.Mass,"mass");
    iom.insert(oo.Length,"length");
    iom.insert(oo.Frequency,"frequency");
    ostr<<iom;
    return ostr;
}
std::istream& operator>> (std::istream& istr, UnitConv& ii)
{
    toolbox::IOMap iom;
    iom.insert(ii.Energy,"energy",UnitConv::Hartree);
    iom.insert(ii.Time,"time",UnitConv::AtomicTime);
    iom.insert(ii.Mass,"mass",UnitConv::AtomicMass);
    iom.insert(ii.Length,"length",UnitConv::AtomicLength);
    iom.insert(ii.Frequency,"frequency",UnitConv::AtomicFrequency);
    istr>>iom;
    return istr;
}
};