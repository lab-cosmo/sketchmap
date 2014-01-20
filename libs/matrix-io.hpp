/* I/O of matrices. Compatible with ioparser
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/

#ifndef __MATRIX_IO_H
#define __MATRIX_IO_H
/*****************************************************************
 defines helper function for matrix i/o, both on streams and 
 for ioparser library.  conditionally includes different sections
 depending on loaded matrix types, so it should be included
 last, after all the matrix-related packages.
*****************************************************************/
namespace toolbox {

//CRS MATRIX
#ifdef __MATRIX_CRS_H
    //iofield for options
    template <class U> 
    class IField<MatrixOptions<CrsMatrix<U> > > : public IFBase
    {
        private:
            MatrixOptions<CrsMatrix<U> >& value;
        public:
	    using IFBase::operator>>;
	    using IFBase::operator<<;

            IField(MatrixOptions<CrsMatrix<U> >& var, std::string nname) : IFBase(nname), value(var) {}
            IField(MatrixOptions<CrsMatrix<U> >& var, std::string nname, const CrsMatrix<U>& defval) : IFBase(nname), value(var) { value=defval; flags |=iff_optional; }

            bool operator>>(std::ostream& ostr) const
            {
                IOMap iom;  std::string atm;
                switch (value.atnorm)
                { 
                    case at_normone: atm="norm_one"; break;
                    case at_norminf: atm="norm_inf"; break;
                    case at_normfrob: atm="norm_frob"; break;
                    default: ERROR("Unsupported norm mode");
                };
                iom.insert(value.atthresh,"at_threshold");
                iom.insert(atm,"at_norm");
                iom.insert(value.atnodiag,"at_nodiag");
                iom.insert(value.atrel,"at_relative");
                ostr << name << " {\n"<<iom<<"\n}\n";
                return true;
            }
    
            bool operator<<(std::istream& istr)
            {
                IOMap iom;  std::string atm;
                iom.insert(value.atthresh,"at_threshold",0.);
                iom.insert(atm,"at_norm",std::string("norm_inf"));
                iom.insert(value.atnodiag,"at_nodiag");
                iom.insert(value.atrel,"at_relative",false);
                
                istr>>iom;
                flags|=iff_set; flags &=~iff_nonvalid;
                if (!istr.good()) { flags|=iff_nonvalid; return false; }
                
                if (atm=="norm_one") value.atnorm=at_normone;
                else if (atm=="norm_inf") value.atnorm=at_norminf;
                else if (atm=="norm_frob") value.atnorm=at_normfrob;
                else { flags |=iff_nonvalid; return false; }
                
                if (value.atthresh<0.) { flags|=iff_nonvalid; return false; }
                
                return true;
        //!!SANITY CHECK MUST BE DONE
            }

            operator MatrixOptions<CrsMatrix<U> >() { return value; }
    };

    template <class U> 
    class IField<const MatrixOptions<CrsMatrix<U> > > : public IFBase
    {
        private:
            const MatrixOptions<CrsMatrix<U> >& value;
        public:
    using IFBase::operator>>;
    using IFBase::operator<<;

            IField(const MatrixOptions<CrsMatrix<U> >& var, std::string nname) : IFBase(nname), value(var) {}

            bool operator>>(std::ostream& ostr) const
            {
                IOMap iom;  std::string atm;
                switch (value.atnorm)
                { 
                    case at_normone: atm="norm_one"; break;
                    case at_norminf: atm="norm_inf"; break;
                    case at_normfrob: atm="norm_frob"; break;
                    default: ERROR("Unsupported norm mode");
                };
                iom.insert(value.atthresh,"at_threshold");
                iom.insert(atm,"at_norm");
                iom.insert(value.atnodiag,"at_nodiag");
                iom.insert(value.atrel,"at_relative");
                ostr << name << " {\n"<<iom<<"\n}\n";
                return true;
            }
            operator MatrixOptions<CrsMatrix<U> >() { return value; }
    };
    
    /***********************************************************
    *             begin IField<CrsMatrix> definitions         *
    ***********************************************************/
    template<typename U> 
    class IField<CrsMatrix<U> > : public IFBase
    {
        private:
            CrsMatrix<U>& value;
        public:
    using IFBase::operator>>;
    using IFBase::operator<<;

            IField(CrsMatrix<U>& var, const std::string& nname) : IFBase(nname), value(var) {}
            IField(CrsMatrix<U>& var, const std::string& nname, const CrsMatrix<U>& defval) : IFBase(nname), value(var) { value=defval; flags |=iff_optional; }

            bool operator>>(std::ostream& ostr) const;
            bool operator<<(std::istream& istr);

            operator CrsMatrix<U>() { return value; }
    };

//!!WE SHOULD THINK ABOUT BINARY IO
    template<typename U> 
    bool IField<CrsMatrix<U> >::operator<<(std::istream& istr) 
    {
        IOMap iom;
        iom.insert(value.pops,"options");
        iom.insert(value.wr,"rows");
        iom.insert(value.wc,"cols");
        iom.insert(value.rpoints,"rowpointers");
        iom.insert(value.indices,"indices");
        iom.insert(value.values,"values");
        iom<<istr;
        flags|=iff_set; flags&=~iff_nonvalid;
    
    //check for input consistency ??MUST BE DONE!!!
        IOMap::iterator it=iom.begin();
        for (; it!=iom.end(); ++it)
        {
            if (it->second->flags & iff_nonvalid) { flags|=iff_nonvalid; return false; }
            if (!(it->second->flags & (iff_set | iff_optional)))  { flags&=~iff_set; return false; }
        }
    
        return true;
    }

    template<typename U> 
    bool IField<CrsMatrix<U> >::operator>>(std::ostream& ostr) const
    {
        IOMap iom;
        iom.insert(value.pops,"options");
        iom.insert(value.wr,"rows");
        iom.insert(value.wc,"cols");
        iom.insert(value.rpoints,"rowpointers");
        iom.insert(value.indices,"indices");
        iom.insert(value.values,"values");
        ostr<<name<<" {\n";
        iom>>ostr;
        ostr<<" }\n";
        return true;
    }

//const version
    template<typename U> 
            class IField<const CrsMatrix<U> > : public IFBase
    {
        private:
            const CrsMatrix<U>& value;
        public:
    using IFBase::operator>>;
    using IFBase::operator<<;

            IField(const CrsMatrix<U>& var, const std::string& nname) : IFBase(nname), value(var) {}
            bool operator>>(std::ostream& ostr) const;
            operator CrsMatrix<U>() { return value; }
    };

    template<typename U> 
    bool IField<const CrsMatrix<U> >::operator>>(std::ostream& ostr) const
    {
        IOMap iom;
        iom.insert(value.pops,"options");
        iom.insert(value.wr,"rows");
        iom.insert(value.wc,"cols");
        iom.insert(value.rpoints,"rowpointers");
        iom.insert(value.indices,"indices");
        iom.insert(value.values,"values");
        ostr<<name<<" {\n";
        iom>>ostr;
        ostr<<" }\n";
        return true;
    }

//in an explosion of laziness, we reuse inputfield specializations for iostream overload
    template<typename U> 
    std::ostream& operator<<(std::ostream& ostr, const CrsMatrix<U>& cm)
    { 
        IOMap iom;
        iom.insert(cm.pops,"options");
        iom.insert(cm.wr,"rows");
        iom.insert(cm.wc,"cols");
        iom.insert(cm.rpoints,"rowpointers");
        iom.insert(cm.indices,"indices");
        iom.insert(cm.values,"values");
        iom>>ostr;
        return ostr;
    }

    template<typename U> std::istream& operator>> (std::istream& istr, CrsMatrix<U>& cm)
    { 
        IField<CrsMatrix<U> > ifcm(cm,"coordmatrix");
        ifcm<<istr;
        return istr;
    }

#endif //ends ifdef __MATRIX_CRS_H

    
//COORD MATRIX
#ifdef __MATRIX_COORD_H
/*****************************************************************
    *   define specialized member functions for IOMap input-output  *
 *****************************************************************/
    template<typename U> 
    class IField<CoordMatrix<U> > : public IFBase
    {
        private:
            CoordMatrix<U>& value;
        public:
    using IFBase::operator>>;
    using IFBase::operator<<;
            IField(CoordMatrix<U>& var, const std::string& nname) : IFBase(nname), value(var) {}
            IField(CoordMatrix<U>& var, const std::string& nname, const CoordMatrix<U>& defval) : IFBase(nname), value(var) { value=defval; flags |=iff_optional; }
    
            bool operator>>(std::ostream& ostr) const;
            bool operator<<(std::istream& istr);
    
            operator CoordMatrix<U>() { return value; }
    };


    template<typename U> 
    bool IField<CoordMatrix<U> >::operator<<(std::istream& istr) 
    {
        IOMap iom; iom.flags=if_warn_all;
        iom.insert(value.wr,"rows");
        iom.insert(value.wc,"cols");
        iom.insert(value.data,"data");
        iom<<istr;
        flags|=iff_set; flags&=~iff_nonvalid;
        if ((iom["rows"].flags | iom["cols"].flags | iom["data"].flags) & iff_nonvalid) { flags|=iff_nonvalid; return false; }
    
    //now we must check the data array, to sort it and to make sure that 
    //there are no elements out of the matrix bounds.
        std::sort(value.data.begin(), value.data.end(), compare<U>);
        value.sz=value.data.size();
        if (value.data[value.sz-1].i>=value.wr)
            ERROR("Index out of row bounds in coord matrix.");
        for (typename CoordMatrix<U>::index_type k=0; k<value.sz; ++k)
        {
            if (value.data[k].j>=value.wc)
                ERROR("Index out of column bounds in coord matrix.");
        }
    
        return true;
    }

    template<typename U> 
    bool IField<CoordMatrix<U> >::operator>>(std::ostream& ostr) const
    {
        IOMap iom;
        iom.insert(value.wr,"rows");
        iom.insert(value.wc,"cols");
        iom.insert(value.data,"data");
        ostr<<name<<" {\n";
        iom>>ostr;
        ostr<<" }\n";
        return true;
    }

//const version
    template<typename U> 
            class IField<const CoordMatrix<U> > : public IFBase
    {
        private:
            const CoordMatrix<U>& value;
        public:
    using IFBase::operator>>;
    using IFBase::operator<<;
            IField(const CoordMatrix<U>& var, std::string nname) : IFBase(nname), value(var) {}
            bool operator>>(std::ostream& ostr) const;
            operator CoordMatrix<U>() { return value; }
    };

    template<typename U> 
            bool IField<const CoordMatrix<U> >::operator>>(std::ostream& ostr) const
    {
        IOMap iom;
        iom.insert(value.wr,"rows");
        iom.insert(value.wc,"cols");
        iom.insert(value.data,"data");
        ostr<<name<<" {\n";
        iom>>ostr;
        ostr<<" }\n";
        return true;
    }

//in an explosion of laziness, we reuse inputfield specializations for iostream overload
    template<typename U> 
    std::ostream& operator<<(std::ostream& ostr, const CoordMatrix<U>& cm)
    { 
        IOMap iom;
        iom.insert(cm.wr,"rows");
        iom.insert(cm.wc,"cols");
        iom.insert(cm.data,"data");
        iom>>ostr;
        return ostr;
    }

    template<typename U> std::istream& operator>> (std::istream& istr, CoordMatrix<U>& cm)
    { 
        IField<CoordMatrix<U> > ifcm(cm, "", cm);
        ifcm<<istr;
        return istr;
    }
#endif //ends #ifdef __MATRIX_COORD_H
    
//FULL MATRIX
#ifdef __MATRIX_FULL_H
/*****************************************************************
    *   define specialized member functions for IOMap input-output  *
 *****************************************************************/
    //iofield for options
template <class U> 
class IField<MatrixOptions<FMatrix<U> > > : public IFBase
{
    private:
        MatrixOptions<FMatrix<U> >& value;
    public:
        IField(MatrixOptions<FMatrix<U> >& var, std::string nname) : IFBase(nname), value(var) {}
        IField(MatrixOptions<FMatrix<U> >& var, std::string nname, const FMatrix<U>& defval) : IFBase(nname), value(var) { value=defval; flags |=iff_optional; }

        bool operator>>(std::ostream& ostr) const
        {
            IOMap iom;  std::string atm;
            switch (value.atnorm)
            { 
                case at_normone: atm="norm_one"; break;
                case at_norminf: atm="norm_inf"; break;
                case at_normfrob: atm="norm_frob"; break;
                default: ERROR("Unsupported norm mode");
            };
            iom.insert(value.atthresh,"at_threshold");
            iom.insert(atm,"at_norm");
            iom.insert(value.atnodiag,"at_nodiag");
            iom.insert(value.atrel,"at_relative");
            ostr << name << " {\n"<<iom<<"\n}\n";
            return true;
        }

        bool operator<<(std::istream& istr)
        {
            IOMap iom;  std::string atm;
            iom.insert(value.atthresh,"at_threshold",0.);
            iom.insert(atm,"at_norm",std::string("norm_inf"));
            iom.insert(value.atnodiag,"at_nodiag");
            iom.insert(value.atrel,"at_relative",false);
            
            istr>>iom;
            flags|=iff_set; flags &=~iff_nonvalid;
            if (!istr.good()) { flags|=iff_nonvalid; return false; }
            
            if (atm=="norm_one") value.atnorm=at_normone;
            else if (atm=="norm_inf") value.atnorm=at_norminf;
            else if (atm=="norm_frob") value.atnorm=at_normfrob;
            else { flags |=iff_nonvalid; return false; }
            
            if (value.atthresh<0.) { flags|=iff_nonvalid; return false; }
            
            return true;
    //!!SANITY CHECK MUST BE DONE
        }

        operator MatrixOptions<FMatrix<U> >() { return value; }
};

template <class U> 
        class IField<const MatrixOptions<FMatrix<U> > > : public IFBase
{
    private:
        const MatrixOptions<FMatrix<U> >& value;
    public:
        IField(const MatrixOptions<FMatrix<U> >& var, std::string nname) : IFBase(nname), value(var) {}

        bool operator>>(std::ostream& ostr) const
        {
            IOMap iom;  std::string atm;
            switch (value.atnorm)
            { 
                case at_normone: atm="norm_one"; break;
                case at_norminf: atm="norm_inf"; break;
                case at_normfrob: atm="norm_frob"; break;
                default: ERROR("Unsupported norm mode");
            };
            iom.insert(value.atthresh,"at_threshold");
            iom.insert(atm,"at_norm");
            iom.insert(value.atnodiag,"at_nodiag");
            iom.insert(value.atrel,"at_relative");
            ostr << name << " {\n"<<iom<<"\n}\n";
            return true;
        }
        operator MatrixOptions<FMatrix<U> >() { return value; }
};

    
    
    template<typename U> 
class IField<FMatrix<U> > : public IFBase
{
private:
    FMatrix<U>& value;
public:
    using IFBase::operator>>;
    using IFBase::operator<<;
    IField(FMatrix<U>& var, const std::string& nname) : IFBase(nname), value(var) {}
    IField(FMatrix<U>& var, const std::string& nname, const FMatrix<U>& defval) : IFBase(nname), value(var) { value=defval; flags |=iff_optional; }

    bool operator>>(std::ostream& ostr) const;
    bool operator<<(std::istream& istr);

    operator FMatrix<U>() { return value; }
};


template<typename U> 
bool IField<FMatrix<U> >::operator<<(std::istream& istr) 
{
    IOMap iom; iom.flags=if_warn_all;
    iom.insert(value.wr,"rows");
    iom.insert(value.wc,"cols");
    iom.insert(value.data,"data");
    iom<<istr;
    flags|=iff_set; flags&=~iff_nonvalid;
    if ((iom["rows"].flags | iom["cols"].flags | iom["data"].flags) & iff_nonvalid) { flags|=iff_nonvalid; return false; }

    return true;
}

template<typename U> 
bool IField<FMatrix<U> >::operator>>(std::ostream& ostr) const
{
    IOMap iom;
    iom.insert(value.wr,"rows");
    iom.insert(value.wc,"cols");
    ostr<<name<<" {\n";
    iom>>ostr;
//we put data afterwards, and even if we print it out just as a vector, we "format" it by rows, to make 
//it more readable
    ostr<<"data "<<value.data.size()<<"\n";
    unsigned long k=0, j=0;
    for (k=0; k<value.data.size();++k)
    {
        ostr<<value.data[k]<<" ";
        ++j; if (j>=value.wc) {j=0; ostr<<"\n";}
    }

    ostr<<" }\n"; //close the field
    return true;
}

//const version
template<typename U> 
class IField<const FMatrix<U> > : public IFBase
{
    private:
        const FMatrix<U>& value;
    public:
        using IFBase::operator>>;
    using IFBase::operator<<;
    IField(const FMatrix<U>& var, const std::string& nname) : IFBase(nname), value(var) {}
        bool operator>>(std::ostream& ostr) const;
        operator FMatrix<U>() { return value; }
};

template<typename U> 
bool IField<const FMatrix<U> >::operator>>(std::ostream& ostr) const
{
    IOMap iom;
    iom.insert(value.wr,"rows");
    iom.insert(value.wc,"cols");
    ostr<<name<<" {\n";
    iom>>ostr;
//we put data afterwards, and even if we print it out just as a vector, we "format" it by rows, to make 
//it more readable
    ostr<<"data "<<value.data.size()<<"\n";
    unsigned long k=0, j=0;
    for (k=0; k<value.data.size();++k)
    {
        ostr<<value.data[k]<<" ";
        ++j; if (j>=value.wc) {j=0; ostr<<"\n";}
    }

    ostr<<" }\n"; //close the field
    return true;
}

//in an explosion of laziness, we reuse inputfield specializations for iostream overload
    template<typename U> 
    std::ostream& operator<<(std::ostream& ostr, const FMatrix<U>& cm)
    {
        IOMap iom;
        iom.insert(cm.wr,"rows");
        iom.insert(cm.wc,"cols");
        iom>>ostr;
    //we put data afterwards, and even if we print it out just as a vector, we "format" it by rows, to make 
    //it more readable
        ostr<<"data "<<cm.data.size()<<"\n";
        unsigned long k=0, j=0;
        for (k=0; k<cm.data.size();++k)
        {
            ostr<<cm.data[k]<<" ";
            ++j; if (j>=cm.wc) {j=0; ostr<<"\n";}
        }
        return ostr;
    }

    template<typename U> std::istream& operator>> (std::istream& istr, FMatrix<U>& cm)
    { 
        IField<FMatrix<U> > ifcm(cm, "", cm);
        ifcm<<istr;
        return istr;
    }
#endif //ends #ifdef __MATRIX_FULL_H

//parallel CRS matrix.
#ifdef __MATRIX_PCRS_H
    //we are just copying the code for CrsMatrix options. probably this can be done in a more elegant fashion
    template <class U> 
    class IField<MatrixOptions<PCrsMatrix<U> > > : public IFBase
    {
        private:
            MatrixOptions<PCrsMatrix<U> >& value;
        public:
    using IFBase::operator>>;
    using IFBase::operator<<;
            IField(MatrixOptions<PCrsMatrix<U> >& var, std::string nname) : IFBase(nname), value(var) {}
            IField(MatrixOptions<PCrsMatrix<U> >& var, std::string nname, const PCrsMatrix<U>& defval) : IFBase(nname), value(var) { value=defval; flags |=iff_optional; }

            bool operator>>(std::ostream& ostr) const
            {
                IOMap iom;  std::string atm;
                switch (value.atnorm)
                { 
                    case at_normone: atm="norm_one"; break;
                    case at_norminf: atm="norm_inf"; break;
                    case at_normfrob: atm="norm_frob"; break;
                    default: ERROR("Unsupported norm mode");
                };
                iom.insert(value.atthresh,"at_threshold");
                iom.insert(atm,"at_norm");
                iom.insert(value.atnodiag,"at_nodiag");
                iom.insert(value.atrel,"at_relative");
                ostr << name << " {\n"<<iom<<"\n}\n";
                return true;
            }
    
            bool operator<<(std::istream& istr)
            {
                IOMap iom;  std::string atm;
                iom.insert(value.atthresh,"at_threshold",0.);
                iom.insert(atm,"at_norm",std::string("norm_inf"));
                iom.insert(value.atnodiag,"at_nodiag");
                iom.insert(value.atrel,"at_relative",false);
                
                istr>>iom;
                flags|=iff_set; flags &=~iff_nonvalid;
                if (!istr.good()) { flags|=iff_nonvalid; return false; }
                
                if (atm=="norm_one") value.atnorm=at_normone;
                else if (atm=="norm_inf") value.atnorm=at_norminf;
                else if (atm=="norm_frob") value.atnorm=at_normfrob;
                else { flags |=iff_nonvalid; return false; }
                
                if (value.atthresh<0.) { flags|=iff_nonvalid; return false; }
                
                return true;
        //!!SANITY CHECK MUST BE DONE
            }

            operator MatrixOptions<PCrsMatrix<U> >() { return value; }
    };
    
    template <class U> 
    class IField<const MatrixOptions<PCrsMatrix<U> > > : public IFBase
    {
        private:
            const MatrixOptions<PCrsMatrix<U> >& value;
        public:
    using IFBase::operator>>;
    using IFBase::operator<<;
            IField(const MatrixOptions<PCrsMatrix<U> >& var, std::string nname) : IFBase(nname), value(var) {}

            bool operator>>(std::ostream& ostr) const
            {
                IOMap iom;  std::string atm;
                switch (value.atnorm)
                { 
                    case at_normone: atm="norm_one"; break;
                    case at_norminf: atm="norm_inf"; break;
                    case at_normfrob: atm="norm_frob"; break;
                    default: ERROR("Unsupported norm mode");
                };
                iom.insert(value.atthresh,"at_threshold");
                iom.insert(atm,"at_norm");
                iom.insert(value.atnodiag,"at_nodiag");
                iom.insert(value.atrel,"at_relative");
                ostr << name << " {\n"<<iom<<"\n}\n";
                return true;
            }
            operator MatrixOptions<PCrsMatrix<U> >() { return value; }
    };
    
    
/*****************************************************************
    *   define specialized member functions for IOMap input-output  *
 *****************************************************************/
    template<typename U> 
    class IField<PCrsMatrix<U> > : public IFBase
    {
        private:
            PCrsMatrix<U>& value;
        public:
    using IFBase::operator>>;
    using IFBase::operator<<;
            IField(PCrsMatrix<U>& var, const std::string& nname) : IFBase(nname), value(var) {}
            IField(PCrsMatrix<U>& var, const std::string& nname, const PCrsMatrix<U>& defval) : IFBase(nname), value(var) { value=defval; flags|=iff_optional; }

            bool operator>>(std::ostream& ostr) const;
            bool operator<<(std::istream& istr);

            operator PCrsMatrix<U>() { return value; }
    };


    template<typename U> 
    bool IField<PCrsMatrix<U> >::operator<<(std::istream& istr) 
    {
        //we just use conversion to CRS matrix: read on node 0 and convert
        CrsMatrix<U> lmat(0,0);
        IField<CrsMatrix<U> > lif(lmat,"");
        
        if (value.myrank==0) {istr >> lif;}  //reads only on node 0
        value=lmat; //calls local CRS to PCRS conversion and hope for the best :-)
        //!we should check that everything was fine.... we'll do that....
        return true;
    }

    template<typename U> 
    bool IField<PCrsMatrix<U> >::operator>>(std::ostream& ostr) const
    {
        CrsMatrix<U> lmat;
        IField<CrsMatrix<U> > lif(lmat,"");
        if (value.myrank==0) lmat.resize(value.wr,value.wc);
        
        lmat=value; //PCRS to local CRS conversion
        if (value.myrank==0) ostr<< lif; //writes only on node 0
        return true;
    }

//const version
    template<typename U> 
    class IField<const PCrsMatrix<U> > : public IFBase
    {
        private:
            const PCrsMatrix<U>& value;
        public:
    using IFBase::operator>>;
    using IFBase::operator<<;
            IField(const PCrsMatrix<U>& var, const std::string& nname) : IFBase(nname), value(var) {}
            bool operator>>(std::ostream& ostr) const;
            operator PCrsMatrix<U>() { return value; }
    };

    template<typename U> 
    bool IField<const PCrsMatrix<U> >::operator>>(std::ostream& ostr) const
    {
        CrsMatrix<U> lmat;
        IField<CrsMatrix<U> > lif(lmat,"");
        if (value.myrank==0) lmat.resize(value.wr,value.wc);
        
        lmat=value; //PCRS to local CRS conversion
        if (value.myrank==0) ostr<< lif; //writes only on node 0
        return true;
    }

//in an explosion of laziness, we reuse inputfield specializations for iostream overload
/*!    template<typename U> 
    std::ostream& operator<<(std::ostream& ostr, const PCrsMatrix<U>& cm)
    {
        IOMap iom;
        iom.insert(cm.wr,"rows");
        iom.insert(cm.wc,"cols");
        iom>>ostr;
    //we put data afterwards, and even if we print it out just as a vector, we "format" it by rows, to make 
    //it more readable
        ostr<<"data "<<cm.data.size()<<"\n";
        unsigned long k=0, j=0;
        for (k=0; k<cm.data.size();++k)
        {
            ostr<<cm.data[k]<<" ";
            ++j; if (j>=cm.wc) {j=0; ostr<<"\n";}
        }
        return ostr;
    }
*/
    template<typename U> std::istream& operator>> (std::istream& istr, PCrsMatrix<U>& cm)
    { 
        IField<PCrsMatrix<U> > ifcm(cm, "", cm);
        ifcm<<istr;
        return istr;
    }
#endif //ends #ifdef __MATRIX_PCRS_H
    
    
};  //ends namespace toolbox


#endif //ends #ifndef __MATRIX_IO_H 
