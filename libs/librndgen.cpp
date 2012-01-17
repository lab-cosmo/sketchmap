#ifndef __RNDGEN_H
#define __RNDGEN_H

#include "rndgen.hpp"

namespace toolbox{
    
/* StdRndUniform */
void StdRndUniform::setstatus(const long& state)
{ srand48(state);  pstate=state; }


/* Mersenne twister */
void MTRndUniformL::seed(const unsigned long& seed)
{
    /* seed the mersenne twister from a long seed */
    unsigned short int i;
    
    unsigned long s=(seed==0?4357:seed); /* the default seed is 4357, to avoid periodicity */
    pstate.mtvec.resize(MT_N);
    std::valarray<unsigned long>& mt=pstate.mtvec;
    
    mt[0]= s & 0xffffffffUL;

    for (i = 1; i < MT_N; i++)
    {
      /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
         Ed. p.106 for multiplier. */

        mt[i] =(1812433253UL * (mt[i-1] ^ (mt[i-1] >> 30)) + i);      
        mt[i] &= 0xffffffffUL;
    }

    pstate.pos = i;
}

} //ends namespace toolbox

#endif //ends __RNDGEN_H
