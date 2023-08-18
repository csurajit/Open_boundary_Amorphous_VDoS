
#include<stdint.h>

    enum{
        A1      = 21,
        A2      = 35,
        A3      = 4,
        D1      = 2685821657736338717ULL,
        PRIME   = 4101842887655102017ULL                        // largest 62-bit prime ?
    };
    static const double INV2TO64  = 5.42101086242752217E-20;    // 1/(2^64)
    static const double INV2TO63  = 1.084202172485504434E-19;   // 1/(2^63)
    static uint64_t val = PRIME;



    // Adapted from Numerical Recipes 3rd Edition, Section 7.1.3
    // The period is 1.8e19

    // The initialisation function
    // It needs to be called once, at the beginning
    // Use time(0) and/or pid as input
    static void initRanq1( uint64_t j ){ val ^= j; };
    // The main RNG function
    // generating 64bit-representable, uniformly distributed integers
    static uint64_t int64Ranq1(){
        val ^= val >> A1;
        val ^= val << A2;
        val ^= val >> A3;
        return val * D1;
    }
    //static uint32_t int32Ranq1(){ return ( uint32_t )int64Ranq1(); }
    // Random numbers between 0.0 and 1.0
    static double dblRanq1(){ return ( INV2TO64 * int64Ranq1() ); };
    // Random numbers between -1.0 and 1.0
    static double dblRanq1Sym(){ return ( ( INV2TO63 * int64Ranq1() ) - 1.0 ); };
