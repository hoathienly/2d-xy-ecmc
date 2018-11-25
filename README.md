# 2d-xy-ecmc
/*  Success on Thursday, November 8th, 2018, Vietnam time
    Posted on Sunday, November 18th, 2018
    Revised in the week after that.
    Full version on Sunday, November 25th, 2018

    This implements the event driven rejection free Monte Carlo for the simple 1D harmonic potential well.

    Ref.: 
    1. (Main) E.A.J.F. Peters and G. de With, Rejection-free Monte-Carlo sampling for general potentials, 2012, arXiv:1112.1263v3
    2. Mersenne Twister libraries, http://www.helsinki.fi/~rummukai/lectures/montecarlo_oulu/prog/index.html

*/
# Compile: 
gcc -o 2d-xy-ecmc 2d-xy-ecmc.c mersenne_inline.c -lm 

#Execute: 
1d-harmonic-ecmc L Tmin Tmax MP s

where   L the lattice size,

        Tmin the minimum temperature,
        
        Tmax the maximum temperature,
        
        MP a multiple of pi, e.g. take a measurement after every MP*pi interval,
        
        s a multiple of the number of samples, the real number of samples is S = s*sS, where sS = 10000
