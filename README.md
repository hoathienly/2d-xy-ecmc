# 2d-xy-ecmc
/*  Success on Thursday, November 15th, 2018, Vietnam time
    Posted on Sunday, November 18th, 2018
    Revised in the week after that.
    Full version posted on Sunday, November 25th, 2018

    This implements the event chain Monte Carlo for the classical 2D XY model.

    Ref.:
    1. (Main paper) Manon Michel, Johannes Mayer, and Werner Krauth, Event-chain Monte Carlo for classical continuous spin models, 2015, arXiv:1508.06541v1
    2. Johannes Mayer, Werner Krauth, Lifting algorithms in the XY model, Internship report, ENS Laboratoire de Physique Statistique, 2015, http://www.lps.ens.fr/~krauth/images/7/72/Stage_Mayer_Johannes_2015.pdf
    3. Manon Michel, Irreversible Markov chains by the factorized Metropolis filter: Algorithms and applications in particle systems and spin models, PhD Thesis, ENS Paris, 2016, https://tel.archives-ouvertes.fr/tel-01394204
    4. Manon Michel, Sebastian C. Kapfer, and Werner Krauth, Generalized event-chain Monte Carlo: Constructing rejection-free global-balance algorithms from infinitesimal steps, 2014, arXiv:1309.7748v2
    5. Yoshihiko Nishikawa, Manon Michel, Werner Krauth, and Koji Hukushima, Event-chain algorithm for the Heisenberg model: Evidence for z ≃ 1 dynamic scaling, 2015, arXiv:1508.05661v2 
    6. Persi Diaconis, Susan Holmes and Radford M. Neal, Analysis of a nonreversible Markov chain sampler, 2000, https://projecteuclid.org/download/pdf_1/euclid.aoap/1019487508
    7. Mersenne Twister libraries, http://www.helsinki.fi/~rummukai/lectures/montecarlo_oulu/prog/index.html

*/

# Compile: 
gcc -o 2d-xy-ecmc 2d-xy-ecmc.c mersenne_inline.c -lm 

# Execute: 
1d-harmonic-ecmc L Tmin Tmax MP s

where  
        L the lattice size,

        Tmin the minimum temperature,
        
        Tmax the maximum temperature,
        
        MP a multiple of pi, e.g. take a measurement after every MP*pi interval,
        
        s a multiple of the number of samples, the real number of samples is S = s*sS, where sS = 10000
        
# Monday, Nov. 26, 2018
I believe the last bug is found today.
