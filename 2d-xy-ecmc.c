/*  Success on Thursday, November 15th, 2018, Vietnam time
    Posted on Sunday, November 18th, 2018

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

#include <studio.h>
#include <stdlib.h>
#include <math.h>

#include "mersenne.h"

struct myst {
    double mean;
    double dev;
};

void prng(long seed, int warm);
void init(double *h);
void flip(int flag);
double get_event();
double pair_event(int i, int j);
void set_index(int *up, int *dn);
void measure(double *lattice, int *up, int *dn, double *m, double *e);
void set_zero();
void sample_collect();
void post_measure();
void post_write(FILE *f);
void draw(double *lattice, int *up, char name[100]);

#define u(i, j) (-J*cos(h[i] - h[j]))
#define du(i, j) (v*J*sin(h[i] - h[j]))

void flip(int flag){
    double vp, temp, tp;

    temp = 0;
    vp = 0;
    tp = 0;

    if(ref<mp){
        tp = get_event();
        vp = tp;
    }

    if(flag == equ){
        h[current] += tp;
        return;
    }

    while(vp + ref < mp){
        h[current] += tp;
        tp = get_event();
        temp = vp;
        vp += tp;
    }

    if (ref > mp){
        h[current] += mp;
        ref -= mp;
        return;
    }

    h[current] += (mp - temp - ref);
    ref = (vp + ref - mp);

    return;
}

double pair_event(int i, int j){
    int np;
    double dinit, init, du, dp, r, temp;

    dinit = du(i, j);

    if(dinit > 0) init = u(i, j);
    else init = -J;

    r = mersenne();
    du = - kB*T*log(r);

    np = du / (2*J);
    du -= (2*J*np);
    du += init;

    dp = 0;
    temp = 0;

    if(du > J){
        du -= (2*J);
        temp = pi;
    }

    if(dinit < 0){
        dp = h[i] - h[j];
        dp = acos(cos(dp));
        if(dp > pi / 2) dp += (pi/2);
        if(dp < pi / 2) dp = dp;
    }

    if(dinit > 0) {
        dp = h[i] - h[j];
        dp = acos(cos(dp));
        if(dp < pi / 2) dp = dp;
        if(dp > pi / 2) dp += (pi/2);
        dp = -dp;
    }

    dp += temp;
    dp += 2*np*pi + acos(-J*du);

    return dp;
}
