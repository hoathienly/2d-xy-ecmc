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

#include <stdio.h>
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

struct myst energy, magnetization;
int current, lift, L, *up, *dn, MP, s;
double *h, e, m, T, S, mp, cv, chi;
long seed = 19880601;
double ref = 0;
double Tmin = 0.5;
double Tstep = 0.01;
double Tmax = 2.0;

#define u(i, j) (-J*cos(h[i] - h[j]))
#define du(i, j) (v*J*sin(h[i] - h[j]))
#define TH 1000000
#define N (L*L)
#define sS 10000
#define kB 1
#define J 1
#define pi 3.141592653589793
#define v 1
#define P 4
#define equ 1
#define upd 2

int main(int argc, char **argv){
    int i, j, k;
    FILE *f;
    char name[100];

    L = atoi(argv[1]);
    Tmin = atof(argv[2]);
    Tmax = atof(argv[3]);
    MP = atoi(argv[4]);
    s = atoi(argv[5]);

    mp = MP * pi;
    S = s * sS;

    h = (double *) malloc (N * sizeof(double));
    j=2*L;

    up = (int *) malloc (j * sizeof(int));
    dn = (int *) malloc (j * sizeof(int));

    sprintf(name, "%2d_%d_%d.txt", L, MP, s);
    f=fopen(name, "w");

    prng(seed, TH);
    init(h);

    set_index(up, dn);

    for(T = Tmin; T < Tmax; T += Tstep) {
        ref = 0;
        i = L * mersenne();
        j = L * mersenne();
        lift = i + L * j;

        for(i = 0; i < TH; i++){
            flip(equ);
        }

        set_zero();

        for(i = 0; i < S; i++){
            flip(upd);
            measure(h, up, dn, &m, &e);
            h[current] += ref;
            sample_collect();
        }
        post_measure();
        post_write(f);
    }   

    fclose(f);
    free(h);
    free(up);
    free(dn);
}

void prng(long seed, int warm){
    int i;

    seed_mersenne(seed);
    for(i = 0; i < warm; i++) mersenne();
}

void init(double *h){
    int i, j;

    for(i = 0; i < L; i++){
        for(j = 0; j < L; j++){
            h[i + L * j] = 1; // 2 * pi * mersenne();
        }
    }
}

void flip(int flag){
    double vp, temp, tp;

    temp = 0;
    vp = 0;
    tp = 0;

    if(ref < mp){
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

double get_event(){
    int j, k;

    current = lift;

    j = current % L;
    k = current / L;

    int i, neighbour[P] = {j + L * up[k], j + L * dn[k], up[j] + L * k, dn[j] + L * k};

    double out, temp;

    out = pair_event(current, neighbour[0]);
    lift = neighbour[0];

    for(i = 0; i < P; i++){
        temp = pair_event(current, neighbour[i]);
        if(out > temp){
            out = temp;
            lift = neighbour[i];
        }
    }
    return out;
}

void set_index(int *up, int *dn){
    int i;

    for(i = 0; i < 2 * L; i++){
        up[i] = (i + 1) % L;
        dn[i] = (i + L - 1) % L;
    }
}

double pair_event(int i, int j){
    int np;
    double dinit, init, du, dp, r, temp;

    dinit = du(i, j);

    if(dinit > 0) init = u(i, j);
    else init = -J;

    r = mersenne();
    while (r < 1e-100) r = mersenne();

    du = - kB * T * log(r);

    np = du / (2 * J);
    du -= (2 * J * np);
    du += init;

    dp = 0;
    temp = 0;

    if(du > J){
        du -= (2 * J);
        temp = pi;
    }

    dp = h[i] - h[j];
    dp = acos(cos(dp));

    if(dinit > 0) dp = -dp;

    dp += temp;
    dp += 2 * np * pi + acos(-J * du);

    return dp;
}

void measure(double *lattice, int *up, int *dn, double *m, double *e){
    double h, m1, m2, spale;
    int i, j;
    h = 0;  m1 = 0; m2 = 0;
    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++) {
            m1 += cos(lattice[j + L * i]);
            m2 += sin(lattice[j + L * i]);
            h += -cos(lattice[j + L * i] - lattice[j + L * up[i]]);
            h += -cos(lattice[j + L * i] - lattice[up[j] + L * i]);
        }
    m1 /= N;
    m2 /= N;
    h /= N;
    spale = m1 * m1 + m2 * m2;
    spale = sqrt(spale);
    *m = spale;
    *e = h;
}

void set_zero(){
	energy.mean = 0;
	energy.dev = 0;
	magnetization.mean = 0;
	magnetization.dev = 0;
}

void sample_collect(){
	energy.mean += e / S;
	energy.dev += e * e / S;
	magnetization.mean += m / S;
	magnetization.dev += m * m / S;
}

void post_measure(){
	energy.dev = energy.dev - energy.mean * energy.mean;
    cv = energy.dev * N / (kB * T * T);
	magnetization.dev = magnetization.dev - magnetization.mean * magnetization.mean;
    chi = magnetization.dev * N / (kB * T);
	energy.dev = sqrt(energy.dev / (S - 1));
	magnetization.dev = sqrt(magnetization.dev / (S - 1));
}

void post_write(FILE * f){
	fprintf(f, "%lf\t", T); // 1
	fprintf(f, "%lf\t%lf\t%.12lf\t", energy.mean, energy.dev, cv); // 2 3
	fprintf(f, "%lf\t%lf\t%.12lf\n", magnetization.mean, magnetization.dev, chi); // 4 5

	printf("%lf\t", T);
	printf("%lf\t%lf\t%.12lf\t", energy.mean, energy.dev, cv); // 2 3
	printf("%lf\t%lf\t%.12lf\n", magnetization.mean, magnetization.dev, chi); // 4 5
}
