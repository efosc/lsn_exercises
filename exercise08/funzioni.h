#ifndef __ES8__
#define __ES8__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

// qui ho una sola osservabile, l'energia
double energy;
double T, V;
//parameters, observables
/*const int m_props=1000;
int n_props, iv, iw, igofr;
double vtail,ptail,bin_size,nbins,sd;*/
//double walker[m_props]; 

// averages
double blk_av,blk_norm,accepted,attempted;
double glob_av,glob_av2;
double stima_ene, err_ene;

//configuration
double x;   //"posizione" della particella (coordinata da cui dipende la psi)

double mu, sigma;

//simulation
std::string volta;
std::string file_name;
int ntimes;
int nstep, nblk;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void Measure(void);
//double Pbc(double);
double Error(double,double,int);
void Equilibrate(void);
double psi2(double);

int conta_accettati(double, int );
double delta_opt(void);

#endif