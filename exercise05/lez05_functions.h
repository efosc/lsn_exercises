#ifndef __Lez05__
#define __Lez05__

//nel file input
double step;  //nell'ottimizzazione è quello che vario
double x_0 , y_0, z_0;
bool orbital, T;
int ntimes; //numero di passi nella simulazione
int nstep; //numero di passi per block
int nblk; //numero di blocchi (nonnè nel file input, si ricava da ntimes/nstep)
int n_acc;  //numero di passi accettati
std::string volta;

// averages -> qui misuro solo la distanza dall'origine 
double blk_av; //media di r nel blocco -> risultato del blocco
double glob_av,glob_av2;
double stima_dist;
double err_dist;

double dist; //distanza dall'origine
Random rnd;

std::string file_name;


void Input(void);
void Reset(int iblk);
void Accumulate();
void Averages(int iblk);

double prob_s(double x, double y, double z);
double prob_p(double x, double y, double z);
int conta_accettati_gauss(double step, double x_old, double y_old, double z_old);
int conta_accettati_unif(double step, double x_old, double y_old, double z_old);  //sbagliata
int conta_accettati_unif2(double step, double x_old, double y_old, double z_old);
void gauss_step();  //move()
void unif_step();

double Error(double sum, double sum2, int iblk);


















#endif