#ifndef __Mol__
#define __Mol__

//parameters, observables
const int m_props=5;
int n_props;  //= m_props
int iv,ik,it,ie, iw;  //indici nell'array di misure
double walker[m_props];

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_pot, stima_ekin, stima_etot, stima_temp, stima_press;  //averages.cpp
//double stima_pot_inst, stima_ekin_inst, stima_etot_inst, stima_temp_inst;  //con_eq3.cpp
double err_pot,err_ekin,err_etot,err_temp, err_press;
  //20000 steps -> 2000 misure; 100 blocchi -> 20 misure per ogni blocco?? ->poche no???
  // ogni blocco è lungo 200 steps (=20 misure)

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];  //positions
double vx[m_part],vy[m_part],vz[m_part];  //speed

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, blk_length, nblk;
double delta;
bool old, equilibration;
std::string filename, dirname; //per cambiare nome al file in cui scrivi senza ricompilare


//functions
void Input(void);
void Read_config0(void);
void Casual_vel_estimate_old(void);
void Read_both(void);
//void estimate_vel(void);
void Move(void);
void ConfFinal(void);
void ConfOld(void);
//void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Rescale(void);
double Error ( double , double , int ); //con valori medi
double Error2 (double, double, int); //con somme
void Read_old(void);
//double Measure_epot();
//double Measure_ekin();
void Reset(int iblk);
void Accumulate(void);
void Averages(int iblk);
void Instant(void);

#endif