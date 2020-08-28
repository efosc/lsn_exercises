/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  if(equi) Equilibration();

  cout << "Equilibration -> Done" << endl;
  equi = false;
  cout << "Setting equi to false" << endl;

  //Simulation
  /*for(int iblk=1; iblk <= nblk; ++iblk) 
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }*/
  ConfFinal();

  //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;
  ifstream ReadConf;
  double s_read;  //to read initial configuration

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations  
// legge la temperatura, n di spin, J, h e se deve ripartire da una configurazione precedente o no
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
  steptot = nblk*nstep;

  ReadInput >> start;

  ReadInput >> equi;

  ReadInput >> filename;

  //ReadInput >> nblk_equi;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Total number of steps = " << steptot << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  if(start){
  	cout << "Reading initial configuration from file 'config.final'"<< endl;
  	ReadConf.open("config.final");
  	for (int i=0; i<nspin; ++i){
  		ReadConf >> s_read;
  		s[i] = s_read;
  	} 
  	ReadConf.close();
	}else if(!start){   //mette a caso ogli spin, quindi è come se fossimo a T infinita
	  for (int i=0; i<nspin; ++i)
	  {
	    if(rnd.Rannyu() >= 0.5) s[i] = 1;
	    else s[i] = -1;
	  }
	}
  
//Evaluate energy etc. of the initial configuration
  Measure();
  //Instant();

//Print initial values for the potential energy and virial //(virial??)
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
  cout << "Initial heat capacity = " << walker[ic] <<endl;
  cout << "Initial magnetization = " << walker[im] <<endl;
  cout << "Initial susceptibility = " << walker[ix] <<endl << endl;
}


void Move(int metro)
{
  int o;
  double p, ene_diff, r, sm;
  //double p, energy_old, energy_new, sm;
  //double energy_up, energy_down;

  for(int i=0; i<nspin; ++i) //-> 1 MC step = flip all the spins at one time
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    r = rnd.Rannyu();
    if(metro==1) //Metropolis
    {
// INCLUDE YOUR CODE HERE
    	sm = -1.*s[o]; //try flip the spin
    	ene_diff = 2.*(Boltzmann(sm, o));
    	p = min(1., exp(-ene_diff/temp));  //acceptance
    	if (r <= p) {
    		s[o] = sm; 
    		accepted ++;
    	}
    	attempted ++;
    }
    else //Gibbs sampling
    {
// INCLUDE YOUR CODE HERE 
//calocolo la probabilità che lo spin sia 1 e se r è <= di quella lo metto a 1
//quindi è come prima ma con sm=1 
    	ene_diff = -2. * Boltzmann(1, o);
    	p = 1./( 1 + exp(- ene_diff/temp));
    	if (r <= p) s[o] = 1.;
    	else s[o] = -1.;
    	accepted ++;
    	attempted ++;
   		if (accepted != attempted) cerr << "Error occured in Gibbs sampling"<<endl;
    }
    
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins  
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);  //somma tutte le energie in quello stato
    m += s[i];  //somma tutti gli spin di quello stato
	}
  walker[iu] = u;  //non dovrebbe essere diviso il numero di spin? sì lo fa dopo non so perché
  walker[ic] = u*u; //energia del sistema ^2 
  walker[im] = m;  //somma degli spin
  walker[ix] = m*m;   
} //per ora walker è come se fosse misura e misura al quad ??

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages  
//fa la somma della somma delle energie e della somma degli spin in un blocco
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
// e aggiorna in glob_av la media aggiungendo il risultato di ogni blocco
{
  string filename1, filename2;
   ofstream Ene, Heat, Mag, Chi;
   ofstream final_ene , final_heat , final_mag, final_chi;
   const int wd=12;

   if(equi){
    filename1 = "equi_";
   }else{
    filename1 = "meas_";
   }

   if(metro){
    filename2 = "metro_";
   }else{
    filename2 = "gibbs_";
   }

    if(equi) cout << "Equilibrating the system" << endl; 
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    if(h==0.0){

      //Ene.open("Outfiles/" + filename1 + filename2 + filename + "_output.ene.1",ios::app);
      stima_u = blk_av[iu]/(double)blk_norm/(double)nspin; //Energy
      glob_av[iu]  += stima_u;
      glob_av2[iu] += stima_u*stima_u;
      err_u=Error(glob_av[iu],glob_av2[iu],iblk);
      //Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      //Ene.close();
      //ora scrivo l'ultimo nel file per i vari valori della temperatura
      if((iblk == nblk) && (!equi)) {
      	final_ene.open("Outfiles/" + filename1 + filename2 + filename + "_final_vals_ene.dat", ios::app);
      	final_ene << iblk << setw(wd) << temp << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      	final_ene.close();
      }

      //Heat.open("Outfiles/" + filename1 + filename2 + filename + "_output.heat.1",ios::app);  //heat capacity (perché non divido per n spin?)
      stima_c = ((blk_av[ic]/(double)blk_norm) - (blk_av[iu]/(double)blk_norm)*(blk_av[iu]/(double)blk_norm))/temp/temp/(double)nspin;
      glob_av[ic]  += stima_c;
      glob_av2[ic] += stima_c*stima_c;
      err_c=Error(glob_av[ic],glob_av2[ic],iblk);
      //Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
      //Heat.close();
      //ora scrivo l'ultimo nel file per i vari valori della temperatura
      if((iblk == nblk) && (!equi)) {
      	final_heat.open("Outfiles/" + filename1 + filename2 + filename + "_final_vals_heat.dat", ios::app);      	
      	final_heat << iblk << setw(wd) << temp << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
      	final_heat.close();
      }

      //Chi.open("Outfiles/" + filename1 + filename2 + filename + "_output.chi.1",ios::app);
      stima_x = blk_av[ix]/blk_norm/(double)nspin/temp; //susceptivity
      glob_av[ix]  += stima_x;
      glob_av2[ix] += stima_x*stima_x;
      err_x=Error(glob_av[ix],glob_av2[ix],iblk);
      //Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      //Chi.close();
      //ora scrivo l'ultimo nel file per i vari valori della temperatura
      if((iblk == nblk) && (!equi)) {
      	final_chi.open("Outfiles/" + filename1 + filename2 + filename + "_final_vals_chi.dat", ios::app);      	      	
      	final_chi << iblk << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      	final_chi.close();
      }
  
    }else if (h==0.02){

      Mag.open("Outfiles/" + filename1 + filename2 + filename + "_output.mag.1",ios::app);
      stima_m = blk_av[im]/blk_norm/(double)nspin; //magnetization
      glob_av[im]  += stima_m;
      glob_av2[im] += stima_m*stima_m;
      err_m=Error(glob_av[im],glob_av2[im],iblk);
      Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      Mag.close();
      //ora scrivo l'ultimo nel file per i vari valori della temperatura
      /*if(iblk == nblk) {
      	final_mag.open("Outfiles/" + filename1 + filename2 + filename + "_final_vals_mag.dat", ios::app);
      	final_mag << iblk << setw(wd) << temp << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      	final_mag.close();
      }*/

    } else {cout << "h should be either 0.0 or 0.02" << endl;}


    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void Equilibration(){  // faccio andare per un certo numero di blocchi prima per equilibrare
	if(!equi) cout << "PROBLEM: equi option should be false" << endl;
  int nstep_equi = 10000;
  int nblk_equi = 100;
  for(int iblk=1; iblk <= nblk_equi; ++iblk) 
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep_equi; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
}

 /*void Instant(){
  const double space = 12.0;
  string filename1, filename2;
  ofstream print;
  if(equi){
    filename1 = "equi_";
  }else{
    filename1 = "meas_";
  }

  if(metro){
    filename2 = "metro_";
  }else{
    filename2 = "gibbs_";
  }
  print.open(filename1 + filename2 + filename + "_instant.dat", ios::app);
  if(h==0){
  	//cerr << walker[iu] << "   " << walker[iu]*walker[iu] << "  " << walker[ic] << endl;
  	//cerr << walker[iu]/(double)nspin << "   " << walker[iu]*walker[iu]/(double)nspin/(double)nspin<< "  " << walker[ic]/(double)nspin << endl;
    stima_u = walker[iu]/(double)nspin; //Energy
    stima_c = (walker[ic] - walker[iu]*walker[iu])/temp/temp/(double)nspin; //heat capacity
    stima_x = walker[ix]/(double)nspin/temp; //susceptivity
    print << stima_u << setw(space) << stima_c << setw(space) << stima_x << endl;
  } else if(h==0.02){
    stima_m = walker[im]/(double)nspin; //magnetization
    print << stima_m << endl;
  } else{ 
    cout << "h should be either 0.0 or 0.02" << endl;
  }

  print.close();

 }*/



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
