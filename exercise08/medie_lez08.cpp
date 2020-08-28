#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "random.h"
#include <string>
#include "funzioni.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization

  ///faccio andare 100000 volte per equilibrare (in realtà non sembra ci sia molta equilibrazione...)
  for (int equi=0; equi<100000; equi++){
    Move();
  }

  //reset 
  accepted = 0;
  attempted = 0;

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }

  return 0;
}


//////////////////////INPUT////////////////////////////////
void Input(){
////////seed///////////////////////
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
      Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  rnd.SaveSeed();

////////////////////////////////////

  ifstream ReadInput;

//legge il file input
  ReadInput.open("input.dat");

  ReadInput >> x;
  ReadInput >> mu;
  ReadInput >> sigma;
  ReadInput >> delta;
  ReadInput >> nblk;   //numero di blocchi
  ReadInput >> nstep;  // numero di step per blocco
  ReadInput >> volta; //per fare nomi file diversi

   ReadInput.close();

/////nome del file in cui scrivere
  file_name = "medie_" + volta + ".dat";

//messaggi sul terminale
  cout << "-----------------------------------------------------" << endl;
  cout << "starting position: " << x << endl;
  cout << "mu: " << mu << endl;
  cout << "sigma: " <<sigma << endl;
  cout << "step length: " << delta << endl;
  cout << "writing into: " << file_name << endl;
  cout << "number of blocks: " << nblk << endl;
  cout << "number of steps per block: " << nstep << endl;
  cout << "-----------------------------------------------------" << endl;

  return;

}

/////////////RESET///////////////

void Reset(int iblk){ //Reset block averages
   
  if(iblk == 1)
   {
    glob_av = 0;
    glob_av2 = 0;
   }
  blk_av = 0;
  accepted = 0;
  attempted = 0;
}

/////////////////PROBABILITY////////////////////////////////////////
double psi2(double x){
  double espo1 = (x-mu)*(x-mu)/2./sigma/sigma;
  double espo2 = (x+mu)*(x+mu)/2./sigma/sigma;
  double psi = exp(- espo1) + exp(-espo2);
  return psi*psi;  //tanto è reale
}


///////////////////MOVE(): PASSO CON PROBABILITA' UNIFORME///////////////////////////////////////
void Move(){
  double A, r, x_new;

 ///////ESTRAGGO UN NUOVO PUNTO//////////////////
      x_new = x + (-1 + rnd.Rannyu()*2)*delta; 

//////////VALUTO ALPHA//////////////////////////////////////////////////////
    A = min(1., psi2(x_new) / psi2(x));

    //accettato?
    r = rnd.Rannyu();
    if (r <= A) {
      x=x_new;
      accepted++;
    }

    attempted ++;

  return;
}

///////////////////////////////////////////////////////


//////////////ACCUMULATE/////////////

///devo misurare l'energia applicata a psi, con
//T = d2(f)/dx2 , V= x^4 - (5/2)*x^2
// e poi fare H*psi/psi

void Accumulate(void){ //Update block averages  
//calcolo l'energia in quel momento e aggiungo alla media del blocco 
  double espo1 = (x-mu)*(x-mu)/2./sigma/sigma;
  double espo2 = (x+mu)*(x+mu)/2./sigma/sigma;
  T = (-1./pow(sigma,2) + (x-mu)*(x-mu)/pow(sigma,4))*exp(-espo1) 
      + ( -1./pow(sigma,2) + (x+mu)*(x+mu)/pow(sigma,4) ) *exp(-espo2);
  T = - 0.5 * T / (exp(- espo1) + exp(-espo2));
  V = (pow(x,4) - 2.5*pow(x,2));
  energy = T + V;
  blk_av = blk_av + energy;
}

////////////////AVERAGES////////////////////
void Averages(int iblk)  //Print results for current block
{
   ofstream Energy;
   const int wd=12;
    
    if((iblk%10)==0){
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << double(accepted)/double(nstep) << endl;
	}
    
    Energy.open(file_name,ios::app);
    stima_ene = blk_av/nstep; //Energy
    glob_av += stima_ene;
    glob_av2 += stima_ene*stima_ene;
    err_ene = Error(glob_av,glob_av2,iblk);
    Energy << setw(wd) << iblk << setw(wd) << accepted << setw(wd) << stima_ene << setw(wd) << glob_av/(double)iblk << setw(wd) << err_ene << endl;
    Energy.close();

    if((iblk%10)==0){
    cout << "<E> in this block: " << stima_ene << endl;
    cout << "----------------------------" << endl << endl;
	}
}

/////////////ERRORS/////////////////////////////////
double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}



