//// bisogna ottimizzare mu e sigma, 
//quindi invece di leggere dal file faccio tipo grid search 

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
  double sigma_v[21];
  double mu_v[21];
  for (int i=0; i<21; i++){  //da 0.5 a 1.5 
    sigma_v[i] = 0.5 + i*0.05;
    mu_v[i] = 0.5 + i*0.05;
  }
  double best_ene = 999., best_sigma=999. , best_mu=999.; //valori a caso alti così nel caso capisco che c'è qualcosa che non va
  ofstream print;

  Input(); //Inizialization without mu and sigma

  print.open(file_name);

  for( int m=0; m<21; m++){
    for(int s=0; s<21; s++){

      sigma = sigma_v[s];
      mu = mu_v[m];
      //trova la delta migliore associata a quei mu e sigma
      delta = delta_opt();
      ///faccio andare 100000 volte per equilibrare (in realtà non sembra ci sia molta equilibrazione...)
      for (int equi=0; equi<100000; equi++){
        Move();
      }

      //reset 
      accepted = 0;
      attempted = 0;

      //cout << "global av vecchia: " << glob_av << endl;

      for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
      {
        Reset(iblk);   //Reset block averages
        //cout << "global av nuova: " << glob_av << endl;
        for(int istep=1; istep <= nstep; ++istep)
        {
          Move();
          Accumulate(); //Update block averages
        }
        Averages(iblk);   //Calculates results for current block and updates global average
      }

      if((glob_av/double(nblk)) <= best_ene){
        best_ene = glob_av/double(nblk);
        best_sigma = sigma_v[s];
        best_mu = mu_v[m];
      }
      cout << "***************************************" << endl;
      cout << mu << "   " << sigma << "   "  <<  glob_av/double(nblk) << "   "  << best_ene << endl;
      cout << "***************************************" << endl;
      print << mu << "   " << sigma << "   "  <<  glob_av/double(nblk) << "   "  << best_ene << endl;
    }
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
  ReadInput.open("input_senzasigmamu.dat");

  ReadInput >> x;
  ReadInput >> nblk;   //numero di blocchi
  ReadInput >> nstep;  // numero di step per blocco
  ReadInput >> volta; //per fare nomi file diversi

   ReadInput.close();

/////nome del file in cui scrivere
  file_name = "opt_sigma_mu_" + volta + ".dat";

//messaggi sul terminale
  cout << "-----------------------------------------------------" << endl;
  cout << "starting position: " << x << endl;
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
//ora mi interessa solo il risultato finale di tutti i blocchi per confrontarlo con quello con altri sigma/mu
//quindi tolgo tutta la parte di stampa su file
void Averages(int iblk)  //Print results for current block
{
    
    if((iblk%10)==0){
    cout << "current mu: " << mu << ", current sigma: " << sigma << ", current delta: " <<delta << endl;
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << double(accepted)/double(nstep) << endl;
	}
    
    stima_ene = blk_av/nstep; //Energy
    glob_av += stima_ene;
    glob_av2 += stima_ene*stima_ene;
    err_ene = Error(glob_av,glob_av2,iblk);

    if((iblk%10)==0){
    cout << "glob_av: " << glob_av/double(iblk) << endl;
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

//////////////////OTTIMIZZAZIONE_DELTA//////////////////
double delta_opt()
{
  int n_acc=0, v = 1000;
  int best=9999;
  double best_delta, delta_loc;

  for (int i=0 ; i<2000; i++){ //vari valore di step
    delta_loc = (double(i) + 1)*0.01;
    n_acc = conta_accettati(delta_loc, v);
    //cerr << n_acc << endl;
    if (abs(best-int(v/2)) >= abs(n_acc-int(v/2))) {
      best=n_acc;
      best_delta = delta_loc;
    }
  }

  cerr <<"best accepted: " << best << ", best delta: " << best_delta << endl;

  return best_delta;

}

////////////////////CONTA_ACCETTATI///////////////////////////
/////data una certa lunghezza dell step e la posizione di partenza della simulazione,
//fa una sim di 1000 passi mi dice quanti sono i passi accettati (ne voglio 500)
// questi sono estratti uniformemente in segmento (x-delta, x+delta)
int conta_accettati(double d, int v){

  int n_acc=0;
  double x_loc = x;
  double A, r, x_new;

  for (int itime=1; itime<=v; itime++){

///////ESTRAGGO UN NUOVO PUNTO//////////////////
    x_new = x_loc + rnd.Rannyu(-d , d);

//////////VALUTO ALPHA//////////////////////////////////////////////////////
    A = min(1., psi2(x_new) / psi2(x_loc));

    //accettato?
    r = rnd.Rannyu();
    if (r <= A) {
      x_loc=x_new;
      n_acc++;
    }

  }

  //cerr << n_acc << "   " << v;

  return n_acc;

}



