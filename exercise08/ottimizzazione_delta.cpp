#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include "random.h"
#include <string>
#include "funzioni.h"

using namespace std;

int main (int argc, char *argv[]){

  Input(); //seed + legge i dati in input
	
  ofstream acc;
  int best=9999;
  double best_delta;
  int n_acc;

  acc.open(file_name);
  for (int i=0 ; i<2000; i++){ //vari valore di step
    delta = (double(i) + 1)*0.01;
    n_acc = conta_accettati(delta, x);
    if (abs(best-int(ntimes/2)) >= abs(n_acc-int(ntimes/2))) {
      best=n_acc;
      best_delta = delta;
    }   
    acc << "delta: " << delta << " accepted delta: " <<n_acc << endl;
    if((i%100) == 0)cout << "delta: " << delta << " accepted steps: " <<n_acc << endl;
  }
 
  cout <<  "best delta: " << best_delta << "  accepted steps: " << best << endl;

  acc.close();

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

  string orbital_name, transition_prob;
  ///char itime;
  ifstream ReadInput;

//legge il file input
  ReadInput.open("input_opt.dat");

  ReadInput >> x;
  ReadInput >> mu;
  ReadInput >> sigma;
  ReadInput >> ntimes;
  ReadInput >> volta; //per fare nomi file diversi

   ReadInput.close();

/////nome del file in cui scrivere
  file_name = "accepted_" + volta + ".dat";

  cout << "start: " << x << endl;
  //cerr << itime << endl;

  cout << "writing into " << file_name << endl;

  return;

}

/////////////////PROBABILITY////////////////////////////////////////
double psi2(double x){
  double espo1 = (x-mu)*(x-mu)/2./sigma/sigma;
  double espo2 = (x+mu)*(x+mu)/2./sigma/sigma;
  double psi = exp(- espo1) + exp(-espo2);
  return psi*psi;  //tanto è reale
}

////////////////////CONTA_ACCETTATI///////////////////////////
/////data una certa lunghezza dell step e la posizione di partenza della simulazione,
//fa una sim di 1000 passi mi dice quanti sono i passi accettati (ne voglio 500)
// questi sono estratti uniformemente in segmento (x-delta, x+delta)
int conta_accettati(double delta, double x){

  int n_acc=0;
  double A, r, x_new;

  for (int itime=1; itime<=ntimes; itime++){

///////ESTRAGGO UN NUOVO PUNTO//////////////////
    x_new = x + rnd.Rannyu(-delta , delta);

//////////VALUTO ALPHA//////////////////////////////////////////////////////
    A = min(1., psi2(x_new) / psi2(x));

    //accettato?
    r = rnd.Rannyu();
    if (r <= A) {
      x=x_new;
      n_acc++;
    }

  }

  return n_acc;

}









