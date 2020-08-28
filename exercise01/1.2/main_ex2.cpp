#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
double error( double, double, int );


int main (int argc, char *argv[]){

   Random rnd;
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

  ofstream out;

////esponenziale, lambda=1./////

  int N[4] = {1,2,10, 100}; //numero di x su cui faccio la media
  int n_volte = 1E4;
  string files[4] = {"exp1_nuovo.dat", "exp2_nuovo.dat" , "exp10_nuovo.dat" , "exp100_nuovo.dat"};
  double medie[4][n_volte]; //faccio una matrice di 4 righe (una per ogni valore di N) e 10^4 cols 
                    //(che contengono le 10^4 ) volte che ho calcolato la media su N x_i
  double somma=0;

  for (int c=0; c<4; c++) {
    for (int i=0; i<n_volte; i++){
      somma = 0;
      for (int j=0; j<N[c]; j++) somma +=  rnd.Exp(1.);  //calcolo la somma su N elementi
      medie[c][i] = somma/double(N[c]);
    }
  out.open(files[c]);
  for (int a=0; a<n_volte; a++) out << medie[c][a] << endl;
  out.close();
  out.clear();
  }

///////Cauchy-Lorentz, gamma = 1., mu=0. //////////////
  string files2[4] = {"cauchy1_nuovo.dat", "cauchy2_nuovo.dat" , "cauchy10_nuovo.dat" , "cauchy100_nuovo.dat"};
  //riazzero tutto, anche se non necessario per medie penso
  for (int j=0; j<4; j++){
    for (int k=0; k<n_volte;k++) medie[j][k]=0;
  }
  somma=0;

  for (int c=0; c<4; c++) {
    for (int i=0; i<n_volte; i++){
      somma = 0;
      for (int j=0; j<N[c]; j++) somma +=  rnd.CauchyLorentz(0., 1.);  //calcolo la somma su N elementi
      medie[c][i] = somma/double(N[c]);
    }
  out.open(files2[c]);
  for (int a=0; a<n_volte; a++) out << medie[c][a] << endl;
  out.close();
  out.clear();
  }

/////////uniforme////////////
  string files3[4] = {"unif1_nuovo.dat", "unif2_nuovo.dat" , "unif10_nuovo.dat" , "unif100_nuovo.dat"};
  //riazzero tutto, anche se non necessario per medie penso
  for (int j=0; j<4; j++){
    for (int k=0; k<n_volte;k++) medie[j][k]=0;
  }
  somma=0;

  for (int c=0; c<4; c++) {
    for (int i=0; i<n_volte; i++){
      somma = 0;
      for (int j=0; j<N[c]; j++) somma +=  rnd.Rannyu();  //calcolo la somma su N elementi
      medie[c][i] = somma/double(N[c]);
    }
  out.open(files3[c]);
  for (int a=0; a<n_volte; a++) out << medie[c][a] << endl;
  out.close();
  out.clear();
  }

  return 0;
}

