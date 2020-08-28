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

  int M=100; //numero di sottointervalli nell'intervallo [0;1]
  double chi2[M]; //vettore che contiene tutti i valori del chi2
  double numero; //numero casuale osservato (vedo in quale intervallo va a finire)
  int n = 1E4; //numero di volte che estraggo un numero casuale
  int conteggi [100];
  double somma;
  ofstream out;

  //for (int b=0; b<(n_tot*M); b++) throws[b]
  for ( int r=0; r < 100; r++){
    for (int b=0; b<M; b++) conteggi[b] =0;
    somma=0;
    for (int j=0; j<n; j++){
      numero = rnd.Rannyu(); //estraggo un numero casuale
      //vedo in che intervallo si trova
      for(int i=0; i<M; i++){
        if( (numero > (i*(1./M))) && (numero <= ((i+1)*(1./M))) ) conteggi[i]++;
      }
    }
    for(int k=0; k<M; k++){
      somma += (conteggi[k] - (double(n)/double(M)))*(conteggi[k] - (double(n)/double(M)));
    }
    chi2[r] = somma /(double(n)/double(M));
  }

  out.open("chi_nuovo.dat");

  for (int a=0; a<M; a++) out << chi2[a] << endl;


  //for(int a=0; a<M; a++) cerr<<conteggi[a] <<endl;

  /*double somma;
  for (int j=0; j<100; j++){
  	somma=0;
		for (int i=1; i<M; i++){
			n_i=0;
			for(int k=0; k<n_tot; k++) {
				x_k = rnd.Rannyu();       //estraggo un numero
				if (x_k < 1./double(M)) n_i++;  //conto quante volte il numero  cade nel primo intervallo
			//if ((k<100) & (j==0) & (i==1)) cerr << x_k << "  " << n_i<< endl;
			}
	  	somma += (n_i - double(n_tot)/double(M)) * (n_i - double(n_tot)/double(M));
	  }
	  chi2[j] = somma/double(n_tot/M);
	}
	out.open("chi2.dat");

	for (int a=0; a<M; a++) out<<chi2[a]<<endl;*/



  return 0;
}