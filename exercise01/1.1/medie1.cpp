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
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
double error( double media , double media_quadrati, int n);


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

   //for(int i=0; i<20; i++){
     // cout << rnd.Rannyu() << endl;
   //}

  int M = 100;  //numero totale di lanci
  int N = 10;   //numero singoli esperimenti
  int L = int(M/N);  //numero di lanci per singolo esperimento
  double throws[M];  // vettore con tutti i lanci
  double av1[N];  //valore medio di r nel singolo esperimento
  double av2[N];  // valore medio ^2 di r nel singolo esperimento
  double cum1[N];  // media delle medie considerando 1 ,2, ... , N esperimenti
  double cum2[N];  // valore medio ^2 delle medie considerando 1, 2, ... , N esperimenti
  double err[N];   //devst delle medie considerando 1, 2, ... , N esperimenti
  ofstream out;

  //cerr << "ciao" << endl;

  //genero il vettore che contiene tutti i lanci casuali

  for (int i=0; i < M; i++){ 
        //cerr << "entro nel for"<<endl;
    throws [i] = rnd.Rannyu();
        //cerr << "fatto numero " << i << endl;
  }

  cerr << "generato vettore throws" << endl;

  //calcolo le medie locali (quindi risultato del singolo esperimento)
  for (int j=0; j<N; j++){
    double media_locale=0;
      for(int k=0; k<L; k++) 
        media_locale = (media_locale * k + throws[(j*L+ k)]) / (k+1) ;
    av1[j] = media_locale;
    av2[j] = media_locale * media_locale;
  }

  cerr << "calcolati av1 e av2" << endl;

  //calcolo la media dei risultati di 1 ,2, ..., N esperimenti (media delle medie)
  // e il valore medio dei quadrati dei risultati degli esperimenti
  for (int n_exp=0; n_exp<N; n_exp++){
    double media = 0;
    double media_quadrati = 0;
    for (int j=0; j<n_exp; j++){
        media = (media * j + av1[j]) / (j+1);
        media_quadrati = (media_quadrati * j + av2[j]) / (j+1);
    }
    cum1[n_exp] = media;
    cum2[n_exp] = media_quadrati;
    err[n_exp] = error ( media , media_quadrati , n_exp+1 );
  }

  cerr << "calcolati cum1 e cum2" << endl;

  out.open("medie.dat");

  for (int a=0; a<N; a++) out<< cum1[a]<<endl;

  out.close();

  out.open("errori.dat");

  for (int a=0; a<N; a++) out<< err[a]<<endl;

   out.close();

  rnd.SaveSeed();
  return 0;
}




double error ( double media , double media_quadrati, int n) {
    if(n==0) return 0;
    else return sqrt((media_quadrati - media*media)/n);
}




/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
