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

  int M = 10E5;  //numero totale di lanci
  int N = 100;   //numero singoli esperimenti
  int L = int(M/N);  //numero di lanci per singolo esperimento

  double throws[M];
  double stdev1[N];      // momento secondo nel singolo esperimento
  double stdev2[N];      // momento secondo ^2 nel singolo esperimento
  double cum1_stdev[N];  // media delle medie considerando 1 ,2, ... , N esperimenti
  double cum2_stdev[N];  // valore medio ^2 delle medie considerando 1, 2, ... , N esperimenti
  double err_stdev[N];
  ofstream out;

  //genera il vettore con tutti i lanci casuali
  for(int i=0; i<M; i++) throws [i] = rnd.Rannyu();
  //singoli esperimenti
  double somma_locale;
  for(int i=0; i<N; i++){
    somma_locale = 0;
    for(int j=0;j<L; j++)
      somma_locale += (throws[i*L + j] - 0.5)*(throws[i*L + j] -0.5);
    stdev1[i] = somma_locale/double(L);
    stdev2[i] = somma_locale/double(L) * somma_locale/double(L);
  }

  //media dei risultati in 1,2,...,N esperimenti
  double somma, somma2;
  for (int i=0; i<N; i++){
    somma=0;
    somma2=0;
    for (int j=0; j<i+1; j++){
      somma += stdev1[j];
      somma2 += stdev2[j];
    }
    cum1_stdev[i] = somma/double(i+1);
    cum2_stdev[i] = somma2/double(i+1);
    err_stdev[i] = error(cum1_stdev[i] , cum2_stdev[i] , i);
  }

  out.open("varianze2.dat");
  for (int a=0; a<N; a++) out<<cum1_stdev[a]<<endl;
  out.close();
  out.clear();

  out.open("errori_var2.dat");
  for (int a=0; a<N; a++) out<<err_stdev[a]<<endl;
  out.close();
  out.clear();


  return 0;

}



double error ( double media , double media_quadrati, int n) {
    if(n==0) return 0;
    else return sqrt((media_quadrati - media*media)/n);
}