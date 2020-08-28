#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double funz (double x);
double error ( double media , double media_quadrati, int n);

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

/////PARTE 1: integrale campionando distribuzione uniforme//////////////

  int M=1E4, N=100, L=int(M/N);  //N numero di eseprimenti, L numero di lanci in ogni esperimento
  double cum[N], errori[N];
  double stima_int=0, stima_int2=0, media=0, media2=0;
  ofstream out;

  for (int j=0; j<N; j++){   //N esperimenti
    stima_int=0;
    stima_int2 =0; //in realtà questo non c'è bisogno penso
    for (int i=0; i<L; i++) stima_int += funz(rnd.Rannyu());   //L prove per ogni esperimento
    stima_int /= L;     
    stima_int2 = stima_int*stima_int;  //risultato singolo esperimento

    media = (media*double(j) + stima_int)/double(j+1);
    media2 = (media2*double(j) + stima_int2)/double(j+1);   //calcolo media di j+1 esperimenti (j=0, ... , N-1)

    cum[j] = media;
    errori[j] = error ( media, media2, j);
  }

  out.open("integrale_unif.dat");
  for (int a=0; a<N; a++) out << cum[a] << endl;
  out.close();
  out.clear();

  out.open("integrale_unif_err.dat");
  for (int a=0; a<N; a++) out << errori[a] << endl;
  out.close();
  out.clear();

//////PARTE 2: importance sampling //////////
  double y;
  for (int a=0; a<N; a++) {cum[a]=0; errori[a]=0;} //inutile
  stima_int=0;
  stima_int2=0;
  media=0;
  media2=0;  //riazzero tutto


  for (int j=0; j<N; j++){        //N esperimenti
    stima_int=0;
    stima_int2=0;
    for (int i=0; i<L; i++) {     //L prove per ogni esperimento
      y = (1 - sqrt(1-rnd.Rannyu()));
      stima_int += funz(y)*0.5/(1-y);
    }
    stima_int /= L;
    stima_int2 = stima_int*stima_int;   //risultato singolo esperimento
    
    media = (media*double(j) + stima_int)/double(j+1);
    media2 = (media2*double(j) + stima_int2)/double(j+1);   //calcolo media di j+1 esperimenti (j=0, ... , N-1)

    cum[j] = media;
    errori[j] = error ( media, media2, j);

  }

  out.open("integrale_imp.dat");
  for (int a=0; a<N; a++) out << cum[a] << endl;
  out.close();
  out.clear();

  out.open("integrale_imp_err.dat");
  for (int a=0; a<N; a++) out << errori[a] << endl;
  out.close();
  out.clear();


  return 0;

}

double funz (double x){
  return M_PI*0.5*cos(x*M_PI*0.5);
}

double error ( double media , double media_quadrati, int n) {
    if(n==0) return 0;
    else return sqrt((media_quadrati - media*media)/n);
}