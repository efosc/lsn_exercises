#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double error ( double media , double media_quadrati, int n);

int main (int argc, char *argv[]){

	//seed
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

// dati
  double S0 = 100;
  double T=1;
  double K = 100;
  double r = 0.1;
  double sigma = 0.25;
// per calcolare le medie
  int M = 1E4;
  int N = 100;  //numero esperimenti
  int L = int(M/N);  //numero prove per ogni esperimento

  ofstream out;

//punto 1: direct sampling S(T)
  double S, somma_call=0, somma_put=0, somma_call2=0, somma_put2=0;
  double media_call=0, media_put=0, media_call2=0, media_put2=0;
  double call_cum[N], put_cum[N], err_call[N], err_put[N];

  for(int j=0; j<N; j++){  // faccio N esperimenti
    somma_call=0;
    somma_put=0;
    somma_call2=0;
    somma_put2=0;
    for (int i=0; i<L; i++){ // ogni esperimento è formato da L prove
      S = S0*exp((r-0.5*sigma*sigma)*T + sigma*rnd.Gauss(0,T));
      somma_call += max(0., S-K)*exp(-r*T);  //sommo all'interno del singolo esperimento
      somma_put += max(0., K-S)*exp(-r*T);   //sommo all'interno del singolo esperimento
    }
    somma_call /= double(L);
    somma_put /= double(L);         //risultato di un singolo esperimento
    somma_call2 = somma_call*somma_call;
    somma_put2 = somma_put*somma_put;      //risultato di un singolo esperimento^2
   	
   	media_call = (media_call*double(j) + somma_call)/double(j+1);   //calcolo la media su j+1 esperimenti (j=0, ... , N-1)
   	media_put = (media_put*double(j) + somma_put)/double(j+1);
   	media_call2 = (media_call2*double(j) + somma_call2)/double(j+1);
   	media_put2 = (media_put2*double(j) + somma_put2)/double(j+1);

   	call_cum[j] = media_call;   //salvo nel vettore
   	put_cum[j] = media_put;
   	err_call[j] = error(media_call , media_call2, j);
   	err_put[j] = error(media_put, media_put2, j);
  }

  //stampo su file
  out.open("direct_call.dat");
  for(int a=0; a<N; a++) out<<call_cum[a]<<endl;
  out.close();
	out.clear();

  out.open("direct_call_err.dat");
  for(int a=0; a<N; a++) out<<err_call[a]<<endl;
  out.close();
	out.clear();

  out.open("direct_put.dat");
  for(int a=0; a<N; a++) out<<put_cum[a]<<endl;
  out.close();
	out.clear();

  out.open("direct_put_err.dat");
  for(int a=0; a<N; a++) out<<err_put[a]<<endl;
  out.close();
	out.clear();

//punto 2: sampling a tempi discreti (ogni 0.01T così dopo 100 volte si arriva a T)
  double z, S_vecchio, deltat = 0.01*T;

  somma_call=0, somma_put=0, somma_call2=0, somma_put2=0;
  media_call=0, media_put=0, media_call2=0, media_put2=0;
  for (int l=0; l<N; l++) {call_cum[l]=0; put_cum[l]=0; err_call[l]=0; err_put[l]=0;}
  	//riazzero per sicurezza, in realtà non dovrebbe esserci bisogno

  for(int j=0; j<N; j++){  // faccio N esperimenti
    somma_call=0;
    somma_put=0;
    somma_call2=0;
    somma_put2=0;
    for (int i=0; i<L; i++){ // ogni esperimento è formato da L prove
      S_vecchio=S0;
      for(int k=0; k<100; k++){  //calcolo S finale 
        z = rnd.Gauss(0,T);
        S = S_vecchio*exp((r-0.5*sigma*sigma)*deltat + sigma*z*sqrt(deltat));
        S_vecchio=S;
      }
      somma_call += max(0., S-K)*exp(-r*T);  //sommo all'interno del singolo esperimento
      somma_put += max(0., K-S)*exp(-r*T);   //sommo all'interno del singolo esperimento
    }
    somma_call /= double(L);
    somma_put /= double(L);         //risultato di un singolo esperimento
    somma_call2 = somma_call*somma_call;
    somma_put2 = somma_put*somma_put;      //risultato di un singolo esperimento^2
   	
   	media_call = (media_call*double(j) + somma_call)/double(j+1);   //calcolo la media su j+1 esperimenti (j=0, ... , N-1)
   	media_put = (media_put*double(j) + somma_put)/double(j+1);
   	media_call2 = (media_call2*double(j) + somma_call2)/double(j+1);
   	media_put2 = (media_put2*double(j) + somma_put2)/double(j+1);

   	call_cum[j] = media_call;   //salvo nel vettore
   	put_cum[j] = media_put;
   	err_call[j] = error(media_call , media_call2, j);
   	err_put[j] = error(media_put, media_put2, j);
  }

	//stampo su file
  out.open("discrete_call.dat");
  for(int a=0; a<N; a++) out<<call_cum[a]<<endl;
  out.close();
	out.clear();

  out.open("discrete_call_err.dat");
  for(int a=0; a<N; a++) out<<err_call[a]<<endl;
  out.close();
	out.clear();

  out.open("discrete_put.dat");
  for(int a=0; a<N; a++) out<<put_cum[a]<<endl;
  out.close();
	out.clear();

  out.open("discrete_put_err.dat");
  for(int a=0; a<N; a++) out<<err_put[a]<<endl;
  out.close();
	out.clear();
  
  return 0;
}

double error ( double media , double media_quadrati, int n) {
    if(n==0) return 0;
    else return sqrt((media_quadrati - media*media)/n);
}
