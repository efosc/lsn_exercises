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

  //int M = 2E7; //numero totale di lanci;
  int N = 100;   //numero singoli esperimenti (=numero di volte che calcolo pi)
  int L = 1E6;  //numero di lanci per singolo esperimento

  //double throws[M];   //vettore con tutti i risultati dei lanci
  double ave[N];      // vettore con risultati dei singoli esperimenti
  double ave2[N];     // vettore con risultati^2 dei singoli esperimenti
  double cum1[N];     // media delle medie considerando 1 ,2, ... , N esperimenti
  double cum2[N];     // valore medio ^2 delle medie considerando 1, 2, ... , N esperimenti
  double err[N];      // errori sulla media delle medie (risultati exp) considerando 1, 2, ... , N esperimenti
  ofstream out;

  double ago=0.4 , d=2.0; //lunghezza ago e distanza tra linee
  int n_hit =0, a=0;
  double angolo, centro, dist, x, y;

  for(int j=0; j<N; j++){
    for (int i=0; i<L; i++){
    	centro=(rnd.Rannyu())*d;
      dist=min(centro, abs(d-centro));
      for (a=0; a<1000; ++a){ //faccio al massimo 1000 prove
        x= rnd.Rannyu();
        y= rnd.Rannyu();
        if(x*x + y*y <= 1) break;
      }
      if(a==1000) cerr << "angolo non trovato"<<endl;
      else angolo = atan(y/x);

      if (dist <= ago*0.5*sin(angolo)) n_hit++;   // se intereseca la linea conto
    }
    cerr << n_hit << "  " <<(2.0*L*ago)/(n_hit*d)<<endl;
    ave[j] = (2.0*L*ago)/(n_hit*d);  //calcolo un valore di pi (risultato di un exp)
    ave2[j] = ave[j]*ave[j];
    n_hit=0;
  }



  double somma, somma2;
  for (int i=0; i<N; i++){
    somma=0;
    somma2=0;
    for (int j=0; j<i+1; j++){
      somma += ave[j];
      somma2 += ave2[j];
    }
    cum1[i] = somma/double(i+1);
    cum2[i] = somma2/double(i+1);
    err[i] = error(cum1[i] , cum2[i] , i);
  }

  out.open("pi_single_exp.dat");   // singoli esperimenti
  for (int b=0; b<N; b++) out << ave[b] << endl;
  out.close();
  out.clear();

  out.open("err_pi_forseok.dat");
  for (int b=0; b<N; b++) out << err[b] << endl;
  out.close();
  out.clear();

  out.open("pi_forseok.dat");
  for (int b=0; b<N; b++) out << cum1[b] << endl;
  out.close();
  out.clear();


  /*int r=0;
  for (int i=0; i<1000; i++){
    centro=(rnd.Rannyu())*d;
    dist=min(centro, abs(d-centro));
    for(r=0; r<1000; r++){   //faccio al massimo 1000 prove
      x=rnd.Gauss(0.0, 1.);
      y=rnd.Gauss(0.0, 1.);  //distribuzione di angolo uniforme in (-pi; pi)
      if((x>0) && (y>0)) break;
    }
    if(r==1000) cerr << "angolo non trovato"<<endl;
    else angolo=atan2(y,x);  //prendo solo quelli in 0;pi/2

    out << angolo / M_PI <<endl;
    }
    //if (d==1000 ) cerr<<"angolo non trovato"<<endl;
    //else

  out.close();*/
  // con gauss non è così uniforme, meglio altro metodo
  return 0;
}


double error ( double media , double media_quadrati, int n) {
    if(n==0) return 0;
    else return sqrt((media_quadrati - media*media)/n);
  }