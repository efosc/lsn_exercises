#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "posizione.h"

using namespace std;


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

  posizione pos(0.,0.,0.);
  int n_max_passi=100;
  int M = 1E4;  
  int N = 100;
  int L = int(M/N);
  double medie[n_max_passi];
  double err[n_max_passi];
  ofstream out;
  double pos_fin=0, pos_fin2=0, media=0, media2=0;

////RETICOLO CUBICO//////

  for (int i=0; i<n_max_passi; i++){
    media=0;
    media2=0;
    for(int w=0; w<N; w++){
      pos_fin=0;
      pos_fin2=0;   //inutile in realtà
      for (int j=0; j<L; j++ ) {
        for(int k=0; k<(i+1); k++) pos.singolo_passo(rnd);  //raggiunge posizione finale
        pos_fin += pos.dist2();  //distanza al quadrato
        pos.Setx(0.);
        pos.Sety(0.);
        pos.Setz(0.);
      }
      pos_fin = sqrt(pos_fin/double(L)); 
      pos_fin2 = pos_fin * pos_fin; //risultato singolo esperimento

      media = (media*double(w) + pos_fin) /double (w+1);
      media2 = (media2*double(w) + pos_fin2) /double (w+1);  //medie degli esperimenti
    }

    medie[i] = media;
    err[i] = error ( media, media2, N-1);   //salvo nel vettore la media finale con i passi
  }

  out.open("rw_cubico.dat");
  for(int a =0; a < n_max_passi; a++) out<<medie[a]<<endl;
  out.close();
  out.clear();

  out.open("rw_cubico_err.dat");
  for(int a =0; a < n_max_passi; a++) out<<err[a]<<endl;
  out.close();
  out.clear();

////RETICOLO CONTINUO 1////////
  //Riazzero tutto (anche se non penso ci sia bisogno)
  for (int b=0; b<n_max_passi; b++) {medie[b]=0; err[b]=0;}
  pos_fin=0; pos_fin2=0; media=0; media2=0;

  for (int i=0; i<n_max_passi; i++){
    media=0;
    media2=0;
    for(int w=0; w<N; w++){
      pos_fin=0;
      pos_fin2=0; //non serve
      for (int j=0; j<L; j++ ) {
        for(int k=0; k<(i+1); k++) pos.singolo_passo_cont(rnd);  //raggiunge posizione finale
        pos_fin += pos.dist2();
        pos.Setx(0.);
        pos.Sety(0.);
        pos.Setz(0.);
      }
      pos_fin = sqrt(pos_fin/double(L)); //risultato singolo esperimento
      pos_fin2 = pos_fin * pos_fin;

      media = (media*double(w) + pos_fin) /double (w+1);
      media2 = (media2*double(w) + pos_fin2) /double (w+1);  //medie degli esperimenti
    }

    medie[i] = media;
    err[i] = error ( media, media2, N-1);   //salvo nel vettore la media finale con i passi
  }


  out.open("rw_cont1.dat");
  for(int a =0; a < n_max_passi; a++) out<<medie[a]<<endl;
  out.close();
  out.clear();

  out.open("rw_cont1_err.dat");
  for(int a =0; a < n_max_passi; a++) out<<err[a]<<endl;
  out.close();
  out.clear();

 ////RETICOLO CONTINUO 2////////
  //Riazzero tutto (anche se non penso ci sia bisogno)
  for (int b=0; b<n_max_passi; b++) {medie[b]=0; err[b]=0;}
  pos_fin=0; pos_fin2=0; media=0; media2=0;

  for (int i=0; i<n_max_passi; i++){
    media=0;
    media2=0;
    for(int w=0; w<N; w++){
      pos_fin=0;
      pos_fin2=0; //inutile
      for (int j=0; j<L; j++ ) {
        for(int k=0; k<(i+1); k++) pos.singolo_passo_cont2(rnd);  //raggiunge posizione finale
        pos_fin += pos.dist2();
        pos.Setx(0.);
        pos.Sety(0.);
        pos.Setz(0.);
      }
      pos_fin = sqrt(pos_fin/double(L)); //risultato singolo esperimento
      pos_fin2 = pos_fin*pos_fin;

      media = (media*double(w) + pos_fin) /double (w+1);
      media2 = (media2*double(w) + pos_fin2) /double (w+1);  //medie degli esperimenti
    }

    medie[i] = media;
    err[i] = error ( media, media2, N-1);   //salvo nel vettore la media finale con i passi
  }


  out.open("rw_cont2.dat");
  for(int a =0; a < n_max_passi; a++) out<<medie[a]<<endl;
  out.close();
  out.clear();

  out.open("rw_cont2_err.dat");
  for(int a =0; a < n_max_passi; a++) out<<err[a]<<endl;
  out.close();
  out.clear();   
    

  return 0;

}

double error ( double media , double media_quadrati, int n) {
    if(n==0) return 0;
    else return sqrt((media_quadrati - media*media)/n);
}
