#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include "random.h"
#include <string>
#include "lez05_functions.h"

using namespace std;

int main (int argc, char *argv[]){

  Input(); //seed + legge i dati in input
  n_acc =0;  //numero di passi accettati
  ofstream writep(file_name);

  ///scrivo posizione iniziale nel file/////////////////////////////
  writep << "0 " << x_0 <<"  " << y_0 <<"  " << z_0 <<"  " 
         << sqrt(x_0*x_0 + y_0*y_0 + z_0*z_0) << endl;

  ///muovo e scrivo nel file la posizione
  for(int i=1; i<=ntimes; i++){
    if(!T) unif_step(); 
    else if (T) gauss_step();
  ///scrivo la posizione attuale (se non viene accettato riscrivo quella vecchia)
    writep << i << "  " << x_0 <<"  " << y_0 <<"  " << z_0<<"  " <<
            sqrt(x_0*x_0 + y_0*y_0 + z_0*z_0) << endl;
    if(i%100==0) cout << "step: " << i << " current position: (" << x_0 <<", " << y_0 <<", " << z_0<<")" <<endl;
  }

  writep.close();
  cout << "final position (" << x_0 <<", " << y_0 <<", " << z_0<<")" <<endl;

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
  ReadInput.open("input.dat");

  ReadInput >> step;
  ReadInput >> x_0;
  ReadInput >> y_0;
  ReadInput >> z_0;
  ReadInput >> orbital;
  ReadInput >> T;
  ReadInput >> ntimes;
  ReadInput >> nstep; //ma per ora non me ne faccio nulla
  ReadInput >> volta;  //per cambiare nome al file se faccio riandare il codice

  cout << nstep << endl;
  cout << volta << endl;

   ReadInput.close();

  if(!orbital) orbital_name="s";
  else if(orbital) orbital_name="p";

  if(!T) {
    transition_prob = "uniform";
  } else if (T){
    transition_prob = "gauss";
  }

/////nome del file in cui scrivere
  file_name = "Optimization/graphic_" + orbital_name + "_" + transition_prob + "_" + volta + ".dat";

/////messaggi sul terminale
  cout << "start: (" << x_0 << ", "<< y_0 <<", "<< z_0 << ") " << endl;
  cout << orbital <<"   " << T<<endl;
  cout << orbital_name <<",  " << transition_prob << " probability for each step" << endl;
  cout << "writing into " << file_name << endl;
  cout << "number of simulation steps: " << ntimes << endl;

  return;

}


///////////////////PASSO CON PROBABILITA' UNIFORME///////////////////////////////////////
void unif_step(){
  double alpha, r, div, x_new, y_new, z_new;

 ///////ESTRAGGO UN NUOVO PUNTO//////////////////
      x_new = x_0 + (-1 + rnd.Rannyu()*2)*step;  //cubo di lato 2*step
      y_new = y_0 + (-1 + rnd.Rannyu()*2)*step;
      z_new = z_0 + (-1 + rnd.Rannyu()*2)*step;

    //cerr << x_1 <<"   " << y_1 <<"  "<< z_1<< endl;

//////////VALUTO ALPHA//////////////////////////////////////////////////////
    if(!orbital) div = (prob_s(x_new , y_new, z_new)/prob_s(x_0 , y_0, z_0));
    else if (orbital) div = (prob_p(x_new , y_new, z_new)/prob_p(x_0 , y_0, z_0));
    alpha = min(1.,div);

    //accettato?
    r = rnd.Rannyu();
    if (r <= alpha) {
      x_0=x_new;
      y_0=y_new;
      z_0=z_new;
      n_acc++;
    }

  return;
}


//////////////PASSO CON PROBABILITA' GAUSSIANA//////////////////////////////////
void gauss_step(){
  double alpha, r, div, x_new, y_new, z_new;

 ///////ESTRAGGO UN NUOVO PUNTO//////////////////
/// con probabilità guassiana su ogni coordinata, con step = sigma e centro il punto vecchio
    x_new = rnd.Gauss( x_0, step);
    y_new = rnd.Gauss( y_0, step);
    z_new = rnd.Gauss( z_0, step);

    //cerr << x_1 <<"   " << y_1 <<"  "<< z_1<< endl;

//////////VALUTO ALPHA//////////////////////////////////////////////////////
    if(!orbital) div = (prob_s(x_new , y_new, z_new)/prob_s(x_0 , y_0, z_0));
    else if (orbital) div = (prob_p(x_new , y_new, z_new)/prob_p(x_0 , y_0, z_0));
    alpha = min(1.,div);

    //accettato?
    r = rnd.Rannyu();
    if (r <= alpha) {
      x_0=x_new;
      y_0=y_new;
      z_0=z_new;
      n_acc++;
    }

  return;
}

double prob_s(double x, double y, double z){
  double r = sqrt(x*x + y*y + z*z);
  return exp(-2*r)/M_PI;
}

double prob_p(double x, double y, double z){
  double r = sqrt(x*x + y*y + z*z);
  return exp(-r)*(z)*(z)/(32*M_PI);
}

