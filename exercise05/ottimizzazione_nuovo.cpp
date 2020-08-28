#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include "random.h"
#include <string>
#include "lez05_functions.h"

using namespace std;

//double error ( double media , double media_quadrati, int n);
//double prob_s(double x, double y, double z);
//double prob_p(double x, double y, double z);
//int conta_accettati_gauss(double step, double x_old, double y_old, double z_old, bool orbital, Random& rnd);
//int conta_accettati_unif(double step, double x_old, double y_old, double z_old, bool orbital, Random& rnd);

int main (int argc, char *argv[]){

  Input(); //seed + legge i dati in input
	
  ofstream acc;
  int best=9999;
  double best_step;
  int n_acc;

  acc.open(file_name);
  for (int i=0 ; i<2000; i++){ //vari valore di step
    step = (double(i) + 1)*0.01;
    if (!T) { n_acc = conta_accettati_unif2(step, x_0, y_0, z_0);}
    else if (T) n_acc = conta_accettati_gauss(step, x_0, y_0, z_0);
    if (abs(best-int(ntimes/2)) >= abs(n_acc-int(ntimes/2))) {
      best=n_acc;
      best_step = step;
    }   
    acc << "step: " << step << " accepted steps: " <<n_acc << endl;
    if((i%100) == 0)cout << "step: " << step << " accepted steps: " <<n_acc << endl;
  }
 
  cout <<  "best step: " << best_step << "  accepted steps: " <<best << endl;

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
  ReadInput.open("input.dat");

  ReadInput >> x_0;
  ReadInput >> y_0;
  ReadInput >> z_0;
  ReadInput >> orbital;
  ReadInput >> T;
  ReadInput >> ntimes;

   ReadInput.close();

  if(!orbital) orbital_name="s";
  else if(orbital) orbital_name="p";

  if(!T) {
    transition_prob = "uniform";
  } else if (T){
    transition_prob = "gauss";
  }

/////nome del file in cui scrivere
  file_name = "Optimization/accepted_" + orbital_name + "_" + transition_prob + "_prova2.dat";

  cout << "start: (" << x_0 << ", "<< y_0 <<", "<< z_0 << ") " << endl;
  cerr << orbital <<"   " << T<<endl;
  //cerr << itime << endl;

  cout << orbital_name <<",  " << transition_prob << " probability for each step" << endl;
  cout << "writing into " << file_name << endl;

  return;

}

/////////////////PROBABILITY////////////////////////////////////////
double prob_s(double x, double y, double z){
  double r = sqrt(x*x + y*y + z*z);
  return exp(-2*r)/M_PI;
}

double prob_p(double x, double y, double z){
  double r = sqrt(x*x + y*y + z*z);
  return exp(-r)*(z)*(z)/(32*M_PI);
}

//////////CONTA_ACCETTATI///////////////////////////////////////////////
/////data una certa lunghezza dell step e la posizione di partenza della simulazione,
//fa una sim di 1000 passi mi dice quanti sono i passi accettati (ne voglio 500)
int conta_accettati_unif(double step, double x_old, double y_old, double z_old){

  int n_acc=0;
  double alpha, r, div, x_new, y_new, z_new;
  double cateto_orizz, forse_cateto_orizz;  //per estrarre uniformemente su una sfera
  double forsex, forsey, forsez;

  for (int itime=1; itime<=ntimes; itime++){

///////ESTRAGGO UN NUOVO PUNTO///////////////////
    int i;
    for ( i=0; i<1000; ++i){   //faccio al massimo 1000 tentativi
      forsex = (-1 + rnd.Rannyu()*2)*step;  //sfera di raggio step  -> cubo di lato 2*step
      forsey = (-1 + rnd.Rannyu()*2)*step;
      forsez = (-1 + rnd.Rannyu()*2)*step;
    //controllo se sono nella sfera di raggio step, se lo sono ho trovato quello che mi interessa ed esco
      if ( sqrt(forsex*forsex + forsey*forsey + forsez*forsez) <= step ) break;
    }
    // se i=1000 sono rimanst* dove ero prima
    if (i==1000) cerr << "nella stessa posizione di prima"<<endl;
    //trasporto sulla sfera per trovare la nuova posizione da proporre
    else{
      forse_cateto_orizz = sqrt(forsex*forsex + forsey*forsey);
      cateto_orizz = forse_cateto_orizz / sqrt(forsex*forsex + forsey*forsey + forsez*forsez); 
      z_new = z_old + forsez / sqrt(forsex*forsex + forsey*forsey + forsez*forsez);
      x_new = x_old + forsex * cateto_orizz / forse_cateto_orizz;
      y_new = y_old + forsey * cateto_orizz / forse_cateto_orizz;
    }

    //cerr << x_1 <<"   " << y_1 <<"  "<< z_1<< endl;

//////////VALUTO ALPHA//////////////////////////////////////////////////////
    if(!orbital) div = (prob_s(x_new , y_new, z_new)/prob_s(x_old , y_old, z_old));
    else if (orbital) div = (prob_p(x_new , y_new, z_new)/prob_p(x_old , y_old, z_old));
    alpha = min(1.,div);

    //accettato?
    r = rnd.Rannyu();
    if (r <= alpha) {
      x_old=x_new;
      y_old=y_new;
      z_old=z_new;
      n_acc++;
    }

  }

  return n_acc;

}

//////////////CONTA_ACCETATI_GAUSS//////////////////////////////////////////////////////
/////data una certa lunghezza dell step e la posizione di partenza della simulazione,
//fa una sim di 1000 passi mi dice quanti sono i passi accettati (ne voglio 500)
int conta_accettati_gauss(double step, double x_old, double y_old, double z_old){
  int n_acc=0;
  double alpha, r, div, x_new, y_new, z_new;

  for (int i=0; i<ntimes; i++){

///////ESTRAGGO UN NUOVO PUNTO/////////////////// 
/// con probabilità guassiana su ogni coordinata, con step = sigma e centro il punto vecchio
    x_new = rnd.Gauss( x_old, step);
    y_new = rnd.Gauss( y_old, step);
    z_new = rnd.Gauss( z_old, step);

    //cerr << x_new <<"   " << y_new <<"  "<< z_new<< endl;

//////////VALUTO ALPHA//////////////////////////////////////////////////////
    if(!orbital) div = (prob_s(x_new , y_new, z_new)/prob_s(x_old , y_old, z_old));
    else if (orbital) div = (prob_p(x_new , y_new, z_new)/prob_p(x_old , y_old, z_old));
    alpha = min(1.,div);

    //accettato?
    r = rnd.Rannyu();

    //cerr << r << "    " << alpha <<endl;
    if (r <= alpha) {
      x_old=x_new;
      y_old=y_new;
      z_old=z_new;
      n_acc++;
    }
  }

  return n_acc;

}

//////////CONTA_ACCETTATI_UNIF"///////////////////////////////////////////////
/////data una certa lunghezza dell step e la posizione di partenza della simulazione,
//fa una sim di 1000 passi mi dice quanti sono i passi accettati (ne voglio 500)
// questi sono estratti uniformemente in un un cubo, nel senso per ogni asse x, y, z
int conta_accettati_unif2(double step, double x_old, double y_old, double z_old){

  int n_acc=0;
  double alpha, r, div, x_new, y_new, z_new;

  for (int itime=1; itime<=ntimes; itime++){

///////ESTRAGGO UN NUOVO PUNTO//////////////////
      x_new = x_old + (-1 + rnd.Rannyu()*2)*step;  //cubo di lato 2*step
      y_new = y_old + (-1 + rnd.Rannyu()*2)*step;
      z_new = z_old + (-1 + rnd.Rannyu()*2)*step;

    //cerr << x_1 <<"   " << y_1 <<"  "<< z_1<< endl;

//////////VALUTO ALPHA//////////////////////////////////////////////////////
    if(!orbital) div = (prob_s(x_new , y_new, z_new)/prob_s(x_old , y_old, z_old));
    else if (orbital) div = (prob_p(x_new , y_new, z_new)/prob_p(x_old , y_old, z_old));
    alpha = min(1.,div);

    //accettato?
    r = rnd.Rannyu();
    if (r <= alpha) {
      x_old=x_new;
      y_old=y_new;
      z_old=z_new;
      n_acc++;
    }

  }

  return n_acc;

}
