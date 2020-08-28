#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "random.h"
#include <string>
#include "lez05_functions.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization

  ///faccio andare 2500 volte per equilibrare
  for (int equi=0; equi<2500; equi++){
    if(!T) unif_step(); 
    else if (T) gauss_step();  //move()
  }

  n_acc = 0;
  //ofstream writepos(file_name + "inst.dat");

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      if(!T) unif_step(); 
      else if (T) gauss_step();  //move()
      //writepos << setw(12) << x_0 << setw(12) << y_0 << setw(12) << z_0 << sqrt(x_0*x_0 + y_0*y_0 + z_0*z_0);
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }

  //writepos.close();

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
  ReadInput >> ntimes;  //step totali
  ReadInput >> nstep; //step per blocco
  ReadInput >> volta;

   ReadInput.close();

   nblk = int(ntimes/nstep); //numero di blocchi

  if(!orbital) orbital_name="s";
  else if(orbital) orbital_name="p";

  if(!T) {
    transition_prob = "uniform";
  } else if (T){
    transition_prob = "gauss";
  }

/////nome del file in cui scrivere
  file_name = "Blocking/averages" + orbital_name + "_" + transition_prob + "_" + volta + ".dat"; //!! aggiungere pezzi .dat!!!

/////messaggi sul terminale
  cout << "-----------------------------------------------------" << endl;
  cout << "start: (" << x_0 << ", "<< y_0 <<", "<< z_0 << ") " << endl;
  cout << orbital <<"   " << T<<endl;
  cout << orbital_name <<",  " << transition_prob << " probability for each step" << endl;
  cout << "step length: " << step << endl;
  cout << "writing into " << file_name << endl;
  cout << "number of simulation steps: " << ntimes << endl;
  cout << "-----------------------------------------------------" << endl;

  return;

}



/////////////RESET///////////////

void Reset(int iblk){ //Reset block averages
   
  if(iblk == 1)
   {
    glob_av = 0;
    glob_av2 = 0;
   }
  blk_av = 0;
  n_acc = 0;
}

//////////////////PASSO CON PROBABILITA' UNIFORME///////////////////////////////////////
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

/////////PROBABILITY///////////////////////////
double prob_s(double x, double y, double z){
  double r = sqrt(x*x + y*y + z*z);
  return exp(-2*r)/M_PI;
}

double prob_p(double x, double y, double z){
  double r = sqrt(x*x + y*y + z*z);
  return exp(-r)*(z)*(z)/(32*M_PI);
}
///////////////////////////////////////////////////////


//////////////ACCUMULATE/////////////

void Accumulate(void){ //Update block averages  
//calcolo r in quel momento e aggiungo alla media del blocco 
  dist = sqrt(x_0*x_0 + y_0*y_0 + z_0*z_0);
  blk_av = blk_av + dist;
}

////////////////AVERAGES////////////////////
void Averages(int iblk)  //Print results for current block
{
   ofstream Dist;
   const int wd=12;
    
    if((iblk%10)==0){
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << double(n_acc)/double(nstep) << endl;
	}
    
    Dist.open(file_name,ios::app);
    stima_dist = blk_av/nstep; //Energy
    glob_av += stima_dist;
    glob_av2 += stima_dist*stima_dist;
    err_dist=Error(glob_av,glob_av2,iblk);
    Dist << setw(wd) << iblk << setw(wd) << n_acc << setw(wd) << stima_dist << setw(wd) << glob_av/(double)iblk << setw(wd) << err_dist << endl;
    Dist.close();

    if((iblk%10)==0){
    cout << "<r> in this block: " << stima_dist << endl;
    cout << "----------------------------" << endl << endl;
	}
}

/////////////ERRORS/////////////////////////////////
double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}



