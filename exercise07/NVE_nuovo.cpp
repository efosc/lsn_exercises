#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "Mol_functions.h"

using namespace std;
int istep;
int iblk;

int main()
{ 
  Input(); //Inizialization (e riscalare le velocità)

  //for(int i=0; i<npart ; ++i) cerr << "POS_DOPO_INPUT: " << x[i] <<"  " << y[i] << "  " << z[i] << endl;
  //for(int i=0; i<npart ; ++i) cerr << "POSOLD_DOPO_INPUT: " << xold[i] <<"  " << yold[i] << "  " << zold[i] << endl;


  for( iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    //for(int i=0; i<npart ; ++i) cerr << iblk <<"  " << "POS_1: " << x[i] <<"  " << y[i] << "  " << z[i] << endl;
    //for(int i=0; i<npart ; ++i) cerr <<  iblk <<"  "<< "POSOLD_1: " << xold[i] <<"  " << yold[i] << "  " << zold[i] << endl;
    for( istep=1; istep <= blk_length; ++istep)
    {
      Move();
      //for(int i=0; i<npart ; ++i) cerr << iblk <<"  " << istep << "POS_2: " << x[i] <<"  " << y[i] << "  " << z[i] << endl;
      //for(int i=0; i<npart ; ++i) cerr <<  iblk <<"  " << istep << "POSOLD_2: " << xold[i] <<"  " << yold[i] << "  " << zold[i] << endl;
      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
        if(istep%10 == 0){
          Measure();     //Properties measurement
         // Instant();
          Accumulate();  //Update block averages
      } 
    }
    Averages(iblk);   //Print results for current block
  }
  ConfOld();
  ConfFinal(); //Write final configuration and old configuration (to restart)

  return 0;
}

/////INPUT/////////////////

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  //double ep, ek, pr, et, vir;  //cosa sono queste variabili?

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator

//legge i dati da input.dat ( temperatura, numero di particelle, densità, rcut,
//timestep, numero di step, ogni quanto calcolare, se prendere la posizione vecchia da old o meno)
//calcola volume del box, lato del box
  
  ReadInput.open("input_mio.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;
  cout << "Target temperature = " << temp << endl<<endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> blk_length;
  ReadInput >> nblk;
  ReadInput >> iprint;
  ReadInput >> old;
  ReadInput >> equilibration;
  ReadInput >> filename; //per cambiare nome al file se lo faccio andare più volte
  ReadInput >> dirname; //per cartelle diverse per fasi diverse

  nstep = nblk*blk_length;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "R_cut = " << rcut << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  if(old){cout << "Old option is on" << endl;} else if (!old){cout << "Old option is off" << endl ;}
  if(equilibration){cout << "Equilibration option is on" << endl;}else if (!equilibration){cout << "Equilibration option is off" << endl << endl;}
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  iw = 4; //pressure
 
  n_props = 5; //Number of observables

  //measurement of g(r)
  igofr = 5;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

// a questo punto ci sono due possibilità: 
// 1. se old è attivo legge old.0 e old.final -> velocità iniziali ricavate da queste due configurazioni
// in realtà Move() ogni volta riscrive le velocità
// ma queste non sono importanti per calcolare la posizione successiva (usa solo grad f) -> 
// -> potrei anche non settare delle velocità iniziali

  if(!old){
  	cout << "Old option: " << old << "-> read initial configuration from file config.0" << endl;
  	Read_config0();
  	cout << "Casual initial velocities" << endl << endl;
  	Casual_vel_estimate_old();
  } else if (old){
  	cout << "Old option: " << old << "-> read both old.0 and old.final" << endl;
  	Read_old();
  	//cout << "Evaluating initial velocities" << endl << endl;
  	//estimate_vel();
  	if(equilibration) Rescale();
  }

  return;

}

/////////////CONFIGURATIONS////////////////////

void Read_config0(){
	if(old) cout << "Read_config0() -> Old option should be false" << endl;
	ifstream ReadConf("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  return;
}
//-------------------------------------------//
void Read_old(){
	if(!old) cout << "Read_old -> old option should be true" << endl;

	ifstream ReadConf;

  ReadConf.open("old.0");
  if(!ReadConf.is_open()) cout << "Problems in opening old.0" <<endl;
  for (int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i]; //prendo dal file
      xold[i] = xold[i] * box; //nel file solo scritti in unità del lato del box
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
  ReadConf.close();

	ReadConf.open("old.final");
  if(!ReadConf.is_open()) cout << "Problems in opening old.final" <<endl;
  for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i]; //prendo dal file
      x[i] = x[i] * box; //nel file solo scritti in unità del lato del box
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
  ReadConf.close();

  return;
}

//////////////INITIAL VELOCITIES////////////////////////////

///NON OLD ----------------------------------------------//
void Casual_vel_estimate_old(){
	if(old) cout << "Casual_vel_estimate_old() -> Old option should be false" << endl;
	double sumv[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<npart; ++i){
    vx[i] = rand()/double(RAND_MAX) - 0.5;
    vy[i] = rand()/double(RAND_MAX) - 0.5;
    vz[i] = rand()/double(RAND_MAX) - 0.5;

    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }
  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
  double sumv2 = 0.0, fs;
  for (int i=0; i<npart; ++i){
    vx[i] = vx[i] - sumv[0];
    vy[i] = vy[i] - sumv[1];
    vz[i] = vz[i] - sumv[2];

    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  sumv2 /= (double)npart;
	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
  for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;
    xold[i] = Pbc(x[i] - vx[i] * delta);
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);
	}

	return;
}


///////////////////////EQUILIBRATION-> RESCALE VELOCITIES///////////////////////////

void Rescale(void){ //Rescale velocities

	if(!equilibration) cout <<"PROBLEM: euilibration option should be true" << endl;
  
  if(equilibration){ // just to be sure
    double t=0., scaling_factor; //kinetic energy and scaling factor
    double initial_temp;
    cout << "Since equilibration option is on -> " <<
    "rescaling initial speeds to match the desired temperature" << endl << endl;

  //compute one step of the Verlet algorithm
    Move();
  //the configuration from the config file is now the old configuration (r(t))
  //and r(t + dt) is the new configuration

  //kinetic energy
    for(int i=0; i<npart; ++i) t+=vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    t = 0.5*t;
    //temperature & rescale factor
    initial_temp = (2./3.)*t/double(npart);
    scaling_factor = temp / initial_temp;  //la temperatura desiderata è quella scritta nell'input file

  //rescaling and estimating a novel r(t) configuration (r(t)=xold, r(t+dt)=x)
    for (int i=0; i<npart; ++i){
      vx[i] *= sqrt(scaling_factor);
      vy[i] *= sqrt(scaling_factor);
      vz[i] *= sqrt(scaling_factor);

      xold[i] = Pbc(x[i] - (vx[i])*delta);
      yold[i] = Pbc(y[i] - (vy[i])*delta);
      zold[i] = Pbc(z[i] - (vz[i])*delta);
    }
  }

   return;
}


/////////////MOVE////////////////////////////////////////////////

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}


/////////////FORCE////////////////////////////////////////////////

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

//////////////////MEASURE//////////////////////////////

void Measure()
{
	int bin;
  double v, t, w, vij, wij;
  double dx, dy, dz, dr;

  v = 0.0; //reset observables
  t = 0.0;
  w = 0.0;

  //reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)
     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

	//update of the histogram of g(r)
     /*for (int k=0; k<nbins; ++k){
      if( (k*bin_size < dr) && ( (k+1)*bin_size < dr) ) walker[k+igofr]+=2; //aggiungo 2 perché è simmetrico
     }*/
     bin=int(dr/bin_size); // bin in cui si trova la paricella
     walker[igofr+bin]+=2; //+ 2 perché è simmetrico

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       wij = 48.0/pow(dr,12) - 24.0/pow(dr,6);
//Potential energy
       v += vij;
       w += wij;
     }
    }          
  }

  v /= double(npart); //energia potenziale per particella

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  t /= double(npart);

  walker[iv] = v; 
  walker[ik] = t;
  walker[ie] = v+t;  
  walker[it] = (2.0/3.0)*t;
  walker[iw] = rho*(2.0/3.0)*t + w/(3.0*vol);
}
//////////////////////////////////////////////////////////////////

/////PRINT INSTANT VALUES ///////////////////////////////////////////////////
void Instant(void){

	ofstream print(dirname + "/instant" + filename + ".dat", ios::app);
	double space = 12.;
  double vol_shell, r ;

	/// epot << ekin << etot << temp << press

	print <<  walker[iv] << setw(space) << walker[ik] << setw(space) << walker[ie] << setw(space)
  << walker[it] << setw(space) << walker[iw] << endl;

  //g(r)
  for(int k=0; k<nbins ; k++){
      r = k*bin_size; //distanza che considero ogni volta
      //"normalize with the volume of the shell of radius dr"
      vol_shell = (4.0/3.0)* pi *(pow((k+1)*bin_size,3) - pow(k*bin_size,3));
      stima_g = walker[k+igofr]/rho/(double)npart/vol_shell; //stima in quel bin
      print << k << setw(space) << r << setw(space) << stima_g << endl;
  }

  print << endl<< endl;

}



//////////////ACCUMULATE//////////////////////////////////////////

void Accumulate(){  //update averages within each block

  for (int i=0; i<m_props; i++){
    blk_av[i] +=  walker[i];
  }

  blk_norm = blk_norm + 1.0; ///conta il numero di misure, non di steps

} 


//////////////RESET/////////////////////////////////////////////

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
}


////////////////////////AVERAGES//////////////////////////////////

void Averages(int iblk) //Print results for current block
// e aggiorna in glob_av la media aggiungendo il risultato di ogni blocco
{
    
   //ofstream Epot, Press;
   //ofstream Ekin, Etot, Temp;
   ofstream Gave, Gofr; 
   const int wd=12;
   double r, vol_shell;

    
    cout << "Block number " << iblk << endl;
    cout << "Step number " << iblk*blk_length << endl;

    //Press.open( dirname + "/press_blocks"+filename+".dat",ios::app);
    stima_press = blk_av[iw]/(double)blk_norm; //Potential energy
    /*glob_av[iw]  = (glob_av[iw]*(double)(iblk-1) + stima_press)/(double)iblk;
    glob_av2[iw] = (glob_av2[iw]*(double)(iblk-1) + stima_press*stima_press)/(double)iblk;
    err_press = Error(glob_av[iw],glob_av2[iw],iblk);
    Press << setw(wd) << iblk <<  setw(wd) << stima_press << setw(wd) << glob_av[iw] << setw(wd) << err_press << endl;
    Press.close();*/
    cout << "Pressure in "<< iblk << ": " << stima_press << endl;

    
    //Epot.open(dirname + "/epot_blocks"+filename+".dat",ios::app);
    stima_pot = blk_av[iv]/(double)blk_norm; //Potential energy
    /*glob_av[iv]  = (glob_av[iv]*(double)(iblk-1) + stima_pot)/(double)iblk;
    glob_av2[iv] = (glob_av2[iv]*(double)(iblk-1) + stima_pot*stima_pot)/(double)iblk;
    err_pot = Error(glob_av[iv],glob_av2[iv],iblk);
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv] << setw(wd) << err_pot << endl;
    Epot.close();*/
    cout << "Potential energy per particle in "<< iblk << ": " << stima_pot << endl;

    /*Ekin.open(dirname + "/ekin_blocks"+filename+".dat",ios::app);
    stima_ekin = blk_av[ik]/(double)blk_norm; //Kinetic energy
    glob_av[ik]  = (glob_av[ik]*(double)(iblk-1) + stima_ekin)/(double)iblk;
    glob_av2[ik] = (glob_av2[ik]*(double)(iblk-1) + stima_ekin*stima_ekin)/(double)iblk;
    err_ekin = Error(glob_av[ik],glob_av2[ik],iblk);
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_ekin << setw(wd) << glob_av[ik] << setw(wd) << err_ekin << endl;
    Ekin.close();
    cout << "Kinetic energy per particle in "<< iblk << ": " << stima_ekin << endl;

    Etot.open(dirname + "/etot_blocks" + filename+ ".dat",ios::app);
    stima_etot = blk_av[ie]/(double)blk_norm; //Total energy
    glob_av[ie]  = (glob_av[ie]*(double)(iblk-1) + stima_etot)/(double)iblk;
    glob_av2[ie] = (glob_av2[ie]*(double)(iblk-1) + stima_etot*stima_etot)/(double)iblk;
    err_etot = Error(glob_av[ie],glob_av2[ie],iblk);
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie] << setw(wd) << err_etot << endl;
    Etot.close();
    cout << "Total energy per particle in "<< iblk << ": " << stima_etot << endl;

    Temp.open(dirname + "/temp_blocks" + filename + ".dat",ios::app);
    stima_temp = blk_av[it]/(double)blk_norm; //Temperature
    //cerr << stima_temp << "   " << glob_av[it] << "    " << glob_av[it]*(double)(iblk-1) << endl;
    glob_av[it]  = (glob_av[it]*(double)(iblk-1) + stima_temp)/(double)iblk;
    glob_av2[it] = (glob_av2[it]*(double)(iblk-1) + stima_temp*stima_temp)/(double)iblk;
    //cerr << stima_temp << "   " << glob_av[it] << "    " << glob_av[it]*(double)(iblk-1) << endl;
    err_temp = Error(glob_av[it],glob_av2[it],iblk);
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it] << setw(wd) << err_temp << endl;
    Temp.close();
    cout << "Temperature in "<< iblk << ": " << stima_temp << endl;*/

    //g(r) 
    Gofr.open(dirname + "/max_blocks" + filename + ".dat",ios::app);
    //Gave.open(dirname + "/gave_blocks" + filename + ".dat",ios::app);

    for(int k=0; k<nbins ; k++){
      r = k*bin_size; //distanza che considero ogni volta
      //"normalize with the volume of the shell of radius dr"
      vol_shell = (4.0/3.0)* pi *(pow((k+1)*bin_size,3) - pow(k*bin_size,3));
      stima_g = blk_av[k+igofr]/blk_norm/rho/(double)npart/vol_shell; //stima in quel bin
      //cerr << glob_av[k+igofr] <<endl;
      glob_av[k+igofr] =  (glob_av[k+igofr]*(double)(iblk-1) + stima_g)/(double)iblk;
      glob_av2[k+igofr] = (glob_av2[k+igofr]*(double)(iblk-1) + stima_g*stima_g)/(double)iblk;
      err_gdir=Error(glob_av[k+igofr],glob_av2[k+igofr],iblk);
      if(k==46) Gofr << setw(wd) << iblk << setw(wd) << k << setw(wd) << r << setw(wd) << stima_g << setw(wd) << glob_av[k+igofr] << setw(wd) << err_gdir << endl;
      // blocco per blocco

    //g(r) finale
    /*if(iblk == nblk){
        Gave << setw(wd) << iblk << setw(wd) << k << setw(wd) << r << setw(wd) << stima_g << setw(wd) << glob_av[k+igofr] << setw(wd) << err_gdir << endl;
      }*/
  	}

    cout << "----------------------------" << endl << endl;
}


//////////////////CONFFINAL///////////////////////////
void ConfOld(void){
  ofstream WriteOld0("old.0");
  ofstream WriteOldF("old.final");
  
  cout << endl <<"Print final configuration to old.final to restart simulation" << endl;
  cout << endl <<"Print (n-1)th configuration to old.0 to restart simulation" << endl;
  
  for (int i=0; i<npart; ++i) WriteOld0 << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  for (int i=0; i<npart; ++i) WriteOldF << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  
  WriteOld0.close();
  WriteOldF.close();
  
  return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

/////////////////PERIODIC BOUNDARY CONDITIONS/////////////////////////////////////////

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);  //rint arrotonda al valore intero più vicino
}

/////
/////////////////////////////////////////////////////
/*double Error(double sum, double sum2, int iblk)
{
    if(iblk==0) {cout << "Problems with iblk! (blk=0)" << endl; return 0.0;}
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}*/


double Error(double ave, double ave2,int iblk){ 
    if(iblk==0) {cout << "Problems with iblk! (blk=0)" << endl; return 0.0;}
    if(iblk==1) return 0.0;
    else return sqrt((ave2-ave*ave)/(double)(iblk-1));
}



