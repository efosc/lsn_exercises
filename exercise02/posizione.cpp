#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "posizione.h"

using namespace std;

posizione :: posizione(): x(0), y(0), z(0) {}
posizione :: posizione(double a , double b, double c): x(a), y(b), z(c){} 

posizione :: ~posizione(){}

void posizione :: singolo_passo( Random & rnd){
  int direzione = int( rnd.Rannyu() / double(1./6.)); //numero intero tra 0 e 5
  if(direzione == 0) x += 1;
  else if (direzione == 1) x += -1;
  else if (direzione == 2) y += 1;
  else if (direzione == 3) y += -1;
  else if (direzione == 4) z += 1;
  else if (direzione == 5) z += -1;
  return;
}

void posizione::singolo_passo_cont(Random & rnd){
	//genera tre candidati x ,y, z
	double cateto_orizz, forse_cateto_orizz;
	double forsex, forsey, forsez;
	int i;
	for ( i=0; i<1000; ++i){   //faccio al massimo 1000 tentativi
		forsex = -1 + rnd.Rannyu()*2;
		forsey = -1 + rnd.Rannyu()*2;
		forsez = -1 + rnd.Rannyu()*2;
	//controllo se sono nella sfera di raggio 1, se lo sono ho trovato quello che mi interessa ed esco
		if ( forsex*forsex + forsey*forsey + forsez*forsez <= 1 ) break;
	}
	//cerr << i << "   " << forsex <<"  " << forsey << "  " << forsez  << "  " << forsex*forsex + forsey*forsey + forsez*forsez<< endl;
	// se i=1000 sono rimanst* dove ero prima
	if (i==1000) cerr << "nella stessa posizione di prima"<<endl;
	//trasporto sulla sfera
	else{
		forse_cateto_orizz = sqrt(forsex*forsex + forsey*forsey);
		cateto_orizz = forse_cateto_orizz / sqrt(forsex*forsex + forsey*forsey + forsez*forsez); 
		z += forsez / sqrt(forsex*forsex + forsey*forsey + forsez*forsez);
		x += forsex * cateto_orizz / forse_cateto_orizz;
		y += forsey * cateto_orizz / forse_cateto_orizz;
	}

	/*cerr << x <<"  " << y << "   "<< z <<"  " << x*x + y*y + z*z << endl;
	cerr<<"forse_cateto_orizz:  "<<forse_cateto_orizz<< "   " << sqrt(forsex*forsex + forsey*forsey)<< endl;
	cerr<<"cateto_orizz:  "<<cateto_orizz<<endl;
	cerr<<"forsex:  "<<forsex<<endl;
	cerr<<"forsey:  "<<forsey<<endl;
	cerr<<"forsez:  "<<forsez<<endl;
	cerr <<"x:  "<<x<<endl;
	cerr <<"y:  "<<y<<endl;
	cerr <<"z:  "<<z<<endl;*/
	return;
}

void posizione::singolo_passo_cont2(Random & rnd){
	//genera tre candidati x ,y, z
	double forsex, forsey, forsez;
	double theta, phi, ro; //angoli che saranno ricavati da forsex, forsey, forsez
	int i;
	for ( i=0; i<1000; ++i){   //faccio al massimo 1000 tentativi
		forsex = -1 + rnd.Rannyu()*2;
		forsey = -1 + rnd.Rannyu()*2;
		forsez = -1 + rnd.Rannyu()*2;
	//controllo se sono nella sfera di raggio 1, se lo sono ho trovato quello che mi interessa ed esco
		if ( forsex*forsex + forsey*forsey + forsez*forsez <= 1 ) break;
	}
	//cerr << i << "   " << forsex <<"  " << forsey << "  " << forsez  << "  " << forsex*forsex + forsey*forsey + forsez*forsez<< endl;
	// se i=1000 sono rimanst* dove ero prima
	if (i==1000) cerr << "nella stessa posizione di prima"<<endl;
	//calcolo angoli e x,y,z veri
	else{
		ro = sqrt(forsex*forsex + forsey*forsey + forsez*forsez);
		theta = acos(forsez/ro);
		phi = atan2(forsey , forsex);
		x += sin(theta)*cos(phi);
		y += sin(theta)*sin(phi);
		z += cos(theta);
	}

	/*cerr << x <<"  " << y << "   "<< z <<"  " << x*x + y*y + z*z << endl;
	//cerr<<"forse_cateto_orizz:  "<<forse_cateto_orizz<< "   " << sqrt(forsex*forsex + forsey*forsey)<< endl;
	//cerr<<"cateto_orizz:  "<<cateto_orizz<<endl;
	cerr<<"forsex:  "<<forsex<<endl;
	cerr<<"forsey:  "<<forsey<<endl;
	cerr<<"forsez:  "<<forsez<<endl;
	cerr <<"x:  "<<x<<endl;
	cerr <<"y:  "<<y<<endl;
	cerr <<"z:  "<<z<<endl;*/
	return;
}


