#include <vector>
#include <string>
#include <iostream>
#include <fstream>  //ifstream
#include <cassert>  //assert
#include <cmath>
#include <cstdlib>  //rand
#include <algorithm>
# include <iomanip> //setw
#include "random.h"
#include "funzioni10.h"

using namespace std;

extern double positions[32][2];
extern int seed[4];
extern Random rnd;

///////////COSE GENERALI///////////
double beta, initial_beta, beta_step1, beta_step2, beta_step3;  //temperatura iniziale 
int n_beta_steps1 , nsteps1, n_beta_steps2 , nsteps2,  n_beta_steps3 , nsteps3;  //di quanto variare la temperatura ogni volta e quanti step per ogni temperatura
int accepted=0, attempted=0;
bool circ; // 1->circ, 0->square
string filename;
Chromosome tour; //ho un percorso possibile di 32 città
Chromosome proposed(tour); //all'inizio uguali

///////////////////////////////////////////////////
int main(){

	const double space = 12.;

	Setting();

	ofstream print( filename + "_temp.dat");

	///proposed evolve per n1, n2 ecc volte
	for(int t=0; t<n_beta_steps1; t++){
		beta = initial_beta + t*beta_step1;
		for(int j=0; j<nsteps1; ++j){
		//muovo per nstep
		Move(beta);
		}
		print << beta << setw(space) << 1.0/beta << setw(space) << accepted << setw(space) << attempted << setw(space) << tour.L2() << endl;
		cerr << beta << setw(space) << 1.0/beta << setw(space) << accepted << setw(space) << attempted << setw(space) << tour.L2() << endl;

		attempted=0;
		accepted = 0;
	}
	for(int t=0; t<n_beta_steps2; t++){
		beta = beta + beta_step2;
		for(int j=0; j<nsteps2; ++j){
		//muovo per nstep
		Move(beta);
		}
		print << beta << setw(space) << 1.0/beta << setw(space) << accepted << setw(space) << attempted << setw(space) << tour.L2() << endl;
		cerr << beta << setw(space) << 1.0/beta << setw(space) << accepted << setw(space) << attempted << setw(space) << tour.L2() << endl;

		attempted=0;
		accepted = 0;
	}
	for(int t=0; t<n_beta_steps3; t++){
		beta = beta + beta_step3;
		for(int j=0; j<nsteps3; ++j){
		//muovo per nstep
		Move(beta);
		}
		print << beta << setw(space) << 1.0/beta << setw(space) << accepted << setw(space) << attempted << setw(space) << tour.L2() << endl;
		cerr << beta << setw(space) << 1.0/beta << setw(space) << accepted << setw(space) << attempted << setw(space) << tour.L2() << endl;

		attempted=0;
		accepted = 0;
	}

	cerr << "Final: " << endl;
	tour.Print_list();
	cerr << "L2 = " << tour.L2() << endl;

	string ciao = filename + "_final_pos.dat";
	tour.Print_positions(ciao);


return 0;

}

/*void Set_annealing_schedule(){
	annealing_schedule[0][0] = 0.3; //beta1
	annealing_schedule[0][1] = 250.;//n1
	annealing_schedule[1][0] = 0.35; //beta1
	annealing_schedule[1][1] = 250.;//n1
	annealing_schedule[2][0] = 0.4; //beta1
	annealing_schedule[2][1] = 250.;//n1
}*/
/////////////

void Setting(){
	ifstream ReadInput("input2.dat");
	accepted = 0;
	attempted = 0;
///printing initial tours
	tour.Print_list();
	proposed.Print_list();
//metto a caso su una circonferenza/quadrato e inizializzo seme ecc
	ReadInput >> circ; 
	ReadInput >> filename; //outfile name
	if(circ){
		Initial_Circ();
		filename = "circ_" + filename;
		cerr << "Cities placed randomly on a circumference (r=1)" << endl;
	} else if(!circ){
		Initial_Square();
		filename = "square_" + filename;
		cerr << "Cities placed randomly within a square (L=2)" << endl;
	}

//setting annealing schedule (beta/steps) reading from input file
	ReadInput >> initial_beta;    //beta iniziale da cui part
	ReadInput >> beta_step1;       //di quanto aumentare beta ogni volta
	ReadInput >> n_beta_steps1;        //quanti step di temperatura
	ReadInput >> nsteps1;         //quanti passi con proposte mutazioni ogni volta
	ReadInput >> beta_step2;       //di quanto aumentare beta ogni volta
	ReadInput >> n_beta_steps2;        //quanti step di temperatura
	ReadInput >> nsteps2;         //quanti passi con proposte mutazioni ogni volta
	ReadInput >> beta_step3;       //di quanto aumentare beta ogni volta
	ReadInput >> n_beta_steps3;        //quanti step di temperatura
	ReadInput >> nsteps3;         //quanti passi con proposte mutazioni ogni volta
	cerr << "Annealing schedule set" << endl;

	ReadInput.close();	

	return;

//

}

////////////////////////////////////////////////////
/*void Set_annealing_schedule()
{
	ifstream ReadInput("input2.dat");

	ReadInput >> initial_beta;    //beta iniziale da cui part
	ReadInput >> beta_step1;       //di quanto aumentare beta ogni volta
	ReadInput >> n_beta_steps1;        //quanti step di temperatura
	ReadInput >> nsteps1;         //quanti passi con proposte mutazioni ogni volta
	ReadInput >> beta_step2;       //di quanto aumentare beta ogni volta
	ReadInput >> n_beta_steps2;        //quanti step di temperatura
	ReadInput >> nsteps2;         //quanti passi con proposte mutazioni ogni volta

	ReadInput.close();
}*/

/////////////MOVE_METROPOLIS///////////////////
void Move(double beta){
	double A, r, L_old, L_new;
	int i , m, which;

	which = rand() % 6;
	if(which == 0){   //swap i-th city and (i+1)th city
			i = (rand() % (proposed.Get_ncities()-1)) +1;
			proposed.pair_permutation(i);
	} else if ((which == 1) || (which == 2)){     //shift all cities by a random number of positions
		i = (rand() % (proposed.Get_ncities()-1));  //quante_pos
		proposed.shift(i); 
	} else if ((which == 3) || (which == 4)){   //swap two cities random
		i = (rand() % (proposed.Get_ncities()-1)) + 1;
		m = (rand() % (proposed.Get_ncities()-1)) + 1;
		proposed.pair_inversion(i,m);
	} else if ((which == 5) || (which == 6)){    //swap group of cities
		i = (rand() % (int(proposed.Get_ncities()/2)));
		proposed.permutation(i);
	}
	proposed.check(); //controlla che non si ripetano città e che la prima sia 1
	//accettato ?
	L_old = tour.L2();
	L_new = proposed.L2();
	//cerr << "L_old = " << L_old << "L_new = " << L_new << endl;
	A = min(1. , exp(-(beta) * (L_new - L_old )));
	r = rnd.Rannyu();
	if(r <= A) { //accepted -> tour = proposed
		for (int k=0; k<tour.Get_ncities(); ++k) tour.Set_element(k, proposed.Get_list()[k]);
		L_old = L_new;
		accepted ++;
	} else { // proposed = tour
		for (int k=0; k<proposed.Get_ncities(); ++k) proposed.Set_element(k, tour.Get_list()[k]);
	}

	attempted ++;
}







