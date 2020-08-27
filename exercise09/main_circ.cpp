#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <cmath>
#include <fstream>
#include "random.h"
#include "funzioni9.h"

using namespace std;

extern double positions[32][2];
extern int seed[4];
extern Random rnd;
extern int ngens, times_per_gen, maxtimes;
extern double sel_pow;
extern string volta;
extern bool circ;

int main(){

	Input();

	string filename="Outfiles/final_pos_" + volta + ".dat";
	int which1, which2, ntimes=0;
	double r, r1, r2, r3, r4;
	double cost_new, cost_old, cost_ave;
	ofstream outl("Outfiles/l2_values_" + volta + ".dat");
	ofstream firstl("Outfiles/l2_first_" + volta + ".dat");

	//genero la prima generazione
	Generation routes(500);
	cout << "starting generation ... " << endl;
	routes.Check_generation();
	//routes.Print_generation();

	cout << "First cities vector at the very beginning: " << endl;
	routes.Get_first().Print_list();
	cout << "Its length is: " << routes.Get_first().L2() << endl;
	cost_old = routes.Get_first().L2();

	firstl << routes.Get_first().L2() << endl;

	//per un tot di volte faccio crossover + mutazioni:
	// -> ordino
	// -> pesco 2 su cui fare crossover 
	// -> aggiungo eventualmente mutazioni
	// -> rifaccio per un tot di volte e ho una nuova generazione

	for (int igen =0; igen<ngens; ++igen){
		if(igen%10==0) cerr << "Generation number: " << igen <<endl;
		for (int itime=0; itime<times_per_gen ; ++itime ){
			//cerr << "volta numero: " << itime << endl;
			routes.my_sort(); //ordina dovrebbe
			r = rnd.Rannyu();
			//cerr << "estraggo un numero " << r << " " ;
			which1 = routes.selector(r, sel_pow);
			//cerr << " e scelgo il cromosoma numero: " << which1 << endl;
			r = rnd.Rannyu();
			//cerr << "estraggo un numero " << r << " " ;
			which2 = routes.selector(r, sel_pow);
			//cerr << " e scelgo il cromosoma numero: " << which2 << endl;			

			//crossover?
			r = rnd.Rannyu();
			if(r<= 0.5) routes.crossover(which1, which2);

			routes.Check_generation();

			//mutations? 
			r1 = rnd.Rannyu();
			r2 = rnd.Rannyu();
			r3 = rnd.Rannyu();
			r4 = rnd.Rannyu();
			routes.single_mutation(routes.Get_nchr() -1 , r1,r2,r3,r4);
			routes.single_mutation(routes.Get_nchr() -2, r1,r2,r3,r4);

			//altre mutazioni

			//int which3 = (rand() % routes.Get_nchr());
			//routes.single_mutation(which3 , 0.01,1.,1.,1.);

		}
		cost_ave = routes.Get_L(); //quella mediata sulla prima metà
		outl << cost_ave << endl;
		cost_new = routes.Get_first().L2(); //solo del primo
		firstl << cost_new << endl;

		cout << "First: " ;
		routes.Get_first().Print_list();

		if(cost_new >= cost_old) ntimes ++;
		else {cost_old = cost_new; ntimes=0;} //mi tengo quello migliore per confrontare con quelli nuovi
		if(ntimes == maxtimes) {
			cout << "Cost didn't improve for " << maxtimes << " generations" <<endl;
			cout << "Stopping at generation number: " << igen << endl;
			break;   //se non miglioro per un tot esco
		}


		routes.Check_generation();
	}

	outl.close();
	firstl.close();

	// dopo n generazioni guardo quello più breve

	routes.my_sort();
	cout << " This may be the shortest route: ";
	routes.Get_first().Print_list();
	cout << "Route length: " << routes.Get_first().L2() <<endl;

	// scrivo il primo vettore (più breve dovrebbe essere?) nel file
	routes.Get_first().Print_positions(filename);


	return 0;
}
