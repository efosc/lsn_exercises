#include <vector>
#include <string>
#include <iostream>
#include <fstream>  //ifstream
#include <cassert>  //assert
#include <cmath>
#include <cstdlib>  //rand
#include <algorithm>
#include "random.h"
#include "funzioni9.h"

using namespace std;
/*Chromosome::Chromosome(int m_cities): n_cities(m_cities) {
		for (int i=1; i<= n_cities; ++i) cities_vec.push_back(i);
	}*/
double positions[32][2];
int seed[4];
Random rnd;
int ngens, times_per_gen, maxtimes;
double sel_pow;
string volta;
bool circ;

Chromosome::Chromosome(): n_cities(32){
		cities_vec.push_back(1);

		int a;
		vector<int>::iterator iter;
		std::vector<int> tutte;
		for (int i=2; i<=n_cities; ++i) tutte.push_back(i);

		for (int i=0; i< n_cities-1; ++i){
			a = (rand() % (tutte.size()));
			iter = tutte.begin() + a;
			cities_vec.push_back(*iter);
			tutte.erase(iter);
			//cout << "city number: " << *iter << endl;
		}
		assert (n_cities == cities_vec.size());
		random_shuffle(cities_vec.begin() +1, cities_vec.end());
		}
void Chromosome::Print_list(){
	cout << "[";
	for (int a=0; a<n_cities; a++) cout<<" " << cities_vec[a] << " ";
	cout << "]"<<endl;
	assert (n_cities == cities_vec.size());
}

void Chromosome::Print_positions(string& filename){
	ofstream out(filename);
	cout << "Printing positions into file: " << filename << endl;
	for (int a=0; a<n_cities; a++){
		out << "città n: " << cities_vec[a] << "  " << positions[cities_vec[a]-1][0] << "  " << positions[cities_vec[a]-1][1] << endl;
	}
	out.close();
	assert (n_cities == cities_vec.size());
}

void Chromosome::Set_element(int npos, int value){
	cities_vec[npos] = value;
	assert (n_cities == cities_vec.size());
}

void Chromosome::check() { //deve controllare che si sia 1 all'inizio e che non ci siano ripetizioni
	assert (n_cities == cities_vec.size());
	if(cities_vec[0] != 1){
		cout << "First city is not number 1!" << endl;
	}
	assert(cities_vec[0] == 1); //decidi quale mettere
	//controllo che ci siano tutte le città e che compaiano una sola volta
	int n_times = 0;
	for (int i=1; i<= n_cities; ++i){
		n_times=0;
		for(int v=0; v < n_cities; ++v) if(cities_vec[v] == i) n_times++;
		if(n_times != 1){
			cout << "City number: " << i << " is visited " << n_times << " times!" << endl;
		}
		assert (n_times == 1);
	}
	assert (n_cities == cities_vec.size());
}

/*double Chromosome::L1(){
	double L1=0.;
	int a,b;
	for (int i=0; i<n_cities-1; ++i){
		a = cities_vec[i];
		b = cities_vec[i + 1];
		//la norma l1 è davvero così no?
		L1 += abs(positions[a][0] - positions[b][0]) + abs(positions[a][1] - positions[b][1]); //deltax + deltay
	}
	a = cities_vec[n_cities-1];
	b = cities_vec[0];
	L1 += abs(positions[a][0] - positions[b][0]) + abs(positions[a][1] - positions[b][1]);
	return L1;
}*/

const double Chromosome::L2(){
	double L2=0.;
	int a,b;
	for (int i=0; i<n_cities-1; ++i){
		a = cities_vec[i];
		b = cities_vec[i + 1];
		//la norma l2 è davvero così no?
		L2 += pow((positions[a-1][0] - positions[b-1][0]),2) + pow((positions[a-1][1] - positions[b-1][1]),2); //deltax + deltay
		//cout << "città cosiderate: " << a << " " << b <<"  " << L2 << endl;
		//cout <<"x di "<<a<<":  "<< positions[a-1][0] << " x di "<<b<<":  " << positions[b-1][0]<<endl<<endl;
	}
	a = cities_vec[n_cities-1];
	b = cities_vec[0];
	L2 += pow((positions[a-1][0] - positions[b-1][0]),2) + pow((positions[a-1][1] - positions[b-1][1]),2);
	//cout << "città cosiderate: " << a << " " << b <<"  " << L2 << endl;
	assert (n_cities == cities_vec.size());
	return L2;
}

void Chromosome::pair_permutation(int i){
	if(i==0 ) cout << "cannot change first city" <<endl;
	assert(i!=0);
	if(i!= 0 && i!=n_cities-1) {swap(cities_vec[i] , cities_vec[i+1] );}
	if(i==n_cities-1 ) {swap(cities_vec[i] , cities_vec[1] );}
	assert (n_cities == cities_vec.size());
}

void Chromosome::pair_inversion(int i, int m){
	if(i==0 || m==0) cout << "cannot change first city" <<endl;
	if(i!= 0 && m!=0) {swap(cities_vec[i] , cities_vec[m] );}
	assert (n_cities == cities_vec.size());
}

void Chromosome::shift(int quante_pos){
	assert(quante_pos < n_cities -1 );
	int nuovo[n_cities -1 ];
	for (int i=0; i< n_cities-1; ++i) nuovo[i] = cities_vec[i+1]; //copio tutto tranne il primo
	for (int i=1; i< n_cities; ++i) {
		cities_vec[i] = nuovo[(i-1 + quante_pos)%(n_cities-1)];
	}
	assert (n_cities == cities_vec.size());
}

void Chromosome::permutation(int m){
	//cout << "valore dentro " << m << "  " << (n_cities/2.) <<endl;
	//cout << "città " << n_cities <<endl;
	assert( m < (n_cities/2.));
	swap_ranges(cities_vec.begin()+1, cities_vec.begin()+m+1 , cities_vec.begin()+m+1);

}

Generation::Generation(int n_individuals): n_chr(n_individuals){
	Chromosome* c;
	int npos;
	for (int i=0; i< n_individuals; ++i){
		//cout << "number " << i << endl;
		c = new Chromosome;
		for(auto a : gen_vec){   //per farli diversi
			if (a.Get_list() == c->Get_list()) {
				npos = rand() % (n_cities-1) +1;
				c->pair_permutation(npos);
			}
		}
		gen_vec.push_back(*c);
	}
}

void Generation::Print_generation(){
	for (auto i : gen_vec) i.Print_list();
}

Chromosome Generation::Get_first(){
	return gen_vec[0];
}

void Generation::Check_generation(){
	for (auto i : gen_vec) i.check();
}

void Generation::my_sort(){
	std::sort(gen_vec.begin() , gen_vec.end(), my_comp);
}

bool my_comp ( Chromosome& a, Chromosome& b){  //usa L2
	bool less=false;
	double prima = a.L2();
	double seconda = b.L2();
	if(prima < seconda) less = true;
	return less;
}
int Generation::selector(double r, double sel_pow){  //deve essere già ordinato // r deve essere un Rannyu()
	//std::sort(gen_vec.begin() , gen_vec.end(), my_comp);
	double a = pow(r,sel_pow);
	return (int(floor(n_chr*a))); 
	//seleziona di più gli ultimi che dovrebbo essere quelli con L2 maggiore
}


void Generation::single_mutation(int pos, double r1, double r2, double r3, double r4){
	int i, m; 
	if(r1 <= 0.1){
		//seleziono un punto in cui fare la mutazione, i = città, numero da 1 a 31 (posizioni nel vettore)
		i = (rand() % (n_cities-1)) +1; 
		gen_vec[pos].pair_permutation(i);
	}
	if(r2 <= 0.1){
		i = (rand() % (n_cities-1));  //quante_pos
		gen_vec[pos].shift(i);
	}
	if(r3 <= 0.1){
		i = (rand() % (n_cities-1)) + 1;
		m = (rand() % (n_cities-1)) + 1;
		//cerr << "pair_inversion" << endl;
		gen_vec[pos].pair_inversion(i,m);
	}
	if(r4 <= 0.1){
		i = (rand() % (int(n_cities/2)));
		//cerr << "valore passato: " << i << endl;
		gen_vec[pos].permutation(i);
	}
}


void Generation::crossover(int a, int b){
	//cerr << "crossover" << endl;
	int cut_pos = (rand() % (n_cities-2)) +1;
	//cerr << "cutpos: " << cut_pos<<endl<<endl;
	std::vector<int> copia1, copia2;
	for(int i=0; i<= cut_pos; i++) {
		copia1.push_back(gen_vec[a].Get_list()[i]);
		copia2.push_back(gen_vec[b].Get_list()[i]);
	}

//a -> 1
//b -> 2

//prendo b e cerco quali non sono già in copia1, aggiungo questi a copia1
	for (int i=0; i<gen_vec[b].Get_list().size(); ++i){
		int l=0;
		for(l=0; l< copia1.size(); ++l ) {
			if(copia1[l] == gen_vec[b].Get_list()[i]) break;
		}
		if( l == copia1.size()) copia1.push_back(gen_vec[b].Get_list()[i]);
	}
//prendo a e cerco quali non sono già in copia2, aggiungo questi a copia2
	for (int i=0; i<gen_vec[a].Get_list().size(); ++i){
		int l=0;
		for(l=0; l< copia2.size(); ++l ) {
			if(copia2[l] == gen_vec[a].Get_list()[i]) break;
		}
		if( l == copia2.size()) copia2.push_back(gen_vec[a].Get_list()[i]);
	}

//sostituisco a->, b->2  //bho provo a sostituirli alla fine così si sostituiscono a quelli "peggiori"
	for(int j=0; j<n_cities; ++j){
		gen_vec[n_chr-1].Set_element(j, copia1[j]);
		gen_vec[n_chr-2].Set_element(j, copia2[j]);
	}

	assert (n_cities == cities_vec.size());

}

void Initial_Circ(){
	//Read seed for random number
	int p1, p2;
   	ifstream Primes;
   	Primes.open("Primes");
   	Primes >> p1 >> p2 ;
   	Primes.close();

   	ifstream input;
   	input.open("seed.in");
   	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   	rnd.SetRandom(seed,p1,p2);
   	input.close();

   	//metto le città a caso su una circonferenza
   	double theta;
	for (int j=0; j<32; ++j){

		theta = -M_PI + rnd.Rannyu()*2*M_PI;

		positions[j][0] = cos(theta);
		positions[j][1] = sin(theta);

			//cerr << positions[j][0] << "  " << positions[j][1] <<endl;

	}
}

void Initial_Square(){
	//Read seed for random number
	int p1, p2;
   	ifstream Primes;
   	Primes.open("Primes");
   	Primes >> p1 >> p2 ;
   	Primes.close();

   	ifstream input;
   	input.open("seed.in");
   	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   	rnd.SetRandom(seed,p1,p2);
   	input.close();

   	//metto le città a caso in n quadrato di lato 2
   	double x, y;
	for (int j=0; j<32; ++j){

		x = -1. + rnd.Rannyu()*2.;
		y = -1. + rnd.Rannyu()*2.;

		positions[j][0] = x;
		positions[j][1] = y;

			//cerr << positions[j][0] << "  " << positions[j][1] <<endl;

	}
}

double Generation::Get_L(){ //trova L2 mediata sulla prima metà della popolazione per ogni generazione
	//ordino
	this->my_sort();
	//calcolo la media di L2 sulla prima metà
	double sum=0.0;
	for (int v=0; v<16; v++) sum += gen_vec[v].L2();
	sum /= 16.0;

	return sum; 
}


void Input(){
	ifstream ReadInput("input.dat");

	ReadInput >> circ; 
	ReadInput >> ngens;
	ReadInput >> times_per_gen;  //quante volte fa crossover + mutazioni in una generazione
	ReadInput >> maxtimes; //dopo quante generazioni si ferma
	ReadInput >> sel_pow; //per selector
	ReadInput >> volta;

	if(circ) Initial_Circ();
	else Initial_Square();

	cout << "32 cities placed " ;
	if(circ) cout << "on a circumference" << endl;
	else cout << "within a square of L=2" << endl;
	cout << "Max munber of generations: " << ngens << endl;
	cout << "Within each generation the program tries mutations and crossover " << times_per_gen << " times" << endl;
	cout << "The program will stop if cost function L2 doesn't improve for " << maxtimes << " generations" << endl;

	cout << volta << endl << endl;

}




