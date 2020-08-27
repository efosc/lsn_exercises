#include <vector>
#include <string>
#include <climits>
#include <iostream>
#include <cassert>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "random.h"
#include "funzioni.h"
#include "mpi.h"   

using namespace std;

int main(int argc, char *argv[]){

	Initial_Square();  //all'inizio metto a caso in un quadrato
   	int ngens = 500, times_per_gen=300, nmigr=50;
	int which1, which2;
	int maxtimes=30, ntimes=0;
	double r, r1, r2, r3, r4;
	double cost_new, cost_old, cost_ave;

///////MPI////////////////////
	int size, rank, scambio, scambio2, scambio3, itag=1;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat[4];
	MPI_Request req, req2;

///ogni continente legge due Primes diversi
	int p1, p2;
	ifstream Primes2("Primes");
	if(rank == 0){
		Primes2 >> p1 >> p2 >> p1 >> p2;
	} else if (rank ==1){
		Primes2 >> p1 >> p2 >> p1 >> p2 >> p1 >> p2;
	}else if (rank ==2){
		Primes2 >> p1 >> p2 >> p1 >> p2 >> p1 >> p2 >> p1 >> p2;
	}else if (rank ==3){
		Primes2 >> p1 >> p2 >> p1 >> p2 >> p1 >> p2 >> p1 >> p2 >> p1 >> p2;
	}	

   	rnd.SetRandom(seed,p1,p2); //ho già letto seed in initial_square()

//genero la prima generazione in ogni continente
	Generation routes(200);
	routes.Check_generation();
	cost_old = routes.Get_first().L2();

///-----------------/////
//per un tot di volte faccio crossover + mutazioni:
	// -> ordino
	// -> pesco 2 su cui fare crossover 
	// -> aggiungo eventualmente mutazioni
	// -> rifaccio per un tot di volte e ho una nuova generazione
		ofstream out;

	for (int igen =0; igen<ngens; ++igen){
		if(igen%10==0) {
			cerr << rank << " Generation number: " << igen <<endl;
			if(rank == 0) out.open("first_vec0d.dat", ios::app);
			else if (rank == 1) out.open("first_vec1d.dat",ios::app);
			else if (rank == 2) out.open("first_vec2d.dat", ios::app);
			else if (rank == 3) out.open("first_vec3d.dat", ios::app);
			out << igen << "      " ;
			for (unsigned int a=0; a<32; a++) {
				out << routes.Get_first().Get_list()[a] << "  ";
			} 
			out << endl;
			out.close();
	
		}
		for (int itime=0; itime<times_per_gen ; ++itime ){
			//cerr << "volta numero: " << itime << endl;
			routes.my_sort(); //ordina dovrebbe
			r = rnd.Rannyu();
			//cerr << "estraggo un numero " << r << " " ;
			which1 = routes.selector(r);
			//cerr << " e scelgo il cromosoma numero: " << which1 << endl;
			r = rnd.Rannyu();
			//cerr << "estraggo un numero " << r << " " ;
			which2 = routes.selector(r);
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
		}

//////migrations////////////	
		if((igen%nmigr)==0 && igen!=0) {
			cout << "Migration" << endl;
			scambio = ((int)rnd.Rannyu(0.,1.)*INT_MAX%3) + 1;
			if(scambio == 1) {scambio2=2; scambio3=3;}
			else if(scambio == 2) {scambio2=1; scambio3=3;}
			else {scambio2=1; scambio3=2;}
			routes.my_sort(); //ordina per prendere il primo
			if(rank == 0){
				MPI_Isend(&routes.Get_first().Get_list()[0], 32, MPI_INTEGER, scambio, itag, MPI_COMM_WORLD, &req);
				MPI_Recv(&routes.Get_first().Get_list()[0], 32, MPI_INTEGER, scambio, itag, MPI_COMM_WORLD, &stat[1]);
			}else if(rank == scambio){
				MPI_Send(&routes.Get_first().Get_list()[0], 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD);
				MPI_Recv(&routes.Get_first().Get_list()[0], 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD, &stat[0]);
			}else if(rank == scambio2){
				MPI_Isend(&routes.Get_first().Get_list()[0], 32, MPI_INTEGER, scambio3, itag, MPI_COMM_WORLD, &req2);
				MPI_Recv(&routes.Get_first().Get_list()[0], 32, MPI_INTEGER, scambio3, itag, MPI_COMM_WORLD, &stat[3]);
			}else if(rank == scambio3){
				MPI_Send(&routes.Get_first().Get_list()[0], 32, MPI_INTEGER, scambio2, itag, MPI_COMM_WORLD);
				MPI_Recv(&routes.Get_first().Get_list()[0], 32, MPI_INTEGER, scambio2, itag, MPI_COMM_WORLD, &stat[2]);
			}
		}

		cost_ave = routes.Get_L(); //quella mediata sulla prima metà
		cost_new = routes.Get_first().L2(); //solo del primo
		Print_lfirst(rank, cost_new);
		Print_lave(rank, cost_ave);
		if(cost_new >= cost_old) ntimes ++;
		else {cost_old = cost_new; ntimes=0;} //mi tengo quello migliore per confrontare con quelli nuovi
		if(ntimes == maxtimes) {
			cout << "Rank: " << rank <<  " Cost didn't improve for " << maxtimes << " generations" <<endl;
			//cout << "Rank: " << rank << "Stopping at generation number: " << igen << endl;
			//break;   //se non miglioro per un tot esco
		}
		
		routes.Check_generation();
	}

	// dopo n generazioni guardo quello più breve
	routes.my_sort();
	if(rank == 0) out.open("first_vec0.dat", ios::app);
	else if (rank == 1) out.open("first_vec1.dat",ios::app);
	else if (rank == 2) out.open("first_vec2.dat", ios::app);
	else if (rank == 3) out.open("first_vec3.dat", ios::app);
	
	for (unsigned int a=0; a<32; a++) {
		int num = routes.Get_first().Get_list()[a];
		out << num << "  " <<positions[num][0] << "  " << positions[num][1] << endl;
	}
	
	MPI_Finalize();

return 0;
}

///////print//////////////////////
void Print_lfirst(int rank, double cost_new){
	ofstream out;
	if(rank == 0) out.open("lfirst0.dat", ios::app);
	else if (rank == 1) out.open("lfirst1.dat",ios::app);
	else if (rank == 2) out.open("lfirst2.dat", ios::app);
	else if (rank == 3) out.open("lfirst3.dat", ios::app);
	
	out << cost_new << endl;
	out.close();
}

void Print_lave(int rank, double cost_ave){
	ofstream out;
	if(rank == 0) out.open("lave0.dat", ios::app);
	else if (rank == 1) out.open("lave1.dat",ios::app);
	else if (rank == 2) out.open("lave2.dat", ios::app);
	else if (rank == 3) out.open("lave3.dat", ios::app);
	
	out << cost_ave << endl;
	out.close();
}

///////////////FUNZIONI.cpp///////////////////////////////////
Chromosome::Chromosome(): n_cities(32){
		cities_vec.push_back(1);

		int a;
		vector<int>::iterator iter;
		std::vector<int> tutte;
		for (unsigned int i=2; i<=n_cities; ++i) tutte.push_back(i);

		for (unsigned int i=0; i< n_cities-1; ++i){
			a = ((int)(rnd.Rannyu(0.,1.)*INT_MAX) % (tutte.size()));
			iter = tutte.begin() + a;
			cities_vec.push_back(*iter);
			tutte.erase(iter);
		}
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
	for (unsigned int i=1; i<= n_cities; ++i){
		n_times=0;
		for(unsigned int v=0; v < n_cities; ++v) if(cities_vec[v] == i) n_times++;
		if(n_times != 1){
			cout << "City number: " << i << " is visited " << n_times << " times!" << endl;
		}
		assert (n_times == 1);
	}
	assert (n_cities == cities_vec.size());
}

void Chromosome::Set_element(int npos, int value){
	cities_vec[npos] = value;
	assert (n_cities == cities_vec.size());
}

const double Chromosome::L2()const {
	double L2=0.;
	int a,b;
	for (unsigned int i=0; i<n_cities-1; ++i){
		a = cities_vec[i];
		b = cities_vec[i + 1];
		L2 += pow((positions[a-1][0] - positions[b-1][0]),2) + pow((positions[a-1][1] - positions[b-1][1]),2); //deltax + deltay
	}
	a = cities_vec[n_cities-1];
	b = cities_vec[0];
	L2 += pow((positions[a-1][0] - positions[b-1][0]),2) + pow((positions[a-1][1] - positions[b-1][1]),2);
	assert (n_cities == cities_vec.size());
	return L2;
}

////////----mutations--------------///////////////
void Chromosome::pair_permutation(unsigned int i){
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

void Chromosome::shift(unsigned int quante_pos){
	assert(quante_pos < n_cities -1 );
	int nuovo[n_cities -1 ];
	for (unsigned int i=0; i< n_cities-1; ++i) nuovo[i] = cities_vec[i+1]; //copio tutto tranne il primo
	for (unsigned int i=1; i< n_cities; ++i) {
		cities_vec[i] = nuovo[(i-1 + quante_pos)%(n_cities-1)];
	}
	assert (n_cities == cities_vec.size());
}

void Chromosome::permutation(int m){
	assert( m < (n_cities/2.));
	swap_ranges(cities_vec.begin()+1, cities_vec.begin()+m+1 , cities_vec.begin()+m+1);

}




///////////////GENERATION///////////////////////////

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

Chromosome& Generation::Get_first(){
	return gen_vec[0];
}

void Generation::Check_generation(){
	for (auto i : gen_vec) i.check();
}

void Generation::my_sort(){
	std::sort(gen_vec.begin() , gen_vec.end(), my_comp);
}

bool my_comp ( const Chromosome& a, const Chromosome& b){  //usa L2
	bool less=false;
	double prima = a.L2();
	double seconda = b.L2();
	if(prima < seconda) less = true;
	return less;
}
int Generation::selector(double r){  //deve essere già ordinato // r deve essere un Rannyu()
	//std::sort(gen_vec.begin() , gen_vec.end(), my_comp);
	double a = pow(r,1.0);
	return (int(floor(n_chr*a))); 
	//seleziona di più gli ultimi che dovrebbo essere quelli con L2 maggiore
}

void Generation::single_mutation(int pos, double r1, double r2, double r3, double r4){
	int i, m;
	if(r1 <= 0.1){
		//seleziono un punto in cui fare la mutazione, i = città, numero da 1 a 31 (posizioni nel vettore)
		i = ((int)(rnd.Rannyu(0.,1.)*INT_MAX) % (n_cities-1)) +1; 
		//dist = pow(positions[(i%n_cities)][0] - positions[((i+1)%n_cities)][0],2) + 
		//cerr << "pair_permutation" << endl;
		gen_vec[pos].pair_permutation(i);
	}
	if(r2 <= 0.1){
		i = ((int)(rnd.Rannyu(0.,1.)*INT_MAX) % (n_cities-1));  //quante_pos
		gen_vec[pos].shift(i);
	}
	if(r3 <= 0.1){
		i = ((int)(rnd.Rannyu(0.,1.)*INT_MAX) % (n_cities-1)) + 1;
		m = ((int)(rnd.Rannyu(0.,1.)*INT_MAX) % (n_cities-1)) + 1;
		//cerr << "pair_inversion" << endl;
		gen_vec[pos].pair_inversion(i,m);
	}
	if(r4 <= 0.1){
		i = ((int)(rnd.Rannyu(0.,1.)*INT_MAX) % (int(n_cities/2)));
		//cerr << "valore passato: " << i << endl;
		gen_vec[pos].permutation(i);
	}
}


void Generation::crossover(int a, int b){
	//cerr << "crossover" << endl;
	int cut_pos = ((int)(rnd.Rannyu(0.,1.)*INT_MAX) % (n_cities-2)) +1;
	//cerr << "cutpos: " << cut_pos<<endl<<endl;
	std::vector<unsigned int> copia1, copia2;
	for(int i=0; i<= cut_pos; i++) {
		copia1.push_back(gen_vec[a].Get_list()[i]);
		copia2.push_back(gen_vec[b].Get_list()[i]);
	}

//a -> 1
//b -> 2

//prendo b e cerco quali non sono già in copia1, aggiungo questi a copia1
	for (unsigned int i=0; i<gen_vec[b].Get_list().size(); ++i){
		unsigned int l=0;
		for(l=0; l< copia1.size(); ++l ) {
			if(copia1[l] == gen_vec[b].Get_list()[i]) break;
		}
		if( l == copia1.size()) copia1.push_back(gen_vec[b].Get_list()[i]);
	}
//prendo a e cerco quali non sono già in copia2, aggiungo questi a copia2
	for (unsigned int i=0; i<gen_vec[a].Get_list().size(); ++i){
		unsigned int l=0;
		for(l=0; l< copia2.size(); ++l ) {
			if(copia2[l] == gen_vec[a].Get_list()[i]) break;
		}
		if( l == copia2.size()) copia2.push_back(gen_vec[a].Get_list()[i]);
	}

	/*cout << "copia1 ";
	for(int i: copia1) cout << i<<  "  ";
	cout << endl << "copia2 ";
	for(int i: copia2) cout << i<<  "  ";
	cout << endl;*/
//sostituisco a->1, b->2  //ora risostituisco nei posti which1 e which2 (ordine non dovrebbe contare)
	for(unsigned int j=0; j<n_cities; ++j){
		gen_vec[n_chr -1 ].Set_element(j, copia1[j]);
		gen_vec[n_chr -2].Set_element(j, copia2[j]);
	}

	assert (n_cities == cities_vec.size());

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

////MPI/////////////////////////////////////////
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


/////////////////////RANDOM////////////////

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}
