#include "random.h"

#ifndef __ES9__
#define __ES9__

class Chromosome {  //per ora faccio vettore con 5 città

protected:
	unsigned int n_cities;
	std::vector<unsigned int> cities_vec;
	//double positions[32][2];

public:
	//Chromosome(int m_cities);
	Chromosome();
	//void Print_list();
	//int Get_first_val();
	std::vector<unsigned int>& Get_list(){return cities_vec;}
	void Set_element(int pos, int value);
	void check(); /// controlla che le città non si ripetano e che la prima sia 1
////////////misura////////////////////////////////
	const double L2() const;
////genetic mutations (within the same chromosome)/////////////////////
	void pair_permutation(unsigned int i);
	void shift(unsigned int m);
	void permutation(int m);
	void pair_inversion(int i, int m);
/////////////////////////////////////////////
	

};

class Generation:Chromosome{

private:
	const int n_chr;
	std::vector<Chromosome> gen_vec; //vettore con tutti i cromosomi di quella generazione


public:
	Generation(int n_individuals);
	void Check_generation(); //controlla che tutti gli individui rispettino i bonds
	void Print_generation();
	Chromosome& Get_first();
	void my_sort();
	int selector(double r); // //ordina sempre L2 spero e seleziona  un cromosoma //r deve essere un Rannyu()
	//void mutations(Random& rnd); // scorre tutti i cromosomi e fa una mutazione con una certa prob
	void crossover(int, int);
	void single_mutation(int,double r1,double r2, double r3, double r4);
	double Get_L();
	int Get_nchr(){return n_chr;}
	// mutazione solo sul cromosoma numero pos 
	// i vari r devono essere Rannyu()
};

//other functions
bool my_comp ( const Chromosome& a, const Chromosome& b);
void Initial_Circ();
void Initial_Square();
void Print_lfirst(int, double);
void Print_lave(int, double);
void Set_gen(int);
void Send_rec_scheme();
//vars
double positions[32][2];
int seed[4], p1, p2;
Random rnd;
int receivers[4] = {9,9,9,9}, senders[4] = {9,9,9,9};
#endif

#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
};

#endif // __Random__

