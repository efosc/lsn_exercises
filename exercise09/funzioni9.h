#include "random.h"

#ifndef __ES9__
#define __ES9__

class Chromosome {  //per ora faccio vettore con 5 città

protected:
	int n_cities;
	std::vector<int> cities_vec;
	//double positions[32][2];

public:
	//Chromosome(int m_cities);
	Chromosome();
	void Print_list();
	std::vector<int> Get_list(){return cities_vec;}
	void Print_positions(std::string& filename);
	void Set_element(int pos, int value);
	void check(); /// controlla che le città non si ripetano e che la prima sia 1
////////////misura/////////////////////////////////
	const double L2();
////genetic mutations (within the same chromosome)/////////////////////
	void pair_permutation(int i);
	void shift(int m);
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
	Chromosome Get_first();
	void my_sort();
	int selector(double r, double sel_pow); // //ordina sempre L2 spero e seleziona  un cromosoma //r deve essere un Rannyu()
	//void mutations(Random& rnd); // scorre tutti i cromosomi e fa una mutazione con una certa prob
	void crossover(int, int);
	void single_mutation(int,double r1,double r2, double r3, double r4);
	double Get_L();
	int Get_nchr(){return n_chr;}
	// mutazione solo sul cromosoma numero pos 
	// i vari r devono essere Rannyu()
};

bool my_comp ( Chromosome& a, Chromosome& b);
//inline double positions[32][2];
//inline int seed[4];
//inline Random rnd;
void Initial_Circ();
void Initial_Square();
void Input();
#endif