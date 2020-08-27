#include "random.h"

#ifndef __ES10__
#define __ES10__

class Chromosome {  //per ora faccio vettore con 32 città

protected:
	int n_cities;
	std::vector<int> cities_vec;
	//double positions[32][2];

public:
	//Chromosome(int m_cities);
	Chromosome();
	Chromosome(Chromosome& altro);
	void Print_list();
	std::vector<int> Get_list(){return cities_vec;}
	int Get_ncities(){return n_cities;}
	void Print_positions(std::string& filename);
	void Set_element(int pos, int value);
	void check(); /// controlla che le città non si ripetano e che la prima sia 1
////////////misura/////////////////////////////////
	//double L1(); //considero solo L2 per ora
	const double L2();
	//bool operator == ( Chromosome& altro);
	//bool operator < ( Chromosome& altro);
	//bool my_comp ( Chromosome& altro);
////genetic mutations (within the same chromosome)/////////////////////
	void pair_permutation(int i);
	void shift(int m);
	void permutation(int m);
	void pair_inversion(int i, int m);
/////////////////////////////////////////////
	

};

bool my_comp ( Chromosome& a, Chromosome& b);
//inline double positions[32][2];
//inline int seed[4];
//inline Random rnd;
void Initial_Circ();
void Initial_Square();
/////GENERALI//////////////
//void Set_annealing_schedule(void);
void Move(double beta);
void Setting(void);

//ora ho un solo cromosoma, quindi questa classe non mi serve giusto?
/*class Generation:Chromosome{

private:
	const int n_chr;
	std::vector<Chromosome> gen_vec; //vettore con tutti i cromosomi di quella generazione


public:
	Generation(int n_individuals);
	void Check_generation(); //controlla che tutti gli individui rispettino i bonds
	void Print_generation();
	Chromosome Get_first();
	void my_sort();
	int selector(double r); // //ordina sempre L2 spero e seleziona  un cromosoma //r deve essere un Rannyu()
	//void mutations(Random& rnd); // scorre tutti i cromosomi e fa una mutazione con una certa prob
	void crossover(int, int);
	void single_mutation(int,double r1,double r2, double r3, double r4);
	double Get_L();
	int Get_nchr(){return n_chr;}
	// mutazione solo sul cromosoma numero pos 
	// i vari r devono essere Rannyu()
};*/

#endif