/*=============================================================================================================
                                                   individual.h
===============================================================================================================

 Definition of the class Indidivual
 
 C++-code accompanying:
 
        authors 
        title
 
 Written by:
        Boris Kramer and G. Sander van Doorn
        Groningen Institute for Evolutionary Life Sciences (Gelifes)
        University of Groningen
        the Netherlands
 
 Program version
        29/03/2016	: first version
 
 =============================================================================================================*/

#ifndef individual_h
#define individual_h

#include<array>
#include<vector>
#include "parameters.h"
#include "corand.h"

class Individual
{
public:
    enum Traits {survival};
    void incrementAge() { ++age; }   // if called just increses the age of the worker or queen by 1 
	virtual bool doesSurvive() const = 0;
    friend std::ostream& operator<< (std::ostream&, const Individual&);
	double getLifeHistoryTrait(const Traits &tr) const { return lifeHistory(static_cast<Matrix::size_type>(tr), age); }
	const Matrix& getLifeHistory() const { return lifeHistory; }
	int getAge() const { return age; }
	
	class Gene
	{
	public:
        Gene() : dx(K, 0.0) { mutate(); }//initial mutations applied to the initial values   // the 0.0 represents the survival offset 
        void mutate() {dx += geneticArchitecture->operator()();}    // i dont get the code here 
        double operator()(const int &cs, const int &tr, const int &ag) {return dx[ag + (tr + cs * nTrait) * ageMax];}
		void getData(std::vector<double>&, std::vector<double>&) const;
		static void initGeneticArchitecture();
	private:
		Vector dx;
		static rnd::MultiNormal const * geneticArchitecture;
        const static int K = nCaste * nTrait * ageMax;
	};
protected:
    Individual() : lifeHistory(nTrait, ageMax), age(0) {}
	//genotype
	typedef std::array<Gene, nGene> Genome;// array with length nGene
	Genome genome;
	
	//phenotype
	int age;
    Matrix lifeHistory;
	virtual void develop() = 0;  // has to be implemented in the derived worker and queen class
	virtual void mutate();    // can be overridden by the derived class
};

class Worker; // forward declaration. 
class Queen: public Individual // derived class
{
public:
	friend class Worker;   // worker access to the private section of the queen class
	Queen();
	Queen(Queen const * const, Queen const * const);
	bool doesSurvive() const;
	void getData(std::vector<double>&, std::vector<double>&) const;
private:
	std::vector<Gene> sperm;
	void develop();
	void mutate();
};

class Worker: public Individual
{
public:
	Worker(Queen const * const mother);
	bool doesSurvive() const;
private:
	Worker() = delete;
	void develop();
};


#endif