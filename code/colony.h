/*=============================================================================================================
                                                   colony.h
===============================================================================================================

 Definition of the class Colony
 
 C++-code accompanying:
 
        authors boris and sander 
        title
 
 Written by:
        Boris Kramer and G. Sander van Doorn
        Groningen Institute for Evolutionary Life Sciences (Gelifes)
        University of Groningen
        the Netherlands
 
 Program version
        29/03/2016	: first version
 
 =============================================================================================================*/

#ifndef colony_h
#define colony_h

#include "individual.h"
#include <vector>

class Colony //defines the class colony
{

	typedef std::vector<Worker*> Workers;

public:

	Colony(Queen * const);
	~Colony(); // destructor? 
	bool growColony();
	Queen const * getQueen() const { return queen; }
	const Workers& getWorkers() const { return workers; }
 	double getFitness() const { return colonyFitness; }
	size_t size() const { return workers.size(); }
	void getData(std::vector<double>&, std::vector<double>&) const;

	///Get the number of dead workers per age class after 'growColony'. Size of the vector equals the number of age classes
	const std::vector<int>& getDeadWorkers() const { return deadWorkers; }

	///Get the number of dead queens per age class after 'growColony'. Size of the vector equals the number of age classes
	const std::vector<int>& getDeadQueens() const { return deadQueens;  }

private:
	double colonyFitness;
	Workers workers;
	Queen* queen;

	///The number of dead workers per age class (after 'growColony')
	std::vector<int> deadWorkers;

	///The number of dead queens per age class (after 'growColony')
	std::vector<int> deadQueens;
};

#endif