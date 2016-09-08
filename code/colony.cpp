/*=============================================================================================================
                                                   colony.cpp
===============================================================================================================

 Implementation of the class Colony
 
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

#include <cassert>
#include "colony.h"
#include "random.h"

///default constructor for Colony
Colony::Colony(Queen * const q)
	: queen(q),
	  colonyFitness(0.0),
	  deadWorkers(std::vector<int>(ageMax, 0)),
	  deadQueens(std::vector<int>(ageMax, 0))
{
	//if (!q) { throw std::invalid_argument("A colony must have a queen"); }
	assert(q);
}

Colony::~Colony() // destructor
{
	delete queen;
	while (!workers.empty()) {
		delete workers.back();
		workers.pop_back();
	}
}

bool Colony::growColony()
//iterates over the generation
{


	double productivity = 0.8;// changed that today from 0 to 1 to see 20160613  --> no changes // in geenral i dont see how productivity rises and where to play around with that  // 20160706 changes from 0 -->0.8
	for (int j = 0; j < ageMax; ++j) {

		//std::cout << workers.size()<<"w " << rnd::uniform()<< "n/"; // check the number of workers produced 
/*
	//survival
	//workers
		for (unsigned i = 0u; i < workers.size(); ) {
			if (!workers[i]->doesSurvive() || rnd::uniform() < extrinsicMortality) { // added || rnd::uniform() <  extrinsicMortality so that worker also die from extrinsic mortality // is this correct? 
			//	std::cout << rnd::uniform() << " < "<<  extrinsicMortality << "d "  << "/n"; // check the number of workers produced 
				delete workers[i];
				workers[i] = workers.back();
				workers.pop_back();
				assert(j >= 0);
				assert(j < static_cast<int>(deadWorkers.size()));
				++deadWorkers[j];
			}
			else ++i;
		} */
		//queen survival
		if (!queen->doesSurvive()|| rnd::uniform() < extrinsicMortality * (1/(2 + workers.size()))){  //rnd::uniform() <  extrinsicMortality * exp(-workers.size()+0.00001)) {   // the extrinsic mortality need to be implemented here  dont know if i can just add it to the parameters 
			assert(j >= 0);
			assert(j < static_cast<int>(deadQueens.size()));
			++deadQueens[j];
			return false;
			}
			

		//accumulate resources  // this should go to up because productivity is used in here? 20160712
		double resources = 0.5 + workers.size();   // worker.size gives the size of the colony   exchangesd from 1 * worker.size()
												   //for (unsigned i = 0u; i < workers.size(); ++i) resources += workers[i]->getLifeHistoryTrait(Individual::productivity); // what is the unsigned i= 0u? 
												   //const double eggs = queen->getLifeHistoryTrait(Individual::productivity);
		//const double eggs = 2.0; // not sure if we should put assumptions on the number off eggs into the modeli guess especially at small colony sized fertility of the queen should not be an issue
		productivity = /*resources;*/ eggs * resources/ (eggs + resources);  // why is that needed? what are the underlying assumptions? leveling off? if eggs is fixed to 1 productivity will always be lower that 1 how does that translate into the production of new workers? 
		//std::cout << productivity << "p ";

		//reproduction
			//produce reproductives
		if(workers.size() > 5 || j == ageMax - 1){   // sets the colony size at which the colony switches to sexual reproduction. in this case the colony gets an increase in colonyFitness which increases the probability that the genes of this colony are drawn for the next generation 
			colonyFitness += productivity; // if there is more that 5 workers then the resources of the colony are invested into the production of sexuals by incresing the colony fitness which in turn increases the genes of the colony to be selected for the next generation
		}
		else 
		{
			//produce workers
			int i = productivity > 0.0 ? rnd::poisson(productivity) : 0; // where are the costs for each individual worker defined? does that setup mean 1 in productivity may produce one worker  (if else statement)  // how does rnd::poisson(productivity) work?  its a poisson distribution with a mean of the productivity so production will be variing around the value of productivity  // ? = conditional ternary operator if prody >0  then  rand(proy) else 0
			while (i)													// as long as there is productivity a random number of workers is produced? 
			{
				workers.push_back(new Worker(queen)); // this is memory allocation stuff? to renew the vector of workers without empty spots? 
				--i;
			}
		}
		//increment ages
		for (unsigned i = 0u; i < workers.size(); ++i) workers[i]->incrementAge(); // this loop runs through all workers in a colony
		queen->incrementAge();



		//survival    20160729 moved down here because otherwise there are never workers in the 1 age class as it inceases before they can die 
		//workers
		for (unsigned i = 0u; i < workers.size(); ) {
			if (!workers[i]->doesSurvive() || rnd::uniform() < extrinsicMortality) { // added || rnd::uniform() <  extrinsicMortality so that worker also die from extrinsic mortality // is this correct? 
																					 //	std::cout << rnd::uniform() << " < "<<  extrinsicMortality << "d "  << "/n"; // check the number of workers produced 
				delete workers[i];
				workers[i] = workers.back();
				workers.pop_back();
				assert(j >= 0);
				assert(j < static_cast<int>(deadWorkers.size()));
				++deadWorkers[j];
			}
			else ++i;
		}
		/*	//queen survival
		if (!queen->doesSurvive()) {//|| rnd::uniform() < extrinsicMortality * (1/(1 + workers.size()))){  //rnd::uniform() <  extrinsicMortality * exp(-workers.size()+0.00001)) {   // the extrinsic mortality need to be implemented here  dont know if i can just add it to the parameters 
			assert(j >= 0);
			assert(j < static_cast<int>(deadQueens.size()));
			++deadQueens[j];
			return false;
		}*/

	}
	return true;
}

		/*
		//accumulate resources 
		double resources = 1.0 * workers.size();   // worker.size gives the size of the colony  
		//for (unsigned i = 0u; i < workers.size(); ++i) resources += workers[i]->getLifeHistoryTrait(Individual::productivity); // what is the unsigned i= 0u? 
		//const double eggs = queen->getLifeHistoryTrait(Individual::productivity);
		const double eggs = 1.0;
		productivity = eggs * resources / (eggs + resources);
		
		//produce new workers
		int i = productivity > 0.0 ? rnd::poisson(productivity) : 0; 
		while (i) {
			workers.push_back(new Worker(queen));
			--i;
		}
		
		//increment ages
		for (unsigned i = 0u; i < workers.size(); ++i) workers[i]->incrementAge();
		queen->incrementAge();     */
	


void Colony::getData(std::vector<double> &sumx, std::vector<double> &sumxx) const
{
	queen->getData(sumx, sumxx);
}