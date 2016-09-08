/*=============================================================================================================
                                                   individual.cpp
===============================================================================================================

 Implementation of the class Indidivual
 
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

#include "individual.h"
#include "random.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

/*=============================================================================================================
                                Implementation of the class Individual::Gene
=============================================================================================================*/

rnd::MultiNormal const * Individual::Gene::geneticArchitecture;

extern long long simulationId;
void Individual::Gene::initGeneticArchitecture()
{
	//create partial correlation matrix
	Matrix P(K,K, 0.0);

	//partial correlations within trait
    for (int cs1 = 0; cs1 < nCaste; ++cs1) {
        const int offset = cs1 * nCaste - ((cs1 + 1) * (cs1 + 2)) / 2;
		for (int tr = 0; tr < nTrait; ++tr)
			for (int i = 0; i < ageMax; ++i) {
				const int k1 = i + (tr + cs1 * nTrait) * ageMax;
				P(k1, k1) = 1.0;
                if (i - 1 >= 0) P(k1, k1 - 1) = -alpha[tr];
                if (i + 1 < ageMax) P(k1, k1 + 1) = -alpha[tr];
                for (int cs2 = cs1 + 1; cs2 < nCaste; ++cs2) {
                    const int k2 = i + (tr + cs2 * nTrait) * ageMax;
                    P(k2, k1) = P(k1, k2) = -beta[tr][offset + cs2] * P(k1, k1);
                    if (i - 1 >= 0) P(k2 - 1, k1) = P(k1, k2 - 1) = beta[tr][offset + cs2] * P(k1, k1 - 1);
                    if (i + 1 < ageMax) P(k2 + 1, k1) = P(k1, k2 + 1) = beta[tr][offset + cs2] * P(k1, k1 + 1);
                }
			}
    }
	/*
	//partial correlations between traits
    for (int tr1 = 0; tr1 < nTrait; ++tr1) {
        const int offset = tr1 * nTrait - ((tr1 + 1) * (tr1 + 2)) / 2;
        for (int cs = 0; cs < nCaste; ++cs)
            for (int i = 0; i < ageMax; ++i) {
                const int k1 = i + (tr1 + cs * nTrait) * ageMax;
                for (int tr2 = tr1 + 1; tr2 < nTrait; ++tr2) {
                    double wgt = 1.0;
                    for (int j = i; j < ageMax; ++j) {
                        const int k2 = j + (tr2 + cs * nTrait) * ageMax;
						P(k1, k2) = P(k2, k1) = - wgt * ggamma[cs][offset + tr2];
                        wgt *= delta;
					}
				}
            }
    }*/

    //write partial correlation matrix to file
	std::ostringstream oss;
	oss << "architecture_" << simulationId << ".csv";
	std::ofstream ofs(oss.str().c_str());
	verify(ofs.is_open());
	ofs.fill(',');
	ofs << "partial correlation matrix\n\n" << P << '\n';

	//calculate the mutational variance/covariance matrix
	P = P.inverse();
    Vector bias(K, eta);
    for(int i = 0; i < K; ++i) {
        const double si = sigma[i / (ageMax * nTrait)][(i / ageMax) % nTrait];
        for(int j = i + 1; j < K; ++j) {
            const double sj = sigma[j / (ageMax * nTrait)][(j / ageMax) % nTrait];
            P(i , j) = P(j , i) *= si * sj / sqrt(P(i, i) * P(j, j));
        }
        P(i, i) = sqr(si);
        bias[i] *= si;
    }
	ofs << "\nmutational variance covariance matrix\n\n" << P << '\n';

	//initialize multivariate Gaussian distribution of mutational effects
	geneticArchitecture = new rnd::MultiNormal(P, bias);
    
    ofs << "\nstandard deviations of the multivariate normal and eigenvectors of the mutational variance covariance matrix\n\n"
        << *geneticArchitecture << '\n';

    ofs << "\nsampled mutations\n\n";
    for (int i = 0; i < K; ++i) ofs << geneticArchitecture->operator()() << '\n';
    
    ofs.close();
}

void Individual::Gene::getData(std::vector<double> &sumx, std::vector<double> &sumxx) const
{
	for (int i = 0; i < K; ++i){
		double tmp = dx[i];
		sumx[i] += tmp;
		sumxx[i] += tmp * tmp;
	}
}

/*=============================================================================================================
                                Implementation of the classes Individual, Worker and Queen
=============================================================================================================*/

Queen::Queen()
//default constructor
{
	sperm = std::vector<Gene>(nGene / 2);
	develop();
}

Worker::Worker(Queen const * const mother)
//constructor for sexual reproduction -> makes a worker
{
	for (int i = 0; i < nGene; i += 2) {
		genome[i] = rnd::uniform() < 0.5 ? mother->genome[i] : mother->genome[i + 1];
		genome[i + 1] = mother->sperm[i / 2];
	}
	mutate();
	develop();
}

Queen::Queen(Queen const * const mother, Queen const * const fathersMother)
//constructor for sexual reproduction -> makes a queen
{
	for (int i = 0; i < nGene; i += 2) {
		assert(mother);
		assert(fathersMother);
		assert(i >= 0);
		assert(i < static_cast<int>(genome.size()));
		assert(i < static_cast<int>(mother->genome.size()));
		assert(i + 1 < static_cast<int>(mother->genome.size()));
		genome[i] = rnd::uniform() < 0.5 ? mother->genome[i] : mother->genome[i + 1];
		assert(i + 1>= 0);
		assert(i + 1 < static_cast<int>(genome.size()));
		genome[i + 1] = mother->sperm[i / 2];
	}
	for (int i = 0; i < nGene; i += 2)
		sperm.push_back(rnd::uniform() < 0.5 ? fathersMother->genome[i] : fathersMother->genome[i + 1]);
	mutate();
	develop();
}

void Individual::mutate()
{
	if (rnd::uniform() < nGene * ageMax * mutationRate) {
		long locus = rnd::integer(nGene);
		genome[locus].mutate();
	}
}

void Queen::mutate()
{
	Individual::mutate();
	if (rnd::uniform() < sperm.size() * ageMax * mutationRate) {
			long locus = rnd::integer(sperm.size());
			sperm[locus].mutate();
	}
}

void Queen::develop()
//determines phenotype from genotype
{
    //accumulate genetic effects
    for (int a = 0; a < ageMax; ++a) {
        //lifeHistory(productivity, a) = fecundityOffset;
        lifeHistory(survival    , a) = survivalOffset;
        for(int tr = 0; tr < nTrait; ++tr)
            for (int i = 0; i < nGene; ++i) lifeHistory(tr, a) += genome[i](0, tr, a);
    }
    
    //calculate fecundity and survival
    for (int a = 0; a < ageMax; ++a) {
        //lifeHistory(productivity, a) = exp(lifeHistory(productivity, a)*2); // how to increase the productivity?  //20160707 original: (     lifeHistory(productivity, a) = exp(lifeHistory(productivity, a)); )
        lifeHistory(survival    , a) = 0.5 * (1.0 + tanh(2.0 * lifeHistory(survival, a)));
    }
}

void Worker::develop()
//determines phenotype from genotype
{
    //accumulate genetic effects
    for (int a = 0; a < ageMax; ++a) {
        //lifeHistory(productivity, a) = productivityOffset;
        lifeHistory(survival    , a) = survivalOffset;
        for(int tr = 0; tr < nTrait; ++tr)
            for (int i = 0; i < nGene; ++i) lifeHistory(tr, a) += genome[i](1, tr, a);
    }
    
    //calculate fecundity and survival /// i dont understand how that works and how the data is stored 
    for (int a = 0; a < ageMax; ++a) {
        //lifeHistory(productivity, a) = exp(lifeHistory(productivity, a));  /// how does this work? ? 
        lifeHistory(survival    , a) = 0.5 * (1.0 + tanh(2.0 * lifeHistory(survival, a))); /// why the tanh ? this looks to me as not diploid but haploid *2
    }
}

std::ostream& operator<< (std::ostream &os, const Individual &obj)
{
    os << obj.lifeHistory << '\n';  // this shoudl be for the architecture file 
    return os;
}


// both of the survival functions need a bit of tweaking
// workers die earlier if there is increased costs of foraging but so far it is just constant this could depend on the number of foragers, size of the colony, different temporal worker castes or the size of the total population in the simulation (density dependency) 
//queens just die in response to the age specific survival probalbility which depends only on age specific (queen)genes 
// in both cases we need to add extrinsic mortality, for worker costs of foraging can represent it but for queen we should have a expontetially dereasing function of extrinsic risk that declines as the colony size increases. 
bool Worker::doesSurvive() const     
{
    /*if(age < ageMax) {
        if (rnd::uniform() < lifeHistory(survival, age) * exp(-costOfForaging * lifeHistory(productivity, age)))  // thrade off between survival and productivity removed 20160708 with sander 
            return true;
    }*/
	const double agespecificsurvival = lifeHistory(survival, age);
	if (age < ageMax) {
		if (rnd::uniform() < agespecificsurvival) {// lifeHistory(survival, age))// needs curly brackets because the if loop is longer than one line because of the out command
		//	std::cout <<"worker " << "," << agespecificsurvival <<"," <<" at age " <<","<< age   << "\n ";
			return true;
												 }
	}
	return false;
}

bool Queen::doesSurvive() const   // to me the does surviva function look the same for queens and workers , how is it possible thet i takes diffeten values fron the genetics? 
{
	const double agespecificsurvival = lifeHistory(survival, age);
	//std::cout << lifeHistory << '\n';
    if(age < ageMax) {
		if (rnd::uniform() < lifeHistory(survival, age)) {
			//std::cout << "queen " << "," << agespecificsurvival << "," << " at age " << "," << age << "\n ";
			return true;  // how do i put the colony size in here ? is there a way to count the workers? exponentially declining function : e^(-workers/2) one with 0 workers and 0 with 10 workers (mathamatica:  Plot[E ^ (-x/2), {x, 0, 10}])
														}
		}
    return false;
}

void Queen::getData(std::vector<double>& sumx, std::vector<double>& sumxx) const    //collect data for the analysis, its calles in the simulation file and also a member function of the class queen? 
{
	std::vector<double> tmpx(nCaste * nTrait * ageMax, 0.0), tmpxx(nCaste * nTrait * ageMax, 0.0); // why always tempx tempss sumx sumxx? 
	for (int i = 0; i < nGene; ++i) genome[i].getData(tmpx, tmpxx);
	for (int i = 0; i < nCaste * nTrait * ageMax; ++i) { // are these multiplications -- > give the lenght of the loop 
		sumx[3 + i] += tmpx[i] / nGene;  /// ?? 
		sumxx[3 + i] += tmpxx[i] / nGene;  // again not entirely sure how this works... 
	}
}