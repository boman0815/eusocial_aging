/*=============================================================================================================
                                                   simulation.cpp
===============================================================================================================

 Model implementation; entry point main() 
 
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
 
 ============================================================================================================*/

#include "simulation.h"
#include "colony.h"
#include "random.h"
#include "utils.h"
#include "parameters.h"
#include <numeric>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iterator>
#include <fstream>
#include <sstream>

//why arend all the files listed here ? are the other dependencies in the other files? foreward declaration?

/*=============================================================================================================
                                                   global variables
=============================================================================================================*/

long long simulationId;   // what does the long long mean? should be the data format? 
Colony * population[nColony]; // light blue indicates classes
std::ofstream ofs; //  storing the data
std::vector<int> deadQueenCounter(ageMax, 0);
std::vector<int> deadWorkerCounter(ageMax, 0);
/*=============================================================================================================
                                                   model implementation
=============================================================================================================*/

void init()
//initializes the population
{
	Individual::Gene::initGeneticArchitecture();
	Queen *q = new Queen();// just this one queen initializes the rest of the starting population 
	//initialize colonies with full sibs
	for (int i = 0; i < nColony; ++i)
		population[i] = new Colony(new Queen(q, q));
	delete q; // why is there a delete q
}


//Sums the elements, e.g. {1,2,3} + {2,3,4} becomes {3,5,7}
std::vector<int> sumElements(const std::vector<int>& v, const std::vector<int>& w)
{
	assert(v.size() == w.size());
	//Copy v as the initial values of z
	std::vector<int> z = v;
	//Add the values of w to z
	const size_t sz = v.size();
	for (size_t i = 0; i != sz; ++i)
	{
		z[i] += w[i];
	}
	return z;
}


bool iterate(
	const int generation, 
	std::ofstream& file, 
	std::vector<double>& sumx, 
	std::vector<double>& sumxx, 
	std::vector<int>& sumDeadWorkers,
	std::vector<int>& sumDeadQueens
)// bool , so just returns true or false 
{
	//get colony fitness build vectors for data collection
	rnd::discrete_distribution fitness(nColony);//from sanders random number generators not sure why this has to be a class?
	double sumFitness = 0.0;
	sumx = sumxx = std::vector<double>(3 + nCaste * nTrait * ageMax, 0.0);// collects the data, why sumx=sumxx for what is that --> just initializes both vectors with zeros
	//grow colonies
	for (int i = 0; i < nColony; ++i) {
		double tmp = 0.0;
		if (population[i]->growColony()) {//it determines if the queen is dead or not
			tmp = 1.0; // how do those 2 tmp doubles work? 
		}
		sumx[0] += tmp;
		sumxx[0] += tmp * tmp;
		sumx[1] += tmp = population[i]->size();
		sumxx[1] += tmp * tmp;
		sumx[2] += tmp = sumFitness += fitness[i] = population[i]->getFitness();
		sumxx[2] += tmp * tmp;  //  i do not fully understand that lines of code this is supposed to crea
		//std::cout << i << '\t' << population[i]->size() << '\t' << population[i]->getFitness() << '\t' << population[i]->getQueen() << '\n';
	}

	//get the life history data : Must be done at the end of the season(generation), when workers are still present 
	if (generation % dataInterval == 0) {  //cout only after each data interval
		saveLifeHistory(generation, population, nColony, file); //saveLifeHistory(generation, population, nColony, ofs);
		}
	
	
	if (sumFitness > 0.0) { // only applied to is extinct labeled colonies why is it not a member of the colony properties? 
		//reproduce
		Colony * tmp[nColony];
		for (int i = 0; i < nColony; ++i) {
			const int sourceColonyQueen = fitness.sample();
			const int sourceColonyDrone = fitness.sample();
			tmp[i] = new Colony(new Queen(population[sourceColonyQueen]->getQueen(), population[sourceColonyDrone]->getQueen())); // here the tmp is used to sample new queens? 
		}

		//Measure then number of dead workers and queens, add these to the current sums
		for (int i = 0; i != nColony; ++i) {
			sumDeadQueens = sumElements(sumDeadQueens, population[i]->getDeadQueens());
			sumDeadWorkers = sumElements(sumDeadWorkers, population[i]->getDeadWorkers());
		}

		//replace old colonies    // how does this select is extinct colonies? 
		for (int i = 0; i < nColony; ++i) {
			delete population[i];
			population[i] = tmp[i];   /// puts new queens? as samplede above for isextinct colonies 
		}
		return true;
	}
	else {
		std::cout << "population extinct!\n"; // if all colonies are extinct , how is the else loop entered? when !isextinct=false for all colonies? 
		return false; // what does the false in comparison of the true 3 lines above mean? 
	}
}

void dataAnalysis(const int& t, std::vector<double>& sumx, std::vector<double>& sumxx)
{
	//collect data
	for (int i = 0; i < nColony; ++i) population[i]->getData(sumx, sumxx);
	
	//process data
	for (int k = 0; k < static_cast<int>(sumx.size()); ++k) {//where is the size of sumx defined?--> in the definition of sumxx(3 + nCaste * nTrait * ageMax, 0.0)
		sumx[k] /= nColony; // probably the same as += but deviding by the new entry? 
		sumxx[k] = sqrt(fabs(sumxx[k] / nColony - sqr(sumx[k])));  // average? 
	}

	//export data
	ofs << t;													//to fill up the text or csv file 
	for (int k = 0; k < 3; ++k)
		ofs << ofs.fill() << ofs.fill() << sumx[k] << ofs.fill() << sumxx[k];
	for (int cs = 0, k = 3; cs < nCaste; ++cs) 
		for (int tr = 0; tr < nTrait; ++tr) {
			ofs << ofs.fill();
			for (int i = 0; i < ageMax; ++i, ++k) ofs << ofs.fill() << sumx[k]; // to me its not clear where stuff if filled or written into.. 
		}
	for (int cs = 0, k = 3; cs < nCaste; ++cs)
		for (int tr = 0; tr < nTrait; ++tr) {
			ofs << ofs.fill();
			for (int i = 0; i < ageMax; ++i, ++k) ofs << ofs.fill() << sumxx[k]; // this gives the ages so it should get the age class specific output in the csv file 
		}

	ofs << '\n';
	ofs.flush();
	std::cout << t << ' ' << sumx[0] << ' ' << sumx[1] << ' ' << sumx[2] << '\n'; // i guess the second entry should be sumx[1] instead of [0] also the output is inly after each 10th simulation round where do i find that? 
}


		////// get data our f the model with help from richel 20160726

		//Get the average Matrix out of a collection of matrices called 'qms'
Matrix getAverage(const std::vector<Matrix>& matrices)
{
	if (matrices.empty())
	{
		throw std::invalid_argument("getAverage: input vector empty");// programm will be terminated if the matrices are empty
	}
	const Matrix first_qm = matrices[0];
	Matrix empty_matrix(first_qm.getRows().size(), first_qm.getCols().size(), 0.0);
	const Matrix sum_qms = std::accumulate( //Sum of the Queen life hostory MatriceS
		matrices.begin(), matrices.end(), empty_matrix
		);
	assert(matrices.size() != 0);
	const Matrix mean_qms = sum_qms / static_cast<double>(matrices.size());
	return mean_qms;
}


//Extract all life history matrices from all queens
std::vector<Matrix> extractLifeHistoryQueens(const Colony * const colonies[], const int sz)
{
	std::vector<Matrix> qms; //Queen Matrices
	qms.reserve(sz);
	for (int i = 0; i != sz; ++i)
	{
		const Colony& c = *colonies[i]; // (this/focal) Colony
		const Matrix& qm = c.getQueen()->getLifeHistory(); // Queen Matrix
		qms.push_back(qm);
	}
	return qms;
}

///Extract all workers from all colonies
std::vector<Worker *> extractWorkers(const Colony * const colonies[], const int sz)
{
	std::vector<Worker *> all_workers;
	for (int i = 0; i != sz; ++i)
	{
		const Colony& c = *colonies[i]; // (this/focal) Colony
		const std::vector<Worker*> this_workers = c.getWorkers(); //Workers in this focal colony
		assert(this_workers.size() == c.size());

		//Append all this colony its workers to the already collected workers
		std::copy(std::begin(this_workers), std::end(this_workers), std::back_inserter(all_workers));

	}
	return all_workers;

}

///Extract all life history matrices from all workers
std::vector<Matrix> extractLifeHistoryWorkers(const Colony * const colonies[], const int sz)
{
	const std::vector<Worker*> all_workers = extractWorkers(colonies, sz);
	const int n_workers = static_cast<int>(all_workers.size());

	std::vector<Matrix> wms; //worker Matrices
	wms.reserve(n_workers);

	for (int i = 0; i != n_workers; ++i)
	{
		const Matrix& wm = all_workers[i]->getLifeHistory(); // Worker Matrix
		wms.push_back(wm);
	}
	return wms;
}

void saveLifeHistoryQueens(const int generation, const Colony * const colonies[], const int sz, std::ofstream& file) // i dont understand why all the ofstream things go into the same file event though its calles offs or file , how do i manage to get two or three different files ? 
{
	const std::vector<Matrix> qms = extractLifeHistoryQueens(colonies, sz);
	const Matrix mean_qms = getAverage(qms);
	file << generation << " , 1 ," << ", queens, " << mean_qms << '\n';
}

void saveLifeHistoryWorkers(const int generation, const Colony * const colonies[], const int sz, std::ofstream& file)
{
	const std::vector<Matrix> wms = extractLifeHistoryWorkers(colonies, sz);
	if (wms.empty())
	{
		file <<  generation <<" , 3 , " << " , worker, " << "NONE" << '\n';
	}
	else
	{
		const Matrix mean_wms = getAverage(wms);
		file << generation << " , 2 , " << " , worker, " << mean_wms << '\n';
	}
}

void saveLifeHistory(const int generation, const Colony * const colonies[], const int sz, std::ofstream& file)
{
	saveLifeHistoryQueens(generation, colonies, sz, file);
	saveLifeHistoryWorkers(generation, colonies, sz, file); //Extract all life history matrices from all workers

}


std::string constructFilename(const long long simulationId)
{
	std::stringstream s;
	s << "data_" << simulationId << ".csv"; // would also be nice to create a nice header line for the files written, one problem here is that  information is printed during the whole simulation also id like two or more seperate files to store the output
	return s.str();
}

///Calculates the dead workers per age class from each population, then sums these per age class
std::vector<int> calcSumDeadWorkers(const Colony * const colonies[], const int sz)
{
	//Create a vector hold the sum per age class
	std::vector<int> sum(ageMax, 0);
	for (int i = 0; i != sz; ++i)
	{
		const Colony& c = *colonies[i];
		const std::vector<int>&  dead_workers = c.getDeadWorkers();
		assert(sum.size() == dead_workers.size());
		//Sum per age class
		for (int age_class_index = 0; age_class_index != ageMax; ++age_class_index)
		{
			assert(age_class_index >= 0);
			assert(age_class_index < static_cast<int>(sum.size()));
			assert(age_class_index < static_cast<int>(dead_workers.size()));
			sum[age_class_index] += dead_workers[age_class_index];
		}
	}
	return sum;
}

///Calculates the dead queens per age class from each population, then sums these per age class
std::vector<int> calcSumDeadQueens(const Colony * const colonies[], const int sz)
{
	//Create a vector hold the sum per age class
	std::vector<int> sum(ageMax, 0);
	for (int i = 0; i != sz; ++i)
	{
		const Colony& c = *colonies[i];
		const std::vector<int>&  dead_queens = c.getDeadQueens();
		assert(sum.size() == dead_queens.size());
		//Sum per age class
		for (int age_class_index = 0; age_class_index != ageMax; ++age_class_index)
		{
			assert(age_class_index >= 0);
			assert(age_class_index < static_cast<int>(sum.size()));
			assert(age_class_index < static_cast<int>(dead_queens.size()));
			sum[age_class_index] += dead_queens[age_class_index];
		}
	}
	return sum;
}

///Converts a vector to a comma-seperated string
std::string vectorToString(const std::vector<int>& v)
{
	std::stringstream s;
	std::copy(std::begin(v), std::end(v), std::ostream_iterator<int>(s, ","));
	return s.str();
}


/*=============================================================================================================
                                                   main()
=============================================================================================================*/
int main()    // the designated start of the program
{
	//preliminaries
	simulationId = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	rnd::set_seed(static_cast<unsigned int>(simulationId));
	const std::string fileName = constructFilename(simulationId);
	ofs = std::ofstream(fileName);
	if (!ofs.is_open()) throw std::runtime_error("Cannot open file");
	ofs.fill(','); // commas between entries? how does c++ know where to write the data and even more interesting can this also be used in a loop to grab the data again.. 
	init();

	int generation = 0;
	std::vector<double> sumx, sumxx; // create the vector to store the data
	std::vector<int> sumDeadQueens(ageMax, 0);
	std::vector<int> sumDeadWorkers(ageMax, 0);

	while (iterate(generation, ofs, sumx, sumxx, sumDeadWorkers, sumDeadQueens) && generation < generationEnd) {  // run the sumulation as often as denined by generationEnd
		if (generation % dataInterval == 0) {  //cout only after each data interval
			dataAnalysis(generation, sumx, sumxx); // dont understand what this is doing
			//saveLifeHistory(generation, population, nColony, ofs);
		//}
		//if (generation == generationEnd - 1) {
			//Calculates the dead workers per age class
			//const std::vector<int> sumDeadWorkers = calcSumDeadWorkers(population, nColony); 
			//Copy that to the file
			ofs << generation << " , 3 , " << " , worker, " << vectorToString(sumDeadWorkers) << '\n';
			//Calculates the dead queens per age class
			//const std::vector<int> sumDeadQueens = calcSumDeadQueens(population, nColony);
			ofs << generation << " , 4 , " << " , queen, " << vectorToString(sumDeadQueens) << '\n';
			//Copy that to the file
		}
		++generation;
	}
	ofs.close(); //if this is on the file will be written automatically 
	wait_for_return(); // from the keyboard to write the files also....  
}

