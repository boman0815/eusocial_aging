#ifndef parameters_h   // checks if parameters are already included
#define parameters_h

class Individual;

// extrinsic risk (boris 20160712) 
const double extrinsicMortality = 0.15;
const double eggs = 2.0;
//double productivity = 0.8;


//colony dynamics
const int nColony = 1000;									//number of colonies
const int generationEnd = 8000;								//simulation time	
const int dataInterval = 10;								// what is the data interval? is it the ages? 	

//life history parameters
const int ageMax = 20;                                      // maximum lifespan of 10 (simulation cycles ?) one gene for each age to determine fertility and survival
const double fecundityOffset        = 1.0;
const double productivityOffset     = 1.0;
const double survivalOffset         = 1.0;
const double costOfForaging         = 1.0;                 //what does this parameter do exactly? 

//genetics
const int nGene = 40;                                       //number of genes; has to be even
const int nCaste = 2;                                       //queens and workers
const int nTrait = 1;                                       //productivity/fecundity and survival
const double mutationRate = 3.0e-5;                         //rate of mutations per gene per generation
const double sigma[nCaste][nTrait] =                        //mutational effect size initial values: {{0.05, 0.05}, {0.05, 0.05}};
    {{0.08}, {0.08}};
const double alpha[nTrait] = {0.0};							//correlation between subsequent age classes within trait initial values: {-0.05, 0.05}
const double beta[nTrait][(nCaste * (nCaste - 1)) / 2] =    //genetic correlation within trait, between castes initial values: {{-0.3},{0.3}}
    {{0.0}};
//const double ggamma[nCaste][(nTrait * (nTrait - 1)) / 2] =  //antagonistic pleiotropy; correlation between trait, within caste initial values:  {{-0.4},{-0.4}}
//    {{-0.1},{-0.1}};
//const double delta = 0.1;                                   //antagonistic pleiotropy delay initial values: 0.5
const double eta = -0.2;                                    //mutation bias (relative to mutational effect size) initial values: -0.2   shifting the mutational effect  if negative mutational effect will be drawn from a distribution with a negative mean, assuming tha mutatons will usually have negative effects that selection has to minimize for the selection mutation balance 
#endif