//
//  main.cpp
//  ConfigurableGA
//
//  Created by Furong Ye on 26/10/2020.
//  Copyright © 2020 Furong Ye. All rights reserved.
//

#include "../algorithm/ga/instance/ea.h"
#include "../algorithm/ga/instance/RLS.h"
#include "../algorithm/ga/instance/fastGA.h"
#include "../algorithm/ga/instance/oneLLEA.h"
#include "../algorithm/ga/instance/twoRateEA.h"
#include "../algorithm/ghc/ghc.h"
#include "../algorithm/eda/estimationOfDistribution.h"
#include "../algorithm/rs/randomsearch.h"

#include "IOHprofiler_string.hpp"
#include "IOHprofiler_PBO_suite.hpp"

#include <vector>

using namespace std;

/**
 * Alrogrithm names:
 *    ea : (1+1)-EA>0 p = 1/n
 *    ea2 : (1+1)-EA>0 p = 2/n
 *    ea23 : (1+1)-EA>0 p = 2/3n
 *    llea  : (1+(lambda,lambda))-EA>_0, initial lambda = 10, p = 1/n
 *    rls : randomized local search
 *    rs : random search
 *    ghc : greedy hill climber
 *    fga : (1+1)-fast GA
 *    2ratega : (1+10)-EA>0 with 2rate self adaptation of mutation rate
 *    umda : univarate marginal distribution algorithm, population size = 50
 */
void runAlgorithm(shared_ptr< IOHprofiler_suite<int> > suite, const string algorithm_name, const string dir, const int budget, const int runs, const unsigned seed)
{
  if (algorithm_name == "ea") {
    staticEA ea(1, 1, 1.0);
    ea.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  } else if (algorithm_name == "ea2") {
    staticEA ea(1, 1, 2.0);
    ea.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  } else if (algorithm_name == "ea23") {
    staticEA ea(1, 1, 1.5);
    ea.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  } else if (algorithm_name == "llea") {
    oneLambdaLambdaEA llEA(1);
    llEA.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  } else if (algorithm_name == "rls") {
    RLS rls;
    rls.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  } else if (algorithm_name == "rs") {
    RandomSearch rs;
    rs.run(dir, algorithm_name, suite, budget,runs, seed);
  } else if (algorithm_name == "ghc") {
    GreedyHillClimber ghc;
    ghc.run(dir, algorithm_name, suite, budget, runs, seed);
  } else if (algorithm_name == "fga") {
    FastGA fGA(1, 1);
    fGA.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  } else if (algorithm_name == "2ratega") {
    TwoRateEA twoRateGA(10);
    twoRateGA.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  } else if (algorithm_name == "umda") {
    EstimationOfDistribution umda(25,50);
    umda.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  } else {
    cout << "Unkownn algorithm : " << algorithm_name << endl;
  }
}

int main(int argc, const char *argv[])
{
  string algorithm_name = "ea"; /// {ea, ea2, ea23, llea, rls, rs, ghc, fga, 2ratega, umda}
  string dir = "./"; /// The folder where to store output
  int runs = 100; /// The number of independent runs
  int budget = 10000; /// The maximal evaluation times for each run
  unsigned seed =  666; /// Random Seed
  
  transform(algorithm_name.begin(),algorithm_name.end(),algorithm_name.begin(),::tolower);
  
  vector<int> problem_id;
  for (int i = 1; i <= 23; i++) {
    problem_id.push_back(i);
  }
  vector<int> instance_id = {1};
  vector<int> dimension = {16,64,100};
  shared_ptr<PBO_suite> suite(new PBO_suite(problem_id, instance_id, dimension));
  runAlgorithm(suite, algorithm_name, dir, budget, runs, seed);
  
  return 0;
}
