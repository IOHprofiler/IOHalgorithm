//
//  main.cpp
//  ConfigurableGA
//
//  Created by Furong Ye on 26/10/2020.
//  Copyright © 2020 Furong Ye. All rights reserved.
//

#include "modularGA.hpp"
#include "ioh.hpp"

#include <vector>

using namespace std;
/// \fn std::string strstrip(std::string s)
/// \brief To erase blanks in the front and end of a string.
///
/// \para s string
/// \return string
static std::string strstrip(std::string s)
{
  if (s.empty()) {
      return s;
  }
  s.erase(0,s.find_first_not_of(' '));
  /// For string line end with '\r' in windows systems.
  s.erase(s.find_last_not_of('\r') + 1);
  s.erase(s.find_last_not_of(' ') + 1);
  return s;
}

/// \fn std::vector<int> get_int_vector_parse_string(std::string input, const int _min, const int _max) {
/// \brief To get a vector of integer values with a range 'n-m'
///
/// Return a vector of integer value. The supported formats are:
///   '-m', [_min, m]
///   'n-', [n, _max]
///   'n-m', [n, m]
///   'n-x-y-z-m', [n,m]
/// \para input string
/// \para _int int
/// \para _max int
/// \return std::vector<int>
static std::vector<int> get_int_vector_parse_string(std::string input, const int _min, const int _max) {
  std::vector<std::string> spiltstring;
  std::string tmp;
  int tmpvalue,tmpvalue1;
  std::vector<int> result;

  size_t n = input.size();
  input = strstrip(input);
  for (size_t i = 0; i < n; ++i) {
    if (input[i] != ',' && input[i] != '-' && !isdigit(input[i])) {
      cerr << "The configuration consists invalid characters.";
    }
  }
  
  std::stringstream raw(input);
  while(getline(raw, tmp, ',')) spiltstring.push_back(tmp);

  n = spiltstring.size();
  for (size_t i = 0; i < n; ++i) {
    //size_t l = spiltstring[i].size();

    if (spiltstring[i][0] == '-') {
      /// The condition beginning with "-m"
      if(i != 0) {
        cerr << "Format error in configuration.";
      } else {
        tmp = spiltstring[i].substr(1);
        if (tmp.find('-') != std::string::npos) {
          cerr << "Format error in configuration.";
        }
        
        tmpvalue = std::stoi(tmp);
        
        if (tmpvalue < _min) {
          cerr << "Input value exceeds lowerbound.";
        }

        for (int value = _min; value <= tmpvalue; ++value) {
          result.push_back(value);
        }
      }
    } else if(spiltstring[i][spiltstring[i].length()-1] == '-') {
      /// The condition endding with "n-"
      if (i != spiltstring.size() - 1) {
        cerr << "Format error in configuration.";
      } else {
        tmp = spiltstring[i].substr(0,spiltstring[i].length() - 1);
        if (tmp.find('-') != std::string::npos) {
          cerr << "Format error in configuration.";
        }
        tmpvalue = std::stoi(tmp);
        if (tmpvalue > _max) {
          cerr << "Input value exceeds upperbound.";
        }
        for (int value = _max; value <= tmpvalue; --value) {
          result.push_back(value);
        }
      }
    } else {
      /// The condition with "n-m,n-x-m"
      std::stringstream tempraw(spiltstring[i]);
      std::vector<std::string> tmpvaluevector;
      while (getline(tempraw, tmp, '-')) {
        tmpvaluevector.push_back(tmp);
      }
      tmpvalue = std::stoi(tmpvaluevector[0]);
      tmpvalue1 = std::stoi(tmpvaluevector[tmpvaluevector.size()-1]);
      if (tmpvalue > tmpvalue1) {
        cerr << "Format error in configuration.";
      }
      if (tmpvalue < _min)  {
        cerr << "Input value exceeds lowerbound.";
      }
      if (tmpvalue1 > _max) {
        cerr << "Input value exceeds upperbound.";
      }
      for(int value = tmpvalue; value <= tmpvalue1; ++value) result.push_back(value);
    } 
  }
  return result;
}
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
 *    sa : simulated annealing
 *    sars : the simulated annealing algorithm with exponential temperature schedule with iterative restarts
 *    fea : (1+1)-EA>0 with frequency fitness assignment
 */
void runAlgorithm(shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite, const string algorithm_name, const string dir, const int budget, const int runs, const unsigned seed)
{
  if (algorithm_name == "ea")
  {
  modularGA::ga::staticEA ea(1, 1, 1.0);
    ea.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  }
  else if (algorithm_name == "ea2")
  {
   modularGA::ga::staticEA ea(1, 1, 2.0);
    ea.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  }
  else if (algorithm_name == "ea23")
  {
    modularGA::ga::staticEA ea(1, 1, 1.5);
    ea.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  }
  else if (algorithm_name == "llea")
  {
    modularGA::ga::oneLambdaLambdaEA llEA(10);
    llEA.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  }
  else if (algorithm_name == "rls")
  {
    modularGA::ga::RLS rls;
    rls.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  }
  else if (algorithm_name == "rs")
  {
    modularGA::rs::RandomSearch rs;
    rs.run(dir, algorithm_name, suite, budget, runs, seed);
  }
  else if (algorithm_name == "ghc")
  {
    modularGA::ghc::GreedyHillClimber ghc;
    ghc.run(dir, algorithm_name, suite, budget, runs, seed);
  }
  else if (algorithm_name == "fga")
  {
    modularGA::ga::FastGA fGA(1, 1);
    fGA.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  }
  else if (algorithm_name == "2ratega")
  {
    modularGA::ga::TwoRateEA twoRateGA(10);
    twoRateGA.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  }
  else if (algorithm_name == "umda")
  {
    modularGA::eda::EstimationOfDistribution umda(25, 50);
    umda.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  }
  else if (algorithm_name == "normea")
  {
    modularGA::ga::NormEA normEA(10);
    normEA.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  }
  else if (algorithm_name == "varea")
  {
    modularGA::ga::VarEA varEA(10);
    varEA.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);
  }
  else if (algorithm_name == "sa")
  {
    modularGA::sa::run_simulated_annealing_exp(dir, suite, budget, runs, seed);
  }
  else if (algorithm_name == "sars")
  {
    modularGA::sars::run_simulated_annealing_exp_rs(dir, suite, budget, runs, seed);
    // }  else if (algorithm_name == "fea") {
    // run_fea1p1(dir, suite, budget, runs, seed);
  }
  else
  {
    cout << "Unkownn algorithm : " << algorithm_name << endl;
  }
}

/**
 * Arguments: suite_name problem_id instance_id dimension
 *  algorithm_name : 'name'
 *  suite_name : 'name', name \in {'PBO','WmodelOneMax', 'WmodelLeadingOnes'}
 *  problem_id : 'start-end' , id \in [1..25] for PBO, [1..product.size()] for Wmodel extensions.
 *  instance_id: 'start-end' , id \in [1..100]
 *  dimension :  'dimension1, dimension2, ...'
 *  dir : 'dir', the directory where the output folder locates.
 *  runs : 'number_of_indepedent_runs', runs > 1
 *  budget : 'budget', the maximum function evaluations. budget > 1
 *  seed : 'seed', a random seed
 *
 * An instance: ./main ea pbo 1-3 1-5 10,100 ./ 10 100 1
 **/

int main(int argc, const char *argv[])
{
  const vector<double> dummy = {0.0, 0.9};
  const vector<int> epistasis = {0, 2, 5, 7};
  const vector<int> neutrality = {1, 5};
  const vector<double> ruggedness = {0, 0.8, 1};
  int number_of_w_problems = dummy.size() * epistasis.size() * neutrality.size() * ruggedness.size();

  string algorithm_name = argv[1];
  string suite_name = argv[2];
  string problem_str = argv[3];
  string instance_str = argv[4];
  string dimension_str = argv[5];
  string dir = argv[6];
  int runs = stoi(argv[7]);
  int budget = stoi(argv[8]);
  unsigned seed =  static_cast<unsigned> (stoi(argv[9]));

  transform(suite_name.begin(),suite_name.end(),suite_name.begin(),::tolower);
  transform(algorithm_name.begin(),algorithm_name.end(),algorithm_name.begin(),::tolower);
  auto problems = get_int_vector_parse_string(problem_str,0,25);
  auto instances = get_int_vector_parse_string(instance_str,0,100);
  auto dimensions = get_int_vector_parse_string(dimension_str,2,20000);

  if (suite_name == "pbo") {
    auto suite = ioh::suite::SuiteRegistry<ioh::problem::IntegerSingleObjective>::instance()
        .create("PBO", problems , instances, dimensions);
    runAlgorithm(suite, algorithm_name, dir, budget, runs, seed);
  } else {
    cout << "Unknown suite : " << suite_name << ", avaliable options are \"PBO\", \"WModelOneMax\", and \"WModelLeadingOnes\"." << endl;
  }
  return 0;
}