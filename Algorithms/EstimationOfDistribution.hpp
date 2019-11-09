/// \file EstimationOfDistribution.hpp
/// \brief implementation of estimation of distribution algorithm
///
/// \author Furong Ye
/// \date 2019-07-26

#ifndef _ESTIMATIONOFDISTRIBUTION_HPP
#define _ESTIMATIONOFDISTRIBUTION_HPP

#include "CommonModule.hpp"

extern int budget;
extern int budget_scale;
extern int budget_power;

extern int DEFAULT_POPULATION_SIZE;
extern double DEFAULT_SELECT_RATE;

/// \fn void estimationOfDistribution(std::shared_ptr<IOHprofiler_problem<int>> problem, std::shared_ptr<IOHprofiler_csv_logger> logger) 
/// Univariate marginal distribution algorithm.
void estimationOfDistribution(std::shared_ptr<IOHprofiler_problem<int>> problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
  int population_size = DEFAULT_POPULATION_SIZE;
  double select_rate = DEFAULT_SELECT_RATE;
  int count = 0;
  int last_found_evaluation = 0;
  int n = problem->IOHprofiler_get_number_of_variables();
  budget = budget_scale * pow((double)(n),budget_power);
  std::vector<std::vector<int> > population(population_size);
  std::vector<double> fitness;
  std::vector<std::vector<double> > distribution;
  std::vector<int> range(n,2);
  std::vector<std::vector<int> > selected_population;

  /***
   * Recording Parameters Preparation
   ***/
  std::shared_ptr<double> best_found_fitness = std::make_shared<double>(-DBL_MAX); /// <  last_found_fitness
  std::shared_ptr<double> update_evaluations = std::make_shared<double>(0.0); /// <  evaluations times since last best fitness has been found
  std::shared_ptr<double> cross_points = std::make_shared<double>(2.0); /// <  evaluations times since last best fitness has been found
  std::vector<std::shared_ptr<double> > parameters;
  parameters.push_back(best_found_fitness);
  parameters.push_back(update_evaluations);
  parameters.push_back(cross_points);
  std::vector<std::string> parameters_name;
  parameters_name.push_back("best_found_fitness");
  parameters_name.push_back("update_evaluations");
  parameters_name.push_back("cross_points");
  logger->set_parameters(parameters,parameters_name);

  for (int i = 0; i != population_size; ++i) {
  initialization(population[i],n,2);
    fitness.push_back(problem->evaluate(population[i]));
    count++;
    if (fitness[i] > *best_found_fitness) {
      *best_found_fitness = fitness[i];
      *update_evaluations = (double)(count - last_found_evaluation);
      last_found_evaluation = count;
    }
    logger->write_line(problem->loggerInfo());
  }

  
  while(count < budget && !problem->IOHprofiler_hit_optimal()) {
    selection(selected_population,population,fitness,select_rate);
    calculateDistribution(distribution,selected_population,range);
    for (int i = 0; i != population_size; ++i) {
      population[i] = sampling(distribution);
      fitness[i] = problem->evaluate(population[i]);
      count++;
      if (fitness[i] > *best_found_fitness) {
        *best_found_fitness = fitness[i];
        *update_evaluations = (double)(count - last_found_evaluation);
        last_found_evaluation = count;
      }
      logger->write_line(problem->loggerInfo());
    }
  }
};
#endif

