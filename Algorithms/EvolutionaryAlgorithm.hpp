/// \file EvolutionaryAlgorithm.hpp
/// \brief implementation of several evolutionary algorithms
///
/// A detail file description. This file includes implementation of evolutionary algorithms
/// with static mutation rate, adaptive mutation rate, adaptive lambda, vanilla genetic algorithm,
/// and genetic algorithm with self adaptive mutation rate and crossover rate.
///
/// \author Furong Ye
/// \date 2019-07-30

#ifndef _EVOLUTIONARYALGORITHM_HPP
#define _EVOLUTIONARYALGORITHM_HPP

#include "CommonModule.hpp"

extern int budget;
extern int budget_scale;
extern int budget_power;

extern int DEFAULT_MU;
extern int DEFAULT_LAMBDA;
extern int DEFAULT_ML;

/// \fn void evolutionaryAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger)
/// A (1+lambda) evolutionary algorithm with static mutation rate 1/n.
void evolutionaryAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
  int n = problem->IOHprofiler_get_number_of_variables(); /// < dimension
  int lambda = DEFAULT_LAMBDA;
  double mutation_rate = (double)DEFAULT_ML / (double)n;
  std::vector<int> parent;
  std::vector<int> x, x_star; /// < temp individuals
  double y; /// < temp fitness
  int count;
  double best_fitness;
  int last_found_evaluation;

  budget = budget_scale * pow((double)(n),budget_power);

  /***
   * Recording Parameters Preparation
   ***/
  std::shared_ptr<double> l = std::make_shared<double>(1.0); /// < mutation length
  std::shared_ptr<double> update_evaluations = std::make_shared<double>(0.0); /// <  evaluations times since last best fitness has been found
  std::vector<std::shared_ptr<double> > parameters;
  parameters.push_back(l);
  parameters.push_back(update_evaluations);;
  std::vector<std::string> parameters_name;
  parameters_name.push_back("l");
  parameters_name.push_back("update_evaluations");
  logger->set_parameters(parameters,parameters_name);

  
  initialization(parent,n,2);
  best_fitness = problem->evaluate(parent);
  logger->write_line(problem->loggerInfo());

  count = 1;
  last_found_evaluation = 1;
  while(count < budget && !problem->IOHprofiler_hit_optimal()) {
    
    x = parent;
    
    for (int i = 0; i != lambda; ++i) {
      *l = (double)sampleLRand(mutation_rate,n);
      x_star = x;
      bitMutation(x_star,(int)(*l+0.5));
      y = problem->evaluate(x_star);
      count++;
      if (y > best_fitness) {
        *update_evaluations = (double)(count - last_found_evaluation);
        last_found_evaluation = count;
      }
      logger->write_line(problem->loggerInfo());

      if (count == budget || problem->IOHprofiler_hit_optimal()) {
        break;
      }

      if (y >= best_fitness) {
        best_fitness = y;
        parent = x_star;
      }

    }
  }
}

/// \fn void twoRateEvolutionaryAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger)
/// (1+lambda) EA with 2-rate self-adaptive mutation rate.
void twoRateEvolutionaryAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
  int n = problem->IOHprofiler_get_number_of_variables(); /// < dimension
  int lambda = DEFAULT_LAMBDA;
  double mutation_rate1, mutation_rate2;
  double r = 2.0; /// < a factor using to adjust mutation rate
  std::vector<int> parent;
  std::vector<int> x, x_star; /// < temp individuals
  double y; /// < temp fitness
  double best_fitness;
  int count;
  double tempb;
  double bestr = 1;
  int last_found_evaluation;

  budget = budget_scale * pow((double)(n),budget_power);
  
  /***
   * Recording Parameters Preparation
   ***/
  std::shared_ptr<double> l = std::make_shared<double>(1.0); /// < mutation length
  std::shared_ptr<double> update_evaluations = std::make_shared<double>(0.0); /// <  evaluations times since last best fitness has been found
  std::vector<std::shared_ptr<double> > parameters;
  parameters.push_back(l);
  parameters.push_back(update_evaluations);;
  std::vector<std::string> parameters_name;
  parameters_name.push_back("l");
  parameters_name.push_back("update_evaluations");
  logger->set_parameters(parameters,parameters_name);

  initialization(parent,n,2);
  best_fitness = problem->evaluate(parent);
  logger->write_line(problem->loggerInfo());
  
  count = 1;
  last_found_evaluation = 1;
  while(count < budget && !problem->IOHprofiler_hit_optimal()) {
    x = parent;
    
    /// Preparing two mutation rates and corresponding probability density for sampling mutation length
    mutation_rate1 = (double)r / 2.0 / (double)n;
    mutation_rate2 = 2.0 * (double)r / (double)n;
    
    tempb = -DBL_MAX;
    for (int i = 0; i < lambda; ++i) {
      
      if (i < lambda / 2) {
        *l = (double)sampleLRand(mutation_rate1,n);
      } else {
        *l = (double)sampleLRand(mutation_rate2,n);
      }
      
      x_star = x;
      bitMutation(x_star,(int)(*l+0.5));
      y = problem->evaluate(x_star);
      count++;
      if (y > best_fitness) {
        *update_evaluations = (double)(count - last_found_evaluation);
        last_found_evaluation = count;
      }
      logger->write_line(problem->loggerInfo());
      
      if (count == budget || problem->IOHprofiler_hit_optimal()) {
         break;
      }

      if(y >= best_fitness) {
        best_fitness = y;
        parent = x_star;
      }

      if(y > tempb) {
        tempb = y;
        if(i < lambda / 2) {
          bestr = r / 2.0;
        } else {
          bestr = r * 2.0;
        }
      }
    }
    
    if (count == budget || problem->IOHprofiler_hit_optimal()) {
      break;
    }

    if(random_generator.IOHprofiler_uniform_rand() > 0.5) {
      r = bestr;
    } else {
      if(random_generator.IOHprofiler_uniform_rand() > 0.5) {
        r = r / 2.0;
      } else {
        r = 2.0 * r;
      }
    }

    r = r < 2.0 ? 2.0 : r;
    r = r > (n / 4.0) ? (n / 4.0) : r;
  }
}

/// \fn void logNromalEvolutionaryAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger)
/// (1+lambda) EA with logNormal self adaptive muataion rate
void logNromalEvolutionaryAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
  int n = problem->IOHprofiler_get_number_of_variables(); /// < dimension
  int lambda = DEFAULT_LAMBDA;
  double mutation_rate = 0.2;
  std::vector<int> parent;
  std::vector<int> x, x_star; /// < temp individuals
  double y; /// < temp fitness
  double tempb = -DBL_MAX;
  double temp_mutation_rate, temp_best_mutation_rate = 0.2;
  int count;
  int last_found_evaluation;

  budget = budget_scale * pow((double)(n),budget_power);

  /***
   * Recording Parameters Preparation
   ***/
  std::shared_ptr<double> l = std::make_shared<double>(1.0); /// < mutation length
  std::shared_ptr<double> update_evaluations = std::make_shared<double>(0.0); /// <  evaluations times since last best fitness has been found
  std::vector<std::shared_ptr<double> > parameters;
  parameters.push_back(l);
  parameters.push_back(update_evaluations);;
  std::vector<std::string> parameters_name;
  parameters_name.push_back("l");
  parameters_name.push_back("update_evaluations");
  logger->set_parameters(parameters,parameters_name);
  
  initialization(parent,n,2);
  double best_fitness = problem->evaluate(parent);
  logger->write_line(problem->loggerInfo());

  count = 1;
  last_found_evaluation = 1;
  while(count < budget && !problem->IOHprofiler_hit_optimal()) {
    tempb = -DBL_MAX;
    mutation_rate = temp_best_mutation_rate;
    x = parent;
    for (int i = 0; i < lambda; ++i) {
      temp_mutation_rate = 1.0 / (1.0 + (((1.0 - mutation_rate) / mutation_rate) * exp(0.22*random_generator.IOHprofiler_normal_rand())));
      *l = (double)sampleLRand(temp_mutation_rate + 1.0 / (double)n,n);
      
      x_star = x;
      bitMutation(x_star,(int)(*l+0.5));
      y = problem->evaluate(x_star);
      count++;
      if (y > best_fitness) {
        *update_evaluations = (double)(count - last_found_evaluation);
        last_found_evaluation = count;
      }
      logger->write_line(problem->loggerInfo());

      if (count == budget || problem->IOHprofiler_hit_optimal()) {
        break;
      }

      if(y >= best_fitness) {
        best_fitness = y;
        parent = x_star;
        temp_best_mutation_rate = temp_mutation_rate;
      }

      if(y > tempb) {
        tempb = y;
        
      }
    }
  }
};

/// \fn void normalizedEvolutionaryAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger)
/// (1+lambda) EA with normalized bit mutation
void normalizedEvolutionaryAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
  int n = problem->IOHprofiler_get_number_of_variables(); /// < dimension
  int lambda = DEFAULT_LAMBDA;
  double mutation_rate = 1.0 / n;
  std::vector<int> parent;
  std::vector<int> x, x_star;
  double y;
  double tempb = -DBL_MAX;
  double bestl = 1, temp_best_l;
  int count;
  int last_found_evaluation;

  budget = budget_scale * pow((double)(n),budget_power);

  /***
   * Recording Parameters Preparation
   ***/
  std::shared_ptr<double> l = std::make_shared<double>(1.0); /// < mutation length
  std::shared_ptr<double> sigma = std::make_shared<double>(sqrt(1 * (1.0- 1.0/(double)n)));
  std::shared_ptr<double> update_evaluations = std::make_shared<double>(0.0); /// <  evaluations times since last best fitness has been found
  std::vector<std::shared_ptr<double> > parameters;
  parameters.push_back(l);
  parameters.push_back(sigma);
  parameters.push_back(update_evaluations);;
  std::vector<std::string> parameters_name;
  parameters_name.push_back("l");
  parameters_name.push_back("sigma");
  parameters_name.push_back("update_evaluations");
  logger->set_parameters(parameters,parameters_name);
  
  
  initialization(parent,n,2);
  double best_fitness = problem->evaluate(parent);
  logger->write_line(problem->loggerInfo());
  
  count = 1;
  last_found_evaluation = 1;
  while(count < budget && !problem->IOHprofiler_hit_optimal()) {
    tempb = -DBL_MAX;
    temp_best_l = bestl;
    x = parent;
    for (int i = 0; i < lambda; ++i) {

      *l = (double)normalRand(bestl,*sigma,n);

      x_star = x;
      bitMutation(x_star,(int)(*l+0.5));

      y = problem->evaluate(x_star);
      count++;
      if (y > best_fitness) {
        *update_evaluations = (double)(count - last_found_evaluation);
        last_found_evaluation = count;
      }
      logger->write_line(problem->loggerInfo());

      if (count == budget || problem->IOHprofiler_hit_optimal()) {
        break;
      }
      
      if(y >= best_fitness) {
        best_fitness = y;
        parent = x_star;
      }

      if(y > tempb) {
        tempb = y;
        temp_best_l = *l;
      }
    }

    if (count == budget || problem->IOHprofiler_hit_optimal()) {
      break;
    }

    bestl = temp_best_l;
    *sigma = sqrt(bestl * (1.0- (double)bestl/(double)n));
  }
};

/// \fn void varianceControlledEvolutionaryAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger)
/// (1+lambda) EA with normalized bit mutation using variance control
void varianceControlledEvolutionaryAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
  int n = problem->IOHprofiler_get_number_of_variables(); /// < dimension
  int lambda = DEFAULT_LAMBDA;
  double mutation_rate = 1.0 / n;
  std::vector<int> parent;
  std::vector<int> x, x_star;
  double y;
  double tempb = -DBL_MAX;
  double bestl = 1, temp_best_l;
  int count;
  int bestl_count = 0;
  int last_found_evaluation;

  budget = budget_scale * pow((double)(n),budget_power);

  /***
   * Recording Parameters Preparation
   ***/
  std::shared_ptr<double> l = std::make_shared<double>(1.0); /// < mutation length
  std::shared_ptr<double> sigma = std::make_shared<double>(sqrt(1 * (1.0- 1.0/(double)n)));
  std::shared_ptr<double> update_evaluations = std::make_shared<double>(0.0); /// <  evaluations times since last best fitness has been found
  std::vector<std::shared_ptr<double> > parameters;
  parameters.push_back(l);
  parameters.push_back(sigma);
  parameters.push_back(update_evaluations);;
  std::vector<std::string> parameters_name;
  parameters_name.push_back("l");
  parameters_name.push_back("sigma");
  parameters_name.push_back("update_evaluations");
  logger->set_parameters(parameters,parameters_name);


  initialization(parent,n,2);
  double best_fitness = problem->evaluate(parent);
  logger->write_line(problem->loggerInfo());
  
  count = 1;
  last_found_evaluation = 1;
  while(count < budget && !problem->IOHprofiler_hit_optimal()) {
    tempb = -DBL_MAX;
    temp_best_l = bestl;
    x = parent;
    for (int i = 0; i < lambda; ++i) {

      *l = (double)normalRand(bestl,*sigma,n);

      x_star = x;
      bitMutation(x_star,(int)(*l+0.5));

      y = problem->evaluate(x_star);
      count++;

      if (y > best_fitness) {
        *update_evaluations = (double)(count - last_found_evaluation);
        last_found_evaluation = count;
      }
      logger->write_line(problem->loggerInfo());
      
      if (count == budget || problem->IOHprofiler_hit_optimal()) {
        break;
      }
      
      if(y >= best_fitness) {
        best_fitness = y;
        parent = x_star;
      }

      if(y > tempb) {
        tempb = y;
        temp_best_l = *l;
      }
    }

    if (count == budget || problem->IOHprofiler_hit_optimal()) {
      break;
    }

    if(count >= lambda && bestl == temp_best_l) {
      bestl_count++;
    } else {
      bestl_count = 0;
    }
    bestl = temp_best_l;
    *sigma = sqrt(bestl * (1.0- (double)bestl/(double)n) * pow(0.98,bestl_count));
  }
};

/// \fn void fastGeneticAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger)
/// (1+lambda) fast GA, which mutation length is sampled from a power law distribution.
void fastGeneticAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
  int n = problem->IOHprofiler_get_number_of_variables(); /// < dimension
  int lambda = DEFAULT_LAMBDA;
  double mutation_rate = 1.0 / n;
  std::vector<int> parent;
  std::vector<int> x, x_star;
  double y;
  int count;
  int last_found_evaluation;
  
  std::vector<double> power_law_distribution;
  powerLawDistribution(power_law_distribution,n/2);

  budget = budget_scale * pow((double)(n),budget_power);

  /***
   * Recording Parameters Preparation
   ***/
  std::shared_ptr<double> l = std::make_shared<double>(1.0); /// < mutation length
  std::shared_ptr<double> update_evaluations = std::make_shared<double>(0.0); /// <  evaluations times since last best fitness has been found
  std::vector<std::shared_ptr<double> > parameters;
  parameters.push_back(l);
  parameters.push_back(update_evaluations);;
  std::vector<std::string> parameters_name;
  parameters_name.push_back("l");
  parameters_name.push_back("update_evaluations");
  logger->set_parameters(parameters,parameters_name);

  initialization(parent,n,2);
  double best_fitness = problem->evaluate(parent);
  logger->write_line(problem->loggerInfo());
  
  count = 1;
  last_found_evaluation = 1;
  while(count < budget && !problem->IOHprofiler_hit_optimal()) {
    x = parent;
    for (int i = 0; i < lambda; ++i) {
      *l = (double)sampleLRand(power_law_distribution);
      x_star = x;
      bitMutation(x_star,(int)(*l+0.5));

      y = problem->evaluate(x_star);
      count++;
      if (y > best_fitness) {
        *update_evaluations = (double)(count - last_found_evaluation);
        last_found_evaluation = count;
      }
      logger->write_line(problem->loggerInfo());

      if(y >= best_fitness) {
        best_fitness = y;
        parent = x_star;
      }
    }
  }
}

/// \fn void oneLambdaLambdaGeneticAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger)
/// (1+(lambda,lambda)) EA with self adaptive lambda, and its mutation rate and crossover rate is 
/// related with lambda.
void oneLambdaLambdaGeneticAlgorithm(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
  int n = problem->IOHprofiler_get_number_of_variables(); /// < dimension
  double lambda = DEFAULT_LAMBDA;
  double mutation_rate, crossover_rate;
  int l;

  budget = budget_scale * pow((double)(n),budget_power);
  double a = pow(1.5,0.25); /// < parameter for adjusting lambda
  double b = 2.0/3.0; /// < parameter for adjusting lambda

  
  std::vector<int> parent;
  std::vector<int> best;
  std::vector<int> mutation_offspring;
  initialization(parent,n,2);
  double best_fitness = problem->evaluate(parent);
  logger->write_line(problem->loggerInfo());
  best = parent;
  double best_value_mutation;

  int update_lambda_flag;
  std::vector<int> x, x_star;
  double y;

  int count = 1;
  
  std::vector<int> flip;
  while(count < budget && !problem->IOHprofiler_hit_optimal()) {

    mutation_rate = (double) lambda / (double) n;
    crossover_rate = 1.0 / (double) lambda;

    x = parent;
    
    update_lambda_flag = 0;
    best_value_mutation = -DBL_MAX;
    l = sampleLRand(mutation_rate,n);
    for (int i = 0; i < (int)(lambda+0.5); ++i) {
      
      x_star = x;
      sampleNFromM(flip,l,n);
      for (int j = 0; j != l; ++j) {
        x_star[flip[j]] = (x_star[flip[j]] + 1) % 2;
      } 

      y = problem->evaluate(x_star);
      logger->write_line(problem->loggerInfo());

      if(y > best_fitness) {
        update_lambda_flag = 1;
      }
      if(y > best_value_mutation) {
        best_value_mutation = y;
        mutation_offspring = x_star;
      }
      if(y >= best_fitness) {
        best_fitness = y;
        best = x_star;
      }

      count++;
      if (count == budget || problem->IOHprofiler_hit_optimal()) {
        break;
      }
    }
    if (count == budget || problem->IOHprofiler_hit_optimal()) {
      break;
    }

    for (int i = 0; i < (int)(lambda + 0.5); ++i) {
      if (crossover(x_star,parent,mutation_offspring,crossover_rate)) {
        if(compareVector(mutation_offspring,x_star)) {
          continue;
        } else {
          y = problem->evaluate(x_star);
          logger->write_line(problem->loggerInfo());
          count++;

          if(y > best_fitness) {
            update_lambda_flag = 1;  
          }

          if(y >= best_fitness) {
            best_fitness = y;
            best = x_star;
          }
        }

        if (count == budget || problem->IOHprofiler_hit_optimal()) {
          break;
        }
      }
    }

    if(update_lambda_flag == 1){
      lambda = (lambda * b) > 1 ? (lambda * b) : 1;
    }
    else{
      lambda =  (lambda * a) < n ? (lambda * a) : n;
    }

    parent = best;
  }
}

/// \fn void randomLocalSearch(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
/// (1+1) randomized local search
void randomLocalSearch(std::shared_ptr<IOHprofiler_problem<int> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
  int n = problem->IOHprofiler_get_number_of_variables();
  int lambda = 1;
  int temp_i;
  std::vector<int> parent;
  std::vector<int> x, x_star;
  double y;
  int count;
  int last_found_evaluation;

  budget = budget_scale * pow((double)(n),budget_power);

   /***
   * Recording Parameters Preparation
   ***/
  std::shared_ptr<double> l = std::make_shared<double>(1.0); /// < mutation length
  std::shared_ptr<double> update_evaluations = std::make_shared<double>(0.0); /// <  evaluations times since last best fitness has been found
  std::vector<std::shared_ptr<double> > parameters;
  parameters.push_back(l);
  parameters.push_back(update_evaluations);;
  std::vector<std::string> parameters_name;
  parameters_name.push_back("l");
  parameters_name.push_back("update_evaluations");
  logger->set_parameters(parameters,parameters_name);


  initialization(parent,n,2);
  double best_fitness = problem->evaluate(parent);
  logger->write_line(problem->loggerInfo());

  count = 1;
  last_found_evaluation = 1;
  while(count < budget && !problem->IOHprofiler_hit_optimal()) {
    x = parent;
    for (int i = 0; i != lambda; ++i) {
      
      //x_star = bitMutation(x,l);
      temp_i = random_generator.IOHprofiler_uniform_rand() * n;
      x_star = x;
      x_star[temp_i] = (x_star[temp_i] + 1)%2;
      y = problem->evaluate(x_star);
      count++;
      if (y > best_fitness) {
        *update_evaluations = (double)(count - last_found_evaluation);
        last_found_evaluation = count;
      }
      logger->write_line(problem->loggerInfo());
      if(y >= best_fitness) {
        best_fitness = y;
        parent = x_star;
      }
    }
  }
}


#endif