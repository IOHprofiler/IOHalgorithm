/// \file estimationOfDistribution.h
/// \brief Header file for class EstimationOfDistribution.
///
/// A class of estimation of distribution algorithm.
///
/// \author Furong Ye
/// \date 2020-11-25

#ifndef _ESTIMATION_OF_DISTRIBUTION_H
#define _ESTIMATION_OF_DISTRIBUTION_H

#include "common.h"

class EstimationOfDistribution {
public:
  EstimationOfDistribution(int mu = 25, int lambda = 50) :
  mu_(mu),
  lambda_(lambda) {}
  
  ~EstimationOfDistribution() {}
  EstimationOfDistribution(const EstimationOfDistribution&) = delete;
  EstimationOfDistribution &operator = (const EstimationOfDistribution&) = delete;
  
  /// \fn virtual Terminate()
  /// \brief Terminate condition of the algorithm.
  ///
  /// You can set up terminate condition in this function. By default, the algorithm terminates when the best fitness value is found or the budget is used out.
  virtual bool Termination();
  
  /// \fn DoEstimationOfDistribution()
  /// \brief Do function of the algorithm.
  ///
  /// The order of processing functions is: Initialization() ->
  /// Loop{SampleIndividuals() -> Selection()-> EstimationDistribution()} -> ~Termination().
  /// The function is virtual, which allows users to implement their own algorithm.
  virtual void DoEstimationOfDistribution();
  
  virtual void Initialization();
  
  void Preparation();
  
  double Evaluate(vector<int> &x);
  
  void EstimateVariablesDistribution();

  void SampleIndividual(const vector< vector <double> > &distributions, vector <int> &individual);

  void SetSeed(unsigned seed);
  
  void AssignProblem(shared_ptr<IOHprofiler_problem<int> > problem_ptr);
  
  void AssignLogger(shared_ptr<IOHprofiler_csv_logger<int> > logger_ptr);

  virtual void Selection();
  
  void set_mu(const int mu);
  void set_lambda(const int lambda);
  void set_evaluation_budget(const size_t evaluation_budget);
  void set_generation_budget(const size_t generation_budget);
  void set_independent_runs(const size_t independent_runs);
  void set_parents_population(const vector< vector<int> > parents_population);
  void set_parents_fitness(const vector<double> &parents_fitness);
  void set_parents_population(const vector<int> &offspring, const size_t index);
  void set_parents_fitness(const double fitness, const size_t index);
  void add_parents_population(const vector<int> & parent);
  void clear_parents_population();
  void add_parents_fitness(const double fitness);
  void clear_parents_fitness();
  void set_offspring_population(const vector< vector<int> > &offspring_population);
  void set_offspring_fitness(const vector<double> &offspring_fitness);
  void set_offspring_population(const vector<int> &offspring, const size_t index);
  void set_offspring_fitness(const double fitness, const size_t index);
  void add_offspring_population(const vector<int> &offspring);
  void clear_offspring_population();
  void add_offspring_fitness(const double fitness);
  void clear_offspring_fitness();
  void set_best_found_fitness(const double best_found_fitness);
  void set_best_individual(const vector <int> best_individual);
  void set_generation(const size_t generation);
  void update_generation();
  void set_variables_distribution(const vector < vector <double> > &variables_distribution);
  void set_variables_distribution(const vector <double> &variables_distribution, size_t index);
  vector < vector <double> > get_variables_distribution();
  vector <double> get_variables_distribution(size_t index);

  int get_mu() const;
  int get_lambda() const;
  int get_dimension() const;
  vector< vector<int> > get_parents_population();
  vector<double> get_parents_fitness();
  vector< vector<int> >  get_offspring_population();
  vector<double> get_offspring_fitness();
  double get_best_found_fitness();
  vector <int> get_best_individual();
  size_t get_generation();
  
  
  /// \fn run()
  /// \brief Run the algorithm once with given parameters.
  void run(string folder_path, string folder_name, shared_ptr<IOHprofiler_suite<int> > suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed);
  void run(shared_ptr<IOHprofiler_suite<int> > suite);

private:
  int mu_; /// < parents population size
  int lambda_; /// < offspring population size

  vector< vector<int> > parents_population_;
  vector< double > parents_fitness_;
  vector< vector<int> > offspring_population_;
  vector< double > offspring_fitness_;
  double best_found_fitness_;
  vector <int> best_individual_;

  vector < vector <double> > variables_distribution_;
  vector < vector <int> > variables_count_;

  size_t evaluation_; /// < evaluation times
  size_t generation_; /// < number of iterations/generations
  size_t evluation_budget_; /// < budget for evaluations
  size_t generation_budget_; /// < budget for generations
  
  size_t independent_runs_; /// < number of independent runs.

  /// TODO: we assume the type of problem are integer only now.
  shared_ptr< IOHprofiler_problem<int> > problem_;
  shared_ptr< IOHprofiler_csv_logger<int> > csv_logger_;
};

#endif // _ESTIMATION_OF_DISTRIBUTION_H
