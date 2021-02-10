/// \file geneticAlgorithm.h
/// \brief Header file for class GeneticAlgorithm.
///
/// A configurable genetic algorithm.
///
/// \author Furong Ye
/// \date 2020-07-02

#ifndef _GENETIC_ALGORITHM_H_
#define _GENETIC_ALGORITHM_H_

#include "common.h"
#include "crossover.h"
#include "mutation.h"
#include "selection.h"

#define DEFAULT_MU_ 1
#define DEFAULT_LAMBDA_ 1
#define DEFAULT_CROSSOVER_PROBABLITY_ 0
#define DEFAULT_CROSSOVER_MUTATION_RELATION_ 0

#define DEFAULT_EVALUATION_BUDGET_ 10000
#define DEFAULT_GENERATION_BUDGET_ 10000

#define DEFAULT_ECDF_TARGET_LENGTH 100
#define DEFAULT_ECDF_BUDGET_LENGTH 100

/// Definition of relation between mutation and crossover.
enum crossover_mutation_relation {
  IND = 1, /// < Always doing mutation, and doing crossover with probability p_c.
  OR = 0 /// < Doing crossover probability p_c, otherwise doing mutation.
};

class GeneticAlgorithm : public Crossover, public Mutation, public Selection {
public:
  GeneticAlgorithm() :
  mu_(DEFAULT_MU_),
  lambda_(DEFAULT_LAMBDA_),
  crossover_probability_(DEFAULT_CROSSOVER_PROBABLITY_),
  crossover_mutation_r_(DEFAULT_CROSSOVER_MUTATION_RELATION_),
  evaluation_(0),
  generation_(0),
  evluation_budget_(DEFAULT_EVALUATION_BUDGET_),
  generation_budget_(DEFAULT_GENERATION_BUDGET_),
  optimum_(numeric_limits<double>::max()) /// < TODO: now we assume doing maximization.
  {}
  
  ~GeneticAlgorithm() {}
  GeneticAlgorithm(const GeneticAlgorithm&) = delete;
  GeneticAlgorithm &operator = (const GeneticAlgorithm&) = delete;
  
  /// \fn virtual AdaptiveStrategy()
  /// \brief A virtual function for adaptive methods.
  ///
  /// If you are about to implement a GA with adaptive paramters, you should implement the adaptive method in this function. This function will be revoked at the end of each iteration/generation.
  virtual void AdaptiveStrategy() {};
  
  /// \fn virtual Terminate()
  /// \brief Terminate condition of the genetic algorithm.
  ///
  /// You can set up terminate condition in this function. By default, the GA terminates when the best fitness value is found or the budget is used out.
  virtual bool Termination();
  
  /// \fn DoGeneticAlgorithm()
  /// \brief Do function of the genetic algorithm.
  ///
  /// The order of processing functions is: Initialization() ->
  /// Loop{ Crossover() -> Mutation() -> Selection() -> ~AdaptiveStrategy() } -> ~Termination().
  /// The function is virtual, which allows users to implement their own algorithm.
  virtual void DoGeneticAlgorithm();
  
  void Initialization();
  
  void Preparation();
  
  void SelectTwoParents();
  
  double Evaluate(vector<int> &x);
  
  void SetSeed(unsigned seed);
  
  /// TODO
  // void Encode();
  
  /// TODO
  // void Decode();
  
  void AssignProblem(shared_ptr<IOHprofiler_problem<int> > problem_ptr);
  
  void AssignLogger(shared_ptr<IOHprofiler_csv_logger<int> > logger_ptr);
  
  void set_mu(const int mu);
  void set_lambda(const int lambda);
  void set_crossover_probability(const double crossover_probablity);
  void set_crossover_mutation_r(const int crossover_mutation_r);
  void set_crossover_mutation_r(string crossover_mutation_r);
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
  
  int get_mu() const;
  int get_lambda() const;
  int get_dimension() const;
  double get_crossover_probability() const;
  int get_crossover_mutation_r() const;
  vector< vector<int> > get_parents_population() const;
  vector<double> get_parents_fitness() const;
  vector< vector<int> >  get_offspring_population() const;
  vector<double> get_offspring_fitness() const;
  double get_best_found_fitness() const;
  vector <int> get_best_individual() const;
  size_t get_generation() const;
  size_t get_independent_runs() const;
  
  /// \fn SetAllParameters()
  /// \brief Set all parameters of the genetic algorithm.
  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
  void SetAllParameters(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para);
  
  /// \fn run()
  /// \brief Run genetic algorithm once with given parameters.
  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
  void run(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para);
  
  /// \fn run_N()
  /// \brief Run genetic algorithm `independent_runs_` times with given parameters.
  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
  void run_N(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para);
  
  /// \fn run()
  /// \brief Run genetic algorithm once with given parameters, on an IOHexperimenter suite.
  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
  /// \param suite a shared pointer of IOHprofiler_suite
  void run(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<IOHprofiler_suite<int> > suite);
  
  /// \fn run_N()
  /// \brief Run genetic algorithm `independent_runs_` times with given parameters, on an IOHexperimenter suite.
  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
  /// \param suite a shared pointer of IOHprofiler_suite
  void run_N(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<IOHprofiler_suite<int> > suite);
  
  /// \fn run_N()
  /// \brief Run genetic algorithm `independent_runs_` times  on an IOHexperimenter suite.
  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
  /// \param suite a shared pointer of IOHprofiler_suite
  void run_N(shared_ptr<IOHprofiler_suite<int> > suite);
  
//  /// \fn EstimateLinearECDF()
//  /// \brief Calculating the area under ECDF (partition budgets and targets with linear scale) of the genetic algorithm with given parameters, on an IOHexperimenter suite with given targets.
//  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
//  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
//  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
//  /// \param suite a shared pointer of IOHprofiler_suite
//  /// \param targets a vector of targets of problems of the suite.
//  double EstimateLinearECDF(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<IOHprofiler_suite<int> > suite, vector<double> targets);
//
//  /// \fn EstimateLogECDF()
//  /// \brief Calculating the area under ECDF (partition budgets and targets with log scale) of the genetic algorithm with given parameters, on an IOHexperimenter suite with given targets.
//  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
//  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
//  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
//  /// \param suite a shared pointer of IOHprofiler_suite
//  /// \param targets a vector of targets of problems of the suite.
//  double EstimateLogECDF(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<IOHprofiler_suite<int> > suite, vector<double> targets);
//
//  /// \fn EstimateLogERT()
//  /// \brief Calculating ERT values of the genetic algorithm with given parameters, on an IOHexperimenter suite with given targets and budgets.
//  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
//  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
//  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
//  /// \param suite a shared pointer of IOHprofiler_suite
//  /// \param targets a vector of targets of problems of the suite.
//  /// \param budgets a vector of budgets of problems.
//  vector<double> EstimateERT(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<IOHprofiler_suite<int> > suite, const vector<double> &targets, const vector<size_t> &budgets);
  
protected:
  int mu_; /// < parents population size
  int lambda_; /// < offspring population size
  double crossover_probability_; /// < probability to do crossover
  int crossover_mutation_r_; /// < a flag for correlation between crossover and mutation
  
  vector< vector<int> > parents_population_;
  vector< double > parents_fitness_;
  vector< vector<int> > offspring_population_;
  vector< double > offspring_fitness_;
  double best_found_fitness_;
  vector <int> best_individual_;
  
  size_t evaluation_; /// < evaluation times
  size_t generation_; /// < number of iterations/generations
  size_t evluation_budget_; /// < budget for evaluations
  size_t generation_budget_; /// < budget for generations
  
  size_t independent_runs_; /// < number of independent runs.
  
  double optimum_;
  
  vector<size_t> selected_parents_;
  
  /// TODO: we assume the type of problem are integer only now.
  shared_ptr< IOHprofiler_problem<int> > problem_;
  shared_ptr< IOHprofiler_csv_logger<int> > csv_logger_;
  //shared_ptr< IOHprofiler_ecdf_logger<int> > ecdf_logger_;
  
  /// For calculating ECDF
  size_t ecdf_sum_;
  double ecdf_ratio_;
  size_t ecdf_budget_width_;
  size_t ecdf_target_width_;
  
  
  /// For calculating ERT
  vector <double> hitting_time_;
  double hiting_target_;
  bool hitting_flag_;
  bool ERT_flag_;
};

#endif // _GENETIC_ALGORITHM_H_
