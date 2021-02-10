/// \file ghc.h
/// \brief Header file for class ghc.
///
/// (1+1) greedy hill climber algorithm goes throught the string from left to right,
/// flipping exactly one bit per iteration, and accepting the offspring if it
/// is as least as good as its parent.
/// 
///
/// \author Furong Ye
/// \date 2020-11-23

#ifndef _GHC_H_
#define _GHC_H_

#include "common.h"

class GreedyHillClimber {
public:
  GreedyHillClimber() {};

  ~GreedyHillClimber() {}
  GreedyHillClimber(const GreedyHillClimber&) = delete;
  GreedyHillClimber &operator = (const GreedyHillClimber&) = delete;

  /// \fn virtual Terminate()
  /// \brief Terminate condition of the genetic algorithm.
  ///
  /// You can set up terminate condition in this function. By default, the GreedyHillClimber terminates when the best fitness value is found or the budget is used out.
  virtual bool Termination();
  
  /// \fn DoGreedyHillClimber()
  /// \brief Do function of the greedy hill climber algorithm.
  ///
  /// The order of processing functions is: Initialization() ->
  /// Loop{ Mutation() -> Selection()} -> ~Termination().
  /// The function is virtual, which allows users to implement their own algorithm.
  virtual void DoGreedyHillClimber();
  
  void Initialization();
  
  void Preparation();
  
  double Evaluate(vector<int> &x);

  void run(shared_ptr<IOHprofiler_suite<int> > suite);

  void run(string folder_path, string folder_name, shared_ptr<IOHprofiler_suite<int> > suite, int eval_budget, int independent_runs, unsigned rand_seed);
  
  void SetSeed(unsigned seed);
  
  void AssignProblem(shared_ptr<IOHprofiler_problem<int> > problem_ptr);
  
  void AssignLogger(shared_ptr<IOHprofiler_csv_logger<int> > logger_ptr);
  
  void set_evaluation_budget(const size_t evaluation_budget);
  void set_independent_runs(const size_t independent_runs);
  void set_parent(const vector< int > &parent);
  void set_parent_fitness(const double parents_fitness);
  void set_offspring(const vector< int > &offspring);
  void set_offspring_fitness(const double offspring_fitness);
  void set_best_found_fitness(const double best_found_fitness);
  void set_best_individual(const vector <int> best_individual);
  void set_generation(const size_t generation);
  void update_generation();
  
  int get_mu() const;
  int get_lambda() const;
  int get_dimension() const;
  vector< int > get_parent();
  double get_parent_fitness();
  vector< int >  get_offspring();
  double get_offspring_fitness();
  double get_best_found_fitness();
  vector< int > get_best_individual();
  size_t get_generation();
  
private:
 
  vector< int > parent_;
  double parent_fitness_;
  vector< int > offspring_;
  double offspring_fitness_;
  double best_found_fitness_;
  vector <int> best_individual_;
  
  size_t evaluation_; /// < evaluation times
  size_t generation_; /// < number of iterations/generations
  size_t evluation_budget_; /// < budget for evaluations
  
  size_t independent_runs_; /// < number of independent runs.
  
  /// TODO: we assume the type of problem are integer only now.
  shared_ptr< IOHprofiler_problem<int> > problem_;
  shared_ptr< IOHprofiler_csv_logger<int> > csv_logger_;
};


#endif // _GHC_H_