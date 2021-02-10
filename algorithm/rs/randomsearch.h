/// \file randomsearch.h
/// \brief Header file for class ghc.
///
/// Random search generates candidates solution randomly.
/// 
///
/// \author Furong Ye
/// \date 2020-11-23

#ifndef _RANDOMSEARCH_H_
#define _RANDOMSEARCH_H_

#include "common.h"

class RandomSearch {
public:
  RandomSearch() {};

  ~RandomSearch() {}
  RandomSearch(const RandomSearch&) = delete;
  RandomSearch &operator = (const RandomSearch&) = delete;

  /// \fn virtual Terminate()
  /// \brief Terminate condition of the algorithm.
  ///
  /// You can set up terminate condition in this function. By default, the RandomSearch terminates when the best fitness value is found or the budget is used out.
  virtual bool Termination();
  
  /// \fn DoRandomSearch()
  /// \brief Do function of the random search algorithm.
  ///
  /// The order of processing functions is:
  /// Loop{ random_sampling } -> ~Termination().
  /// The function is virtual, which allows users to implement their own algorithm.
  virtual void DoRandomSearch();
    
  void Preparation();
  
  double Evaluate(vector<int> &x);

  void run(shared_ptr<IOHprofiler_suite<int> > suite);

  void run(string folder_path, string folder_name, shared_ptr<IOHprofiler_suite<int> > suite, int eval_budget, int independent_runs, unsigned rand_seed);
  
  void SetSeed(unsigned seed);
  
  void AssignProblem(shared_ptr<IOHprofiler_problem<int> > problem_ptr);
  
  void AssignLogger(shared_ptr<IOHprofiler_csv_logger<int> > logger_ptr);
  
  void set_evaluation_budget(const size_t evaluation_budget);
  void set_independent_runs(const size_t independent_runs);
  void set_solution(const vector< int > &solution);
  void set_solution_fitness(const double solution_fitness);
  void set_best_found_fitness(const double best_found_fitness);
  void set_best_individual(const vector <int> best_individual);
  void set_generation(const size_t generation);
  void update_generation();
  
  int get_mu() const;
  int get_lambda() const;
  int get_dimension() const;
  vector< int > get_solution();
  double get_solution_fitness();
  double get_best_found_fitness();
  vector< int > get_best_individual();
  size_t get_generation();
  
private:
 
  vector< int > solution_;
  double solution_fitness_;
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


#endif // _RANDOMSEARCH_H_