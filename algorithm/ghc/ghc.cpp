#include "ghc.h"

#include "geneticAlgorithm.h"

bool GreedyHillClimber::Termination() {
  if (!this->problem_->IOHprofiler_hit_optimal() && this->evaluation_ < this->evluation_budget_) {
    return false;
  } else {
    return true;
  }
}

void GreedyHillClimber::Preparation() {
  this->problem_->reset_problem();
  if (this->csv_logger_ != nullptr) {
    this->csv_logger_->track_problem(*this->problem_);
  }

  this->parent_.clear();
  this->offspring_.clear();
  this->evaluation_ = 0;
  this->generation_ = 0;
  this->best_individual_ = vector<int>(this->problem_->IOHprofiler_get_number_of_variables());
  if (Opt == optimizationType::MAXIMIZATION) {
    this->best_found_fitness_ = numeric_limits<double>::lowest();
  } else {
    this->best_found_fitness_ = numeric_limits<double>::max();
  }
}

void GreedyHillClimber::AssignProblem(shared_ptr<IOHprofiler_problem<int> > problem_ptr) {
  this->problem_ = problem_ptr;
}

void GreedyHillClimber::AssignLogger(shared_ptr<IOHprofiler_csv_logger<int> > logger_ptr) {
  this->csv_logger_ = logger_ptr;
}

double GreedyHillClimber::Evaluate(vector<int> &x) {
  double result;
  result = this->problem_->evaluate(x);
  if (this->csv_logger_ != nullptr) {
    this->csv_logger_->do_log(this->problem_->loggerInfo()); /// TODO: we assume only using PBO suite now.
  }

  ++this->evaluation_;
  
  if (Opt == optimizationType::MAXIMIZATION) {
    if (result > this->best_found_fitness_) {
      this->best_found_fitness_ = result;
      this->best_individual_ = x;
    }
  } else {
    if (result < this->best_found_fitness_) {
      this->best_found_fitness_ = result;
      this->best_individual_ = x;
    }
  }
  return result;
}

void GreedyHillClimber::SetSeed(unsigned seed) {
  random_gen.seed(seed);
}

void GreedyHillClimber::set_evaluation_budget(const size_t evaluation_budget) {
  this->evluation_budget_ = evaluation_budget;
}


void GreedyHillClimber::set_independent_runs(const size_t independent_runs) {
  this->independent_runs_ = independent_runs;
}

void GreedyHillClimber::set_parent(const vector< int > &parent) {
  this->parent_ = parent;
}

void GreedyHillClimber::set_parent_fitness(const double parent_fitness) {
  this->parent_fitness_ = parent_fitness;
}

void GreedyHillClimber::set_offspring(const vector< int > &offspring) {
  this->offspring_ = offspring;
}

void GreedyHillClimber::set_offspring_fitness(const double offspring_fitness) {
  this->offspring_fitness_ = offspring_fitness;
}

void GreedyHillClimber::set_best_found_fitness(const double best_found_fitness) {
  this->best_found_fitness_ = best_found_fitness;
}

void GreedyHillClimber::set_best_individual(const vector <int> best_individual) {
  this->best_individual_ = best_individual;
}

void GreedyHillClimber::set_generation(const size_t generation) {
  this->generation_ = generation;
}

int GreedyHillClimber::get_mu() const {
  return 1;
}

int GreedyHillClimber::get_lambda() const {
  return 1;
}

int GreedyHillClimber::get_dimension() const {
  return this->problem_->IOHprofiler_get_number_of_variables();
}

vector< int > GreedyHillClimber::get_parent() {
  return this->parent_;
}

double GreedyHillClimber::get_parent_fitness() {
  return this->parent_fitness_;
}

vector< int >  GreedyHillClimber::get_offspring() {
  return this->offspring_;
}

double GreedyHillClimber::get_offspring_fitness() {
  return this->offspring_fitness_;
}

double GreedyHillClimber::get_best_found_fitness() {
  return this->best_found_fitness_;
}

vector <int> GreedyHillClimber::get_best_individual() {
  return this->best_individual_;
}

size_t GreedyHillClimber::get_generation() {
  return this->generation_;
}

void GreedyHillClimber::DoGreedyHillClimber() {
  this->Preparation();
  
  this->Initialization();
  size_t flip_index = 0;
  while (!this->Termination()) {
    ++this->generation_;
    
    this->offspring_.clear();
    this->offspring_ = this->parent_;
  
    // Filp the bit at flip_index.
    this->offspring_[flip_index] = (this->offspring_[flip_index] + 1) % 2;
    if(++flip_index >= this->get_dimension()) {
      flip_index = 0;
    }
    
    this->offspring_fitness_ = this->Evaluate(this->offspring_);

    if (this->offspring_fitness_ >= this->parent_fitness_) {
      this->parent_ = this->offspring_;
      this->parent_fitness_ = this->offspring_fitness_;
    }
  }
}

void GreedyHillClimber::Initialization() {
  int n = this->problem_->IOHprofiler_get_number_of_variables();
 
  this->parent_ = vector<int> (n,0);
  for (int i = 0; i != n; ++i) {
    if (uniform_random() < 0.5) {
      this->parent_[i] = 1;
    }
  }
  this->parent_fitness_ = this->Evaluate(this->parent_);
}

void GreedyHillClimber::run(shared_ptr<IOHprofiler_suite<int> > suite) {
  while ((this->problem_ = suite->get_next_problem()) != nullptr) {
    size_t i = 0;
    while(++i <= this->independent_runs_) {
      size_t i = 0;
      this->DoGreedyHillClimber();
    }
  }
}

void GreedyHillClimber::run(string folder_path, string folder_name, shared_ptr<IOHprofiler_suite<int> > suite, int eval_budget, int independent_runs, unsigned rand_seed) {
  string algorithm_name = "gHC";
  std::shared_ptr<IOHprofiler_csv_logger<int>> logger(new IOHprofiler_csv_logger<int>(folder_path,folder_name,algorithm_name,algorithm_name) );
  logger->activate_logger();
  this->AssignLogger(logger);
    

  this->set_evaluation_budget(eval_budget);
  this->set_independent_runs(independent_runs);
  this->SetSeed(rand_seed);
  
  this->run(suite);
  logger->clear_logger();
}
  
