#include "randomsearch.h"

bool RandomSearch::Termination() {
  if (!this->problem_->IOHprofiler_hit_optimal() && this->evaluation_ < this->evluation_budget_ ) {
    return false;
  } else {
    return true;
  }
}

void RandomSearch::Preparation() {
  this->problem_->reset_problem();
  if (this->csv_logger_ != nullptr) {
    this->csv_logger_->track_problem(*this->problem_);
  }

  this->solution_ = vector <int> (this->get_dimension());
  this->evaluation_ = 0;
  this->generation_ = 0;
  this->best_individual_ = vector<int>(this->get_dimension());
  if (Opt == optimizationType::MAXIMIZATION) {
    this->best_found_fitness_ = numeric_limits<double>::lowest();
  } else {
    this->best_found_fitness_ = numeric_limits<double>::max();
  }
}

void RandomSearch::AssignProblem(shared_ptr<IOHprofiler_problem<int> > problem_ptr) {
  this->problem_ = problem_ptr;
}

void RandomSearch::AssignLogger(shared_ptr<IOHprofiler_csv_logger<int> > logger_ptr) {
  this->csv_logger_ = logger_ptr;
}

double RandomSearch::Evaluate(vector<int> &x) {
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

void RandomSearch::SetSeed(unsigned seed) {
  random_gen.seed(seed);
}

void RandomSearch::set_evaluation_budget(const size_t evaluation_budget) {
  this->evluation_budget_ = evaluation_budget;
}

void RandomSearch::set_independent_runs(const size_t independent_runs) {
  this->independent_runs_ = independent_runs;
}

void RandomSearch::set_solution(const vector< int > &solution) {
  this->solution_ = solution;
}

void RandomSearch::set_solution_fitness(const double solution_fitness) {
  this->solution_fitness_ = solution_fitness;
}

void RandomSearch::set_best_found_fitness(const double best_found_fitness) {
  this->best_found_fitness_ = best_found_fitness;
}

void RandomSearch::set_best_individual(const vector <int> best_individual) {
  this->best_individual_ = best_individual;
}

void RandomSearch::set_generation(const size_t generation) {
  this->generation_ = generation;
}

int RandomSearch::get_mu() const {
  return 1;
}

int RandomSearch::get_lambda() const {
  return 1;
}

int RandomSearch::get_dimension() const {
  return this->problem_->IOHprofiler_get_number_of_variables();
}

vector< int > RandomSearch::get_solution() {
  return this->solution_;
}

double RandomSearch::get_solution_fitness() {
  return this->solution_fitness_;
}

double RandomSearch::get_best_found_fitness() {
  return this->best_found_fitness_;
}

vector <int> RandomSearch::get_best_individual() {
  return this->best_individual_;
}

size_t RandomSearch::get_generation() {
  return this->generation_;
}

void RandomSearch::DoRandomSearch() {
  this->Preparation();
  
  while (!this->Termination()) {
    ++this->generation_;
    
    for (int i = 0; i != this->get_dimension(); ++i) {
      if (uniform_random() < 0.5) {
        this->solution_[i] = 1;
      } else {
        this->solution_[i] = 0;
      }
    } 

    this->solution_fitness_ = this->Evaluate(this->solution_);
  }
}

void RandomSearch::run(shared_ptr<IOHprofiler_suite<int> > suite) {
  while ((this->problem_ = suite->get_next_problem()) != nullptr) {
    size_t i = 0;
    while(++i <= this->independent_runs_) {
      size_t i = 0;
      this->DoRandomSearch();
    }
  }
}

void RandomSearch::run(string folder_path, string folder_name, shared_ptr<IOHprofiler_suite<int> > suite, int eval_budget, int independent_runs, unsigned rand_seed) {
  string algorithm_name = "random search";
  std::shared_ptr<IOHprofiler_csv_logger<int>> logger(new IOHprofiler_csv_logger<int>(folder_path,folder_name,algorithm_name,algorithm_name) );
  logger->activate_logger();
  this->AssignLogger(logger);
    

  this->set_evaluation_budget(eval_budget);
  this->set_independent_runs(independent_runs);
  this->SetSeed(rand_seed);
  
  this->run(suite);
  logger->clear_logger();
}
  
