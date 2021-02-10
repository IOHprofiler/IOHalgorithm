#include "estimationOfDistribution.h"

bool EstimationOfDistribution::Termination() {
  if (!this->problem_->IOHprofiler_hit_optimal() && this->evaluation_ < this->evluation_budget_ && this->generation_ <= this->generation_budget_) {
    return false;
  } else {
    return true;
  }
}

void EstimationOfDistribution::run(string folder_path, string folder_name, shared_ptr<IOHprofiler_suite<int> > suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed) {
  string algorithm_name = "UMDA";
  std::shared_ptr<IOHprofiler_csv_logger<int>> logger(new IOHprofiler_csv_logger<int>(folder_path,folder_name,algorithm_name,algorithm_name) );
  logger->activate_logger();
  this->AssignLogger(logger);
  

  this->set_evaluation_budget(eval_budget);
  this->set_generation_budget(gene_budget);
  this->set_independent_runs(independent_runs);
  this->SetSeed(rand_seed);
  
  this->run(suite);
  logger->clear_logger();
}


void EstimationOfDistribution::run(shared_ptr<IOHprofiler_suite<int> > suite) {
  while ((this->problem_ = suite->get_next_problem()) != nullptr) {
    size_t i = 0;
    while(++i <= this->independent_runs_) {
      size_t i = 0;
      this->DoEstimationOfDistribution();
    }
  }
}

void EstimationOfDistribution::DoEstimationOfDistribution() {
  double rand;
  this->Preparation();
  this->Initialization();
  while (!this->Termination()) {
    ++this->generation_;
    if (this->offspring_population_.size() != this->lambda_) {
      this->offspring_population_.clear();
      this->offspring_fitness_.clear();
      this->offspring_population_ = vector < vector <int> > (this->lambda_, vector < int > (this->get_dimension(), 0));
      this->offspring_fitness_ = vector <double> (this->lambda_, 0.0);
    }
    for (size_t i = 0; i < this->lambda_; ++i) {
      this->SampleIndividual(this->variables_distribution_, this->offspring_population_[i]);
      this->offspring_fitness_[i] = this->Evaluate(this->offspring_population_[i]);
      if (this->Termination()) break;
    }
  
    if (this->Termination()) break;
    this->Selection();
    this->EstimateVariablesDistribution();
  }
}

void EstimationOfDistribution::Initialization() {
  vector< double > default_distribution (2,0.5);
  vector< int > default_count (2,0);
  for (size_t i = 0; i != this->get_dimension(); ++i) {
    this->variables_distribution_.push_back(default_distribution);
    this->variables_count_.push_back(default_count);
  }
}

void EstimationOfDistribution::Selection() {
  if (this->parents_population_.size() != this->mu_) {
    this->parents_population_.clear();
    this->parents_fitness_.clear();
    this->parents_population_ = vector < vector <int> > (this->mu_, vector< int > (this->get_dimension(), 0) );
    this->parents_fitness_ = vector <double > (this->mu_, 0);
  }

  vector< int > index(this->lambda_, 0);
  assert(this->lambda_ ==  this->offspring_fitness_.size());
  for (int i = 0; i != this->lambda_; ++i) {
    index[i] = i;
  }

  assert (this->mu_ <= this->lambda_);
  if (Opt == MAXIMIZATION) {
    partial_sort(index.begin(), index.begin() + this->mu_, index.end(),
    [&](const int a, const int b){return this->offspring_fitness_[a] >= this->offspring_fitness_[b];});
  } else {
    partial_sort(index.begin(), index.begin() + this->mu_, index.end(),
    [&](const int a, const int b){return this->offspring_fitness_[a] < this->offspring_fitness_[b];});
  }

  for (int i = 0; i != this->mu_; ++i) {
    this->parents_population_[i] = this->offspring_population_[index[i]];
    this->parents_fitness_[i] = this->offspring_fitness_[index[i]];
  }
}

void EstimationOfDistribution::Preparation() {
  this->problem_->reset_problem();
  if (this->csv_logger_ != nullptr) {
    this->csv_logger_->track_problem(*this->problem_);
  }

  this->parents_fitness_.clear();
  this->parents_population_.clear();
  this->offspring_fitness_.clear();
  this->offspring_population_.clear();
  this->evaluation_ = 0;
  this->generation_ = 0;
  this->best_individual_ = vector<int>(this->problem_->IOHprofiler_get_number_of_variables());
  if (Opt == optimizationType::MAXIMIZATION) {
    this->best_found_fitness_ = numeric_limits<double>::lowest();
  } else {
    this->best_found_fitness_ = numeric_limits<double>::max();
  }

  this->variables_distribution_.clear();
  this->variables_count_.clear();
}

void EstimationOfDistribution::AssignProblem(shared_ptr<IOHprofiler_problem<int> > problem_ptr) {
  this->problem_ = problem_ptr;
}

void EstimationOfDistribution::AssignLogger(shared_ptr<IOHprofiler_csv_logger<int> > logger_ptr) {
  this->csv_logger_ = logger_ptr;
}

double EstimationOfDistribution::Evaluate(vector<int> &x) {
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


void EstimationOfDistribution::EstimateVariablesDistribution() {
  for (size_t i = 0; i != this->get_dimension(); ++i) {
    for (size_t j = 0; j != this->variables_count_[i].size(); ++j) {
      this->variables_count_[i][j] = 0;
    }

    for (size_t j = 0; j != this->mu_; ++j) {
      this->variables_count_[i][this->parents_population_[j][i]]++;
    }

    for (size_t j = 0; j != this->variables_distribution_[i].size(); ++j) {
      this->variables_distribution_[i][j] = static_cast<double>(this->variables_count_[i][j]) / static_cast<double>(this->mu_);
      // Restriction for p
      this->variables_distribution_[i][j] = this->variables_distribution_[i][j] > (1.0 - 1.0 / static_cast<double>(this->get_dimension())) ?  
                                            (1.0 - 1.0 / static_cast<double>(this->get_dimension())) : this->variables_distribution_[i][j];
      this->variables_distribution_[i][j] = this->variables_distribution_[i][j] < 1.0 / static_cast<double>(this->get_dimension()) ? 
                                            (1.0 / static_cast<double>(this->get_dimension())) : this->variables_distribution_[i][j]; 
    }
  }
}

void EstimationOfDistribution::SampleIndividual(const vector< vector <double> > &distributions, vector <int> &individual) {
  if (individual.size() != this->get_dimension()) {
    individual.clear();
    individual = vector <int> (this->get_dimension(), 0);
  }

  size_t d = static_cast<size_t>( this->get_dimension() );
  assert(distributions.size() == d);
  int j = 0;
  double r,p;
  for (size_t i = 0; i != d; ++i) {
    r = uniform_random();
    j = 0;
    p = 0;
    while (j < distributions[i].size()) {
      if (r > p && r < p + distributions[i][j]) {
        break;
      }

      p += distributions[i][j];
      ++j;
    }
    individual[i] = j;
  }
}

void EstimationOfDistribution::SetSeed(unsigned seed) {
  random_gen.seed(seed);
}

void EstimationOfDistribution::set_mu(const int mu) {
  this->mu_ = mu;
}

void EstimationOfDistribution::set_lambda(const int lambda) {
  this->lambda_ = lambda;
}

void EstimationOfDistribution::set_evaluation_budget(const size_t evaluation_budget) {
  this->evluation_budget_ = evaluation_budget;
}

void EstimationOfDistribution::set_generation_budget(const size_t generation_budget) {
  this->generation_budget_ = generation_budget;
}

void EstimationOfDistribution::set_independent_runs(const size_t independent_runs) {
  this->independent_runs_ = independent_runs;
}

void EstimationOfDistribution::set_parents_population(const vector< vector<int> > parents_population) {
  this->parents_population_ = parents_population;
}

void EstimationOfDistribution::set_parents_fitness(const vector<double> &parents_fitness) {
  this->parents_fitness_ = parents_fitness;
}

void EstimationOfDistribution::set_parents_population(const vector<int> &parent, const size_t index) {
  assert(this->parents_population_.size() > index);
  this->parents_population_[index] = parent;
}

void EstimationOfDistribution::set_parents_fitness(const double fitness, const size_t index) {
  assert(this->parents_fitness_.size() > index);
  this->parents_fitness_[index] = fitness;
}

void EstimationOfDistribution::add_parents_population(const vector<int> & parent) {
  this->parents_population_.push_back(parent);
  
}
void EstimationOfDistribution::clear_parents_population() {
  this->parents_population_.clear();
}

void EstimationOfDistribution::add_parents_fitness(const double fitness) {
  this->parents_fitness_.push_back(fitness);
}

void EstimationOfDistribution::clear_parents_fitness() {
  this->parents_fitness_.clear();
}

void EstimationOfDistribution::set_offspring_population(const vector< vector<int> > &offspring_population) {
  this->offspring_population_ = offspring_population;
}

void EstimationOfDistribution::set_offspring_fitness(const vector<double> &offspring_fitness) {
  this->offspring_fitness_ = offspring_fitness;
}

void EstimationOfDistribution::set_offspring_population(const vector<int> &offspring, const size_t index) {
  assert(this->offspring_population_.size() > index);
  this->offspring_population_[index] = offspring;
}

void EstimationOfDistribution::set_offspring_fitness(const double fitness, const size_t index) {
  assert(this->offspring_fitness_.size() > index);
  this->offspring_fitness_[index] = fitness;
}

void EstimationOfDistribution::set_variables_distribution(const vector < vector <double> > &variables_distribution) {
  this->variables_distribution_ = variables_distribution;
}

void EstimationOfDistribution::set_variables_distribution(const vector <double> &variables_distribution, size_t index) {
  this->variables_distribution_[index] = variables_distribution;
}

void EstimationOfDistribution::add_offspring_population(const vector<int> & parent) {
  this->offspring_population_.push_back(parent);
  
}
void EstimationOfDistribution::clear_offspring_population() {
  this->offspring_population_.clear();
}

void EstimationOfDistribution::add_offspring_fitness(const double fitness) {
  this->offspring_fitness_.push_back(fitness);
}

void EstimationOfDistribution::clear_offspring_fitness() {
  this->offspring_fitness_.clear();
}

void EstimationOfDistribution::set_best_found_fitness(const double best_found_fitness) {
  this->best_found_fitness_ = best_found_fitness;
}

void EstimationOfDistribution::set_best_individual(const vector <int> best_individual) {
  this->best_individual_ = best_individual;
}

void EstimationOfDistribution::set_generation(const size_t generation) {
  this->generation_ = generation;
}

void EstimationOfDistribution::update_generation() {
  this->generation_++;
}

int EstimationOfDistribution::get_mu() const {
  return this->mu_;
}

int EstimationOfDistribution::get_lambda() const {
  return this->lambda_;
}

int EstimationOfDistribution::get_dimension() const {
  return this->problem_->IOHprofiler_get_number_of_variables();
}

vector< vector<int> > EstimationOfDistribution::get_parents_population() {
  return this->parents_population_;
}

vector<double> EstimationOfDistribution::get_parents_fitness() {
  return this->parents_fitness_;
}

vector< vector<int> >  EstimationOfDistribution::get_offspring_population() {
  return this->offspring_population_;
}

vector<double> EstimationOfDistribution::get_offspring_fitness() {
  return this->offspring_fitness_;
}

double EstimationOfDistribution::get_best_found_fitness() {
  return this->best_found_fitness_;
}

vector <int> EstimationOfDistribution::get_best_individual() {
  return this->best_individual_;
}

size_t EstimationOfDistribution::get_generation() {
  return this->generation_;
}

vector < vector <double> > EstimationOfDistribution::get_variables_distribution() {
  return this->variables_distribution_;
}

vector <double> EstimationOfDistribution::get_variables_distribution(size_t index) {
  return this->variables_distribution_[index];
}