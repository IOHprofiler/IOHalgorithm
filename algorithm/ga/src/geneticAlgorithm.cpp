#include "geneticAlgorithm.h"

bool GeneticAlgorithm::Termination() {
  if (!this->problem_->IOHprofiler_hit_optimal() && this->evaluation_ < this->evluation_budget_ && this->generation_ <= this->generation_budget_) {
    return false;
  } else {
    return true;
  }
}

void GeneticAlgorithm::SetAllParameters(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para) {
  assert(integer_para.size() == 4);
  this->set_mu(integer_para[0]);
  this->set_lambda(integer_para[1]);
  this->set_l(integer_para[2]);
  this->set_tournament_k(integer_para[3]);
  
  assert(continuous_para.size() == 6);
  this->set_crossover_probability(continuous_para[0]);
  this->set_p_u(continuous_para[1]);
  this->set_mutation_rate(continuous_para[2]);
  this->set_r_n(continuous_para[3]);
  this->set_sigma_n(continuous_para[4]);
  this->set_beta_f(continuous_para[5]);
  
  assert(category_para.size() == 4);
  this->set_crossover_mutation_r(category_para[0]);
  this->set_crossover_operator(category_para[1]);
  this->set_mutation_operator(category_para[2]);
  this->set_selection_operator(category_para[3]);
}

void GeneticAlgorithm::run(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para) {
  this->SetAllParameters(integer_para, continuous_para, category_para);
  this->DoGeneticAlgorithm();
}

void GeneticAlgorithm::run(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para,  shared_ptr<IOHprofiler_suite<int> > suite) {
  this->SetAllParameters(integer_para, continuous_para, category_para);
  
  while ((this->problem_ = suite->get_next_problem()) != nullptr) {
    this->DoGeneticAlgorithm();
  }
}

void GeneticAlgorithm::run_N(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para) {
  this->SetAllParameters(integer_para, continuous_para, category_para);
  
  size_t r = 0;
  this->ecdf_sum_ = 0;
  while (r < this->independent_runs_) {
    this->DoGeneticAlgorithm();
    r++;
  }
  
  this->ecdf_sum_ = this->ecdf_sum_ / this->independent_runs_;
  this->ecdf_ratio_ = (double)this->ecdf_sum_ / (double)this->ecdf_budget_width_ / (double)this->ecdf_target_width_ / this->independent_runs_;
}

void GeneticAlgorithm::run_N(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<IOHprofiler_suite<int> > suite) {
  this->SetAllParameters(integer_para, continuous_para, category_para);
  while ((this->problem_ = suite->get_next_problem()) != nullptr) {
    this->run_N(integer_para,continuous_para,category_para);
  }
}

void GeneticAlgorithm::run_N(shared_ptr<IOHprofiler_suite<int> > suite) {
  while ((this->problem_ = suite->get_next_problem()) != nullptr) {
    size_t r = 0;
    this->ecdf_sum_ = 0;
    while (r < this->independent_runs_) {
      this->DoGeneticAlgorithm();
      r++;
    }
    
    this->ecdf_sum_ = this->ecdf_sum_ / this->independent_runs_;
    this->ecdf_ratio_ = (double)this->ecdf_sum_ / (double)this->ecdf_budget_width_ / (double)this->ecdf_target_width_ / this->independent_runs_;
  }
}

void GeneticAlgorithm::DoGeneticAlgorithm() {
  double rand;
  this->Preparation();
  
  this->Initialization();
  while (!this->Termination()) {
    ++this->generation_;
    
    this->offspring_population_.clear();
    this->offspring_fitness_.clear();
    for (size_t i = 0; i < this->lambda_; ++i) {
      this->SelectTwoParents();
      this->offspring_population_.push_back(this->parents_population_[this->selected_parents_[0]]);
      
      rand = uniform_random();
      this->c_flipped_index.clear();
      this->m_flipped_index.clear();
      if (rand < this->crossover_probability_) {
        this->DoCrossover(this->offspring_population_[i], this->parents_population_[this->selected_parents_[0]], this->parents_population_[this->selected_parents_[1]]);
      }
      
      if (this->crossover_mutation_r_) {
        this->DoMutation(this->offspring_population_[i]);
      } else if (rand >= this->crossover_probability_) {
        this->DoMutation(this->offspring_population_[i]);
      }
      
      if (this->c_flipped_index == this->m_flipped_index) { /// If the flipping indexes of crossover and mutation are identical, the individual remains the same.
        this->offspring_fitness_.push_back(this->parents_fitness_[this->selected_parents_[0]]);
      } else if (rand < this->crossover_probability_ && this->offspring_population_[i] == this->parents_population_[this->selected_parents_[1]]) { /// If the offspring is identical with the second parent.
        /// TODO: Do something to save time for this comparison.
        this->offspring_fitness_.push_back(this->parents_fitness_[this->selected_parents_[1]]);
      } else { /// otherwise evaluate.
        this->offspring_fitness_.push_back(this->Evaluate(this->offspring_population_[i]));
      }
      
      if (this->Termination()) break;
    }
    
    if (this->Termination()) break;
    
    this->DoSelection(this->parents_population_, this->parents_fitness_, this->offspring_population_, this->offspring_fitness_);
    this->AdaptiveStrategy();
  }
}

void GeneticAlgorithm::Initialization() {
  int n = this->problem_->IOHprofiler_get_number_of_variables();
  for (int i = 0; i != this->mu_; ++i) {
    vector<int> tmp(n,0);
    for (int j = 0; j != n; ++j) {
      if (uniform_random() < 0.5) {
        tmp[j] = 1;
      }
    }
    this->parents_population_.push_back(tmp);
    this->parents_fitness_.push_back(this->Evaluate(tmp));
  }
}

void GeneticAlgorithm::Preparation() {
  this->problem_->reset_problem();
  if (this->csv_logger_ != nullptr) {
    this->csv_logger_->track_problem(*this->problem_);
  }
//  if (this->ecdf_logger_ != nullptr) {
//    this->ecdf_logger_->track_problem(*this->problem_);
//  }
  this->selected_parents_ = vector<size_t>(2);
  this->optimum_ = this->problem_->IOHprofiler_get_optimal()[0];
  this->PowerLawDistribution(this->problem_->IOHprofiler_get_number_of_variables());
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
  this->hitting_flag_ = false;
}

void GeneticAlgorithm::SelectTwoParents() {
  this->selected_parents_[0] = static_cast<size_t>( floor(uniform_random() * this->mu_) );
  if (this->mu_ >= 2) {
    do {
      this->selected_parents_[1] = static_cast<size_t>( floor(uniform_random() * this->mu_) );
    } while(this->selected_parents_[0] == this->selected_parents_[1]);
  }
}

void GeneticAlgorithm::AssignProblem(shared_ptr<IOHprofiler_problem<int> > problem_ptr) {
  this->problem_ = problem_ptr;
}

void GeneticAlgorithm::AssignLogger(shared_ptr<IOHprofiler_csv_logger<int> > logger_ptr) {
  this->csv_logger_ = logger_ptr;
}

double GeneticAlgorithm::Evaluate(vector<int> &x) {
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
  
  // This is only used to calculate ERT
  if (this->ERT_flag_) {
    if (Opt == optimizationType::MAXIMIZATION) {
      if (!this->hitting_flag_ && this->best_found_fitness_ >= this->hiting_target_) {
        this->hitting_time_.push_back(this->evaluation_);
        this->hitting_flag_ = true;
      }
    } else {
      if (!this->hitting_flag_ && this->best_found_fitness_ <= this->hiting_target_) {
        this->hitting_time_.push_back(this->evaluation_);
        this->hitting_flag_ = true;
      }
    }
  }
  return result;
}

void GeneticAlgorithm::SetSeed(unsigned seed) {
  random_gen.seed(seed);
}

void GeneticAlgorithm::set_mu(const int mu) {
  this->mu_ = mu;
}

void GeneticAlgorithm::set_lambda(const int lambda) {
  this->lambda_ = lambda;
}

void GeneticAlgorithm::set_crossover_probability(const double crossover_probability) {
  this->crossover_probability_ = crossover_probability;
}

void GeneticAlgorithm::set_crossover_mutation_r(const int crossover_mutation_r) {
  this->crossover_mutation_r_ = crossover_mutation_r;
}

void GeneticAlgorithm::set_crossover_mutation_r(string crossover_mutation_r) {
  transform(crossover_mutation_r.begin(), crossover_mutation_r.end(), crossover_mutation_r.begin(), ::toupper);
  if (crossover_mutation_r == "IND") {
    this->set_crossover_mutation_r(IND);
  } else if (crossover_mutation_r == "OR") {
    this->set_crossover_mutation_r(OR);
  } else {
    cerr << "invalid value for set_crossover_mutation_r";
    assert(false);
  }
}

void GeneticAlgorithm::set_evaluation_budget(const size_t evaluation_budget) {
  this->evluation_budget_ = evaluation_budget;
}

void GeneticAlgorithm::set_generation_budget(const size_t generation_budget) {
  this->generation_budget_ = generation_budget;
}

void GeneticAlgorithm::set_independent_runs(const size_t independent_runs) {
  this->independent_runs_ = independent_runs;
}

void GeneticAlgorithm::set_parents_population(const vector< vector<int> > parents_population) {
  this->parents_population_ = parents_population;
}

void GeneticAlgorithm::set_parents_fitness(const vector<double> &parents_fitness) {
  this->parents_fitness_ = parents_fitness;
}

void GeneticAlgorithm::set_parents_population(const vector<int> &parent, const size_t index) {
  assert(this->parents_population_.size() > index);
  this->parents_population_[index] = parent;
}

void GeneticAlgorithm::set_parents_fitness(const double fitness, const size_t index) {
  assert(this->parents_fitness_.size() > index);
  this->parents_fitness_[index] = fitness;
}

void GeneticAlgorithm::add_parents_population(const vector<int> & parent) {
  this->parents_population_.push_back(parent);
  
}
void GeneticAlgorithm::clear_parents_population() {
  this->parents_population_.clear();
}

void GeneticAlgorithm::add_parents_fitness(const double fitness) {
  this->parents_fitness_.push_back(fitness);
}

void GeneticAlgorithm::clear_parents_fitness() {
  this->parents_fitness_.clear();
}

void GeneticAlgorithm::set_offspring_population(const vector< vector<int> > &offspring_population) {
  this->offspring_population_ = offspring_population;
}

void GeneticAlgorithm::set_offspring_fitness(const vector<double> &offspring_fitness) {
  this->offspring_fitness_ = offspring_fitness;
}

void GeneticAlgorithm::set_offspring_population(const vector<int> &offspring, const size_t index) {
  assert(this->offspring_population_.size() > index);
  this->offspring_population_[index] = offspring;
}

void GeneticAlgorithm::set_offspring_fitness(const double fitness, const size_t index) {
  assert(this->offspring_fitness_.size() > index);
  this->offspring_fitness_[index] = fitness;
}

void GeneticAlgorithm::add_offspring_population(const vector<int> & parent) {
  this->offspring_population_.push_back(parent);
  
}
void GeneticAlgorithm::clear_offspring_population() {
  this->offspring_population_.clear();
}

void GeneticAlgorithm::add_offspring_fitness(const double fitness) {
  this->offspring_fitness_.push_back(fitness);
}

void GeneticAlgorithm::clear_offspring_fitness() {
  this->offspring_fitness_.clear();
}

void GeneticAlgorithm::set_best_found_fitness(const double best_found_fitness) {
  this->best_found_fitness_ = best_found_fitness;
}

void GeneticAlgorithm::set_best_individual(const vector <int> best_individual) {
  this->best_individual_ = best_individual;
}

void GeneticAlgorithm::set_generation(const size_t generation) {
  this->generation_ = generation;
}

void GeneticAlgorithm::update_generation() {
  this->generation_++;
}

int GeneticAlgorithm::get_mu() const {
  return this->mu_;
}

int GeneticAlgorithm::get_lambda() const {
  return this->lambda_;
}

int GeneticAlgorithm::get_dimension() const {
  return this->problem_->IOHprofiler_get_number_of_variables();
}

double GeneticAlgorithm::get_crossover_probability() const {
  return this->crossover_probability_;
}

int GeneticAlgorithm::get_crossover_mutation_r() const {
  return this->crossover_mutation_r_;
}

vector< vector<int> > GeneticAlgorithm::get_parents_population() const{
  return this->parents_population_;
}

vector<double> GeneticAlgorithm::get_parents_fitness() const{
  return this->parents_fitness_;
}

vector< vector<int> >  GeneticAlgorithm::get_offspring_population() const{
  return this->offspring_population_;
}

vector<double> GeneticAlgorithm::get_offspring_fitness() const{
  return this->offspring_fitness_;
}

double GeneticAlgorithm::get_best_found_fitness() const{
  return this->best_found_fitness_;
}

vector <int> GeneticAlgorithm::get_best_individual() const{
  return this->best_individual_;
}

size_t GeneticAlgorithm::get_generation() const{
  return this->generation_;
}

size_t GeneticAlgorithm::get_independent_runs() const {
  return this->independent_runs_;
}
