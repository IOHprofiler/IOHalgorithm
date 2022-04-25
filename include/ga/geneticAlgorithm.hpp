#pragma once

/// \file geneticAlgorithm.h
/// \brief Header file for class GeneticAlgorithm.
///
/// A configurable genetic algorithm.
///
/// \author Furong Ye
/// \date 2020-07-02

#include "common.hpp"
#include "crossover.hpp"
#include "mutation.hpp"
#include "selection.hpp"

#define DEFAULT_MU_ 1
#define DEFAULT_LAMBDA_ 1
#define DEFAULT_CROSSOVER_PROBABLITY_ 0
#define DEFAULT_CROSSOVER_MUTATION_RELATION_ 0

#define DEFAULT_EVALUATION_BUDGET_ 10000
#define DEFAULT_GENERATION_BUDGET_ 10000

#define DEFAULT_ECDF_TARGET_LENGTH 100
#define DEFAULT_ECDF_BUDGET_LENGTH 100

namespace modularGA {


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
      GeneticAlgorithm& operator = (const GeneticAlgorithm&) = delete;

      /// \fn virtual AdaptiveStrategy()
      /// \brief A virtual function for adaptive methods.
      ///
      /// If you are about to implement a GA with adaptive paramters, you should implement the adaptive method in this function. This function will be revoked at the end of each iteration/generation.
      virtual void AdaptiveStrategy() {};

      /// \fn virtual Terminate()
      /// \brief Terminate condition of the genetic algorithm.
      ///
      /// You can set up terminate condition in this function. By default, the GA terminates when the best fitness value is found or the budget is used out.
      bool Termination() {
        if (!this->problem_->state().optimum_found && this->problem_->state().evaluations < this->evluation_budget_ && this->generation_ <= this->generation_budget_) {
          return false;
        }
        else {
          return true;
        }
      }

      /// \fn SetAllParameters()
            /// \brief Set all parameters of the genetic algorithm.
            /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
            /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
            /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
      void SetAllParameters(const vector<int>& integer_para, const vector<double>& continuous_para, const vector<string>& category_para) {
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

      /// \fn run()
            /// \brief Run genetic algorithm once with given parameters.
            /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
            /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
            /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}

      virtual void run(const vector<int>& integer_para, const vector<double>& continuous_para, const vector<string>& category_para) {
        this->SetAllParameters(integer_para, continuous_para, category_para);
        this->DoGeneticAlgorithm();
      }

      /// \fn run()
            /// \brief Run genetic algorithm once with given parameters, on an IOHexperimenter suite.
            /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
            /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
            /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
            /// \param suite a shared pointer of IOHprofiler_suite

      virtual void run(const vector<int>& integer_para, const vector<double>& continuous_para, const vector<string>& category_para, shared_ptr<ioh::suite::Suite<ioh::problem::Integer> > suite) {
        this->SetAllParameters(integer_para, continuous_para, category_para);

        for (const auto &p : *suite) {
          this->problem_ = p;
          this->DoGeneticAlgorithm();
        }
      }


      /// \fn run_N()
            /// \brief Run genetic algorithm `independent_runs_` times with given parameters.
            /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
            /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
            /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}

      void run_N(const vector<int>& integer_para, const vector<double>& continuous_para, const vector<string>& category_para) {
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

      /// \fn run_N()
            /// \brief Run genetic algorithm `independent_runs_` times with given parameters, on an IOHexperimenter suite.
            /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
            /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
            /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
            /// \param suite a shared pointer of IOHprofiler_suite
      virtual void run_N(const vector<int>& integer_para, const vector<double>& continuous_para, const vector<string>& category_para, shared_ptr<ioh::suite::Suite<ioh::problem::Integer> > suite) {
        this->SetAllParameters(integer_para, continuous_para, category_para);
        for (const auto &p : *suite) {
          this->problem_ = p;
          this->run_N(integer_para, continuous_para, category_para);
        }
      }

      /// \fn run_N()
            /// \brief Run genetic algorithm `independent_runs_` times  on an IOHexperimenter suite.
            /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
            /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
            /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
            /// \param suite a shared pointer of IOHprofiler_suite

      virtual void run_N(shared_ptr<ioh::suite::Suite<ioh::problem::Integer> > suite) {
         for (const auto &p : *suite) {
          this->problem_ = p;
          size_t r = 0;
          this->ecdf_sum_ = 0;
          while (r < this->independent_runs_) {
            this->DoGeneticAlgorithm();
            r++;
          }
          this->problem_->reset();
          this->ecdf_sum_ = this->ecdf_sum_ / this->independent_runs_;
          this->ecdf_ratio_ = (double)this->ecdf_sum_ / (double)this->ecdf_budget_width_ / (double)this->ecdf_target_width_ / this->independent_runs_;
        }
      }

      virtual void DoGeneticAlgorithm() {
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
            }
            else if (rand >= this->crossover_probability_) {
              this->DoMutation(this->offspring_population_[i]);
            }

            if (this->c_flipped_index == this->m_flipped_index) { /// If the flipping indexes of crossover and mutation are identical, the individual remains the same.
              this->offspring_fitness_.push_back(this->parents_fitness_[this->selected_parents_[0]]);
            }
            else if (rand < this->crossover_probability_ && this->offspring_population_[i] == this->parents_population_[this->selected_parents_[1]]) { /// If the offspring is identical with the second parent.
           /// TODO: Do something to save time for this comparison.
              this->offspring_fitness_.push_back(this->parents_fitness_[this->selected_parents_[1]]);
            }
            else { /// otherwise evaluate.
              this->offspring_fitness_.push_back(this->Evaluate(this->offspring_population_[i]));
            }

            if (this->Termination()) break;
          }

          if (this->Termination()) break;

          this->DoSelection(this->parents_population_, this->parents_fitness_, this->offspring_population_, this->offspring_fitness_);
          this->AdaptiveStrategy();
        }
      }

      virtual void Initialization() {
        int n = this->problem_->meta_data().n_variables;
        for (int i = 0; i != this->mu_; ++i) {
          vector<int> tmp(n, 0);
          for (int j = 0; j != n; ++j) {
            if (uniform_random() < 0.5) {
              tmp[j] = 1;
            }
          }
          this->parents_population_.push_back(tmp);
          this->parents_fitness_.push_back(this->Evaluate(tmp));
        }
      }

      void Preparation() {
        this->problem_->reset();
        if (this->csv_logger_ != nullptr) {
          this->problem_->attach_logger(*this->csv_logger_);
        }
        //  if (this->ecdf_logger_ != nullptr) {
        //    this->ecdf_logger_->track_problem(*this->problem_);
        //  }
        this->selected_parents_ = vector<size_t>(2);
        this->optimum_ = this->problem_->objective().y;
        this->PowerLawDistribution(this->problem_->meta_data().n_variables);
        this->parents_fitness_.clear();
        this->parents_population_.clear();
        this->offspring_fitness_.clear();
        this->offspring_population_.clear();
        this->evaluation_ = 0;
        this->generation_ = 0;
        this->best_individual_ = vector<int>(this->problem_->meta_data().n_variables);
        if (Opt == optimizationType::MAXIMIZATION) {
          this->best_found_fitness_ = numeric_limits<double>::lowest();
        }
        else {
          this->best_found_fitness_ = numeric_limits<double>::max();
        }
        this->hitting_flag_ = false;
      }

      void SelectTwoParents() {
        this->selected_parents_[0] = static_cast<size_t>(floor(uniform_random() * this->mu_));
        if (this->mu_ >= 2) {
          do {
            this->selected_parents_[1] = static_cast<size_t>(floor(uniform_random() * this->mu_));
          } while (this->selected_parents_[0] == this->selected_parents_[1]);
        }
      }

      void AssignProblem(shared_ptr<ioh::problem::Integer> problem_ptr) {
        this->problem_ = problem_ptr;
      }

      void AssignLogger(shared_ptr<ioh::logger::Analyzer > logger_ptr) {
        this->csv_logger_ = logger_ptr;
      }

      double Evaluate(vector<int>& x) {
        double result;
        result = (*this->problem_)(x);
        // if (this->csv_logger_ != nullptr) {
        //   this->csv_logger_->do_log(this->problem_->loggerInfo()); /// TODO: we assume only using PBO suite now.
        // }
        ++this->evaluation_;

        if (Opt == optimizationType::MAXIMIZATION) {
          if (result > this->best_found_fitness_) {
            this->best_found_fitness_ = result;
            this->best_individual_ = x;
          }
        }
        else {
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
          }
          else {
            if (!this->hitting_flag_ && this->best_found_fitness_ <= this->hiting_target_) {
              this->hitting_time_.push_back(this->evaluation_);
              this->hitting_flag_ = true;
            }
          }
        }
        return result;
      }

      void SetSeed(unsigned seed) {
        random_gen.seed(seed);
      }

      void set_mu(const int mu) {
        this->mu_ = mu;
      }

      void set_lambda(const int lambda) {
        this->lambda_ = lambda;
      }

      void set_crossover_probability(const double crossover_probability) {
        this->crossover_probability_ = crossover_probability;
      }

      void set_crossover_mutation_r(const int crossover_mutation_r) {
        this->crossover_mutation_r_ = crossover_mutation_r;
      }

      void set_crossover_mutation_r(string crossover_mutation_r) {
        transform(crossover_mutation_r.begin(), crossover_mutation_r.end(), crossover_mutation_r.begin(), ::toupper);
        if (crossover_mutation_r == "IND") {
          this->set_crossover_mutation_r(IND);
        }
        else if (crossover_mutation_r == "OR") {
          this->set_crossover_mutation_r(OR);
        }
        else {
          cerr << "invalid value for set_crossover_mutation_r";
          assert(false);
        }
      }

      void set_evaluation_budget(const size_t evaluation_budget) {
        this->evluation_budget_ = evaluation_budget;
      }

      void set_generation_budget(const size_t generation_budget) {
        this->generation_budget_ = generation_budget;
      }

      void set_independent_runs(const size_t independent_runs) {
        this->independent_runs_ = independent_runs;
      }

      void set_parents_population(const vector< vector<int> > parents_population) {
        this->parents_population_ = parents_population;
      }

      void set_parents_fitness(const vector<double>& parents_fitness) {
        this->parents_fitness_ = parents_fitness;
      }

      void set_parents_population(const vector<int>& parent, const size_t index) {
        assert(this->parents_population_.size() > index);
        this->parents_population_[index] = parent;
      }

      void set_parents_fitness(const double fitness, const size_t index) {
        assert(this->parents_fitness_.size() > index);
        this->parents_fitness_[index] = fitness;
      }

      void add_parents_population(const vector<int>& parent) {
        this->parents_population_.push_back(parent);

      }
      void clear_parents_population() {
        this->parents_population_.clear();
      }

      void add_parents_fitness(const double fitness) {
        this->parents_fitness_.push_back(fitness);
      }

      void clear_parents_fitness() {
        this->parents_fitness_.clear();
      }

      void set_offspring_population(const vector< vector<int> >& offspring_population) {
        this->offspring_population_ = offspring_population;
      }

      void set_offspring_fitness(const vector<double>& offspring_fitness) {
        this->offspring_fitness_ = offspring_fitness;
      }

      void set_offspring_population(const vector<int>& offspring, const size_t index) {
        assert(this->offspring_population_.size() > index);
        this->offspring_population_[index] = offspring;
      }

      void set_offspring_fitness(const double fitness, const size_t index) {
        assert(this->offspring_fitness_.size() > index);
        this->offspring_fitness_[index] = fitness;
      }

      void add_offspring_population(const vector<int>& parent) {
        this->offspring_population_.push_back(parent);

      }
      void clear_offspring_population() {
        this->offspring_population_.clear();
      }

      void add_offspring_fitness(const double fitness) {
        this->offspring_fitness_.push_back(fitness);
      }

      void clear_offspring_fitness() {
        this->offspring_fitness_.clear();
      }

      void set_best_found_fitness(const double best_found_fitness) {
        this->best_found_fitness_ = best_found_fitness;
      }

      void set_best_individual(const vector <int> best_individual) {
        this->best_individual_ = best_individual;
      }

      void set_generation(const size_t generation) {
        this->generation_ = generation;
      }

      void update_generation() {
        this->generation_++;
      }

      int get_mu() const {
        return this->mu_;
      }

      int get_lambda() const {
        return this->lambda_;
      }

      int get_dimension() const {
        return this->problem_->meta_data().n_variables;
      }

      double get_crossover_probability() const {
        return this->crossover_probability_;
      }

      int get_crossover_mutation_r() const {
        return this->crossover_mutation_r_;
      }

      vector< vector<int> > get_parents_population() const {
        return this->parents_population_;
      }

      vector<double> get_parents_fitness() const {
        return this->parents_fitness_;
      }

      vector< vector<int> >  get_offspring_population() const {
        return this->offspring_population_;
      }

      vector<double> get_offspring_fitness() const {
        return this->offspring_fitness_;
      }

      double get_best_found_fitness() const {
        return this->best_found_fitness_;
      }

      vector <int> get_best_individual() const {
        return this->best_individual_;
      }

      size_t get_generation() const {
        return this->generation_;
      }

      size_t get_independent_runs() const {
        return this->independent_runs_;
      }

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

    private:
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
      std::shared_ptr<ioh::problem::Integer> problem_;
      std::shared_ptr<ioh::logger::Analyzer> csv_logger_;
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

} // namespace modularGA
