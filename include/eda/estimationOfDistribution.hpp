/// \file estimationOfDistribution.h
/// \brief Header file for class EstimationOfDistribution.
///
/// A class of estimation of distribution algorithm.
///
/// \author Furong Ye
/// \date 2020-11-25

#ifndef _ESTIMATION_OF_DISTRIBUTION_H
#define _ESTIMATION_OF_DISTRIBUTION_H

#include "utils/common.hpp"
namespace modularGA
{

  namespace eda
  {
    class EstimationOfDistribution
    {
    public:
      EstimationOfDistribution(int mu = 25, int lambda = 50) : mu_(mu),
                                                               lambda_(lambda) {}

      ~EstimationOfDistribution() {}
      EstimationOfDistribution(const EstimationOfDistribution &) = delete;
      EstimationOfDistribution &operator=(const EstimationOfDistribution &) = delete;

      /// \fn virtual Terminate()
      /// \brief Terminate condition of the algorithm.
      ///
      /// You can set up terminate condition in this function. By default, the algorithm terminates when the best fitness value is found or the budget is used out.
      virtual bool Termination()
      {
        if (!this->problem_->state().optimum_found && this->evaluation_ < this->evluation_budget_ && this->generation_ <= this->generation_budget_)
        {
          return false;
        }
        else
        {
          return true;
        }
      };

      /// \fn DoEstimationOfDistribution()
      /// \brief Do function of the algorithm.
      ///
      /// The order of processing functions is: Initialization() ->
      /// Loop{SampleIndividuals() -> Selection()-> EstimationDistribution()} -> ~Termination().
      /// The function is virtual, which allows users to implement their own algorithm.
      virtual void DoEstimationOfDistribution()
      {
        double rand;
        this->Preparation();
        this->Initialization();
        while (!this->Termination())
        {
          ++this->generation_;
          if (this->offspring_population_.size() != this->lambda_)
          {
            this->offspring_population_.clear();
            this->offspring_fitness_.clear();
            this->offspring_population_ = vector<vector<int>>(this->lambda_, vector<int>(this->get_dimension(), 0));
            this->offspring_fitness_ = vector<double>(this->lambda_, 0.0);
          }
          for (size_t i = 0; i < this->lambda_; ++i)
          {
            this->SampleIndividual(this->variables_distribution_, this->offspring_population_[i]);
            this->offspring_fitness_[i] = this->Evaluate(this->offspring_population_[i]);
            if (this->Termination())
              break;
          }

          if (this->Termination())
            break;
          this->Selection();
          this->EstimateVariablesDistribution();
        }
      };

      virtual void Initialization()
      {
        vector<double> default_distribution(2, 0.5);
        vector<int> default_count(2, 0);
        for (size_t i = 0; i != this->get_dimension(); ++i)
        {
          this->variables_distribution_.push_back(default_distribution);
          this->variables_count_.push_back(default_count);
        }
      };

      void Preparation()
      {
        this->problem_->reset();
        if (this->csv_logger_ != nullptr)
        {
          this->problem_->attach_logger(*this->csv_logger_);
        }

        this->parents_fitness_.clear();
        this->parents_population_.clear();
        this->offspring_fitness_.clear();
        this->offspring_population_.clear();
        this->evaluation_ = 0;
        this->generation_ = 0;
        this->best_individual_ = vector<int>(this->problem_->meta_data().n_variables);
        if (Opt == optimizationType::MAXIMIZATION)
        {
          this->best_found_fitness_ = numeric_limits<double>::lowest();
        }
        else
        {
          this->best_found_fitness_ = numeric_limits<double>::max();
        }

        this->variables_distribution_.clear();
        this->variables_count_.clear();
      };

      double Evaluate(vector<int> &x)
      {
        double result;
        result = (*this->problem_)(x);

        ++this->evaluation_;

        if (Opt == optimizationType::MAXIMIZATION)
        {
          if (result > this->best_found_fitness_)
          {
            this->best_found_fitness_ = result;
            this->best_individual_ = x;
          }
        }
        else
        {
          if (result < this->best_found_fitness_)
          {
            this->best_found_fitness_ = result;
            this->best_individual_ = x;
          }
        }
        return result;
      };

      void EstimateVariablesDistribution()
      {
        for (size_t i = 0; i != this->get_dimension(); ++i)
        {
          for (size_t j = 0; j != this->variables_count_[i].size(); ++j)
          {
            this->variables_count_[i][j] = 0;
          }

          for (size_t j = 0; j != this->mu_; ++j)
          {
            this->variables_count_[i][this->parents_population_[j][i]]++;
          }

          for (size_t j = 0; j != this->variables_distribution_[i].size(); ++j)
          {
            this->variables_distribution_[i][j] = static_cast<double>(this->variables_count_[i][j]) / static_cast<double>(this->mu_);
            // Restriction for p
            this->variables_distribution_[i][j] = this->variables_distribution_[i][j] > (1.0 - 1.0 / static_cast<double>(this->get_dimension())) ? (1.0 - 1.0 / static_cast<double>(this->get_dimension())) : this->variables_distribution_[i][j];
            this->variables_distribution_[i][j] = this->variables_distribution_[i][j] < 1.0 / static_cast<double>(this->get_dimension()) ? (1.0 / static_cast<double>(this->get_dimension())) : this->variables_distribution_[i][j];
          }
        }
      };

      void SampleIndividual(const vector<vector<double>> &distributions, vector<int> &individual)
      {
        if (individual.size() != this->get_dimension())
        {
          individual.clear();
          individual = vector<int>(this->get_dimension(), 0);
        }

        size_t d = static_cast<size_t>(this->get_dimension());
        assert(distributions.size() == d);
        int j = 0;
        double r, p;
        for (size_t i = 0; i != d; ++i)
        {
          r = uniform_random();
          j = 0;
          p = 0;
          while (j < distributions[i].size())
          {
            if (r > p && r < p + distributions[i][j])
            {
              break;
            }

            p += distributions[i][j];
            ++j;
          }
          individual[i] = j;
        }
      };
      void SetSeed(unsigned seed)
      {
        random_gen.seed(seed);
      };

      void AssignProblem(shared_ptr<ioh::problem::IntegerSingleObjective> problem_ptr)
      {
        this->problem_ = problem_ptr;
      };

      void AssignLogger(shared_ptr<ioh::logger::Analyzer> logger_ptr)
      {
        this->csv_logger_ = logger_ptr;
      };

      virtual void Selection()
      {
        if (this->parents_population_.size() != this->mu_)
        {
          this->parents_population_.clear();
          this->parents_fitness_.clear();
          this->parents_population_ = vector<vector<int>>(this->mu_, vector<int>(this->get_dimension(), 0));
          this->parents_fitness_ = vector<double>(this->mu_, 0);
        }

        vector<int> index(this->lambda_, 0);
        assert(this->lambda_ == this->offspring_fitness_.size());
        for (int i = 0; i != this->lambda_; ++i)
        {
          index[i] = i;
        }

        assert(this->mu_ <= this->lambda_);
        if (Opt == MAXIMIZATION)
        {
          partial_sort(index.begin(), index.begin() + this->mu_, index.end(),
                       [&](const int a, const int b)
                       { return this->offspring_fitness_[a] >= this->offspring_fitness_[b]; });
        }
        else
        {
          partial_sort(index.begin(), index.begin() + this->mu_, index.end(),
                       [&](const int a, const int b)
                       { return this->offspring_fitness_[a] < this->offspring_fitness_[b]; });
        }

        for (int i = 0; i != this->mu_; ++i)
        {
          this->parents_population_[i] = this->offspring_population_[index[i]];
          this->parents_fitness_[i] = this->offspring_fitness_[index[i]];
        }
      };

      void set_mu(const int mu)
      {
        this->mu_ = mu;
      }

      void set_lambda(const int lambda)
      {
        this->lambda_ = lambda;
      }

      void set_evaluation_budget(const size_t evaluation_budget)
      {
        this->evluation_budget_ = evaluation_budget;
      }

      void set_generation_budget(const size_t generation_budget)
      {
        this->generation_budget_ = generation_budget;
      }

      void set_independent_runs(const size_t independent_runs)
      {
        this->independent_runs_ = independent_runs;
      }

      void set_parents_population(const vector<vector<int>> parents_population)
      {
        this->parents_population_ = parents_population;
      }

      void set_parents_fitness(const vector<double> &parents_fitness)
      {
        this->parents_fitness_ = parents_fitness;
      }

      void set_parents_population(const vector<int> &parent, const size_t index)
      {
        assert(this->parents_population_.size() > index);
        this->parents_population_[index] = parent;
      }

      void set_parents_fitness(const double fitness, const size_t index)
      {
        assert(this->parents_fitness_.size() > index);
        this->parents_fitness_[index] = fitness;
      }

      void add_parents_population(const vector<int> &parent)
      {
        this->parents_population_.push_back(parent);
      }
      void clear_parents_population()
      {
        this->parents_population_.clear();
      }

      void add_parents_fitness(const double fitness)
      {
        this->parents_fitness_.push_back(fitness);
      }

      void clear_parents_fitness()
      {
        this->parents_fitness_.clear();
      }

      void set_offspring_population(const vector<vector<int>> &offspring_population)
      {
        this->offspring_population_ = offspring_population;
      }

      void set_offspring_fitness(const vector<double> &offspring_fitness)
      {
        this->offspring_fitness_ = offspring_fitness;
      }

      void set_offspring_population(const vector<int> &offspring, const size_t index)
      {
        assert(this->offspring_population_.size() > index);
        this->offspring_population_[index] = offspring;
      }

      void set_offspring_fitness(const double fitness, const size_t index)
      {
        assert(this->offspring_fitness_.size() > index);
        this->offspring_fitness_[index] = fitness;
      }

      void set_variables_distribution(const vector<vector<double>> &variables_distribution)
      {
        this->variables_distribution_ = variables_distribution;
      }

      void set_variables_distribution(const vector<double> &variables_distribution, size_t index)
      {
        this->variables_distribution_[index] = variables_distribution;
      }

      void add_offspring_population(const vector<int> &parent)
      {
        this->offspring_population_.push_back(parent);
      }
      void clear_offspring_population()
      {
        this->offspring_population_.clear();
      }

      void add_offspring_fitness(const double fitness)
      {
        this->offspring_fitness_.push_back(fitness);
      }

      void clear_offspring_fitness()
      {
        this->offspring_fitness_.clear();
      }

      void set_best_found_fitness(const double best_found_fitness)
      {
        this->best_found_fitness_ = best_found_fitness;
      }

      void set_best_individual(const vector<int> best_individual)
      {
        this->best_individual_ = best_individual;
      }

      void set_generation(const size_t generation)
      {
        this->generation_ = generation;
      }

      void update_generation()
      {
        this->generation_++;
      }

      int get_mu() const
      {
        return this->mu_;
      }

      int get_lambda() const
      {
        return this->lambda_;
      }

      int get_dimension() const
      {
        return this->problem_->meta_data().n_variables;
      }

      vector<vector<int>> get_parents_population()
      {
        return this->parents_population_;
      }

      vector<double> get_parents_fitness()
      {
        return this->parents_fitness_;
      }

      vector<vector<int>> get_offspring_population()
      {
        return this->offspring_population_;
      }

      vector<double> get_offspring_fitness()
      {
        return this->offspring_fitness_;
      }

      double get_best_found_fitness()
      {
        return this->best_found_fitness_;
      }

      vector<int> get_best_individual()
      {
        return this->best_individual_;
      }

      size_t get_generation()
      {
        return this->generation_;
      }

      vector<vector<double>> get_variables_distribution()
      {
        return this->variables_distribution_;
      }

      vector<double> get_variables_distribution(size_t index)
      {
        return this->variables_distribution_[index];
      }

      /// \fn run()
      /// \brief Run the algorithm once with given parameters.
      void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
      {
        string algorithm_name = "UMDA";
        std::shared_ptr<ioh::logger::Analyzer> logger(new ioh::logger::Analyzer(
            {ioh::trigger::on_improvement}, // trigger when the objective value improves
            {},                             // no additional properties
            folder_path,                    // path to store data
            folder_name,                    // name of the folder in path, which will be newly created
            algorithm_name,                 // name of the algoritm
            algorithm_name,                 // additional info about the algorithm
            false                           // where to store x positions in the data files
            ));
        this->AssignLogger(logger);

        this->set_evaluation_budget(eval_budget);
        this->set_generation_budget(gene_budget);
        this->set_independent_runs(independent_runs);
        this->SetSeed(rand_seed);

        this->run(suite);
      };

      void run(shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite)
      {
        for (auto &p : *suite)
        {
          this->problem_ = p;
          size_t i = 0;
          while (++i <= this->independent_runs_)
          {
            size_t i = 0;
            this->DoEstimationOfDistribution();
          }
          this->problem_->reset();
        }
      };

    private:
      int mu_;     /// < parents population size
      int lambda_; /// < offspring population size

      vector<vector<int>> parents_population_;
      vector<double> parents_fitness_;
      vector<vector<int>> offspring_population_;
      vector<double> offspring_fitness_;
      double best_found_fitness_;
      vector<int> best_individual_;

      vector<vector<double>> variables_distribution_;
      vector<vector<int>> variables_count_;

      size_t evaluation_;        /// < evaluation times
      size_t generation_;        /// < number of iterations/generations
      size_t evluation_budget_;  /// < budget for evaluations
      size_t generation_budget_; /// < budget for generations

      size_t independent_runs_; /// < number of independent runs.

      /// TODO: we assume the type of problem are integer only now.
      shared_ptr<ioh::problem::IntegerSingleObjective> problem_;
      shared_ptr<ioh::logger::Analyzer> csv_logger_;
    };
  }
}

#endif // _ESTIMATION_OF_DISTRIBUTION_H
