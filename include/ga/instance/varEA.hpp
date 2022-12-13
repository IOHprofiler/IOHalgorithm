#pragma once
#include "ga/geneticAlgorithm.hpp"
namespace modularGA
{

  namespace ga
  {
    class VarEA : public GeneticAlgorithm
    {

    public:
      VarEA(int lambda, double init_r_n_ = 1.0)
      {
        this->set_mu(1);
        this->set_lambda(lambda);
        this->set_mutation_operator("NORMALSAMPLE");
        this->set_r_n(init_r_n_);
      }

      void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
      {
        string algorithm_name = "(1+" + to_string(this->get_lambda()) + ")-varEA>0";
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

        this->run_N(suite);
      }

      void DoGeneticAlgorithm()
      {
        vector<int> best_ind, x;
        double r_n = this->get_r_n();
        double best_f_x, f_x;
        int flip_n, best_r_n;
        double F = init_F;
        int c = 0;
        this->Preparation();
        this->Initialization();


        while (!this->Termination())
        {
          this->update_generation();
          best_r_n = r_n;
          best_f_x = std::numeric_limits<int>::min();
          this->set_r_n(r_n);
          this->set_sigma_n(sqrt(pow(F, c) *
                                 static_cast<double>(r_n) *
                                 (1.0 - r_n / static_cast<double>(this->get_dimension()))));
          for (size_t i = 0; i < this->get_lambda(); ++i)
          {
            x = this->get_parents_population()[0];
            flip_n = this->DoMutation(x);
            f_x = this->Evaluate(x);
            if (f_x > best_f_x)
            {
              best_f_x = f_x;
              r_n = flip_n;
              best_ind = x;
            }

            if (this->Termination())
              break;
          }

          if (this->Termination())
            break;

          if (best_r_n == r_n)
          {
            c++;
          }
          else
          {
            c = 0;
          }

          if (best_f_x >= this->get_parents_fitness()[0])
          {
            this->set_parents_population(best_ind, 0);
            this->set_parents_fitness(best_f_x, 0);
          }
        }
      }
      double init_F = 0.98;
    };
  }
}