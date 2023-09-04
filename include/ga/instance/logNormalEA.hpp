#pragma once
#include "algorithms/ga/geneticAlgorithm.hpp"
namespace ioha
{
  namespace alg
  {
    class LogNormalEA : public GeneticAlgorithm
    {

    public:

      LogNormalEA(int lambda, double init_mutation_rate = 0.2)
      {
        this->set_mu(1);
        this->set_lambda(lambda);
        this->set_crossover_mutation_r("OR");
        this->set_crossover_probability(0);
        this->set_mutation_operator("BINOMIALSAMPLE");
        this->set_selection_operator("BESTPLUS");
        this->init_mutation_rate_ = init_mutation_rate;
      }

      void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
      {
        string algorithm_name = "(1+" + to_string(this->get_lambda()) + ")-lognormal-EA>0";
        ioh::logger::Analyzer logger(
            {ioh::trigger::on_improvement}, // trigger when the objective value improves
            {},                             // no additional properties
            folder_path,                    // path to store data
            folder_name,                    // name of the folder in path, which will be newly created
            algorithm_name,                 // name of the algoritm
            algorithm_name,                 // additional info about the algorithm
            false                           // where to store x positions in the data files
            );
        suite->attach_logger(logger);
        this->set_evaluation_budget(eval_budget);
        this->set_generation_budget(gene_budget);
        this->set_independent_runs(independent_runs);
        this->set_seed(rand_seed);

        this->run_N(suite);
      }

      void opt()
      {
        vector<int> best_ind, x;
        double tmp_mr, tmp_best_mr;
        double f_x, best_f_x;

        this->preparation();
        this->initialization();
        
        double mutation_rate = this->init_mutation_rate_;
        while (!this->termination())
        {
          this->update_generation();

          best_f_x = std::numeric_limits<int>::min();
          for (size_t i = 0; i < this->get_lambda(); ++i)
          {
            x = this->get_parents_population()[0];
            
            tmp_mr = 1.0 / (1.0 + (((1.0 - mutation_rate) / mutation_rate) * exp(0.22 * operators::normal_random())));
            this->set_mutation_rate(tmp_mr);
            this->do_mutation(x);
           
            

            f_x = this->evaluate(x);
            if (f_x > best_f_x)
            {
              best_f_x = f_x;
              best_ind = x;
              tmp_best_mr = tmp_mr;
            }

            if (this->termination())
              break;
          }

          if (this->termination())
            break;


          if (best_f_x >= this->get_parents_fitness()[0])
          {
            this->set_parents_population(best_ind, 0);
            this->set_parents_fitness(best_f_x, 0);
          }
          mutation_rate = tmp_best_mr;
        }
      }

      double init_mutation_rate_;
    };
  }
}
