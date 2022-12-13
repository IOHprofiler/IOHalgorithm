#pragma once

#include "ga/geneticAlgorithm.hpp"
namespace modularGA
{

  namespace ga
  {
    class FastGA : public GeneticAlgorithm
    {

    public:
      FastGA(int mu, int lambda)
      {
        this->set_mu(mu);
        this->set_lambda(lambda);
        this->set_crossover_mutation_r("OR");
        this->set_crossover_probability(0);
        this->set_mutation_operator("POWERLAWSAMPLE");
        this->set_beta_f(1.5);
        this->set_selection_operator("BESTPLUS");
      }

      void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
      {
        string algorithm_name = "(" + to_string(this->get_mu()) + "+" + to_string(this->get_lambda()) + ") fast GA";
         std::shared_ptr<ioh::logger::Analyzer > logger(new ioh::logger::Analyzer(
                                                      {ioh::trigger::on_improvement}, // trigger when the objective value improves
                                                      {},                   // no additional properties 
                                                      folder_path,        // path to store data
                                                      folder_name,               // name of the folder in path, which will be newly created
                                                      algorithm_name,                     // name of the algoritm 
                                                      algorithm_name,                   // additional info about the algorithm              
                                                      false            // where to store x positions in the data files 
                                                    ));                                                  
          this->AssignLogger(logger);

        this->set_evaluation_budget(eval_budget);
        this->set_generation_budget(gene_budget);
        this->set_independent_runs(independent_runs);
        this->SetSeed(rand_seed);

        this->run_N(suite);
      }
    };
  }
}
