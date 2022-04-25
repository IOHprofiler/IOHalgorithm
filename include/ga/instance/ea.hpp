#pragma once
#include "ga/geneticAlgorithm.hpp"

namespace modularGA
{

  namespace ga
  {
    class staticEA : public GeneticAlgorithm
    {

    public:
      staticEA(int mu, int lambda, double mutation_rate_scale = 1)
      {
        this->set_mu(mu);
        this->set_lambda(lambda);
        this->set_crossover_mutation_r("OR");
        this->set_crossover_probability(0);
        this->set_mutation_operator("BINOMIALSAMPLE");
        this->set_selection_operator("BESTPLUS");
        mutation_rate_scale_ = mutation_rate_scale;
      }

      void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::Integer>> suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
      {
        string algorithm_name = "(" + to_string(this->get_mu()) + "+" + to_string(this->get_lambda()) + ")>_0 EA";
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
        for (auto &problem_ptr : *suite)
        {
          this->AssignProblem(problem_ptr);
          size_t r = 0;
          this->set_mutation_rate(this->mutation_rate_scale_ / static_cast<double>(this->get_dimension()));
          while (r < this->get_independent_runs())
          {
            this->DoGeneticAlgorithm();
            r++;
          }
          problem_ptr->reset();
        }
      }

      double mutation_rate_scale_;
    };
  }
}
