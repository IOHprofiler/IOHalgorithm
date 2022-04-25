#pragma once
#include "ga/geneticAlgorithm.hpp"
namespace modularGA
{

  namespace ga
  {

      class oneLambdaLambdaEA : public GeneticAlgorithm
      {

      public:
        oneLambdaLambdaEA(int lambda)
        {
          this->set_mu(1);
          this->set_lambda(lambda);
          this->set_mutation_operator("BINOMIALSAMPLE");
          this->set_crossover_operator("UNIFORMCROSSOVER");
          this->set_selection_operator("BESTPLUS");
        }

        void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::Integer> > suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
        {
          string algorithm_name = "(1+(" + to_string(this->get_lambda()) + "," + to_string(this->get_lambda()) + "))>_0 EA";
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

        double a = pow(1.5, 0.25); /// < parameter for adjusting lambda
        double b = 2.0 / 3.0;      /// < parameter for adjusting lambda

        void DoGeneticAlgorithm()
        {
          double rand, best_f = numeric_limits<double>::lowest(), best_mutation_f;
          vector<int> mutation_offspring, best_ind, x;
          bool update_lambda_flag;
          int lambda, dimension;
          double mutation_rate, p_u;
          int mutation_strength;

          this->Preparation();

          lambda = this->get_lambda();
          dimension = this->get_dimension();

          this->Initialization();

          while (!this->Termination())
          {

            this->update_generation();
            update_lambda_flag = false;
            best_mutation_f = numeric_limits<double>::lowest();
            this->clear_offspring_population();
            this->clear_offspring_fitness();

            mutation_rate = static_cast<double>(lambda) / static_cast<double>(dimension);
            this->set_p_u(1.0 / static_cast<double>(lambda));

            /**Mutation stage
             */
            mutation_strength = this->SampleConditionalBinomial(mutation_rate, dimension);
            for (size_t i = 0; i < this->get_lambda(); ++i)
            {
              x = this->get_parents_population()[0];

              this->Flip(x, mutation_strength);
              this->add_offspring_fitness(this->Evaluate(x));

              if (this->get_offspring_fitness()[i] >= best_f)
              {
                best_ind = x;
                best_f = this->get_offspring_fitness()[i];
                if (this->get_offspring_fitness()[i] > best_f)
                {
                  update_lambda_flag = true;
                }
              }

              if (this->get_offspring_fitness()[i] > best_mutation_f)
              {
                best_mutation_f = this->get_offspring_fitness()[i];
                mutation_offspring = x;
              }

              if (this->Termination())
                break;
            }

            /** Crossover Stage
             */
            for (size_t i = 0; i < this->get_lambda(); ++i)
            {

              this->DoCrossover(x, this->get_parents_population()[0], mutation_offspring);
              if (this->c_flipped_index.size() == 0)
              {
                this->set_offspring_fitness(this->get_parents_fitness()[0], i);
              }
              else if (x == mutation_offspring)
              {
                this->set_offspring_fitness(best_mutation_f, i);
              }
              else
              {
                this->set_offspring_fitness(this->Evaluate(x), i);
              }

              if (this->get_offspring_fitness()[i] >= best_f)
              {
                best_ind = x;
                best_f = this->get_offspring_fitness()[i];
                if (this->get_offspring_fitness()[i] > best_f)
                {
                  update_lambda_flag = true;
                }
              }

              if (this->Termination())
                break;
            }

            if (this->Termination())
              break;

            this->set_parents_population(best_ind, 0);
            this->set_parents_fitness(best_f, 0);
            if (update_lambda_flag == true)
            {
              lambda = (lambda * b) > 1 ? (lambda * b) : 1;
            }
            else
            {
              lambda = (lambda * a) < dimension ? (lambda * a) : dimension;
            }
          }
        }
      };
    
  }
}