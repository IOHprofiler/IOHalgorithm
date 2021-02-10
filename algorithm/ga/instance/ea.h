#include "geneticAlgorithm.h"

class staticEA : public GeneticAlgorithm {
  
public:
  staticEA(int mu, int lambda, double mutation_rate_scale = 1) {
    this->set_mu(mu);
    this->set_lambda(lambda);
    this->set_crossover_mutation_r("OR");
    this->set_crossover_probability(0);
    this->set_mutation_operator("BINOMIALSAMPLE");
    this->set_selection_operator("BESTPLUS");
    mutation_rate_scale_ = mutation_rate_scale;
  }
  
  void run(string folder_path, string folder_name, shared_ptr<IOHprofiler_suite<int> > suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed){
    string algorithm_name = "(" + to_string(this->get_mu()) + "+" + to_string(this->get_lambda()) + ")>_0 EA";
    std::shared_ptr<IOHprofiler_csv_logger<int>> logger(new IOHprofiler_csv_logger<int>(folder_path,folder_name,algorithm_name,algorithm_name) );
    logger->activate_logger();
    this->AssignLogger(logger);

    this->set_evaluation_budget(eval_budget);
    this->set_generation_budget(gene_budget);
    this->set_independent_runs(independent_runs);
    this->SetSeed(rand_seed);
    shared_ptr<IOHprofiler_problem<int> > problem_ptr;
    while ((problem_ptr = suite->get_next_problem()) != nullptr) {
      this->AssignProblem(problem_ptr);
      size_t r = 0;
      this->set_mutation_rate(this->mutation_rate_scale_ / static_cast<double> (this->get_dimension()));
      while (r < this->get_independent_runs()) {
        this->DoGeneticAlgorithm();
        r++;
      }  
    }
    logger->clear_logger();
  }

  double mutation_rate_scale_;
};
