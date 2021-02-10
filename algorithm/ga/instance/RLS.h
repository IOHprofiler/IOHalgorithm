#include "geneticAlgorithm.h"

class RLS : public GeneticAlgorithm {
  
public:
  RLS() {
    this->set_mu(1);
    this->set_lambda(1);
    this->set_crossover_mutation_r("IND");
    this->set_crossover_probability(0);
    this->set_mutation_operator("STATICSAMPLE");
    this->set_l(1);
    this->set_selection_operator("BESTPLUS");
  }
  
  void run(string folder_path, string folder_name, shared_ptr<IOHprofiler_suite<int> > suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
  {
    string algorithm_name = "RLS";
    std::shared_ptr<IOHprofiler_csv_logger<int>> logger(new IOHprofiler_csv_logger<int>(folder_path,folder_name,algorithm_name,algorithm_name) );
    logger->activate_logger();
    this->AssignLogger(logger);
    

    this->set_evaluation_budget(eval_budget);
    this->set_generation_budget(gene_budget);
    this->set_independent_runs(independent_runs);
    this->SetSeed(rand_seed);
    
    this->run_N(suite);
    logger->clear_logger();
  }
};
