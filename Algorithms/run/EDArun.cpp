#include "../../src/Template/Experiments/IOHprofiler_experimenter.hpp"
#include "../EstimationOfDistribution.hpp"

int DEFAULT_POPULATION_SIZE = 50;
double DEFAULT_SELECT_RATE = 0.5;

int budget;
int budget_scale = 5;
int budget_power = 2;

IOHprofiler_random random_generator(1);

void _run_experiment() {
  std::string configName = "./EDAconfiguration.ini";
  IOHprofiler_experimenter<int> experimenter(configName,estimationOfDistribution);
  experimenter._set_independent_runs(11);
  experimenter._run();
}

int main(){
  _run_experiment();
  return 0;
}