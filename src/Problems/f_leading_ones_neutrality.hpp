/// \file f_leading_ones_neutrality.hpp
/// \brief cpp file for class f_leading_ones_neutrality.
///
/// This file implements a LeadingOnes problem with neutrality transformation method from w-model.
/// The parameter mu is chosen as 3.
///
/// \author Furong Ye
/// \date 2019-06-27
#ifndef _F_LEADING_ONES_NEUTRALITY_H
#define _F_LEADING_ONES_NEUTRALITY_H

#include "../Template/IOHprofiler_problem.hpp"
#include "common_used_functions/wmodels.hpp"

class LeadingOnes_Neutrality : public IOHprofiler_problem<int> {
public:
   LeadingOnes_Neutrality() {

    IOHprofiler_set_problem_name("LeadingOnes_Neutrality");
    IOHprofiler_set_problem_type("pseudo_Boolean_problem");
    IOHprofiler_set_number_of_objectives(1);
    IOHprofiler_set_lowerbound(0);
    IOHprofiler_set_upperbound(1);
    IOHprofiler_set_best_variables(1);
  }
  LeadingOnes_Neutrality(int instance_id, int dimension) {

    IOHprofiler_set_instance_id(instance_id);
    IOHprofiler_set_problem_name("LeadingOnes_Neutrality");
    IOHprofiler_set_problem_type("pseudo_Boolean_problem");
    IOHprofiler_set_number_of_objectives(1);
    IOHprofiler_set_lowerbound(0);
    IOHprofiler_set_upperbound(1);
    IOHprofiler_set_best_variables(1);
    Initilize_problem(dimension);
  }

  ~LeadingOnes_Neutrality() {};

  void Initilize_problem(int dimension) {
    IOHprofiler_set_number_of_variables(dimension);

  };

  double internal_evaluate(const std::vector<int> &x) {
    std::vector<int> new_variables = neutrality(x,3);
    int n = new_variables.size();
    int result = 0;
    for (int i = 0; i != n; ++i) {
      if (new_variables[i] == 1) {
        result = i + 1;
      } else {
        break;
      }
    }
    return (double)result;
  };

  static LeadingOnes_Neutrality * createInstance() {
    return new LeadingOnes_Neutrality();
  };

  static LeadingOnes_Neutrality * createInstance(int instance_id, int dimension) {
    return new LeadingOnes_Neutrality(instance_id, dimension);
  };
};

#endif