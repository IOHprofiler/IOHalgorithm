/// \file crossover.h
/// \brief Header file for class Crossover.
///
/// Implementation of 3 crossover operators: uniform crossover, one-point crossover, and two-point crossover.
/// Each operator only returns an individual.
/// Uniform crossover allows to set the probability that switches points from the other individual.
///
/// \author Furong Ye
/// \date 2020-07-02

#ifndef _CROSSOVER_H_
#define _CROSSOVER_H_

#include "common.h"

/// Definition of crossover operator id.
enum crossover_operator {
  UNIFORMCROSSOVER = 1, /// < uniform crossover.
  ONEPOINTCROSSOVER = 2, /// < one-point crossover.
  TWOPOINTCROSSOVER = 3 /// < two-point crossover.
};

class Crossover {
public:
  Crossover() :
  crossover_operator_(1),
  p_u_(0.5) {}
  
  ~Crossover() {}
  Crossover (const Crossover&) = delete;
  Crossover &operator = (const Crossover&) = delete;
  
  void DoCrossover(vector<int> &y, const vector<int> &x1, const vector<int> &x2);
  
  void UniformCrossover(vector<int> &y, const vector<int> &x1, const vector<int> &x2);
  void OnePointCrossover(vector<int> &y, const vector<int> &x1, const vector<int> &x2);
  void TwoPointCrossover(vector<int> &y, const vector<int> &x1, const vector<int> &x2);
  
  void set_crossover_operator(const int c);
  void set_crossover_operator(string c);
  void set_p_u(const double p_u);
  
  int get_crossover_operator() const;
  double get_p_u() const;
  
  vector<size_t> c_flipped_index; /// < recording indexes where bits flip.
  
private:
  int crossover_operator_; /// < a flag for operator to be used for crossover. 1: uniform crossover; 2: one-point crossover, 3: two-point crossover.
  double p_u_; /// < probability that a bit to be replaced by the bit of the other parent in uniform crossover
};

#endif // _CROSSOVER_H_
