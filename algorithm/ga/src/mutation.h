/// \file mutation.h
/// \brief Header file for class Mutation.
///
/// It contains functions of mutation operators.
/// 
/// \author Furong Ye
/// \date 2020-07-02

#ifndef _MUTATION_H_
#define _MUTATION_H_

#include "common.h"

#define DEFAULT_MUTATION_STRENGTH_ 1
#define DEFAULT_MUTATION_RATE_ 0.01
#define DEFAULT_NORMALIZED_MUTATION_STRENGTH_MEAN_ 1
#define DEFAULT_NORMALIZED_MUTATION_STRENGTH_SD_ 0.01
#define DEFAULT_FAST_MUTATION_BETA_ 1.5


/// Definition of mutation operator id.
enum mutation_operator {
  STATICSAMPLE = 1, /// < mutation strength is a fixed value.
  BINOMIALSAMPLE = 2, /// < sampling mutation strength from a binomial distribution.
  NORMALSAMPLE = 3, /// < sampling mutation strength from a normal distribution.
  POWERLAWSAMPLE = 4 /// < sampling mutation strength from a power-law distribution.
};

class Mutation {
public:
  Mutation() :
  mutation_operator_(2),
  l_(DEFAULT_MUTATION_STRENGTH_),
  mutation_rate_(DEFAULT_MUTATION_RATE_),
  r_n_(DEFAULT_NORMALIZED_MUTATION_STRENGTH_MEAN_),
  sigma_n_(DEFAULT_NORMALIZED_MUTATION_STRENGTH_SD_),
  beta_f_(DEFAULT_FAST_MUTATION_BETA_) {}
  
  ~Mutation() {}
  Mutation(const Mutation&) = delete;
  Mutation &operator = (const Mutation&) = delete;
  
  void DoMutation(vector<int> &y);
  
  void Flip(vector<int> &y, const int l);
  int SampleL(const int N);
  int SampleBinomial(const double p, const int n);
  int SampleConditionalBinomial(const double p, const int n);
  int SampleNormal(double mu, double sigma);
  int SampleConditionalNormal(double mu, double sigma, int upperbound);
  int SampleLogNormal(double mu, double sigma);
  int SampleConditionalLogNormal(double mu, double sigma, int upperbound);
  int SampleFromDistribution(const std::vector<double> &distribution);
  void PowerLawDistribution(int N);
  
  void set_mutation_operator(const int m);
  void set_mutation_operator(string m);
  void set_l(const int l);
  void set_mutation_rate(const double mutation_rate);
  void set_r_n(const double r_n);
  void set_sigma_n(const double sigma_n);
  void set_beta_f(const double beta_f);
  
  int get_mutation_operator() const;
  int get_l() const;
  double get_mutation_rate() const;
  double get_r_n() const;
  double get_sigma_n() const;
  double get_beta_f() const;
  
  vector <size_t> m_flipped_index;
  
private:
  std::vector<double> power_law_distribution_;
  
  int mutation_operator_; /// < a flag for operator to be used for mutation
  int l_; /// < a static mutation strength
  double mutation_rate_; /// < mutation rate for standard bit mutation
  double r_n_; /// < mean value for normalized bit mutation
  double sigma_n_; /// < sd for normalized bit mutation
  double beta_f_; /// < beta for fast mutation
};

#endif // _MUTATION_H_
