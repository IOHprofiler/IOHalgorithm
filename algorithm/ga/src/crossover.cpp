#include "crossover.h"

void Crossover::DoCrossover(vector<int> &y, const vector<int> &x1, const vector<int> &x2) {
  switch (this->crossover_operator_) {
    case 1:
      this->UniformCrossover(y,x1,x2);
      break;
    case 2:
      this->OnePointCrossover(y,x1,x2);
      break;
    case 3:
      this->TwoPointCrossover(y,x1,x2);
      break;
    default:
      cerr << "unknown crossover operator" << endl;
      assert(false);
  }
}

void Crossover::UniformCrossover(vector<int> &y, const vector<int> &x1, const vector<int> &x2) {
  assert(x1.size() == x2.size());
  
  size_t n = x1.size();
  if (y.size() != 0) {
    y.clear();
    y.reserve(n);
  }
  this->c_flipped_index.clear();
  
  for (size_t i = 0; i != n; ++i) {
    if (uniform_random() < this->p_u_) {
      y.push_back(x2[i]);
      if (x2[i] != x1[i]) {
        this->c_flipped_index.push_back(i);
      }
    } else {
      y.push_back(x1[i]);
    }
  }
}

void Crossover::OnePointCrossover(vector<int> &y, const vector<int> &x1, const vector<int> &x2) {
  assert(x1.size() == x2.size());
  
  size_t n = x1.size();
  if (y.size() != 0) {
    y.clear();
    y.reserve(n);
  }
  this->c_flipped_index.clear();
  
  size_t point = (size_t)(uniform_random() * n);
  for (size_t i = 0; i != n; ++i) {
    if (i < point) {
      y.push_back(x1[i]);
    } else {
      y.push_back(x2[i]);
      if (x2[i] != x1[i]) {
        this->c_flipped_index.push_back(i);
      }
    }
  }
}

void Crossover::TwoPointCrossover(vector<int> &y, const vector<int> &x1, const vector<int> &x2) {
  assert(x1.size() == x2.size());
  
  size_t n = x1.size();
  if (y.size() != 0) {
    y.clear();
    y.reserve(n);
  }
  this->c_flipped_index.clear();
  
  vector<size_t> twopoints;
  twopoints.push_back((size_t)(uniform_random() * n));
  twopoints.push_back((size_t)(uniform_random() * n));
  while(twopoints[0] == twopoints[1]) {
    twopoints[1] = (size_t)(uniform_random() * n);
  }
  if (twopoints[0] > twopoints[1]) {
    size_t temp = twopoints[1];
    twopoints[1] = twopoints[0];
    twopoints[0] = temp;
  }
  
  for (size_t i = 0; i != n; ++i) {
    if (i >= twopoints[0]  && i <= twopoints[1]) {
      y.push_back(x2[i]);
      if (x2[i] != x1[i]) {
        this->c_flipped_index.push_back(i);
      }
    } else {
      y.push_back(x1[i]);
    }
  }
}

void Crossover::set_crossover_operator(const int c) {
  assert( (c >= 1) && (c <= 3) );
  this->crossover_operator_ = c;
}

void Crossover::set_crossover_operator(string c) {
  transform(c.begin(), c.end(), c.begin(), ::toupper);
  if (c == "UNIFORMCROSSOVER") {
    this->set_crossover_operator(UNIFORMCROSSOVER);
  } else if (c == "ONEPOINTCROSSOVER") {
    this->set_crossover_operator(ONEPOINTCROSSOVER);
  } else if (c == "TWOPOINTCROSSOVER") {
    this->set_crossover_operator(TWOPOINTCROSSOVER);
  } else {
    cerr << "invalid name for crossover operator" << endl;
    assert(false);
  }
}

void Crossover::set_p_u(const double p_u) {
  this->p_u_ = p_u;
}

int Crossover::get_crossover_operator() const {
  return this->crossover_operator_;
}

double Crossover::get_p_u() const {
  return this->p_u_;
}
