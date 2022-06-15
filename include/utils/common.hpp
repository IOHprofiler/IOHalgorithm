#pragma once
/// \file common.h
/// \brief Functions commonly used by the projects, including random generators (from c++ random class), and sampling functions.
///
/// A configurable genetic algorithm.
///
/// \author Furong Ye
/// \date 2020-07-02

#ifndef _COMMON_H_
#define _COMMON_H_

#include <iostream>
#include <string>
#include <cstddef>
#include <vector>
#include <memory>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>
#include <assert.h>

#include "ioh.hpp"

using namespace std;

static default_random_engine random_gen(1);
static normal_distribution<double> normal_dis(0,1);
static uniform_real_distribution<double> uniform_dis(0.0,1.0);

enum optimizationType {
  MINIMIZATION = 0,
  MAXIMIZATION = 1
};

static optimizationType Opt = MAXIMIZATION;

static double normal_random() {
  return normal_dis(random_gen);
}

static double uniform_random() {
  return uniform_dis(random_gen);
}

/// \fn sampleNFromM
/// \brief Sampling n different indexes from length m
static void sampleNFromM(vector<size_t> &sampled_number, size_t n, size_t m) {
  if (sampled_number.size() != 0) {
    sampled_number.clear();
  }
  
  if (n == 0) {
    std::clog <<  "sampled zero number" << std::endl;
  }
  
  size_t randPos;
  sampled_number.reserve(n);
  
  if (n > m/2) { /// If n is larger than m/2, we sample random indexes by reordering a permutation.
    vector<size_t> population;
    population.reserve(m);
    for (size_t i = 0; i < m; ++i) {
      population.push_back(i);
    }
    
    int temp;
    for (size_t i = m-1; i > 0; --i) {
      randPos = static_cast<size_t>( floor(uniform_random() * (i+1)) );
      temp = population[i];
      population[i] = population[randPos];
      population[randPos] = temp;
      sampled_number.push_back(population[i]);
      if(m-i-1 == n-1) {
        break;
      }
    }
    if(n == m) {
      sampled_number.push_back(population[0]);
    }
  } else { /// If n is smaller than m/2, we sample indexes repeatly until getting different values.
    bool resample = false;
    for (size_t i = 0; i != n; ++i) {
      do {
        resample = false;
        randPos =  static_cast<size_t>( floor(uniform_random() * m) );
        for (size_t j = 0; j != i; ++j) {
          if(randPos == sampled_number[j]) {
            resample = true;
            break;
          }
        }
      } while (resample);
      sampled_number.push_back(randPos);
    }
  }
};

#endif // _COMMON_H_
