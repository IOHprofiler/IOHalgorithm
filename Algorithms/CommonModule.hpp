#ifndef _COMMONMODULE_HPP
#define _COMMONMODULE_HPP

#include "../src/Template/Experiments/IOHprofiler_experimenter.hpp"
#include <vector>
#include <set>
#include <iostream>

extern IOHprofiler_random random_generator;

enum optimizationType{
  MAXIMIZATION,
  MINIMIZATION
};

optimizationType Opt = MAXIMIZATION;

void initialization(std::vector<int> &individual, int dimension, int range);
void eaInitialiazation(std::vector<std::vector<int> > p, const int mu, const int dimension, const int range);
int proportionalSelection(std::vector<double> fitness);
void selection(std::vector<std::vector<int> > &selected_population, const std::vector<std::vector<int> > &population, const std::vector<double>& fitness, const double select_rate);
void selection(std::vector<std::vector<int> > &selected_population, const std::vector<std::vector<int> > &population,const std::vector<double> &fitness, const int select_size);
void eaSelection_plus(std::vector<std::vector<int> > &parents, std::vector<double> &parents_fitness, const std::vector<std::vector<int> > &offsprings, const std::vector<double> &offsprings_fitness);
void eaSelection_comma(std::vector<std::vector<int> > &parents, std::vector<double> &parents_fitness, const std::vector<std::vector<int> > &offsprings,const std::vector<double> &offsprings_fitness);
void bitMutation(std::vector<int> &x, const std::vector<double> mutation_rate);
void bitMutation(std::vector<int> &x, const double mutation_rate);
void bitMutation(std::vector<int> &x, const int mutation_length);
bool crossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2, const double pc);
bool twopointsCrossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2);
bool uniformCrossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2);
void recombination (std::vector<std::vector<int> > &new_individuals, const std::vector<std::vector<int> > &parents, const std::vector<double> &fitness, const double &pc);
void sampleNFromM(std::vector<int> &sampled_number, int n, int m);
std::vector<int> sampling(const std::vector<std::vector<double> > &distribution);
int sampling(const std::vector<double> &distribution);
int sampleLRand(const std::vector<double> &distribution);
int sampleLRand(const double mutation_rate, const int dimension);
int normalRand(int mu, double sigma, int upperbound);
int normalRand_int(double mu, double sigma);
double normalRand(double mu, double sigma);
int logNormal_random_int(double mu, double sigma);
double sigmod(const double v);
void calculateDistribution(std::vector<std::vector<double> > &distribution, const std::vector<std::vector<int> > &selected_population,const std::vector<int> &range);
void binomial(const int N, const double p, std::vector<double> &result);
void powerLawDistribution(std::vector<double> &p, int N);


/***********************************************************************************************
 * 
 * INITIALIZATION SESSION
 * 
 * *********************************************************************************************/


/// \fn initialization(std::vector<int> &individual, int dimension, int range)
///
/// It creates a vector of uniformly random integer value in [0,range[
/// \param individual & the generated vector
/// \param dimension
/// \param range
void initialization(std::vector<int> &individual, int dimension, int range) {
  if (individual.size() != 0) {
    individual.clear();
  }

  individual.reserve(dimension);
  for (int i = 0; i != dimension; ++i) {
    individual.push_back( (int) ( random_generator.IOHprofiler_uniform_rand() * range ) );
  }
};

/// \fn void eaInitialiazation(std::vector<std::vector<int> > p, const int mu, const int dimension, const int range)
/// \brief Initialization of evolutionary algorithms.
///
/// It creates a vector of individuals, the size of the vector is mu, dimension is the length of 
/// individuals, and range presents the range of valuse of individuals' bits.
/// \param p & the generated population
/// \param mu the size of population
/// \param dimension
/// \param range the range of bit values are [0,range[
void eaInitialiazation(std::vector<std::vector<int> > p, const int mu, const int dimension, const int range) {
  if (p.size() != 0) {
    p.clear();
  }

  p = std::vector<std::vector<int> >(mu);
  for (int i = 0; i != mu; ++i) {
    initialization(p[i], dimension,range);
  }
}

/***********************************************************************************************
 * 
 * SELECTION SESSION
 * 
 * *********************************************************************************************/
/// \fn proportionalSelection(std::vector<double> fitness)
/// \brief Proportional selection based on a vector of fitness.
///
/// It returns the index of the selected fitness
/// \param fitness a vector of fitness values
int proportionalSelection(std::vector<double> fitness) {
  int n = fitness.size();
  std::vector<double> cumsum_fitness(n);
  
  double min_fitness = DBL_MAX;
  double max_fitness = -DBL_MAX;
  for (int i = 0; i != n; ++i) {
    if (fitness[i] < min_fitness) {
      min_fitness = fitness[i];
    }
    if (fitness[i] > max_fitness) {
      max_fitness = fitness[i];
    }
  }
  if (max_fitness == min_fitness) {
    return (int)floor(random_generator.IOHprofiler_uniform_rand() * n);
  }

  if (min_fitness < 0) {
    for (int i = 0; i != n; ++i) {
      fitness[i] += (-min_fitness);
    }
  }

  cumsum_fitness[0] = fitness[0];
  
  for (int i = 1; i != n; ++i) {
    cumsum_fitness[i] = fitness[i] + cumsum_fitness[i-1];
  }

  if(Opt == MINIMIZATION) {  
    std::vector<double> transfered_fitness(n);
    for (int i = 0; i != n; ++i) {
      transfered_fitness[i] = cumsum_fitness[n-1] - fitness[i] / cumsum_fitness[n-1];
    }
    cumsum_fitness[0] = 0;
    for (int i = 0; i != n; ++i) {
      cumsum_fitness[i] = transfered_fitness[i] + cumsum_fitness[i-1];
    }
  }

  double r = random_generator.IOHprofiler_uniform_rand() * cumsum_fitness[n-1];
  int index = 0;
  while (cumsum_fitness[index] < r && index < n - 1) {
    index++;
  }

  return index;
};


/// \fn void selection(std::vector<std::vector<int> > &selected_population, const std::vector<std::vector<int> > &population, const std::vector<double>& fitness, const double select_rate)
///
/// It returns a vector of selected individuals based on the ranking of fitness values.
/// \param selected_population
/// \param population
/// \param fitness
/// \param select_rate
void selection(std::vector<std::vector<int> > &selected_population, const std::vector<std::vector<int> > &population, const std::vector<double>& fitness, const double select_rate) {
  int select_size = (int)floor(population.size() * select_rate);
  
  if (selected_population.size() != 0) {
    selected_population.clear();
  }
  
  std::vector<int> index( fitness.size() );
  for (int i = 0; i != fitness.size(); ++i) {
    index[i] = i;
  }

  if(Opt == MAXIMIZATION) {
    std::partial_sort( index.begin(), index.begin() + select_size, index.end(), 
      [&](const int a, const int b) {return fitness[a] >= fitness[b];});
  } else if(Opt == MINIMIZATION) {
    std::partial_sort( index.begin(), index.begin() + select_size, index.end(), 
      [&](const int a, const int b) {return fitness[a] < fitness[b];});
  }

  selected_population = std::vector<std::vector<int> >(select_size);  
  for (int i = 0; i != select_size; ++i) {
    selected_population[i] = population[index[i]];
  }
};

/// \fn void selection(std::vector<std::vector<int> > &selected_population, std::vector<std::vector<int> > population, std::vector<double> fitness, const int select_size)
///
/// It creates a vector of selected individuals based on the ranking of fitness values. The best select_size individuals are selected.
/// \param selected_population
/// \param population
/// \param fitness
/// \param select_size
void selection(std::vector<std::vector<int> > &selected_population, const std::vector<std::vector<int> > &population,const std::vector<double> &fitness, const int select_size) {
  if (selected_population.size() != 0) {
    selected_population.clear();
  }
  
  std::vector<int> index( fitness.size() );
  for (int i = 0; i != fitness.size(); ++i) {
    index[i] = i;
  }

  if(Opt == MAXIMIZATION) {
    std::partial_sort( index.begin(), index.begin() + select_size, index.end(), 
      [&](const int a, const int b) {return fitness[a] >= fitness[b];});
  } else if(Opt == MINIMIZATION) {
    std::partial_sort( index.begin(), index.begin() + select_size, index.end(), 
      [&](const int a, const int b) {return fitness[a] < fitness[b];});
  }

  selected_population = std::vector<std::vector<int> >(select_size);  
  for (int i = 0; i != select_size; ++i) {
    selected_population[i] = population[index[i]];
  }
};

/// \fn void eaSelection(std::vector<std::vector<int> > &parents, std::vector<double> &parents_fitness, std::vector<std::vector<int> > offsprings, std::vector<double> offsprings_fitness)
/// \brief Selection operator.
///
/// It updateds paretns vector and parents fitness vector by a vector of selected individuals and
/// corresponding fitness values, using mu+lambda strategy.
/// \param parents
/// \param parents_fitness
/// \param offspring
/// \param offspring_fitness 
void eaSelection_plus(std::vector<std::vector<int> > &parents, std::vector<double> &parents_fitness, const std::vector<std::vector<int> > &offsprings, const std::vector<double> &offsprings_fitness) {
  int mu = parents_fitness.size();
  int lambda = offsprings_fitness.size();
  int n = parents[0].size();
  std::vector<int> index(mu+lambda);
  std::vector<std::vector<int> > population(mu+lambda, std::vector<int> (n));
  std::vector<double> fitness(mu+lambda);

  for (int i = 0; i != mu; ++i) {
    index[i] = i;
    population[i] = parents[i];
    fitness[i] = parents_fitness[i];
  }
  for (int i = 0; i != lambda; ++i) {
    index[i + mu] = i + mu;
    population[i + mu] = offsprings[i];
    fitness[i + mu] = offsprings_fitness[i];
  }

  if(Opt == MAXIMIZATION) {
    std::partial_sort(index.begin(), index.begin() + mu, index.end(), 
      [&](const int a, const int b) {return fitness[a] >= fitness[b];});
  } else if(Opt == MINIMIZATION) {
    std::partial_sort(index.begin(), index.begin() + mu, index.end(), 
      [&](const int a, const int b) {return fitness[a] < fitness[b];});
  }

  for (int i = 0; i != mu; ++i) {
    parents[i] = population[index[i]];
    parents_fitness[i] = fitness[index[i]];
  }
};

/// \fn void eaSelection(std::vector<std::vector<int> > &parents, std::vector<double> &parents_fitness, std::vector<std::vector<int> > offsprings, std::vector<double> offsprings_fitness)
/// \brief Selection operator.
///
/// It updateds paretns vector and parents fitness vector by a vector of selected individuals and
/// corresponding fitness values, using (mu,lambda) strategy.
/// \param parents
/// \param parents_fitness
/// \param offspring
/// \param offspring_fitness 
void eaSelection_comma(std::vector<std::vector<int> > &parents, std::vector<double> &parents_fitness, const std::vector<std::vector<int> > &offsprings,const std::vector<double> &offsprings_fitness) {
  int mu = parents_fitness.size();
  int lambda = offsprings_fitness.size();
  if (mu > lambda) {
    std::cout << "the size of offsprings is smaller than the size of parents" << std::endl;
    exit(1);
  }
  int n = offsprings[0].size();
  std::vector<int> index(lambda);

  for (int i = 0; i != lambda; ++i) {
    index[i] = i;
  }

  if(Opt == MAXIMIZATION) {
    std::partial_sort(index.begin(), index.begin() + mu, index.end(), 
      [&](const int a, const int b) {return offsprings_fitness[a] >= offsprings_fitness[b];});
  } else if(Opt == MINIMIZATION) {
    std::partial_sort(index.begin(), index.begin() + mu, index.end(), 
      [&](const int a, const int b) {return offsprings_fitness[a] < offsprings_fitness[b];});
  }

  for (int i = 0; i != mu; ++i) {
    parents[i] = offsprings[index[i]];
    parents_fitness[i] = offsprings_fitness[index[i]];
  }
};


/***********************************************************************************************
 * 
 * MUTATION SESSION
 * 
 * *********************************************************************************************/
/// void bitMutation(std::vector<int> &x, const std::vector<double> mutation_rate)
/// \brief bitMutation
///
/// It creates a new individual by flipping random bits of x, each bit of x has its own mutation rate.
/// \param y &
/// \param x
/// \param mutation_length
void bitMutation(std::vector<int> &x, const std::vector<double> mutation_rate) {
  if ( x.size() != mutation_rate.size() ) {
    std::cout << "bit_mutation, mutation rate size does not match" << std::endl;
    exit(1);
  }

  for (int i = 0; i != x.size(); ++i) {
    if ( random_generator.IOHprofiler_uniform_rand() < mutation_rate[i] ) {
      x[i] = (x[i] + 1) % 2;
    }
  }
};


/// void bitMutation(std::vector<int> &x, const double mutation_rate) 
/// \brief standard bit mutation
///
/// It returns a new individual after mutation.
/// The mutation operator works as flipping each bit with a probability = mutation rate.
/// \param x
/// \param mutation_rate
void bitMutation(std::vector<int> &x, const double mutation_rate) {
  for (int i = 0; i != x.size(); ++i) {
    if (random_generator.IOHprofiler_uniform_rand() < mutation_rate) {
      x[i] = (x[i] + 1) % 2;
    }
  }
};

/// std::vector<int> bitMutation(std::vector<int> x, const int mutation_length)
/// \brief standard bit mutation
///
/// It returns a new individual after mutation.
/// The mutation operator works as flipping l randomly selecteed bits, where l is mutation length.
/// \param x
/// \param mutation_length
void bitMutation(std::vector<int> &x, const int mutation_length) {
  // std::vector<int> l = sampleNFromM(mutation_length,x.size());
  
  std::vector<int> l;
  int flag, temp;
  for(int i = 0; i < mutation_length; ++i) {
    while(1) {
      flag = 0;
      temp = (int)(random_generator.IOHprofiler_uniform_rand() * x.size());
      for (int j = 0; j < i; ++j) {
        if (temp == l[j]) {
          flag = 1;
          break;
        }
      }
      if (flag == 0) {
        break;
      }
    }
    l.push_back(temp);
  }

  if (l.size() != mutation_length) {
    std::cout << "bit_mutation, l and mutation length do not match" << l.size() << "&" << mutation_length<<std::endl;
    exit(1);
  }
  
  for (int i = 0; i != l.size(); ++i) {
    x[l[i]] = (x[l[i]] + 1) % 2;
  }
};

/***********************************************************************************************
 * 
 * CROSSOVER SESSION
 * 
 * *********************************************************************************************/

/// \fn bool crossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2, const double &pc)
/// \brief Crossover operator.
///
/// It creates a new individual with uniform crossover.
/// \param ind1 one individual
/// \param ind2 one individual
/// \param pc crossover rate
bool crossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2, const double pc) {
  if (ind1.size() != ind2.size()) {
    std::cout << "the size of two individuals are different" << std::endl;
    exit(1);
  }

  int n = ind1.size();
  bool flag = false;
  new_individual = ind1;


  int l = 0;
  
  for (int i = 0; i != n; ++i) {
    if (random_generator.IOHprofiler_uniform_rand() < pc) {
      l++;
    }
  }
  
  if(l == 0 || l == n) {
    return false;
  }
  std::vector<int> flip;
  sampleNFromM(flip,l,n);
  for (int i = 0; i != l; ++i) {
    if (new_individual[flip[i]] != ind2[flip[i]]) {
      flag = true;
    }
    new_individual[flip[i]] = ind2[flip[i]];
  }

  return flag;
}


/// \fn bool twopointsCrossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2)
/// \brief Crossover operator.
///
/// It creates a new individual with two-bits crossover.
/// \param ind1 one individual
/// \param ind2 one individual
bool twopointsCrossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2) {
  if (ind1.size() != ind2.size()) {
    std::cout << "the size of two individuals are different" << std::endl;
    exit(1);
  }

  bool is_changed = false;
  int n = ind1.size();
  new_individual = ind1;
  
  std::vector<int> twopoints;
  sampleNFromM(twopoints,2,n);
  if (twopoints[0] > twopoints[1]) {
    double temp = twopoints[1];
    twopoints[1] = twopoints[0];
    twopoints[0] = temp;
  }
  for (int i = 0; i != n; ++i) {
    if (i >= twopoints[0]  && i < twopoints[1]) {
      new_individual[i] = ind2[i];
      if (ind1[i] != ind2[i]) {
        is_changed = true;
      }
    }
  }

  return is_changed;
}

/// \fn bool twopointsCrossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2)
/// \brief Crossover operator.
///
/// It creates a new individual with one-bit crossover.
/// \param ind1 one individual
/// \param ind2 one individual
bool onepointsCrossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2) {
  if (ind1.size() != ind2.size()) {
    std::cout << "the size of two individuals are different" << std::endl;
    exit(1);
  }

  bool is_changed = false;
  int n = ind1.size();
  new_individual = ind1;
  
  int point = (int)(random_generator.IOHprofiler_uniform_rand() * n);
  for (int i = 0; i != n; ++i) {
    if (i >= point) {
      new_individual[i] = ind2[i];
      if (ind1[i] != ind2[i]) {
        is_changed = true;
      }
    }
  }

  return is_changed;
}

/// \fn bool uniformCrossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2)
/// \brief Uniform crossover operator.
///
/// It creates a new individual with uniform crossover.
/// \param ind1 one individual
/// \param ind2 one individual
bool uniformCrossover(std::vector<int> &new_individual, const std::vector<int> &ind1, const std::vector<int> &ind2) {
  if (ind1.size() != ind2.size()) {
    std::cout << "the size of two individuals are different" << std::endl;
    exit(1);
  }

  int is_changed = 0;
  int n = ind1.size();
  new_individual = ind1;
  
  for (int i = 0; i != n; ++i) {
    if (random_generator.IOHprofiler_uniform_rand() < 0.5) {
      new_individual[i] = ind2[i];
      if (ind2[i] != ind1[i]) {
        is_changed = true;
      }
    }
  }
  return is_changed;
}


/// \fn std::vector<std::vector<int> > recombination (std::vector<std::vector<int> > parents, std::vector<double> fitness, double pc)
/// \brief Recombination operator.
///
/// It returns a vector of new individuals which size is two.
/// \param parents
/// \param fitness a vector of fitness values corresponding to parents
/// \param pc recombination rate
void recombination (std::vector<std::vector<int> > &new_individuals, const std::vector<std::vector<int> > &parents, const std::vector<double> &fitness, const double &pc) {
  if (parents.size() < 2) {
    std::cout << "the size of parents is smaller than 2" << std::endl;
    exit(1);
  }

  if (new_individuals.size() != 0) {
    new_individuals.clear();
  }

  int n = parents[0].size();
  int idx1 = proportionalSelection(fitness);
  int idx2 = proportionalSelection(fitness);
  

  while(idx1 == idx2) {
    idx2 = proportionalSelection(fitness);
  }

  new_individuals.push_back(parents[idx1]);
  new_individuals.push_back(parents[idx2]);
  
  if(random_generator.IOHprofiler_uniform_rand() < pc) {
    int cross_point = random_generator.IOHprofiler_uniform_rand() * (n - 1) + 1;
    for (int i = 0; i != n; ++i) {
      if (i  < cross_point) {
        new_individuals[0][i] = parents[idx2][i];
      } else {
        new_individuals[1][i] = parents[idx1][i];
      }
    }
  }
}

/***********************************************************************************************
 * 
 * SAMPLING METHODS SESSION
 * 
 * *********************************************************************************************/

/// \fn void sampleNFromM(std::vector<int> &sampled_number, int n, int m)
///
/// It creates a vector of n different ramdom values selected from [0,m[.
void sampleNFromM(std::vector<int> &sampled_number, int n, int m) {
  if (sampled_number.size() != 0) {
    sampled_number.clear();
  }

  if (n == 0) {
    std::cout <<  "sampled zero number" << std::endl;
  }

  std::vector<int> population;
  population.reserve(m);
  for (int i = 0; i < m; ++i) {
    population.push_back(i);
  }

  int temp,randPos;
  sampled_number.reserve(n);
  for (int i = m-1; i > 0; --i) {
    randPos = (int)floor(random_generator.IOHprofiler_uniform_rand() * (i+1));
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
};

/// \fn std::vector<int> sampling(std::vector<std::vector<double> > distribution)
/// It returns a vector of integer values sampling from the input distribution
///
/// \param distribution a vector of probablity desity function of each bit, elements
/// of the vector are probability of the corresponding variable (value of the index).
std::vector<int> sampling(const std::vector<std::vector<double> > &distribution) {
  std::vector<int> individual;
  for (int i = 0; i < distribution.size(); ++i) {
    double p = 0;
    double r = random_generator.IOHprofiler_uniform_rand();
    int j = 0;
    p = 0;
    while (j < distribution[i].size()) {
      if(r > p && r < p + distribution[i][j]) {
        break;
      }
      p +=  distribution[i][j];
      ++j;
    }
    individual.push_back(j);
  }
  return individual;
};

/// \fn int sampling(std::vector<double> distribution)
/// It returns an integer value sampling from the input distribution
///
/// \param distribution a probablity desity function, elements
/// of the vector are probability of the corresponding variable (value of the index).
int sampling(const std::vector<double> &distribution) {
  double p = 0.0;
  double r = random_generator.IOHprofiler_uniform_rand();
  int j = 0;
  while (j < distribution.size()) {
    if(r > p && r < p + distribution[j]) {
      break;
    }
    p +=  distribution[j];
    ++j;
  }
  return j;
}

/// \fn int sampleRand(std::vector<double> distribution)
/// \brief Sampling a number from a distribution.
///
/// It returns a random number from a given distribution. In algorithms, it is used to sample the number of
/// bits to be flipped in mutation or to be swapped in crossover.
/// \param distribution a vector of double values, which presents probability density.
int sampleLRand(const std::vector<double> &distribution) {
  int l = 0;
  while (l == 0) {
    double r = random_generator.IOHprofiler_uniform_rand();
    l = 0;
    double p = 0;
    while (l < distribution.size()) {
      if(r > p && r < p + distribution[l]) {
        break;
      }
      p +=  distribution[l];
      ++l;
    }
  }
  return l;
}


/// \fn int sampleLRand(const double mutation_rate)
/// \brief Sampling a number from a distribution.
///
/// It returns a mutation length based on mutation rate and dimension.s
int sampleLRand(const double mutation_rate, const int dimension) {
  int l = 0;
  while (l == 0) {
    for (int i = 0; i != dimension; ++i) {
      if (random_generator.IOHprofiler_uniform_rand() < mutation_rate) {
        l++;
      }
    }
  }
  return l;
}

/***********************************************************************************************
 * 
 * RANDOM METHODS SESSION
 * 
 * *********************************************************************************************/

/// \fn int normalRand(int mu, double sigma, int upperbound)
/// \brief Samping a positive integer random number ]0,upperbound]Wfrom N(mu,sigma) that is smaller than upperbound.
///
/// \param mu mean
/// \param sigma variance
int normalRand(int mu, double sigma, int upperbound) {
  int l = 0;
  while(l <= 0 || l > upperbound) {
    l = mu + (int)(random_generator.IOHprofiler_normal_rand() * sigma);
  }
  return l;
}

/// \fn int normalRand_int(double mu, double sigma)
/// \brief Samping a positive integer random number from N(mu,sigma) that is smaller than upperbound.
///
/// \param mu mean
/// \param sigma variance
int normalRand_int(double mu, double sigma) {
  double randN = random_generator.IOHprofiler_normal_rand()* sigma + mu;
  return (int)(randN + 0.5);

}


/// \fn int normalRand(int mu, double sigma)
/// \brief Samping a integer random number from N(mu,sigma).
///
/// \param mu mean
/// \param sigma variance
double normalRand(double mu, double sigma) {
  double result = mu + (random_generator.IOHprofiler_normal_rand() * sigma);
  return result;
}

/// \fn int logNormal_random_int(double mu, double sigma, int upperbound)
/// \brief Samping a positive integer random number from N(mu,sigma) that is smaller than upperbound.
///
/// \param mu mean
/// \param sigma variance
int logNormal_random_int(double mu, double sigma) {
  double randN = random_generator.IOHprofiler_normal_rand()* sigma + mu;
  return (int)(exp(randN) + 0.5);
}

/***********************************************************************************************
 * 
 * OTHER METHODS SESSION
 * 
 * *********************************************************************************************/


/// \fn double sigmod(const double v)
/// It returns 1.0 / (1 + exp(-v))
double sigmod(const double v) {
  return (1.0 / (1 + exp(-v)));
};

/// \fn std::vector<std::vector<double> > calculateDistribution(std::vector<std::vector<int> > selected_population, std::vector<int> range)
///
/// It returns a vector of double values (distribution)
/// \param selection_population
/// \param range defines the range of variables as [0,range-1]
void calculateDistribution(std::vector<std::vector<double> > &distribution, const std::vector<std::vector<int> > &selected_population,const std::vector<int> &range) {
  if (distribution.size() != 0) {
    distribution.clear();
  }

  int n = selected_population[0].size();
  distribution.reserve(n);
  for (int i = 0; i != n; ++i) {
    std::vector<int> count(range[i],0);
    std::vector<double> p(range[i],0.0);
    for (int j = 0; j != selected_population.size(); ++j) {    
      count[selected_population[j][i]]++;
    }
    for (int j = 0; j != range[i]; ++j) {
      p[j] = (double)count[j] / (double)(selected_population.size());

      /**
       * restriction for p
       */
      p[j] = p[j] > (1.0 - 1.0 / (double)n) ? (1.0 - 1.0 / (double)n) : p[j];
      p[j] = p[j] < 1.0 / (double)n ? 1.0/ (double)n : p[j];
    }
    distribution.push_back(p);
  }
};

/// \fn std::vector<double> binomial(const int N, const double p) 
/// \brief a binomial distribution Bin(N,p).
///
/// It returns a vector of probablity that Pr(x = n) that n in [0,N] based on binomial disritbuion Bin(N,p)
void binomial(const int N, const double p, std::vector<double> &result) {
  if (result.size() != 0) {
    result.clear();
  }
  double **binomalP = new double*[N + 1];
  for (int i = 0; i != N + 1; ++i) {
    binomalP[i] = new double[N + 1];
  }
    
  binomalP[0][0] = 1.0;

  for (int i = 1; i != N + 1; ++i) {
    binomalP[i][0] = (1.0 - p) * binomalP[i - 1][0];
  }
  
  for (int j = 1; j != N + 1; ++j) {
    binomalP[0][j] = 0.0;
  }
    
  for (int i = 1; i != N + 1; ++i) {
    for (int j = 1; j != N + 1; ++j) {
      binomalP[i][j] = (1.0 - p) * binomalP[i - 1][j] + p * binomalP[i - 1][j - 1];
    }
  }
  
  for (int j = 0; j != N+1; ++j) {
    result.push_back(binomalP[N][j]);
  }

  for (int i = 0; i != N + 1; ++i) {
    delete [] binomalP[i];
  }
  delete [] binomalP;
}

/// \fn void powerLawDistribution(std::vector<double> &p, int N)
/// \brief a power law distribution.
void powerLawDistribution(std::vector<double> &p, int N) {
  p = std::vector<double>(N+1);
  double C,beta;
  size_t i;

  beta = 1.5;
  C = 0.0;
  for(i = 0; i < N; ++i){
    C += pow(i+1,-beta);
  }
  for(i = 0; i < N; ++i){
    p[i+1] = 1/C * pow(i+1,-beta);
  }

  p[0] = 0.0;
}

#endif