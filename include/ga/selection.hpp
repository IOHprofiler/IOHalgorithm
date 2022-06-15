/// \file selection.h
/// \brief Header file for class Selection.
///
/// It contains functions of selection operators.
///
/// \author Furong Ye
/// \date 2020-07-02

#ifndef _SELECTION_H_
#define _SELECTION_H_

#include "utils/common.hpp"

namespace modularGA
{


    /// Definition of selection operator id.
    enum selection_operator
    {
      BESTPLUS = 1,         /// < (mu+lambda) using elitist selection.
      BESTCOMMA = 2,        /// < (mu,lambda)  using elitist selection.
      TOURNAMENTPLUS = 3,   /// < (mu+lambda) using tournament selection.
      TOURNAMENTCOMMA = 4,  /// < (mu+lambda) using tournament selection.
      PROPORTIONALPLUS = 5, /// < (mu+lambda) using proportional selection.
      PROPORTIONALCOMMA = 6 /// < (mu+lambda) using proportional selection.
    };

    class Selection
    {
    public:
      Selection() {}
      ~Selection() {}
      Selection(const Selection &) = delete;
      Selection &operator=(const Selection &) = delete;
      void DoSelection(vector<vector<int>> &parents, vector<double> &parents_fitness, const vector<vector<int>> &offspring, const vector<double> &offspring_fitness)
      {
        switch (this->selection_operator_)
        {
        case 1:
          this->BestPlusStrategy(parents, parents_fitness, offspring, offspring_fitness);
          break;
        case 2:
          this->BestCommaStrategy(parents, parents_fitness, offspring, offspring_fitness);
          break;
        case 3:
          this->TournamentPlusStrategy(parents, parents_fitness, offspring, offspring_fitness);
          break;
        case 4:
          this->TournamentCommaStrategy(parents, parents_fitness, offspring, offspring_fitness);
          break;
        case 5:
          this->ProportionalPlusStrategy(parents, parents_fitness, offspring, offspring_fitness);
          break;
        case 6:
          this->ProportionalCommaStrategy(parents, parents_fitness, offspring, offspring_fitness);
          break;
        default:
          cerr << "unknown selection operator" << endl;
          assert(false);
          break;
        }
      }

      void BestCommaStrategy(vector<vector<int>> &parents, vector<double> &parents_fitness, const vector<vector<int>> &offspring, const vector<double> &offspring_fitness)
      {
        size_t mu = parents_fitness.size();
        size_t lambda = offspring_fitness.size();

        assert(mu <= lambda); /// < because it is best comma strategy
        assert(parents.size() == parents_fitness.size());
        assert(offspring.size() == offspring_fitness.size());

        vector<size_t> index(lambda);
        for (size_t i = 0; i != lambda; ++i)
        {
          index[i] = i;
        }

        if (Opt == optimizationType::MAXIMIZATION)
        {
          partial_sort(index.begin(), index.begin() + mu, index.end(),
                       [&](const size_t a, const size_t b)
                       { return offspring_fitness[a] >= offspring_fitness[b]; });
        }
        else if (Opt == optimizationType::MINIMIZATION)
        {
          partial_sort(index.begin(), index.begin() + mu, index.end(),
                       [&](const size_t a, const size_t b)
                       { return offspring_fitness[a] <= offspring_fitness[b]; });
        }

        for (size_t i = 0; i != mu; ++i)
        {
          parents[i] = offspring[index[i]];
          parents_fitness[i] = offspring_fitness[index[i]];
        }
      }

      void BestPlusStrategy(vector<vector<int>> &parents, vector<double> &parents_fitness, const vector<vector<int>> &offspring, const vector<double> &offspring_fitness)
      {
        assert(parents.size() == parents_fitness.size());
        assert(offspring.size() == offspring_fitness.size());

        size_t mu = parents_fitness.size();
        size_t lambda = offspring_fitness.size();

        vector<size_t> index_parents(mu);
        vector<size_t> index_offspring(lambda);
        for (size_t i = 0; i != mu; ++i)
        {
          index_parents[i] = i;
        }
        for (size_t i = 0; i != lambda; ++i)
        {
          index_offspring[i] = i;
        }

        if (Opt == optimizationType::MAXIMIZATION)
        {
          sort(index_parents.begin(), index_parents.end(),
               [&](const size_t a, const size_t b)
               { return parents_fitness[a] < parents_fitness[b]; });
          sort(index_offspring.begin(), index_offspring.end(),
               [&](const size_t a, const size_t b)
               { return offspring_fitness[a] > offspring_fitness[b]; });

          size_t j = mu - 1, i = 0, replace_flag = 0, pick_count = 0;
          while (pick_count < mu && i < lambda)
          {
            if (offspring_fitness[index_offspring[i]] >= parents_fitness[index_parents[j]])
            {
              replace_flag++;
              i++;
            }
            else
            {
              j--;
            }
            pick_count++;
          }
          for (size_t i = 0; i != replace_flag; ++i)
          {
            parents[index_parents[i]] = offspring[index_offspring[i]];
            parents_fitness[index_parents[i]] = offspring_fitness[index_offspring[i]];
          }
        }
        else if (Opt == optimizationType::MINIMIZATION)
        {
          std::sort(index_parents.begin(), index_parents.end(),
                    [&](const size_t a, const size_t b)
                    { return parents_fitness[a] > parents_fitness[b]; });
          std::sort(index_offspring.begin(), index_offspring.end(),
                    [&](const size_t a, const size_t b)
                    { return offspring_fitness[a] < offspring_fitness[b]; });

          size_t j = mu - 1, i = 0, replace_flag = 0, pick_count = 0;
          while (pick_count < mu && i < lambda)
          {
            if (offspring_fitness[index_offspring[i]] <= parents_fitness[index_parents[j]])
            {
              replace_flag++;
              i++;
            }
            else
            {
              j--;
            }
            pick_count++;
          }

          for (size_t i = 0; i != replace_flag; ++i)
          {
            parents[index_parents[i]] = offspring[index_offspring[i]];
            parents_fitness[index_parents[i]] = offspring_fitness[index_offspring[i]];
          }
        }
      }

      void TournamentCommaStrategy(vector<vector<int>> &parents, vector<double> &parents_fitness, const vector<vector<int>> &offspring, const vector<double> &offspring_fitness)
      {
        size_t mu = parents_fitness.size();
        size_t lambda = offspring_fitness.size();

        assert(this->tournament_k_ <= lambda);

        double tmp_best;
        vector<size_t> sample_k;
        vector<size_t> select_index(mu);
        for (size_t i = 0; i != mu; ++i)
        {
          sampleNFromM(sample_k, this->tournament_k_, lambda);
          if (Opt == optimizationType::MAXIMIZATION)
          {
            tmp_best = numeric_limits<double>::lowest();
            for (size_t j = 0; j != this->tournament_k_; ++j)
            {
              if (offspring_fitness[sample_k[j]] > tmp_best)
              {
                tmp_best = offspring_fitness[sample_k[j]];
                select_index[i] = sample_k[j];
              }
            }
          }
          else
          {
            tmp_best = numeric_limits<double>::max();
            for (size_t j = 0; j != this->tournament_k_; ++j)
            {
              if (offspring_fitness[sample_k[j]] < tmp_best)
              {
                tmp_best = offspring_fitness[sample_k[j]];
                select_index[i] = sample_k[j];
              }
            }
          }
        }

        for (size_t i = 0; i != mu; ++i)
        {
          parents[i] = offspring[select_index[i]];
          parents_fitness[i] = offspring_fitness[select_index[i]];
        }
      }

      void TournamentPlusStrategy(vector<vector<int>> &parents, vector<double> &parents_fitness, const vector<vector<int>> &offspring, const vector<double> &offspring_fitness)
      {
        vector<vector<int>> backup_parents = parents;
        vector<double> backup_parents_fitness = parents_fitness;
        size_t mu = parents_fitness.size();
        size_t lambda = offspring_fitness.size();

        assert(this->tournament_k_ <= lambda + mu);

        double tmp_best;
        vector<size_t> sample_k;
        vector<size_t> select_index(mu);
        for (int i = 0; i != mu; ++i)
        {
          sampleNFromM(sample_k, this->tournament_k_, mu + lambda);
          if (Opt == optimizationType::MAXIMIZATION)
          {
            tmp_best = numeric_limits<double>::lowest();
            for (size_t j = 0; j != this->tournament_k_; ++j)
            {
              if (sample_k[j] >= mu)
              {
                if (offspring_fitness[sample_k[j] - mu] > tmp_best)
                {
                  tmp_best = offspring_fitness[sample_k[j] - mu];
                  select_index[i] = sample_k[j];
                }
              }
              else
              {
                if (parents_fitness[sample_k[j]] > tmp_best)
                {
                  tmp_best = parents_fitness[sample_k[j]];
                  select_index[i] = sample_k[j];
                }
              }
            }
          }
          else
          {
            tmp_best = numeric_limits<double>::max();
            for (size_t j = 0; j != this->tournament_k_; ++j)
            {
              if (sample_k[j] >= mu)
              {
                if (offspring_fitness[sample_k[j] - mu] < tmp_best)
                {
                  tmp_best = offspring_fitness[sample_k[j] - mu];
                  select_index[i] = sample_k[j];
                }
              }
              else
              {
                if (parents_fitness[sample_k[j]] < tmp_best)
                {
                  tmp_best = parents_fitness[sample_k[j]];
                  select_index[i] = sample_k[j];
                }
              }
            }
          }
        }

        for (size_t i = 0; i != mu; ++i)
        {
          if (select_index[i] >= mu)
          {
            parents[i] = offspring[select_index[i] - mu];
            parents_fitness[i] = offspring_fitness[select_index[i] - mu];
          }
          else
          {
            parents[i] = backup_parents[select_index[i]];
            parents_fitness[i] = backup_parents_fitness[select_index[i]];
          }
        }
      }

      void ProportionalCommaStrategy(vector<vector<int>> &parents, vector<double> &parents_fitness, const vector<vector<int>> &offspring, const vector<double> &offspring_fitness)
      {
        size_t mu = parents_fitness.size();
        size_t lambda = offspring_fitness.size();

        vector<double> proportional_f(lambda);
        double fitness_sum = 0;

        for (size_t i = 0; i != lambda; ++i)
        {
          fitness_sum += offspring_fitness[i];
        }

        if (Opt == optimizationType::MAXIMIZATION)
        {
          for (size_t i = 0; i != lambda; ++i)
          {
            proportional_f[i] = offspring_fitness[i] / fitness_sum;
          }
        }
        else
        {
          double max_fitness = numeric_limits<double>::lowest();
          for (size_t i = 0; i != lambda; ++i)
          {
            if (offspring_fitness[i] > max_fitness)
            {
              max_fitness = offspring_fitness[i];
            }
          }
          fitness_sum = max_fitness * lambda - fitness_sum;
          for (size_t i = 0; i != lambda; ++i)
          {
            proportional_f[i] = (max_fitness - offspring_fitness[i]) / fitness_sum;
          }
        }

        for (size_t i = 1; i != lambda; ++i)
        {
          proportional_f[i] += proportional_f[i - 1];
        }

        double r;
        vector<size_t> select_index(mu);
        for (int i = 0; i != mu; ++i)
        {
          r = uniform_random();
          size_t j = 0;
          while (r > proportional_f[j] && r < mu)
          {
            ++j;
          }
          select_index[i] = j;
        }

        for (size_t i = 0; i != mu; ++i)
        {
          parents[i] = offspring[select_index[i]];
          parents_fitness[i] = offspring_fitness[select_index[i]];
        }
      }

      void ProportionalPlusStrategy(vector<vector<int>> &parents, vector<double> &parents_fitness, const vector<vector<int>> &offspring, const vector<double> &offspring_fitness)
      {
        vector<vector<int>> backup_parents = parents;
        vector<double> backup_parents_fitness = parents_fitness;
        size_t mu = parents_fitness.size();
        size_t lambda = offspring_fitness.size();

        vector<double> proportional_f(mu + lambda);
        double fitness_sum = 0;

        for (size_t i = 0; i != mu; ++i)
        {
          fitness_sum += parents_fitness[i];
        }
        for (size_t i = 0; i != lambda; ++i)
        {
          fitness_sum += offspring_fitness[i];
        }

        if (Opt == optimizationType::MAXIMIZATION)
        {
          for (size_t i = 0; i != mu; ++i)
          {
            proportional_f[i] = parents_fitness[i] / fitness_sum;
          }
          for (size_t i = mu; i != mu + lambda; ++i)
          {
            proportional_f[i] = offspring_fitness[i - mu] / fitness_sum;
          }
        }
        else
        {
          double max_fitness = numeric_limits<double>::lowest();
          for (size_t i = 0; i != mu; ++i)
          {
            if (parents_fitness[i] > max_fitness)
            {
              max_fitness = parents_fitness[i];
            }
          }
          for (size_t i = 0; i != lambda; ++i)
          {
            if (offspring_fitness[i] > max_fitness)
            {
              max_fitness = offspring_fitness[i];
            }
          }
          fitness_sum = max_fitness * (mu + lambda) - fitness_sum;
          for (size_t i = 0; i != mu; ++i)
          {
            proportional_f[i] = (max_fitness - parents_fitness[i]) / fitness_sum;
          }
          for (size_t i = mu; i != mu + lambda; ++i)
          {
            proportional_f[i] = (max_fitness - offspring_fitness[i - mu]) / fitness_sum;
          }
        }

        for (size_t i = 1; i != mu + lambda; ++i)
        {
          proportional_f[i] += proportional_f[i - 1];
        }

        double r;
        vector<size_t> select_index(mu);
        for (size_t i = 0; i != mu; ++i)
        {
          r = uniform_random();
          int j = 0;
          while (r > proportional_f[j] && j < mu + lambda)
          {
            ++j;
          }
          select_index[i] = j;
        }

        for (size_t i = 0; i != mu; ++i)
        {
          if (select_index[i] >= mu)
          {
            parents[i] = offspring[select_index[i] - mu];
            parents_fitness[i] = offspring_fitness[select_index[i] - mu];
          }
          else
          {
            parents[i] = backup_parents[select_index[i]];
            parents_fitness[i] = backup_parents_fitness[select_index[i]];
          }
        }
      }

      void set_selection_operator(const int s)
      {
        this->selection_operator_ = s;
      }

      void set_selection_operator(string s)
      {
        transform(s.begin(), s.end(), s.begin(), ::toupper);
        if (s == "BESTPLUS")
        {
          this->set_selection_operator(BESTPLUS);
        }
        else if (s == "BESTCOMMA")
        {
          this->set_selection_operator(BESTCOMMA);
        }
        else if (s == "TOURNAMENTPLUS")
        {
          this->set_selection_operator(TOURNAMENTPLUS);
        }
        else if (s == "TOURNAMENTCOMMA")
        {
          this->set_selection_operator(TOURNAMENTCOMMA);
        }
        else if (s == "PROPORTIONALPLUS")
        {
          this->set_selection_operator(PROPORTIONALPLUS);
        }
        else if (s == "PROPORTIONALCOMMA")
        {
          this->set_selection_operator(PROPORTIONALCOMMA);
        }
        else
        {
          cerr << "invalid name for selection operator" << endl;
          assert(false);
        }
      }

      void set_tournament_k(const int k)
      {
        this->tournament_k_ = k;
      }

      int get_selection_operator() const
      {
        return this->selection_operator_;
      }

      int get_tournament_k() const
      {
        return this->tournament_k_;
      }

    private:
      int selection_operator_;
      int tournament_k_;
    };
}
#endif // _SELECTION_H_
