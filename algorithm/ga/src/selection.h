/// \file selection.h
/// \brief Header file for class Selection.
///
/// It contains functions of selection operators.
/// 
/// \author Furong Ye
/// \date 2020-07-02

#ifndef _SELECTION_H_
#define _SELECTION_H_

#include "common.h"

/// Definition of selection operator id.
enum selection_operator {
  BESTPLUS = 1, /// < (mu+lambda) using elitist selection.
  BESTCOMMA = 2, /// < (mu,lambda)  using elitist selection.
  TOURNAMENTPLUS = 3, /// < (mu+lambda) using tournament selection.
  TOURNAMENTCOMMA = 4, /// < (mu+lambda) using tournament selection.
  PROPORTIONALPLUS = 5, /// < (mu+lambda) using proportional selection.
  PROPORTIONALCOMMA = 6 /// < (mu+lambda) using proportional selection.
};

class Selection {
public:
  Selection(){}
  ~Selection() {}
  Selection(const Selection&) = delete;
  Selection &operator = (const Selection&) = delete;
  
  void DoSelection(vector< vector<int> > &parents, vector <double> &parents_fitness, const vector< vector<int> > &offspring, const vector<double> &offspring_fitness);
  
  void BestCommaStrategy(vector< vector<int> > &parents, vector <double> &parents_fitness, const vector< vector<int> > &offspring, const vector<double> &offspring_fitness);
  void BestPlusStrategy(vector< vector<int> > &parents, vector <double> &parents_fitness, const vector< vector<int> > &offspring, const vector<double> &offspring_fitness);
  void TournamentCommaStrategy(vector< vector<int> > &parents, vector <double> &parents_fitness, const vector< vector<int> > &offspring, const vector<double> &offspring_fitness);
  void TournamentPlusStrategy(vector< vector<int> > &parents, vector <double> &parents_fitness, const vector< vector<int> > &offspring, const vector<double> &offspring_fitness);
  void ProportionalCommaStrategy(vector< vector<int> > &parents, vector <double> &parents_fitness, const vector< vector<int> > &offspring, const vector<double> &offspring_fitness);
  void ProportionalPlusStrategy(vector< vector<int> > &parents, vector <double> &parents_fitness, const vector< vector<int> > &offspring, const vector<double> &offspring_fitness);
  
  void set_selection_operator(const int s);
  void set_selection_operator(string s);
  void set_tournament_k(const int k);
  
  int get_selection_operator() const;
  int get_tournament_k() const;
  
private:
  int selection_operator_;
  int tournament_k_;
};

#endif // _SELECTION_H_
