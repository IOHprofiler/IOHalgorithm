/// \file w_model_OneMax_suite.h
/// \brief Hpp file for class w_model_OneMax_suite.
///
/// A suite for test w-model extensions of OneMax.
///
/// \author Furong Ye
/// \date 2019-12-05

#ifndef _W_MODEL_ONEMAX_SUITE_HPP
#define _W_MODEL_ONEMAX_SUITE_HPP

#include <sstream>
#include <iomanip>
#include "f_one_max.hpp"
#include "f_w_model_one_max.hpp"
#include "IOHprofiler_suite.h"

using namespace std;

vector<double> default_om_dummy = {0.0, 0.1, 0.25, 0.9};
vector<int> default_om_epistasis = {0, 2, 3, 4, 5, 7};
vector<int> default_om_neutrality = {1, 3, 5};
vector<double> default_om_ruggedness = {0, 0.2, 0.4, 0.6, 0.8, 1};

class W_Model_OneMax_suite : public IOHprofiler_suite<int>
{
public:
  W_Model_OneMax_suite(vector<double> dummy_para = default_om_dummy,
                       vector<int> epistasis_para = default_om_epistasis,
                       vector<int> neturality_para = default_om_neutrality,
                       vector<double> ruggedness_para = default_om_ruggedness)
      : dummy_para_(dummy_para),
        epistasis_para_(epistasis_para),
        neturality_para_(neturality_para),
        ruggedness_para_(ruggedness_para)
  {

    for (size_t du = 0; du != this->dummy_para_.size(); ++du)
    {
      for (size_t ep = 0; ep != this->epistasis_para_.size(); ++ep)
      {
        for (size_t ne = 0; ne != this->neturality_para_.size(); ++ne)
        {
          for (size_t ru = 0; ru != this->ruggedness_para_.size(); ++ru)
          {
            this->para_product_.push_back({du, ep, ne, ru});
          }
        }
      }
    }

    vector<int> problem_id;
    vector<int> instance_id = {1};
    vector<int> dimension = {20, 30};
    for (int i = 0; i < para_product_.size(); ++i)
    {
      problem_id.push_back(i + 1);
    }

    IOHprofiler_set_suite_problem_id(problem_id);
    IOHprofiler_set_suite_instance_id(instance_id);
    IOHprofiler_set_suite_dimension(dimension);
    IOHprofiler_set_suite_name("W_Model_OneMax_suite");
    this->loadProblem();
  }

  W_Model_OneMax_suite(std::vector<int> problem_id, std::vector<int> instance_id, std::vector<int> dimension,
                       vector<double> dummy_para = default_om_dummy,
                       vector<int> epistasis_para = default_om_epistasis,
                       vector<int> neturality_para = default_om_neutrality,
                       vector<double> ruggedness_para = default_om_ruggedness)
      : dummy_para_(dummy_para),
        epistasis_para_(epistasis_para),
        neturality_para_(neturality_para),
        ruggedness_para_(ruggedness_para)
  {

    for (size_t du = 0; du != this->dummy_para_.size(); ++du)
    {
      for (size_t ep = 0; ep != this->epistasis_para_.size(); ++ep)
      {
        for (size_t ne = 0; ne != this->neturality_para_.size(); ++ne)
        {
          for (size_t ru = 0; ru != this->ruggedness_para_.size(); ++ru)
          {
            this->para_product_.push_back({du, ep, ne, ru});
          }
        }
      }
    }

    for (size_t i = 0; i < problem_id.size(); ++i)
    {
      if (problem_id[i] < 0 || problem_id[i] > this->para_product_.size())
      {
        IOH_error("problem_id " + std::to_string(problem_id[i]) + " is not in W_Model_OneMax_suite");
      }
    }

    for (size_t i = 0; i < instance_id.size(); ++i)
    {
      if (instance_id[i] < 0 || instance_id[i] > 100)
      {
        IOH_error("instance_id " + std::to_string(instance_id[i]) + " is not in W_Model_OneMax_suite");
      }
    }

    for (size_t i = 0; i < dimension.size(); ++i)
    {
      if (dimension[i] < 0 || dimension[i] > 20000)
      {
        IOH_error("dimension " + std::to_string(dimension[i]) + " is not in W_Model_OneMax_suite");
      }
    }

    IOHprofiler_set_suite_problem_id(problem_id);
    IOHprofiler_set_suite_instance_id(instance_id);
    IOHprofiler_set_suite_dimension(dimension);
    IOHprofiler_set_suite_name("W_Model_OneMax_suite");
    this->loadProblem();
  }

  void loadProblem()
  {

    if (this->size() != 0)
    {
      this->clear();
    }
    this->IOHprofiler_set_size_of_problem_list(this->IOHprofiler_suite_get_number_of_problems() *
                                               this->IOHprofiler_suite_get_number_of_instances() *
                                               this->IOHprofiler_suite_get_number_of_dimensions());

    vector<int> p_id = this->IOHprofiler_suite_get_problem_id();
    vector<int> i_id = this->IOHprofiler_suite_get_instance_id();
    vector<int> d = this->IOHprofiler_suite_get_dimension();
    for (int i = 0; i != p_id.size(); ++i)
    {
      for (int j = 0; j != d.size(); ++j)
      {
        for (int h = 0; h != i_id.size(); ++h)
        {
          shared_ptr<W_Model_OneMax> p(new W_Model_OneMax());
          p->set_w_setting(this->dummy_para_[this->para_product_[p_id[i] - 1][0]],
                           this->epistasis_para_[this->para_product_[p_id[i] - 1][1]],
                           this->neturality_para_[this->para_product_[p_id[i] - 1][2]],
                           static_cast<int>(floor(this->IOHprofiler_suite_get_dimension()[j] * this->ruggedness_para_[this->para_product_[p_id[i] - 1][3]])));

          // Set the problem name.
          string problem_name = "Onemax";
          std::stringstream dss;
          dss << std::setprecision(3) << this->dummy_para_[this->para_product_[p_id[i] - 1][0]];
          problem_name += "_D" + dss.str();
          problem_name += "_E" + std::to_string(this->epistasis_para_[this->para_product_[p_id[i] - 1][1]]);
          problem_name += "_N" + std::to_string(this->neturality_para_[this->para_product_[p_id[i] - 1][2]]);
          std::stringstream rss;
          rss << std::setprecision(3) << this->ruggedness_para_[this->para_product_[p_id[i] - 1][3]];
          problem_name += "_R" + rss.str();
          p->IOHprofiler_set_problem_name(problem_name);
          p->IOHprofiler_set_problem_id(p_id[i]);
          p->IOHprofiler_set_number_of_variables(d[j]);
          p->IOHprofiler_set_instance_id(i_id[h]);
          this->push_back(p);
          mapIDTOName(p_id[i], problem_name);
        }
      }
    }

    assert(this->IOHprofiler_get_size_of_problem_list() == this->size());
    this->IOHprofiler_set_get_problem_flag(false);
    this->IOHprofiler_set_problem_list_index(0);
    this->IOHprofiler_set_load_problem_flag(true);
  }

  vector<double> get_dummy_para()
  {
    return this->dummy_para_;
  }

  vector<int> get_epistasis_para()
  {
    return this->epistasis_para_;
  }

  vector<int> get_neturality_para()
  {
    return this->neturality_para_;
  }

  vector<double> get_ruggedness_para()
  {
    return this->ruggedness_para_;
  }

  void set_dummy_para(const vector<double> &dummy_para)
  {
    this->dummy_para_ = dummy_para;
  }

  void set_epistasis_para(const vector<int> &epistasis_para)
  {
    this->epistasis_para_ = epistasis_para;
  }

  void set_neutrality_para(const vector<int> &neutrality_para)
  {
    this->neturality_para_ = neutrality_para;
  }

  void set_ruggedness_para(const vector<double> &ruggedness_para)
  {
    this->ruggedness_para_ = ruggedness_para;
  }

  static W_Model_OneMax_suite *createInstance()
  {
    return new W_Model_OneMax_suite();
  }

  static W_Model_OneMax_suite *createInstance(std::vector<int> problem_id, std::vector<int> instance_id, std::vector<int> dimension)
  {
    return new W_Model_OneMax_suite(problem_id, instance_id, dimension);
  }

private:
  vector<double> dummy_para_;
  vector<int> epistasis_para_;
  vector<int> neturality_para_;
  vector<double> ruggedness_para_;
  vector<vector<size_t>> para_product_;
};

#endif //_W_MODEL_ONEMAX_SUITE_HPP
