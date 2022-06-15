#pragma once
/**
 * Here we implement the Simulated Annealing algorithm with
 * exponential temperature schedule.
 *
 * The problem with simulated annealing is that it has several
 * parameters that need to be configured well in order for it
 * to work.
 * The proper configuration of the parameters depends on the
 * scale of the objective function value as well as on the
 * available computation budget, and even stuff like, e.g.,
 * the expected differences in objective values of local optima.
 *
 * This makes it very complicated to properly use it for
 * black-box settings where absolutely nothing is known
 * about the objective function.
 *
 * Still, to be able to do some experiments, we propose the
 * following automatic setup based on the number 'n' of
 * variables of the problem:
 *
 * - the temperature schedule is exponential
 * - the start temperature is such that it time 1 the
 *   probability to accept a solution which is n/4 worse
 *   than the current solution would be 0.1
 * - the end temperature is such that the probability to
 *   accept a solution which is 1 worse than the current
 *   solution is 1/sqrt(n)
 * - the value of the epsilon parameter of the schedule
 *   is computed accordingly
 *
 * This will have several implications:
 *
 * - Problems that can be solved easily (OneMax,
 *   LeadingOnes) will be _still_ solved late in the runs.
 *   This is normal for simulated annealing, as the goal
 *   is to transcend from random walk to hill climbing
 *   behavior.
 * - Problems whose objective value range differs largely
 *   from n (e.g., LARS) will probably not be solved well.
 *
 * Author: Thomas Weise
 *         Institute of Applied Optimization
 *         Hefei University
 *         Hefei, Anhui, China
 * Email: tweise@hfuu.edu.cn, tweise@ustc.edu.cn
 */

#ifndef _SA_H_
#define _SA_H_
#include "utils/common.hpp"
#include <math.h>
namespace modularGA
{

    namespace sa
    {
        /**
         * Here we implement the Simulated Annealing algorithm with
         * exponential temperature schedule.
         *
         * The problem with simulated annealing is that it has several
         * parameters that need to be configured well in order for it
         * to work.
         * The proper configuration of the parameters depends on the
         * scale of the objective function value as well as on the
         * available computation budget, and even stuff like, e.g.,
         * the expected differences in objective values of local optima.
         *
         * This makes it very complicated to properly use it for
         * black-box settings where absolutely nothing is known
         * about the objective function.
         *
         * Still, to be able to do some experiments, we propose the
         * following automatic setup based on the number 'n' of
         * variables of the problem:
         *
         * - the temperature schedule is exponential
         * - the start temperature is such that it time 1 the
         *   probability to accept a solution which is n/4 worse
         *   than the current solution would be 0.1
         * - the end temperature is such that the probability to
         *   accept a solution which is 1 worse than the current
         *   solution is 1/sqrt(n)
         * - the value of the epsilon parameter of the schedule
         *   is computed accordingly
         *
         * This will have several implications:
         *
         * - Problems that can be solved easily (OneMax,
         *   LeadingOnes) will be _still_ solved late in the runs.
         *   This is normal for simulated annealing, as the goal
         *   is to transcend from random walk to hill climbing
         *   behavior.
         * - Problems whose objective value range differs largely
         *   from n (e.g., LARS) will probably not be solved well.
         *
         * Author: Thomas Weise
         *         Institute of Applied Optimization
         *         Hefei University
         *         Hefei, Anhui, China
         * Email: tweise@hfuu.edu.cn, tweise@ustc.edu.cn
         */

        // Compute the acceptance probability from DeltaE and the temperature
        // DeltaE must be positive
        inline static double p_accept(const double DeltaE, const double temperature)
        {
            return exp(-DeltaE / temperature);
        }

        // Compute the temperature at a given step
        // Tstart = start temperatre
        // epsilon = temperature decrease parameter
        // step = step index; must be >= 1
        inline static double temperature(const double Tstart, const double epsilon,
                                         const unsigned long long step)
        {
            return Tstart * pow(1.0 - epsilon, step - 1);
        }

        // compute a temperature from a given DeltaE and acceptance probability P
        inline static double T_from_DeltaE_and_P(const double DeltaE, const double P)
        {
            if ((!isfinite(DeltaE)) || (DeltaE <= 0.0))
                throw "DeltaE must be positive and finite.";
            if ((!isfinite(P)) || (P <= 0.0) || (P >= 1.0))
                throw "P must be finite and from (0,1).";
            return -DeltaE / log(P);
        }

        inline static double epsilon_from_T_and_step(const double Tstart,
                                                     const double Tstep, const unsigned long long step)
        {
            if ((!isfinite(Tstart)) || (Tstart <= 0.0))
                throw "Tstart must be positive and finite.";
            if ((!isfinite(Tstep)) || (Tstep >= Tstart))
                throw "Tstep must finite and > Tstart";
            const double epsilon = 1.0 - pow(Tstep / Tstart, 1.0 / (step - 1));
            if ((!isfinite(epsilon)) || (epsilon <= 0.0) || (epsilon >= 1.0))
                throw "epsilon must be positive and from (0,1)";
            return epsilon;
        }

        void simulated_annealing_exp(shared_ptr<ioh::problem::Integer> problem,
                                     shared_ptr<ioh::logger::Analyzer> logger,
                                     const unsigned long long eval_budget)
        {
            // check input variables
            if (eval_budget <= 1)
                throw "eval_budget must be > 1";
            if (!problem)
                throw "problem must not be null";
            if (!logger)
                throw "logger must not be null";

            // n be the number of variables
            const int n = problem->meta_data().n_variables;
            if (n <= 0)
                throw "number of variables must be positive";
            // the bit flip probability
            const double p = 1.0 / ((double)n);
            if ((!isfinite(p)) || (p <= 0.0) || (p >= 1.0))
                throw "p must be from (0,1)";

            // perform the auto-configuration
            const double Tstart = T_from_DeltaE_and_P(max(1.0, n / 4.0), 0.1);
            if ((!isfinite(Tstart)) || (Tstart < 1))
                throw "Tstart must be >= 1";
            const double Tend = T_from_DeltaE_and_P(1.0, 1.0 / sqrt(eval_budget));
            if ((!isfinite(Tend)) || (Tend >= Tstart))
                throw "Tend must be < Tstart";
            const double epsilon = epsilon_from_T_and_step(Tstart, Tend, eval_budget);
            if ((!isfinite(epsilon)) || (epsilon <= 0) || (epsilon >= 1))
                throw "epsilon must be in (0,1)";

            // store all the parameters
            logger->watch(ioh::watch::address("Tstart",&Tstart));
            logger->watch(ioh::watch::address("Tend",&Tend));
            logger->watch(ioh::watch::address("epsilon",&epsilon));
            logger->watch(ioh::watch::address("p",&p));


            // xcur is the current best candidate solution (based on frequency fitness)
            std::vector<int> xcur;
            // xnew is the new candidate solution generated in each step
            std::vector<int> xnew;
            // ycur is the objective value of the current solution
            double ycur = std::numeric_limits<double>::infinity();
            // ynew is the objective value of the new solution
            double ynew = std::numeric_limits<double>::infinity();

            // first we generate the random initial solution
            xcur.reserve(n);
            for (int i = 0; i < n; i++)
            {
                xcur.push_back((int)(2 * uniform_random()));
            }
            // we evaluate the random initial solution
            ycur = (*problem)(xcur);
            // we perform iterations until either the optimum is discovered or the budget has been exhausted
            unsigned long long int step = 1;
            while (((++step) <= eval_budget) && (!problem->state().optimum_found))
            {

                // copy the current solution to the new solution
                xnew = xcur;
                bool unchanged = true;
                // until the solution changes, repeat
                do
                {
                    // flip each bit with the independent probability of 1/n
                    for (int i = n; (--i) >= 0;)
                    {
                        if (uniform_random() < p)
                        {
                            unchanged = false; // there was a change
                            xnew[i] ^= 1;      // flip the bit
                        }
                    }
                } while (unchanged); // repeat until at least one change

                // evaluate the new candidate solution
                ynew = (*problem)(xnew);

                // if new solution is at least as good as current one, accept it
                // otherwise: if check if it is acceptable at the current temperature
                if ((ynew >= ycur) || (uniform_random() < p_accept(ycur - ynew, temperature(Tstart, epsilon, step))))
                {
                    ycur = ynew;
                    xcur = xnew;
                }
            }
        }

        // run the simulated annealing algorithm with automatic configuration
        void run_simulated_annealing_exp(const string folder_path,
                                         shared_ptr<ioh::suite::Suite<ioh::problem::Integer>> suite,
                                         const unsigned long long eval_budget,
                                         const unsigned long long independent_runs,
                                         const unsigned long long rand_seed)
        {
            if (folder_path.empty())
                throw "folder path cannot be empty";
            if (!suite)
                throw "suite cannot be null";
            if (eval_budget <= 1)
                throw "eval_budget must be > 1";
            if (independent_runs < 1)
                throw "independent_runs must be > 0";

            const string algorithm_name = "sa_auto";
            std::shared_ptr<ioh::logger::Analyzer> logger(new ioh::logger::Analyzer(
                {ioh::trigger::on_improvement}, // trigger when the objective value improves
                {},                             // no additional properties
                folder_path,                    // path to store data
                algorithm_name,                 // name of the folder in path, which will be newly created
                algorithm_name,                 // name of the algoritm
                algorithm_name,                 // additional info about the algorithm
                false                           // where to store x positions in the data files
                ));

            random_gen.seed(rand_seed);

            for (auto &problem : *suite)
            {
                problem->attach_logger(*logger);
                for (unsigned long long i = 0; i < independent_runs; i++)
                {
                    simulated_annealing_exp(problem, logger, eval_budget);
                    problem->reset();
                }
            }
        }
    }
}

#endif
