#pragma once
/**
 * Here we implement the Simulated Annealing algorithm with
 * exponential temperature schedule with iterative restarts.
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
 * The dependency on the budget is reduced here by performing
 * internal restarts. The first restart has a budget of 1024
 * FEs, the second one of 2048 FEs, and so on: the "inner"
 * runs always double in budget. Thus, the advantages of SA
 * can be obtained faster by wasting about half of the
 * computational budget.
 *
 * Still, to be able to do some experiments, we need to also
 * configure the temperatures. We propose the following
 * automatic setup based on the number 'n' of variables of
 * the problem and the budgets of the single internal runs:
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
#ifndef _SARS_H_
#define _SARS_H_
#include "utils/common.hpp"
#include <math.h>
namespace modularGA
{

    namespace sars
    {
        /**
         * Here we implement the Simulated Annealing algorithm with
         * exponential temperature schedule with iterative restarts.
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
         * The dependency on the budget is reduced here by performing
         * internal restarts. The first restart has a budget of 1024
         * FEs, the second one of 2048 FEs, and so on: the "inner"
         * runs always double in budget. Thus, the advantages of SA
         * can be obtained faster by wasting about half of the
         * computational budget.
         *
         * Still, to be able to do some experiments, we need to also
         * configure the temperatures. We propose the following
         * automatic setup based on the number 'n' of variables of
         * the problem and the budgets of the single internal runs:
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

        void simulated_annealing_exp_rs(shared_ptr<ioh::problem::IntegerSingleObjective> problem,
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

            // xcur is the current best candidate solution (based on frequency fitness)
            std::vector<int> xcur;
            xcur.reserve(n);
            for (int i = n; (--i) >= 0;)
            {
                xcur.push_back(0);
            }

            // xnew is the new candidate solution generated in each step
            std::vector<int> xnew;
            // ycur is the objective value of the current solution
            double ycur = std::numeric_limits<double>::infinity();
            // ynew is the objective value of the new solution
            double ynew = std::numeric_limits<double>::infinity();

            const double Tstart = T_from_DeltaE_and_P(max(1.0, n / 4.0), 0.1);
            if ((!isfinite(Tstart)) || (Tstart < 1))
                throw "Tstart must be >= 1";

            double Tend = 0;
            double epsilon = 0;
            double innerBudgetParam = 0;

            // store all the parameters
            logger->watch(ioh::watch::address("Tstart",&Tstart));
            logger->watch(ioh::watch::address("Tend",&Tend));
            logger->watch(ioh::watch::address("epsilon",&epsilon));
            logger->watch(ioh::watch::address("innerBudget",&innerBudgetParam));
            logger->watch(ioh::watch::address("p",&p));

            unsigned long long int innerBudget = 512;
            unsigned long long int stepMain = 1;

            // the main loop containing the inner, independent runs
            while (((++stepMain) <= eval_budget) && (!problem->state().optimum_found))
            {
                // set the budget for the next inner run
                innerBudget += innerBudget;
                innerBudgetParam = (double)innerBudget;
                ++stepMain;

                // perform the auto-configuration for the inner runs
                Tend = T_from_DeltaE_and_P(1.0, 1.0 / sqrt(innerBudget));
                if ((!isfinite(Tend)) || (Tend >= Tstart))
                    throw "Tend must be < Tstart";
                epsilon = epsilon_from_T_and_step(Tstart, Tend, innerBudget);
                if ((!isfinite(epsilon)) || (epsilon <= 0) || (epsilon >= 1))
                    throw "epsilon must be in (0,1)";

                ycur = std::numeric_limits<double>::infinity();
                ynew = std::numeric_limits<double>::infinity();

                // first we generate the random initial solution
                for (int i = 0; i < n; i++)
                {
                    xcur[i] = (int)(2 * uniform_random());
                }
                // we evaluate the random initial solution
                ycur = (*problem)(xcur);

                // we perform iterations until either the optimum is discovered or the budget has been exhausted
                unsigned long long int step = 1;
                while (((++stepMain) <= eval_budget) && ((++step) <= innerBudget) && (!problem->state().optimum_found))
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
        }

        // run the simulated annealing algorithm with automatic configuration and restarts
        void run_simulated_annealing_exp_rs(const string folder_path,
                                            shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite,
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

            const string algorithm_name = "sars_auto";
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
                    simulated_annealing_exp_rs(problem, logger, eval_budget);
                    problem->reset();
                }
            }
        }
    }
}

#endif
