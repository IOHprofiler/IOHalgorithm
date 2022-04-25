/**
 * Here we implement the (1+1)-FEA, i.e., the (1+1)-EA>0 with frequency
 * fitness assignment
 *
 * The (1+1)-FEA is based on the (1+1)-EA:
 * The (1+1)-EA>0, in turn, is a slight extension of
 * the off-the-shelf evolutionary algorithm with one parent
 * (mu=1) that generates one single offspring solution in each
 * step (lambda=1) by flipping each bit of the parent solution
 * with the same probability of 1/n.
 * The difference of the (1+1)-EA<0 and the (1+1)-EA is that
 * in the (1+1)-EA<0, it is ensured that  always at least one
 * bit is flipped.
 * In the (1+1)-FEA, we use the encounter frequency of objective
 * values as fitness. This fitness is subject to minimization.
 *
 * The (1+1)-FEA is introduced in:
 * 1. Thomas Weise, Zhize Wu, Xinlu Li, and Yan Chen. Frequency Fitness
 *    Assignment: Making Optimization Algorithms Invariant under Bijective
 *    Transformations of the Objective Function Value. IEEE Transactions
 *    on Evolutionary Computation. doi:10.1109/TEVC.2020.3032090
 * 2. The preprint of the above paper at arXiv:2001.01416v5 [cs.NE] 15 Oct
 *    2020. http://arxiv.org/abs/2001.01416
 *
 * Author: Thomas Weise
 *         Institute of Applied Optimization
 *         Hefei University
 *         Hefei, Anhui, China
 * Email: tweise@hfuu.edu.cn, tweise@ustc.edu.cn
 */

#ifndef _OPOFEA_H_
#define _OPOFEA_H_
#include "common.hpp"
#include <tr1/unordered_map>
namespace modularGA
{
    namespace fea
    {

        // void run_fea1p1(const string folder_path,
        //                 shared_ptr<IOHprofiler_suite<int>> suite,
        //                 const unsigned long long eval_budget,
        //                 const unsigned long long independent_runs,
        //                 const unsigned long long rand_seed);

        /**
         * Here we implement the (1+1)-FEA, i.e., the (1+1)-EA>0 with frequency
         * fitness assignment
         *
         * The (1+1)-FEA is based on the (1+1)-EA:
         * The (1+1)-EA>0, in turn, is a slight extension of
         * the off-the-shelf evolutionary algorithm with one parent
         * (mu=1) that generates one single offspring solution in each
         * step (lambda=1) by flipping each bit of the parent solution
         * with the same probability of 1/n.
         * The difference of the (1+1)-EA<0 and the (1+1)-EA is that
         * in the (1+1)-EA<0, it is ensured that  always at least one
         * bit is flipped.
         * In the (1+1)-FEA, we use the encounter frequency of objective
         * values as fitness. This fitness is subject to minimization.
         *
         * The (1+1)-FEA is introduced in:
         * 1. Thomas Weise, Zhize Wu, Xinlu Li, and Yan Chen. Frequency Fitness
         *    Assignment: Making Optimization Algorithms Invariant under Bijective
         *    Transformations of the Objective Function Value. IEEE Transactions
         *    on Evolutionary Computation. doi:10.1109/TEVC.2020.3032090
         * 2. The preprint of the above paper at arXiv:2001.01416v5 [cs.NE] 15 Oct
         *    2020. http://arxiv.org/abs/2001.01416
         *
         * Author: Thomas Weise
         *         Institute of Applied Optimization
         *         Hefei University
         *         Hefei, Anhui, China
         * Email: tweise@hfuu.edu.cn, tweise@ustc.edu.cn
         */

        void fea1p1(shared_ptr<ioh::problem::Integer> problem,
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
            // xnew is the new candidate solution generated in each step
            std::vector<int> xnew;
            // ycur is the objective value of the current solution
            double ycur = std::numeric_limits<double>::infinity();
            // ynew is the objective value of the new solution
            double ynew = std::numeric_limits<double>::infinity();

            // H is the frequency table, storing the encounter frequency of objective values
            std::tr1::unordered_map<double, unsigned long long> H;

            // first we generate the random initial solution
            xcur.reserve(n);
            for (int i = 0; i < n; i++)
            {
                xcur.push_back((int)(2 * uniform_random()));
            }

            // we evaluate the random initial solution
            ycur = (*problem)(xcur);
            H[ycur] = 0; // and initialize its frequency

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

                // is it completely new?
                if (H.find(ynew) == H.end())
                {
                    H[ynew] = 0; // store frequency 0
                }
                // update the frequencies of both solutions
                ++H[ycur];
                ++H[ynew];
                // if the new solution has a lower or equal frequency, take it
                if (H[ynew] <= H[ycur])
                {
                    ycur = ynew;
                    xcur = xnew;
                }
            }
        }

        // run the ea with frequency fitness assignment
        void run_fea1p1(const string folder_path,
                        shared_ptr<ioh::suite::Suite<ioh::problem::Integer> > suite,
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

            const string algorithm_name = "opofea";
            std::shared_ptr<ioh::logger::Analyzer> logger(
                new ioh::logger::Analyzer(
                                                      {ioh::trigger::on_improvement}, // trigger when the objective value improves
                                                      {},                   // no additional properties 
                                                      folder_path,        // path to store data
                                                      algorithm_name,               // name of the folder in path, which will be newly created
                                                      algorithm_name,                     // name of the algoritm 
                                                      algorithm_name,                   // additional info about the algorithm              
                                                      false            // where to store x positions in the data files 
                                                    ));

            random_gen.seed(rand_seed);
            shared_ptr<ioh::problem::Integer> problem;
            for (auto &problem : *suite) 
            {
                problem->attach_logger(*logger);
                for (unsigned long long i = 0; i < independent_runs; i++)
                {
                    fea1p1(problem, logger, eval_budget);
                    problem->reset();
                }
            }
        }
    }
}
#endif
