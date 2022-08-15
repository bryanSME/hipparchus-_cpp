#pragma once
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

//package org.hipparchus.optim.nonlinear.scalar.noderiv;

//import java.util.Array_list;
//import java.util.Arrays;
//import java.util.List;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.linear.Eigen_Decomposition;
//import org.hipparchus.linear.Matrix_Utils;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.optim.Point_valuePair;
//import org.hipparchus.optim.nonlinear.scalar.Goal_Type;
//import org.hipparchus.optim.nonlinear.scalar.Multivariate_Optimizer;
//import org.hipparchus.random.Random_Generator;
//import org.hipparchus.util.FastMath;
#include "../../../../core/linear/MatrixUtils.h"

/**
 * An implementation of the active Covariance Matrix Adaptation Evolution Strategy (CMA-ES)
 * for non-linear, non-convex, non-smooth, global function minimization.
 * <p>
 * The CMA-Evolution Strategy (CMA-ES) is a reliable stochastic optimization method
 * which should be applied if derivative-based methods, e.g. quasi-Newton BFGS or
 * conjugate gradient, fail due to a rugged search landscape (e.g. noise, local
 * optima, outlier, etc.) of the objective function. Like a
 * quasi-Newton method, the CMA-ES learns and applies a variable metric
 * on the underlying search space. Unlike a quasi-Newton method, the
 * CMA-ES neither estimates nor uses gradients, making it considerably more
 * reliable in terms of finding a good, or even close to optimal, solution.
 * <p>
 * In general, on smooth objective functions the CMA-ES is roughly ten times
 * slower than BFGS (counting objective function evaluations, no gradients provided).
 * For up to <math>N=10</math> variables also the derivative-free simplex
 * direct search method (Nelder and Mead) can be faster, but it is
 * far less reliable than CMA-ES.
 * <p>
 * The CMA-ES is particularly well suited for non-separable
 * and/or badly conditioned problems. To observe the advantage of CMA compared
 * to a conventional evolution strategy, it will usually take about
 * <math>30 N</math> function evaluations. On difficult problems the complete
 * optimization (a single run) is expected to take <em>roughly</em> between
 * <math>30 N</math> and <math>300 N<sup>2</sup></math>
 * function evaluations.
 * <p>
 * This implementation is translated and adapted from the Matlab version
 * of the CMA-ES algorithm as implemented in module {@code cmaes.m} version 3.51.
 * <p>
 * For more information, please refer to the following links:
 * <ul>
 *  <li><a href="http://www.lri.fr/~hansen/cmaes.m">Matlab code</a></li>
 *  <li><a href="http://www.lri.fr/~hansen/cmaesintro.html">Introduction to CMA-ES</a></li>
 *  <li><a href="http://en.wikipedia.org/wiki/CMA-ES">Wikipedia</a></li>
 * </ul>
 *
 */
class CMAES_Optimizer
    extends Multivariate_Optimizer 
    {
    // global search parameters
    /**
     * Population size, offspring number. The primary strategy parameter to play
     * with, which can be increased from its default value. Increasing the
     * population size improves global search properties in exchange to speed.
     * Speed decreases, as a rule, at most linearly with increasing population
     * size. It is advisable to begin with the default small population size.
     */
    private int lambda; // population size
    /**
     * Covariance update mechanism, default is active CMA. is_active_c_m_a = true
     * turns on "active CMA" with a negative update of the covariance matrix and
     * checks for positive definiteness. OPTS.CMA.active = 2 does not check for
     * pos. def. and is numerically faster. Active CMA usually speeds up the
     * adaptation.
     */
    private const bool is_active_c_m_a;
    /**
     * Determines how often a random offspring is generated in case it is
     * not feasible / beyond the defined limits, default is 0.
     */
    private const int check_feasable_count;
    /**
     * @see Sigma
     */
    private std::vector<double> input_sigma;
    /** Number of objective variables/problem dimension */
    private int dimension;
    /**
     * Defines the number of initial iterations, where the covariance matrix
     * remains diagonal and the algorithm has internally linear time complexity.
     * diagonal_only = 1 means keeping the covariance matrix always diagonal and
     * this setting also exhibits linear space complexity. This can be
     * particularly useful for dimension > 100.
     * @see <a href="http://hal.archives-ouvertes.fr/inria-00287367/en">A Simple Modification in CMA-ES</a>
     */
    private int diagonal_only;
    /** Number of objective variables/problem dimension */
    private bool is_minimize = true;
    /** Indicates whether statistic data is collected. */
    private const bool generate_statistics;

    // termination criteria
    /** Maximal number of iterations allowed. */
    private const int max_iterations;
    /** Limit for fitness value. */
    private const double stop_fitness;
    /** Stop if x-changes larger stop_tol_up_x. */
    private double stop_tol_up_x;
    /** Stop if x-change smaller stop_tol_x. */
    private double stop_tol_x;
    /** Stop if fun-changes smaller stop_tol_fun. */
    private double stop_tol_fun;
    /** Stop if back fun-changes smaller stop_tol_hist_fun. */
    private double stop_tol_hist_fun;

    // selection strategy parameters
    /** Number of parents/points for recombination. */
    private int mu; //
    /** log(mu + 0.5), stored for efficiency. */
    private double log_mu2;   // NOPMD - using a field here is for performance reasons
    /** Array for weighted recombination. */
    private Real_Matrix weights;
    /** Variance-effectiveness of sum w_i x_i. */
    private double mueff; //

    // dynamic strategy parameters and constants
    /** Overall standard deviation - search volume. */
    private double sigma;
    /** Cumulation constant. */
    private double cc;
    /** Cumulation constant for step-size. */
    private double cs;
    /** Damping for step-size. */
    private double damps;
    /** Learning rate for rank-one update. */
    private double ccov1;
    /** Learning rate for rank-mu update' */
    private double ccovmu;
    /** Expectation of ||N(0,I)|| == norm(randn(N,1)). */
    private double chi_n;
    /** Learning rate for rank-one update - diagonal_only */
    private double ccov1_sep;
    /** Learning rate for rank-mu update - diagonal_only */
    private double ccovmu_sep;

    // CMA internal values - updated each generation
    /** Objective variables. */
    private Real_Matrix xmean;
    /** Evolution path. */
    private Real_Matrix pc;
    /** Evolution path for sigma. */
    private Real_Matrix ps;
    /** Norm of ps, stored for efficiency. */
    private double normps;
    /** Coordinate system. */
    private Real_Matrix B;
    /** Scaling. */
    private Real_Matrix D;
    /** B*D, stored for efficiency. */
    private Real_Matrix BD;
    /** Diagonal of sqrt(D), stored for efficiency. */
    private Real_Matrix diag_d;
    /** Covariance matrix. */
    private Real_Matrix C;
    /** Diagonal of C, used for diagonal_only. */
    private Real_Matrix diag_c;
    /** Number of iterations already performed. */
    private int iterations;

    /** History queue of best values. */
    private std::vector<double> fitness_history;

    /** Random generator. */
    private const Random_Generator random;

    /** History of sigma values. */
    private const List<Double> statistics_sigma_history;
    /** History of mean matrix. */
    private const List<Real_Matrix> statistics_mean_history;
    /** History of fitness values. */
    private const List<Double> statistics_fitness_history;
    /** History of D matrix. */
    private const List<Real_Matrix> statistics_d_history;

    /**
     * @param max_iterations Maximal number of iterations.
     * @param stop_fitness Whether to stop if objective function value is smaller than
     * {@code stop_fitness}.
     * @param is_active_c_m_a Chooses the covariance matrix update method.
     * @param diagonal_only Number of initial iterations, where the covariance matrix
     * remains diagonal.
     * @param check_feasable_count Determines how often random objective variables are
     * generated in case they are out of bounds.
     * @param random Random generator.
     * @param generate_statistics Whether statistic data is collected.
     * @param checker Convergence checker.
     *
     */
    public CMAES_Optimizer(const int& max_iterations, double stop_fitness, bool is_active_c_m_a, int diagonal_only, int check_feasable_count, Random_Generator random, bool generate_statistics, Convergence_Checker<Point_valuePair> checker) 
    {
        super(checker);
        this.max_iterations = max_iterations;
        this.stop_fitness = stop_fitness;
        this.is_active_c_m_a = is_active_c_m_a;
        this.diagonal_only = diagonal_only;
        this.check_feasable_count = check_feasable_count;
        this.random = random;
        this.generate_statistics = generate_statistics;
        this.statistics_sigma_history = Array_list<>();
        this.statistics_mean_history = Array_list<>();
        this.statistics_fitness_history = Array_list<>();
        this.statistics_d_history = Array_list<>();
    }

    /**
     * @return History of sigma values.
     */
    public List<Double> get_statistics_sigma_history() 
    {
        return statistics_sigma_history;
    }

    /**
     * @return History of mean matrix.
     */
    public List<Real_Matrix> get_statistics_mean_history() 
    {
        return statistics_mean_history;
    }

    /**
     * @return History of fitness values.
     */
    public List<Double> get_statistics_fitness_history() 
    {
        return statistics_fitness_history;
    }

    /**
     * @return History of D matrix.
     */
    public List<Real_Matrix> get_statistics_d_history() 
    {
        return statistics_d_history;
    }

    /**
     * Input sigma values.
     * They define the initial coordinate-wise standard deviations for
     * sampling search points around the initial guess.
     * It is suggested to set them to the estimated distance from the
     * initial to the desired optimum.
     * Small values induce the search to be more local (and very small
     * values are more likely to find a local optimum close to the initial
     * guess).
     * Too small values might however lead to early termination.
     */
    public static class Sigma : Optimization_data 
    {
        /** Sigma values. */
        private const std::vector<double> s;

        /**
         * @param s Sigma values.
         * @ if any of the array entries is smaller
         * than zero.
         */
        public Sigma(std::vector<double> s)
             
            {
            for (int i{}; i < s.size(); i++) 
            {
                if (s[i] < 0) 
                {
                    throw (Localized_Core_Formats.NUMBER_TOO_SMALL, s[i], 0);
                }
            }

            this.s = s.clone();
        }

        /**
         * @return the sigma values.
         */
        public std::vector<double> get_sigma() 
        {
            return s.clone();
        }
    }

    /**
     * Population size.
     * The number of offspring is the primary strategy parameter.
     * In the absence of better clues, a good default could be an
     * integer close to {@code 4 + 3 ln(n)}, where {@code n} is the
     * number of optimized parameters.
     * Increasing the population size improves global search properties
     * at the expense of speed (which in general decreases at most
     * linearly with increasing population size).
     */
    public static class Population_Size : Optimization_data 
    {
        /** Population size. */
        private const int lambda;

        /**
         * @param size Population size.
         * @ if {@code size <= 0}.
         */
        public Population_Size(const int& size)
             
            {
            if (size <= 0) 
            {
                throw (Localized_Core_Formats.NUMBER_TOO_SMALL_BOUND_EXCLUDED, size, 0);
            }
            lambda = size;
        }

        /**
         * @return the population size.
         */
        public int get_population_size() 
        {
            return lambda;
        }
    }

    /**
     * {@inherit_doc}
     *
     * @param opt_data Optimization data. In addition to those documented in
     * {@link Multivariate_Optimizer#parse_optimization_data(Optimization_data[])
     * Multivariate_Optimizer}, this method will register the following data:
     * <ul>
     *  <li>{@link Sigma}</li>
     *  <li>{@link Population_Size}</li>
     * </ul>
     * @return {@inherit_doc}
     * @Math_Illegal_State_Exception if the maximal number of
     * evaluations is exceeded.
     * @ if the initial guess, target, and weight
     * arguments have inconsistent dimensions.
     */
    //override
    public Point_valuePair optimize(Optimization_data... opt_data)
        , Math_Illegal_State_Exception 
        {
        // Set up base class and perform computation.
        return super.optimize(opt_data);
    }

    /** {@inherit_doc} */
    //override
    protected Point_valuePair do_optimize() 
    {
         // -------------------- Initialization --------------------------------
        is_minimize = get_goal_type().equals(Goal_Type.MINIMIZE);
        const Fitness_Function fitfun = Fitness_Function();
        const std::vector<double> guess = get_start_point();
        // number of objective variables/problem dimension
        dimension = guess.size();
        initialize_c_m_a(guess);
        iterations = 0;
        Value_Penalty_Pair value_penalty = fitfun.value(guess);
        double best_value = value_penalty.value+value_penalty.penalty;
        push(fitness_history, best_value);
        Point_valuePair optimum
            = Point_valuePair(get_start_point(), is_minimize ? best_value : -best_value);
        Point_valuePair last_result = NULL;

        // -------------------- Generation Loop --------------------------------

        generation_loop:
        for (iterations = 1; iterations <= max_iterations; iterations++) 
        {
            increment_iteration_count();

            // Generate and evaluate lambda offspring
            const Real_Matrix& arz = randn1(dimension, lambda);
            const Real_Matrix& arx = zeros(dimension, lambda);
            const std::vector<double> fitness = std::vector<double>(lambda];
            const Value_Penalty_Pair[] value_penalty_pairs = Value_Penalty_Pair[lambda];
            // generate random offspring
            for (int k{}; k < lambda; k++) 
            {
                Real_Matrix arxk = NULL;
                for (int i{}; i < check_feasable_count + 1; i++) 
                {
                    if (diagonal_only <= 0) 
                    {
                        arxk = xmean.add(BD.multiply(arz.get_column_matrix(k))
                                         .scalar_multiply(sigma)); // m + sig * Normal(0,C)
                    }
else 
                    {
                        arxk = xmean.add(times(diag_d,arz.get_column_matrix(k))
                                         .scalar_multiply(sigma));
                    }
                    if (i >= check_feasable_count ||
                        fitfun.is_feasible(arxk.get_column(0))) 
                        {
                        break;
                    }
                    // regenerate random arguments for row
                    arz.set_column(k, randn(dimension));
                }
                copy_column(arxk, 0, arx, k);
                try 
                {
                    value_penalty_pairs[k] = fitfun.value(arx.get_column(k)); // compute fitness
                }
catch (Math_Illegal_State_Exception e) 
                {
                    break generation_loop;
                }
            }

            // Compute fitnesses by adding value and penalty after scaling by value range.
            double value_range = value_range(value_penalty_pairs);
            for (const int& i_value=0;i_value<value_penalty_pairs.size();i_value++) 
            {
                 fitness[i_value] = value_penalty_pairs[i_value].value + value_penalty_pairs[i_value].penalty*value_range;
            }

            // Sort by fitness and compute weighted mean into xmean
            const std::vector<int> arindex = sorted_indices(fitness);
            // Calculate xmean, this is selection and recombination
            const Real_Matrix xold = xmean; // for speed up of Eq. (2) and (3)
            const Real_Matrix best_arx = select_columns(arx, Arrays.copy_of(arindex, mu));
            xmean = best_arx.multiply(weights);
            const Real_Matrix best_arz = select_columns(arz, Arrays.copy_of(arindex, mu));
            const Real_Matrix zmean = best_arz.multiply(weights);
            const bool hsig = update_evolution_paths(zmean, xold);
            if (diagonal_only <= 0) 
            {
                update_covariance(hsig, best_arx, arz, arindex, xold);
            }
else 
            {
                update_covariance_diagonal_only(hsig, best_arz);
            }
            // Adapt step size sigma - Eq. (5)
            sigma *= std::exp(std::min(1, (normps/chi_n - 1) * cs / damps));
            const double best_fitness = fitness[arindex[0]];
            const double worst_fitness = fitness[arindex[arindex.size() - 1]];
            if (best_value > best_fitness) 
            {
                best_value = best_fitness;
                last_result = optimum;
                optimum = Point_valuePair(fitfun.repair(best_arx.get_column(0)), is_minimize ? best_fitness : -best_fitness);
                if (get_convergence_checker() != NULL && last_result != NULL &&
                    get_convergence_checker().converged(iterations, optimum, last_result)) 
                    {
                    break generation_loop;
                }
            }
            // handle termination criteria
            // Break, if fitness is good enough
            if (stop_fitness != 0 && best_fitness < (is_minimize ? stop_fitness : -stop_fitness)) 
            {
                break generation_loop;
            }
            const std::vector<double> sqrt_diagC = sqrt(diag_c).get_column(0);
            const std::vector<double> pc_col = pc.get_column(0);
            for (int i{}; i < dimension; i++) 
            {
                if (sigma * std::max(std::abs(pc_col[i]), sqrt_diagC[i]) > stop_tol_x) 
                {
                    break;
                }
                if (i >= dimension - 1) 
                {
                    break generation_loop;
                }
            }
            for (int i{}; i < dimension; i++) 
            {
                if (sigma * sqrt_diagC[i] > stop_tol_up_x) 
                {
                    break generation_loop;
                }
            }
            const double history_best = min(fitness_history);
            const double history_worst = max(fitness_history);
            if (iterations > 2 &&
                std::max(history_worst, worst_fitness) -
                std::min(history_best, best_fitness) < stop_tol_fun) 
                {
                break generation_loop;
            }
            if (iterations > fitness_history.size() &&
                history_worst - history_best < stop_tol_hist_fun) 
                {
                break generation_loop;
            }
            // condition number of the covariance matrix exceeds 1e14
            if (max(diag_d) / min(diag_d) > 1e7) 
            {
                break generation_loop;
            }
            // user defined termination
            if (get_convergence_checker() != NULL) 
            {
                const Point_valuePair current
                    = Point_valuePair(best_arx.get_column(0), is_minimize ? best_fitness : -best_fitness);
                if (last_result != NULL &&
                    get_convergence_checker().converged(iterations, current, last_result)) 
                    {
                    break generation_loop;
                    }
                last_result = current;
            }
            // Adjust step size in case of equal function values (flat fitness)
            if (best_value == fitness[arindex[static_cast<int>((0.1+lambda/4.)]]) 
            {
                sigma *= std::exp(0.2 + cs / damps);
            }
            if (iterations > 2 && std::max(history_worst, best_fitness) -
                std::min(history_best, best_fitness) == 0) 
                {
                sigma *= std::exp(0.2 + cs / damps);
            }
            // store best in history
            push(fitness_history,best_fitness);
            if (generate_statistics) 
            {
                statistics_sigma_history.add(sigma);
                statistics_fitness_history.add(best_fitness);
                statistics_mean_history.add(xmean.transpose());
                statistics_d_history.add(diag_d.transpose().scalar_multiply(1E5));
            }
        }
        return optimum;
    }

    /**
     * Scans the list of (required and optional) optimization data that
     * characterize the problem.
     *
     * @param opt_data Optimization data. The following data will be looked for:
     * <ul>
     *  <li>{@link Sigma}</li>
     *  <li>{@link Population_Size}</li>
     * </ul>
     */
    //override
    protected void parse_optimization_data(Optimization_data... opt_data) 
    {
        // Allow base class to register its own data.
        super.parse_optimization_data(opt_data);

        // The existing values (as set by the previous call) are reused if
        // not provided in the argument list.
        for (Optimization_data data : opt_data) 
        {
            if (data instanceof Sigma) 
            {
                input_sigma = ((Sigma) data).get_sigma();
                continue;
            }
            if (data instanceof Population_Size) 
            {
                lambda = ((Population_Size) data).get_population_size();
                continue;
            }
        }

        check_parameters();
    }

    /**
     * Checks dimensions and values of boundaries and input_sigma if defined.
     */
    private void check_parameters() 
    {
        if (input_sigma != NULL) 
        {
            const std::vector<double> init = get_start_point();

            if (input_sigma.size() != init.size()) 
            {
                throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, input_sigma.size(), init.size());
            }

            const std::vector<double> lB = get_lower_bound();
            const std::vector<double> uB = get_upper_bound();

            for (int i{}; i < init.size(); i++) 
            {
                if (input_sigma[i] > uB[i] - lB[i]) 
                {
                    throw (Localized_Core_Formats.OUT_OF_RANGE_SIMPLE, input_sigma[i], 0, uB[i] - lB[i]);
                }
            }
        }
    }

    /**
     * Initialization of the dynamic search parameters
     *
     * @param guess Initial guess for the arguments of the fitness function.
     */
    private void initialize_c_m_a(std::vector<double> guess) 
    {
        if (lambda <= 0) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL_BOUND_EXCLUDED, lambda, 0);
        }
        // initialize sigma
        const std::vector<std::vector<double>> sigma_array = std::vector<double>(guess.size()][1];
        for (int i{}; i < guess.size(); i++) 
        {
            sigma_array[i][0] = input_sigma[i];
        }
        const Real_Matrix insigma = Array_2D_Row_Real_Matrix(sigma_array, false);
        sigma = max(insigma); // overall standard deviation

        // initialize termination criteria
        stop_tol_up_x = 1e3 * max(insigma);
        stop_tol_x = 1e-11 * max(insigma);
        stop_tol_fun = 1e-12;
        stop_tol_hist_fun = 1e-13;

        // initialize selection strategy parameters
        mu = lambda / 2; // number of parents/points for recombination
        log_mu2 = std::log(mu + 0.5);
        weights = log(sequence(1, mu, 1)).scalar_multiply(-1).scalar_add(log_mu2);
        double sumw = 0;
        double sumwq = 0;
        for (int i{}; i < mu; i++) 
        {
            double w = weights.get_entry(i, 0);
            sumw += w;
            sumwq += w * w;
        }
        weights = weights.scalar_multiply(1 / sumw);
        mueff = sumw * sumw / sumwq; // variance-effectiveness of sum w_i x_i

        // initialize dynamic strategy parameters and constants
        cc = (4 + mueff / dimension) /
                (dimension + 4 + 2 * mueff / dimension);
        cs = (mueff + 2) / (dimension + mueff + 3.);
        damps = (1 + 2 * std::max(0, std::sqrt((mueff - 1) /
                                                       (dimension + 1)) - 1)) *
            std::max(0.3, 1 - dimension / (1e-6 + max_iterations)) + cs; // minor increment
        ccov1 = 2 / ((dimension + 1.3) * (dimension + 1.3) + mueff);
        ccovmu = std::min(1 - ccov1, 2 * (mueff - 2 + 1 / mueff) /
                              ((dimension + 2) * (dimension + 2) + mueff));
        ccov1_sep = std::min(1, ccov1 * (dimension + 1.5) / 3);
        ccovmu_sep = std::min(1 - ccov1, ccovmu * (dimension + 1.5) / 3);
        chi_n = std::sqrt(dimension) *
                (1 - 1 / (static_cast<double>( 4 * dimension) + 1 / (static_cast<double>( 21 * dimension * dimension));
        // intialize CMA internal values - updated each generation
        xmean = Matrix_Utils.create_column_real__matrix(guess); // objective variables
        diag_d = insigma.scalar_multiply(1 / sigma);
        diag_c = square(diag_d);
        pc = zeros(dimension, 1); // evolution paths for C and sigma
        ps = zeros(dimension, 1); // B defines the coordinate system
        normps = ps.get_frobenius_norm();

        B = eye(dimension, dimension);
        D = ones(dimension, 1); // diagonal D defines the scaling
        BD = times(B, repmat(diag_d.transpose(), dimension, 1));
        C = B.multiply(diag(square(D)).multiply(B.transpose())); // covariance
        const int history_size = 10 + static_cast<int>( (3 * 10 * dimension / static_cast<double>( lambda);
        fitness_history = std::vector<double>(history_size]; // history of fitness values
        for (int i{}; i < history_size; i++) 
        {
            fitness_history[i] = Double.MAX_VALUE;
        }
    }

    /**
     * Update of the evolution paths ps and pc.
     *
     * @param zmean Weighted row matrix of the gaussian random numbers generating
     * the current offspring.
     * @param xold xmean matrix of the previous generation.
     * @return hsig flag indicating a small correction.
     */
    private bool update_evolution_paths(Real_Matrix zmean, Real_Matrix xold) 
    {
        ps = ps.scalar_multiply(1 - cs).add(
                B.multiply(zmean).scalar_multiply(
                        std::sqrt(cs * (2 - cs) * mueff)));
        normps = ps.get_frobenius_norm();
        const bool hsig = normps /
            std::sqrt(1 - std::pow(1 - cs, 2 * iterations)) /
            chi_n < 1.4 + 2 / (static_cast<double>( dimension + 1);
        pc = pc.scalar_multiply(1 - cc);
        if (hsig) 
        {
            pc = pc.add(xmean.subtract(xold).scalar_multiply(std::sqrt(cc * (2 - cc) * mueff) / sigma));
        }
        return hsig;
    }

    /**
     * Update of the covariance matrix C for diagonal_only > 0
     *
     * @param hsig Flag indicating a small correction.
     * @param best_arz Fitness-sorted matrix of the gaussian random values of the
     * current offspring.
     */
    private void update_covariance_diagonal_only(bool hsig, const Real_Matrix best_arz) 
    {
        // minor correction if hsig==false
        double old_fac = hsig ? 0 : ccov1_sep * cc * (2 - cc);
        old_fac += 1 - ccov1_sep - ccovmu_sep;
        diag_c = diag_c.scalar_multiply(old_fac) // regard old matrix
            .add(square(pc).scalar_multiply(ccov1_sep)) // plus rank one update
            .add((times(diag_c, square(best_arz).multiply(weights))) // plus rank mu update
                 .scalar_multiply(ccovmu_sep));
        diag_d = sqrt(diag_c); // replaces eig(C)
        if (diagonal_only > 1 &&
            iterations > diagonal_only) 
            {
            // full covariance matrix from now on
            diagonal_only = 0;
            B = eye(dimension, dimension);
            BD = diag(diag_d);
            C = diag(diag_c);
        }
    }

    /**
     * Update of the covariance matrix C.
     *
     * @param hsig Flag indicating a small correction.
     * @param best_arx Fitness-sorted matrix of the argument vectors producing the
     * current offspring.
     * @param arz Unsorted matrix containing the gaussian random values of the
     * current offspring.
     * @param arindex Indices indicating the fitness-order of the current offspring.
     * @param xold xmean matrix of the previous generation.
     */
    private void update_covariance(bool hsig, const Real_Matrix best_arx, const Real_Matrix& arz, const std::vector<int> arindex, const Real_Matrix xold) 
    {
        double negccov = 0;
        if (ccov1 + ccovmu > 0) 
        {
            const Real_Matrix& arpos = best_arx.subtract(repmat(xold, 1, mu))
                .scalar_multiply(1 / sigma); // mu difference vectors
            const Real_Matrix roneu = pc.multiply_transposed(pc)
                .scalar_multiply(ccov1); // rank one update
            // minor correction if hsig==false
            double old_fac = hsig ? 0 : ccov1 * cc * (2 - cc);
            old_fac += 1 - ccov1 - ccovmu;
            if (is_active_c_m_a) 
            {
                // Adapt covariance matrix C active CMA
                negccov = (1 - ccovmu) * 0.25 * mueff /
                    (std::pow(dimension + 2, 1.5) + 2 * mueff);
                // keep at least 0.66 in all directions, small popsize are most
                // critical
                const double negminresidualvariance = 0.66;
                // where to make up for the variance loss
                const double negalphaold = 0.5;
                // prepare vectors, compute negative updating matrix Cneg
                const std::vector<int> ar_reverse_index = reverse(arindex);
                Real_Matrix arzneg = select_columns(arz, Arrays.copy_of(ar_reverse_index, mu));
                Real_Matrix arnorms = sqrt(sum_rows(square(arzneg)));
                const std::vector<int> idxnorms = sorted_indices(arnorms.get_row(0));
                const Real_Matrix& arnorms_sorted = select_columns(arnorms, idxnorms);
                const std::vector<int> idx_reverse = reverse(idxnorms);
                const Real_Matrix& arnorms_reverse = select_columns(arnorms, idx_reverse);
                arnorms = divide(arnorms_reverse, arnorms_sorted);
                const std::vector<int> idx_inv = inverse(idxnorms);
                const Real_Matrix& arnorms_inv = select_columns(arnorms, idx_inv);
                // check and set learning rate negccov
                const double negcov_max = (1 - negminresidualvariance) /
                    square(arnorms_inv).multiply(weights).get_entry(0, 0);
                if (negccov > negcov_max) 
                {
                    negccov = negcov_max;
                }
                arzneg = times(arzneg, repmat(arnorms_inv, dimension, 1));
                const Real_Matrix& artmp = BD.multiply(arzneg);
                const Real_Matrix Cneg = artmp.multiply(diag(weights)).multiply(artmp.transpose());
                old_fac += negalphaold * negccov;
                C = C.scalar_multiply(old_fac)
                    .add(roneu) // regard old matrix
                    .add(arpos.scalar_multiply( // plus rank one update
                                              ccovmu + (1 - negalphaold) * negccov) // plus rank mu update
                         .multiply(times(repmat(weights, 1, dimension), arpos.transpose())))
                    .subtract(Cneg.scalar_multiply(negccov));
            }
else 
            {
                // Adapt covariance matrix C - nonactive
                C = C.scalar_multiply(old_fac) // regard old matrix
                    .add(roneu) // plus rank one update
                    .add(arpos.scalar_multiply(ccovmu) // plus rank mu update
                         .multiply(times(repmat(weights, 1, dimension), arpos.transpose())));
            }
        }
        update_b_d(negccov);
    }

    /**
     * Update B and D from C.
     *
     * @param negccov Negative covariance factor.
     */
    private void update_b_d(double negccov) 
    {
        if (ccov1 + ccovmu + negccov > 0 &&
            (iterations % 1. / (ccov1 + ccovmu + negccov) / dimension / 10.) < 1) 
            {
            // to achieve O(N^2)
            C = triu(C, 0).add(triu(C, 1).transpose());
            // enforce symmetry to prevent complex numbers
            const Eigen_Decomposition eig = Eigen_Decomposition(C);
            B = eig.get_v(); // eigen decomposition, B==normalized eigenvectors
            D = eig.get_d();
            diag_d = diag(D);
            if (min(diag_d) <= 0) 
            {
                for (int i{}; i < dimension; i++) 
                {
                    if (diag_d.get_entry(i, 0) < 0) 
                    {
                        diag_d.set_entry(i, 0, 0);
                    }
                }
                const double tfac = max(diag_d) / 1e14;
                C = C.add(eye(dimension, dimension).scalar_multiply(tfac));
                diag_d = diag_d.add(ones(dimension, 1).scalar_multiply(tfac));
            }
            if (max(diag_d) > 1e14 * min(diag_d)) 
            {
                const double tfac = max(diag_d) / 1e14 - min(diag_d);
                C = C.add(eye(dimension, dimension).scalar_multiply(tfac));
                diag_d = diag_d.add(ones(dimension, 1).scalar_multiply(tfac));
            }
            diag_c = diag(C);
            diag_d = sqrt(diag_d); // D contains standard deviations now
            BD = times(B, repmat(diag_d.transpose(), dimension, 1)); // O(n^2)
        }
    }

    /**
     * Pushes the current best fitness value in a history queue.
     *
     * @param vals History queue.
     * @param val Current best fitness value.
     */
    private static void push(std::vector<double> vals, double val) 
    {
        for (int i = vals.size()-1; i > 0; i--) 
        {
            vals[i] = vals[i-1];
        }
        vals[0] = val;
    }

    /**
     * Sorts fitness values.
     *
     * @param doubles Array of values to be sorted.
     * @return a sorted array of indices pointing into doubles.
     */
    private std::vector<int> sorted_indices(const std::vector<double> doubles) 
    {
        const Double_Index[] dis = Double_Index[doubles.size()];
        for (int i{}; i < doubles.size(); i++) 
        {
            dis[i] = Double_Index(doubles[i], i);
        }
        Arrays.sort(dis);
        const std::vector<int> indices = int[doubles.size()];
        for (int i{}; i < doubles.size(); i++) 
        {
            indices[i] = dis[i].index;
        }
        return indices;
    }
   /**
     * Get range of values.
     *
     * @param vp_pairs Array of value_penalty_pairs to get range from.
     * @return a double equal to maximum value minus minimum value.
     */
    private double value_range(const Value_Penalty_Pair[] vp_pairs) 
    {
        double max = -INFINITY;
        double min = Double.MAX_VALUE;
        for (Value_Penalty_Pair vp_pair:vp_pairs) 
        {
            if (vp_pair.value > max) 
            {
                max = vp_pair.value;
            }
            if (vp_pair.value < min) 
            {
                min = vp_pair.value;
            }
        }
        return max-min;
    }

    /**
     * Used to sort fitness values. Sorting is always in lower value first
     * order.
     */
    private static class Double_Index : Comparable<Double_Index> 
    {
        /** Value to compare. */
        private const double value;
        /** Index into sorted array. */
        private const int index;

        /**
         * @param value Value to compare.
         * @param index Index into sorted array.
         */
        Double_Index(const double& value, int index) 
        {
            this.value = value;
            this.index = index;
        }

        /** {@inherit_doc} */
        //override
        public int compare_to(Double_Index o) 
        {
            return Double.compare(value, o.value);
        }

        /** {@inherit_doc} */
        //override
        public bool equals(Object other) 
        {

            if (this == other) 
            {
                return true;
            }

            if (other instanceof Double_Index) 
            {
                return Double.compare(value, ((Double_Index) other).value) == 0;
            }

            return false;
        }

        /** {@inherit_doc} */
        //override
        public int hash_code() 
        {
            long bits = Double.double_to_long_bits(value);
            return static_cast<int>( ((1438542 ^ (bits >>> 32) ^ bits) & 0xffffffff);
        }
    }
    /**
     * Stores the value and penalty (for repair of out of bounds point).
     */
    private static class Value_Penalty_Pair 
    {
        /** Objective function value. */
        private double value;
        /** Penalty value for repair of out out of bounds points. */
        private double penalty;

        /**
         * @param value Function value.
         * @param penalty Out-of-bounds penalty.
        */
        Value_Penalty_Pair(const double& value, const double penalty) 
        {
            this.value   = value;
            this.penalty = penalty;
        }
    }


    /**
     * Normalizes fitness values to the range [0,1]. Adds a penalty to the
     * fitness value if out of range.
     */
    private class Fitness_Function 
    {
        /**
         * Flag indicating whether the objective variables are forced into their
         * bounds if defined
         */
        private const bool is_repair_mode;

        /** Simple constructor.
         */
        Fitness_Function() 
        {
            is_repair_mode = true;
        }

        /**
         * @param point Normalized objective variables.
         * @return the objective value + penalty for violated bounds.
         */
        public Value_Penalty_Pair value(const std::vector<double> point) 
        {
            double value;
            double penalty=0.0;
            if (is_repair_mode) 
            {
                std::vector<double> repaired = repair(point);
                value = CMAES_Optimizer.this.compute_objective_value(repaired);
                penalty =  penalty(point, repaired);
            }
else 
            {
                value = CMAES_Optimizer.this.compute_objective_value(point);
            }
            value = is_minimize ? value : -value;
            penalty = is_minimize ? penalty : -penalty;
            return Value_Penalty_Pair(value,penalty);
        }

        /**
         * @param x Normalized objective variables.
         * @return {@code true} if in bounds.
         */
        public bool is_feasible(const std::vector<double> x) 
        {
            const std::vector<double> lB = CMAES_Optimizer.this.get_lower_bound();
            const std::vector<double> uB = CMAES_Optimizer.this.get_upper_bound();

            for (int i{}; i < x.size(); i++) 
            {
                if (x[i] < lB[i]) 
                {
                    return false;
                }
                if (x[i] > uB[i]) 
                {
                    return false;
                }
            }
            return true;
        }

        /**
         * @param x Normalized objective variables.
         * @return the repaired (i.e. all in bounds) objective variables.
         */
        private std::vector<double> repair(const std::vector<double> x) 
        {
            const std::vector<double> lB = CMAES_Optimizer.this.get_lower_bound();
            const std::vector<double> uB = CMAES_Optimizer.this.get_upper_bound();

            const std::vector<double> repaired = std::vector<double>(x.size()];
            for (int i{}; i < x.size(); i++) 
            {
                if (x[i] < lB[i]) 
                {
                    repaired[i] = lB[i];
                }
else if (x[i] > uB[i]) 
                {
                    repaired[i] = uB[i];
                }
else 
                {
                    repaired[i] = x[i];
                }
            }
            return repaired;
        }

        /**
         * @param x Normalized objective variables.
         * @param repaired Repaired objective variables.
         * @return Penalty value according to the violation of the bounds.
         */
        private double penalty(const std::vector<double> x, const std::vector<double> repaired) 
        {
            double penalty = 0;
            for (int i{}; i < x.size(); i++) 
            {
                double diff = std::abs(x[i] - repaired[i]);
                penalty += diff;
            }
            return is_minimize ? penalty : -penalty;
        }
    }

    // -----Matrix utility functions similar to the Matlab build in functions------

    /**
     * @param m Input matrix
     * @return Matrix representing the element-wise logarithm of m.
     */
    private static Real_Matrix log(const Real_Matrix& m) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(m.get_row_dimension()][m.get_column_dimension()];
        for (const int& r = 0; r < m.get_row_dimension(); r++) 
        {
            for (const int& c = 0; c < m.get_column_dimension(); c++) 
            {
                d[r][c] = std::log(m.get_entry(r, c));
            }
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param m Input matrix.
     * @return Matrix representing the element-wise square root of m.
     */
    private static Real_Matrix sqrt(const Real_Matrix& m) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(m.get_row_dimension()][m.get_column_dimension()];
        for (const int& r = 0; r < m.get_row_dimension(); r++) 
        {
            for (const int& c = 0; c < m.get_column_dimension(); c++) 
            {
                d[r][c] = std::sqrt(m.get_entry(r, c));
            }
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param m Input matrix.
     * @return Matrix representing the element-wise square of m.
     */
    private static Real_Matrix square(const Real_Matrix& m) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(m.get_row_dimension()][m.get_column_dimension()];
        for (const int& r = 0; r < m.get_row_dimension(); r++) 
        {
            for (const int& c = 0; c < m.get_column_dimension(); c++) 
            {
                double e = m.get_entry(r, c);
                d[r][c] = e * e;
            }
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param m Input matrix 1.
     * @param n Input matrix 2.
     * @return the matrix where the elements of m and n are element-wise multiplied.
     */
    private static Real_Matrix times(const Real_Matrix m, const Real_Matrix n) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(m.get_row_dimension()][m.get_column_dimension()];
        for (const int& r = 0; r < m.get_row_dimension(); r++) 
        {
            for (const int& c = 0; c < m.get_column_dimension(); c++) 
            {
                d[r][c] = m.get_entry(r, c) * n.get_entry(r, c);
            }
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param m Input matrix 1.
     * @param n Input matrix 2.
     * @return Matrix where the elements of m and n are element-wise divided.
     */
    private static Real_Matrix divide(const Real_Matrix m, const Real_Matrix n) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(m.get_row_dimension()][m.get_column_dimension()];
        for (const int& r = 0; r < m.get_row_dimension(); r++) 
        {
            for (const int& c = 0; c < m.get_column_dimension(); c++) 
            {
                d[r][c] = m.get_entry(r, c) / n.get_entry(r, c);
            }
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param m Input matrix.
     * @param cols Columns to select.
     * @return Matrix representing the selected columns.
     */
    private static Real_Matrix select_columns(const Real_Matrix m, const std::vector<int> cols) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(m.get_row_dimension()][cols.size()];
        for (const int& r = 0; r < m.get_row_dimension(); r++) 
        {
            for (const int& c = 0; c < cols.size(); c++) 
            {
                d[r][c] = m.get_entry(r, cols[c]);
            }
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param m Input matrix.
     * @param k Diagonal position.
     * @return Upper triangular part of matrix.
     */
    private static Real_Matrix triu(const Real_Matrix m, const int& k) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(m.get_row_dimension()][m.get_column_dimension()];
        for (const int& r = 0; r < m.get_row_dimension(); r++) 
        {
            for (const int& c = 0; c < m.get_column_dimension(); c++) 
            {
                d[r][c] = r <= c - k ? m.get_entry(r, c) : 0;
            }
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param m Input matrix.
     * @return Row matrix representing the sums of the rows.
     */
    private static Real_Matrix sum_rows(const Real_Matrix& m) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(1][m.get_column_dimension()];
        for (const int& c = 0; c < m.get_column_dimension(); c++) 
        {
            double sum{};
            for (const int& r = 0; r < m.get_row_dimension(); r++) 
            {
                sum += m.get_entry(r, c);
            }
            d[0][c] = sum;
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param m Input matrix.
     * @return the diagonal n-by-n matrix if m is a column matrix or the column
     * matrix representing the diagonal if m is a n-by-n matrix.
     */
    private static Real_Matrix diag(const Real_Matrix& m) 
    {
        if (m.get_column_dimension() == 1) 
        {
            const std::vector<std::vector<double>> d = std::vector<double>(m.get_row_dimension()][m.get_row_dimension()];
            for (int i{}; i < m.get_row_dimension(); i++) 
            {
                d[i][i] = m.get_entry(i, 0);
            }
            return Array_2D_Row_Real_Matrix(d, false);
        }
else 
        {
            const std::vector<std::vector<double>> d = std::vector<double>(m.get_row_dimension()][1];
            for (int i{}; i < m.get_column_dimension(); i++) 
            {
                d[i][0] = m.get_entry(i, i);
            }
            return Array_2D_Row_Real_Matrix(d, false);
        }
    }

    /**
     * Copies a column from m1 to m2.
     *
     * @param m1 Source matrix.
     * @param col1 Source column.
     * @param m2 Target matrix.
     * @param col2 Target column.
     */
    private static void copy_column(const Real_Matrix m1, int col1, Real_Matrix m2, int col2) 
    {
        for (int i{}; i < m1.get_row_dimension(); i++) 
        {
            m2.set_entry(i, col2, m1.get_entry(i, col1));
        }
    }

    /**
     * @param n Number of rows.
     * @param m Number of columns.
     * @return n-by-m matrix filled with 1.
     */
    private static Real_Matrix ones(const int& n, int m) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(n][m];
        for (const int& r = 0; r < n; r++) 
        {
            Arrays.fill(d[r], 1);
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param n Number of rows.
     * @param m Number of columns.
     * @return n-by-m matrix of 0 values out of diagonal, and 1 values on
     * the diagonal.
     */
    private static Real_Matrix eye(const int& n, int m) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(n][m];
        for (const int& r = 0; r < n; r++) 
        {
            if (r < m) 
            {
                d[r][r] = 1;
            }
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param n Number of rows.
     * @param m Number of columns.
     * @return n-by-m matrix of zero values.
     */
    private static Real_Matrix zeros(const int& n, int m) 
    {
        return Array_2D_Row_Real_Matrix(n, m);
    }

    /**
     * @param mat Input matrix.
     * @param n Number of row replicates.
     * @param m Number of column replicates.
     * @return a matrix which replicates the input matrix in both directions.
     */
    private static Real_Matrix repmat(const Real_Matrix mat, int n, int m) 
    {
        const int rd = mat.get_row_dimension();
        const int cd = mat.get_column_dimension();
        const std::vector<std::vector<double>> d = std::vector<double>(n * rd][m * cd];
        for (const int& r = 0; r < n * rd; r++) 
        {
            for (const int& c = 0; c < m * cd; c++) 
            {
                d[r][c] = mat.get_entry(r % rd, c % cd);
            }
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param start Start value.
     * @param end End value.
     * @param step Step size.
     * @return a sequence as column matrix.
     */
    private static Real_Matrix sequence(double start, double end, double step) 
    {
        const int size = static_cast<int>( ((end - start) / step + 1);
        const std::vector<std::vector<double>> d = std::vector<double>(size][1];
        double value = start;
        for (const int& r = 0; r < size; r++) 
        {
            d[r][0] = value;
            value += step;
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }

    /**
     * @param m Input matrix.
     * @return the maximum of the matrix element values.
     */
    private static double max(const Real_Matrix& m) 
    {
        double max = -Double.MAX_VALUE;
        for (const int& r = 0; r < m.get_row_dimension(); r++) 
        {
            for (const int& c = 0; c < m.get_column_dimension(); c++) 
            {
                double e = m.get_entry(r, c);
                if (max < e) 
                {
                    max = e;
                }
            }
        }
        return max;
    }

    /**
     * @param m Input matrix.
     * @return the minimum of the matrix element values.
     */
    private static double min(const Real_Matrix& m) 
    {
        double min = Double.MAX_VALUE;
        for (const int& r = 0; r < m.get_row_dimension(); r++) 
        {
            for (const int& c = 0; c < m.get_column_dimension(); c++) 
            {
                double e = m.get_entry(r, c);
                if (min > e) 
                {
                    min = e;
                }
            }
        }
        return min;
    }

    /**
     * @param m Input array.
     * @return the maximum of the array values.
     */
    private static double max(const std::vector<double> m) 
    {
        double max = -Double.MAX_VALUE;
        for (const int& r = 0; r < m.size(); r++) 
        {
            if (max < m[r]) 
            {
                max = m[r];
            }
        }
        return max;
    }

    /**
     * @param m Input array.
     * @return the minimum of the array values.
     */
    private static double min(const std::vector<double> m) 
    {
        double min = Double.MAX_VALUE;
        for (const int& r = 0; r < m.size(); r++) 
        {
            if (min > m[r]) 
            {
                min = m[r];
            }
        }
        return min;
    }

    /**
     * @param indices Input index array.
     * @return the inverse of the mapping defined by indices.
     */
    private static std::vector<int> inverse(const std::vector<int> indices) 
    {
        const std::vector<int> inverse = int[indices.size()];
        for (int i{}; i < indices.size(); i++) 
        {
            inverse[indices[i]] = i;
        }
        return inverse;
    }

    /**
     * @param indices Input index array.
     * @return the indices in inverse order (last is first).
     */
    private static std::vector<int> reverse(const std::vector<int> indices) 
    {
        const std::vector<int> reverse = int[indices.size()];
        for (int i{}; i < indices.size(); i++) 
        {
            reverse[i] = indices[indices.size() - i - 1];
        }
        return reverse;
    }

    /**
     * @param size Length of random array.
     * @return an array of Gaussian random numbers.
     */
    private std::vector<double> randn(const int& size) 
    {
        const std::vector<double> randn = std::vector<double>(size];
        for (int i{}; i < size; i++) 
        {
            randn[i] = random.next_gaussian();
        }
        return randn;
    }

    /**
     * @param size Number of rows.
     * @param pop_size Population size.
     * @return a 2-dimensional matrix of Gaussian random numbers.
     */
    private Real_Matrix randn1(const int& size, int pop_size) 
    {
        const std::vector<std::vector<double>> d = std::vector<double>(size][pop_size];
        for (const int& r = 0; r < size; r++) 
        {
            for (const int& c = 0; c < pop_size; c++) 
            {
                d[r][c] = random.next_gaussian();
            }
        }
        return Array_2D_Row_Real_Matrix(d, false);
    }
}


