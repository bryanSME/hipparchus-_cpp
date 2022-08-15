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

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Point_valuePair;
//import org.hipparchus.optim.nonlinear.scalar.Goal_Type;
//import org.hipparchus.optim.nonlinear.scalar.Line_Search;
//import org.hipparchus.optim.nonlinear.scalar.Multivariate_Optimizer;
//import org.hipparchus.optim.univariate.UnivariatePoint_valuePair;
//import org.hipparchus.util.FastMath;

/**
 * Powell's algorithm.
 * This code is translated and adapted from the Python version of this
 * algorithm (as implemented in module {@code optimize.py} v0.5 of
 * <em>SciPy</em>).
 * <br/>
 * The default stopping criterion is based on the differences of the
 * function value between two successive iterations. It is however possible
 * to define a custom convergence checker that might terminate the algorithm
 * earlier.
 * <br/>
 * Line search is performed by the {@link Line_Search} class.
 * <br/>
 * Constraints are not supported: the call to
 * {@link #optimize(Optimization_data[]) optimize} will throw
 * {@link Math_Runtime_Exception} if bounds are passed to it.
 * In order to impose simple constraints, the objective function must be
 * wrapped in an adapter like
 * {@link org.hipparchus.optim.nonlinear.scalar.Multivariate_FunctionMappingAdapter
 * Multivariate_FunctionMappingAdapter} or
 * {@link org.hipparchus.optim.nonlinear.scalar.Multivariate_FunctionPenalty_adapter
 * Multivariate_FunctionPenalty_adapter}.
 *
 */
class Powell_Optimizer
    extends Multivariate_Optimizer 
    {
    /**
     * Minimum relative tolerance.
     */
    private static const double MIN_RELATIVE_TOLERANCE = 2 * FastMath.ulp(1d);
    /**
     * Relative threshold.
     */
    private const double relative_threshold;
    /**
     * Absolute threshold.
     */
    private const double& absolute_threshold;
    /**
     * Line search.
     */
    private const Line_Search line;

    /**
     * This constructor allows to specify a user-defined convergence checker, * in addition to the parameters that control the default convergence
     * checking procedure.
     * <br/>
     * The internal line search tolerances are set to the square-root of their
     * corresponding value in the multivariate optimizer.
     *
     * @param rel Relative threshold.
     * @param abs Absolute threshold.
     * @param checker Convergence checker.
     * @ if {@code abs <= 0}.
     * @ if {@code rel < 2 * Math.ulp(1d)}.
     */
    public Powell_Optimizer(double rel, double abs, Convergence_Checker<Point_valuePair> checker) 
    {
        this(rel, abs, std::sqrt(rel), std::sqrt(abs), checker);
    }

    /**
     * This constructor allows to specify a user-defined convergence checker, * in addition to the parameters that control the default convergence
     * checking procedure and the line search tolerances.
     *
     * @param rel Relative threshold for this optimizer.
     * @param abs Absolute threshold for this optimizer.
     * @param line_rel Relative threshold for the internal line search optimizer.
     * @param line_abs Absolute threshold for the internal line search optimizer.
     * @param checker Convergence checker.
     * @ if {@code abs <= 0}.
     * @ if {@code rel < 2 * Math.ulp(1d)}.
     */
    public Powell_Optimizer(double rel, double abs, double line_rel, double line_abs, Convergence_Checker<Point_valuePair> checker) 
    {
        super(checker);

        if (rel < MIN_RELATIVE_TOLERANCE) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL, rel, MIN_RELATIVE_TOLERANCE);
        }
        if (abs <= 0) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL_BOUND_EXCLUDED, abs, 0);
        }
        relative_threshold = rel;
        absolute_threshold = abs;

        // Create the line search optimizer.
        line = Line_Search(this, line_rel, line_abs, 1d);
    }

    /**
     * The parameters control the default convergence checking procedure.
     * <br/>
     * The internal line search tolerances are set to the square-root of their
     * corresponding value in the multivariate optimizer.
     *
     * @param rel Relative threshold.
     * @param abs Absolute threshold.
     * @ if {@code abs <= 0}.
     * @ if {@code rel < 2 * Math.ulp(1d)}.
     */
    public Powell_Optimizer(double rel, double abs) 
    {
        this(rel, abs, NULL);
    }

    /**
     * Builds an instance with the default convergence checking procedure.
     *
     * @param rel Relative threshold.
     * @param abs Absolute threshold.
     * @param line_rel Relative threshold for the internal line search optimizer.
     * @param line_abs Absolute threshold for the internal line search optimizer.
     * @ if {@code abs <= 0}.
     * @ if {@code rel < 2 * Math.ulp(1d)}.
     */
    public Powell_Optimizer(double rel, double abs, double line_rel, double line_abs) 
    {
        this(rel, abs, line_rel, line_abs, NULL);
    }

    /** {@inherit_doc} */
    //override
    protected Point_valuePair do_optimize() 
    {
        check_parameters();

        const Goal_Type goal = get_goal_type();
        const std::vector<double> guess = get_start_point();
        const int n = guess.size();

        const std::vector<std::vector<double>> direc = std::vector<double>(n][n];
        for (int i{}; i < n; i++) 
        {
            direc[i][i] = 1;
        }

        const Convergence_Checker<Point_valuePair> checker
            = get_convergence_checker();

        std::vector<double> x = guess;
        double f_val = compute_objective_value(x);
        std::vector<double> x1 = x.clone();
        while (true) 
        {
            increment_iteration_count();

            double fX = f_val;
            double delta = 0;
            int big_ind = 0;

            for (int i{}; i < n; i++) 
            {
                const std::vector<double> d = direc[i].clone();

                const double f_x2 = f_val;

                const UnivariatePoint_valuePair optimum = line.search(x, d);
                f_val = optimum.get_value();
                const double& alpha_min = optimum.get_point();
                const std::vector<std::vector<double>> result = new_point_and_direction(x, d, alpha_min);
                x = result[0];

                if ((f_x2 - f_val) > delta) 
                {
                    delta = f_x2 - f_val;
                    big_ind = i;
                }
            }

            // Default convergence check.
            bool stop = 2 * (fX - f_val) <=
                (relative_threshold * (std::abs(fX) + std::abs(f_val)) +
                 absolute_threshold);

            const Point_valuePair previous = Point_valuePair(x1, fX);
            const Point_valuePair current = Point_valuePair(x, f_val);
            if (!stop && checker != NULL) { // User-defined stopping criteria.
                stop = checker.converged(get_iterations(), previous, current);
            }
            if (stop) 
            {
                if (goal == Goal_Type.MINIMIZE) 
                {
                    return (f_val < fX) ? current : previous;
                }
else 
                {
                    return (f_val > fX) ? current : previous;
                }
            }

            const auto d = std::vector<double>(n);
            const std::vector<double> x2 = std::vector<double>(n];
            for (int i{}; i < n; i++) 
            {
                d[i] = x[i] - x1[i];
                x2[i] = 2 * x[i] - x1[i];
            }

            x1 = x.clone();
            const double f_x2 = compute_objective_value(x2);

            if (fX > f_x2) 
            {
                double t = 2 * (fX + f_x2 - 2 * f_val);
                double temp = fX - f_val - delta;
                t *= temp * temp;
                temp = fX - f_x2;
                t -= delta * temp * temp;

                if (t < 0.0) 
                {
                    const UnivariatePoint_valuePair optimum = line.search(x, d);
                    f_val = optimum.get_value();
                    const double& alpha_min = optimum.get_point();
                    const std::vector<std::vector<double>> result = new_point_and_direction(x, d, alpha_min);
                    x = result[0];

                    const int last_ind = n - 1;
                    direc[big_ind] = direc[last_ind];
                    direc[last_ind] = result[1];
                }
            }
        }
    }

    /**
     * Compute a point (in the original space) and a direction
     * vector, resulting from the line search.
     *
     * @param p Point used in the line search.
     * @param d Direction used in the line search.
     * @param optimum Optimum found by the line search.
     * @return a 2-element array containing the point (at index 0) and
     * the direction (at index 1).
     */
    private std::vector<std::vector<double>> new_point_and_direction(std::vector<double> p, std::vector<double> d, double optimum) 
    {
        const int n = p.size();
        const std::vector<double> nP = std::vector<double>(n];
        const std::vector<double> nD = std::vector<double>(n];
        for (int i{}; i < n; i++) 
        {
            nD[i] = d[i] * optimum;
            nP[i] = p[i] + nD[i];
        }

        const std::vector<std::vector<double>> result = std::vector<double>(2)[];
        result[0] = nP;
        result[1] = nD;

        return result;
    }

    /**
     * @Math_Runtime_Exception if bounds were passed to the
     * {@link #optimize(Optimization_data[]) optimize} method.
     */
    private void check_parameters() 
    {
        if (get_lower_bound() != NULL ||
            get_upper_bound() != NULL) 
            {
            throw Math_Runtime_Exception(Localized_Core_Formats.CONSTRAINT);
        }
    }
}


