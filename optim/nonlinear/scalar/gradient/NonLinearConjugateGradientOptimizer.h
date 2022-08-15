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

//package org.hipparchus.optim.nonlinear.scalar.gradient;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.optim.Point_valuePair;
//import org.hipparchus.optim.nonlinear.scalar.Goal_Type;
//import org.hipparchus.optim.nonlinear.scalar.GradientMultivariate_Optimizer;
//import org.hipparchus.optim.nonlinear.scalar.Line_Search;


/**
 * Non-linear conjugate gradient optimizer.
 * <br/>
 * This class supports both the Fletcher-_Reeves and the Polak-Ribière
 * update formulas for the conjugate search directions.
 * It also supports optional preconditioning.
 * <br/>
 * Constraints are not supported: the call to
 * {@link #optimize(Optimization_data[]) optimize} will throw
 * {@link Math_Runtime_Exception} if bounds are passed to it.
 *
 */
class NonLinearConjugate_GradientOptimizer
    extends GradientMultivariate_Optimizer 
    {
    /** Update formula for the beta parameter. */
    private const Formula update_formula;
    /** Preconditioner (may be NULL). */
    private const Preconditioner preconditioner;
    /** Line search algorithm. */
    private const Line_Search line;

    /**
     * Available choices of update formulas for the updating the parameter
     * that is used to compute the successive conjugate search directions.
     * For non-linear conjugate gradients, there are
     * two formulas:
     * <ul>
     *   <li>Fletcher-_Reeves formula</li>
     *   <li>Polak-Ribière formula</li>
     * </ul>
     *
     * On the one hand, the Fletcher-_Reeves formula is guaranteed to converge
     * if the start point is close enough of the optimum whether the
     * Polak-Ribière formula may not converge in rare cases. On the
     * other hand, the Polak-Ribière formula is often faster when it
     * does converge. Polak-Ribière is often used.
     *
     */
    public enum Formula 
    {
        /** Fletcher-_Reeves formula. */
        FLETCHER_REEVES, /** Polak-Ribière formula. */
        POLAK_RIBIERE
    }

    /**
     * Constructor with default tolerances for the line search (1e-8) and
     * {@link Identity_Preconditioner preconditioner}.
     *
     * @param update_formula formula to use for updating the &beta; parameter, * must be one of {@link Formula#FLETCHER_REEVES} or
     * {@link Formula#POLAK_RIBIERE}.
     * @param checker Convergence checker.
     */
    public NonLinearConjugate_GradientOptimizer(const Formula update_formula, Convergence_Checker<Point_valuePair> checker) 
    {
        this(update_formula, checker, 1e-8, 1e-8, 1e-8, Identity_Preconditioner());
    }

    /**
     * Constructor with default {@link Identity_Preconditioner preconditioner}.
     *
     * @param update_formula formula to use for updating the &beta; parameter, * must be one of {@link Formula#FLETCHER_REEVES} or
     * {@link Formula#POLAK_RIBIERE}.
     * @param checker Convergence checker.
     * @param relative_tolerance Relative threshold for line search.
     * @param absolute_tolerance Absolute threshold for line search.
     * @param initial_bracketing_range Extent of the initial interval used to
     * find an interval that brackets the optimum in order to perform the
     * line search.
     *
     * @see Line_Search#Line_Search(org.hipparchus.optim.nonlinear.scalar.Multivariate_Optimizer,double,double,double)
     */
    public NonLinearConjugate_GradientOptimizer(const Formula update_formula, Convergence_Checker<Point_valuePair> checker, double relative_tolerance, double absolute_tolerance, double initial_bracketing_range) 
    {
        this(update_formula, checker, relative_tolerance, absolute_tolerance, initial_bracketing_range, Identity_Preconditioner());
    }

    /**
     * @param update_formula formula to use for updating the &beta; parameter, * must be one of {@link Formula#FLETCHER_REEVES} or
     * {@link Formula#POLAK_RIBIERE}.
     * @param checker Convergence checker.
     * @param preconditioner Preconditioner.
     * @param relative_tolerance Relative threshold for line search.
     * @param absolute_tolerance Absolute threshold for line search.
     * @param initial_bracketing_range Extent of the initial interval used to
     * find an interval that brackets the optimum in order to perform the
     * line search.
     *
     * @see Line_Search#Line_Search(org.hipparchus.optim.nonlinear.scalar.Multivariate_Optimizer,double,double,double)
     */
    public NonLinearConjugate_GradientOptimizer(const Formula update_formula, Convergence_Checker<Point_valuePair> checker, double relative_tolerance, double absolute_tolerance, double initial_bracketing_range, const Preconditioner preconditioner) 
    {
        super(checker);

        this.update_formula = update_formula;
        this.preconditioner = preconditioner;
        line = Line_Search(this, relative_tolerance, absolute_tolerance, initial_bracketing_range);
    }

    /**
     * {@inherit_doc}
     */
    //override
    public Point_valuePair optimize(Optimization_data... opt_data)
        Math_Illegal_State_Exception 
        {
        // Set up base class and perform computation.
        return super.optimize(opt_data);
    }

    /** {@inherit_doc} */
    //override
    protected Point_valuePair do_optimize() 
    {
        const Convergence_Checker<Point_valuePair> checker = get_convergence_checker();
        const std::vector<double> point = get_start_point();
        const Goal_Type goal = get_goal_type();
        const int n = point.size();
        std::vector<double> r = compute_objective_gradient(point);
        if (goal == Goal_Type.MINIMIZE) 
        {
            for (int i{}; i < n; i++) 
            {
                r[i] = -r[i];
            }
        }

        // Initial search direction.
        std::vector<double> steepest_descent = preconditioner.precondition(point, r);
        std::vector<double> search_direction = steepest_descent.clone();

        double delta = 0;
        for (int i{}; i < n; ++i) 
        {
            delta += r[i] * search_direction[i];
        }

        Point_valuePair current = NULL;
        while (true) 
        {
            increment_iteration_count();

            const double objective = compute_objective_value(point);
            Point_valuePair previous = current;
            current = Point_valuePair(point, objective);
            if (previous != NULL && checker.converged(get_iterations(), previous, current)) 
            {
                // We have found an optimum.
                return current;
            }

            const double step = line.search(point, search_direction).get_point();

            // Validate point.
            for (int i{}; i < point.size(); ++i) 
            {
                point[i] += step * search_direction[i];
            }

            r = compute_objective_gradient(point);
            if (goal == Goal_Type.MINIMIZE) 
            {
                for (int i{}; i < n; ++i) 
                {
                    r[i] = -r[i];
                }
            }

            // Compute beta.
            const double delta_old = delta;
            const std::vector<double> new_steepest_descent = preconditioner.precondition(point, r);
            delta = 0;
            for (int i{}; i < n; ++i) 
            {
                delta += r[i] * new_steepest_descent[i];
            }

            const double beta;
            switch (update_formula) 
            {
            case FLETCHER_REEVES:
                beta = delta / delta_old;
                break;
            case POLAK_RIBIERE:
                double delta_mid = 0;
                for (int i{}; i < r.size(); ++i) 
                {
                    delta_mid += r[i] * steepest_descent[i];
                }
                beta = (delta - delta_mid) / delta_old;
                break;
            default:
                // Should never happen.
                throw Math_Runtime_Exception.create_internal_error();
            }
            steepest_descent = new_steepest_descent;

            // Compute conjugate search direction.
            if (get_iterations() % n == 0 ||
                beta < 0) 
                {
                // Break conjugation: reset search direction.
                search_direction = steepest_descent.clone();
            }
else 
            {
                // Compute conjugate search direction.
                for (int i{}; i < n; ++i) 
                {
                    search_direction[i] = steepest_descent[i] + beta * search_direction[i];
                }
            }
        }
    }

    /**
     * {@inherit_doc}
     */
    //override
    protected void parse_optimization_data(Optimization_data... opt_data) 
    {
        // Allow base class to register its own data.
        super.parse_optimization_data(opt_data);

        check_parameters();
    }

    /** Default identity preconditioner. */
    public static class Identity_Preconditioner : Preconditioner 
    {
        /** {@inherit_doc} */
        //override
        public std::vector<double> precondition(std::vector<double> variables, std::vector<double> r) 
        {
            return r.clone();
        }
    }

    // Class is not used anymore (cf. MATH-1092). However, it might
    // be interesting to create a class similar to "Line_Search", but
    // that will take advantage that the model's gradient is available.
//     /**
//      * Internal class for line search.
//      * <p>
//      * The function represented by this class is the dot product of
//      * the objective function gradient and the search direction. Its
//      * value is zero when the gradient is orthogonal to the search
//      * direction, i.e. when the objective function value is a local
//      * extremum along the search direction.
//      * </p>
//      */
//     private class Line_SearchFunction : Univariate_Function 
{
//         /** Current point. */
//         private const std::vector<double> current_point;
//         /** Search direction. */
//         private const std::vector<double> search_direction;

//         /**
//          * @param point Current point.
//          * @param direction Search direction.
//          */
//         public Line_SearchFunction(std::vector<double> point, //                                   std::vector<double> direction) 
{
//             current_point = point.clone();
//             search_direction = direction.clone();
//         }

//         /** {@inherit_doc} */
//         public double value(double x) 
{
//             // current point in the search direction
//             const std::vector<double> shifted_point = current_point.clone();
//             for (int i{}; i < shifted_point.size(); ++i) 
{
//                 shifted_point[i] += x * search_direction[i];
//             }

//             // gradient of the objective function
//             const std::vector<double> gradient = compute_objective_gradient(shifted_point);

//             // dot product with the search direction
//             double dot_product = 0;
//             for (int i{}; i < gradient.size(); ++i) 
{
//                 dot_product += gradient[i] * search_direction[i];
//             }

//             return dot_product;
//         }
//     }

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


