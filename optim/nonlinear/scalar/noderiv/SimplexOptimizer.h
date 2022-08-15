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

//import java.util.Comparator;

//import org.hipparchus.analysis.Multivariate_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.optim.Point_valuePair;
//import org.hipparchus.optim.Simple_Value_Checker;
//import org.hipparchus.optim.nonlinear.scalar.Goal_Type;
//import org.hipparchus.optim.nonlinear.scalar.Multivariate_Optimizer;

/**
 * This class : simplex-based direct search optimization.
 *
 * <p>
 *  Direct search methods only use objective function values, they do
 *  not need derivatives and don't either try to compute approximation
 *  of the derivatives. According to a 1996 paper by Margaret H. Wright
 *  (<a href="http://cm.bell-labs.com/cm/cs/doc/96/4-02.ps.gz">Direct
 *  Search Methods: Once Scorned, Now Respectable</a>), they are used
 *  when either the computation of the derivative is impossible (noisy
 *  functions, unpredictable discontinuities) or difficult (complexity, *  computation cost). In the first cases, rather than an optimum, a
 *  <em>not too bad</em> point is desired. In the latter cases, an
 *  optimum is desired but cannot be reasonably found. In all cases
 *  direct search methods can be useful.
 * </p>
 * <p>
 *  Simplex-based direct search methods are based on comparison of
 *  the objective function values at the vertices of a simplex (which is a
 *  set of n+1 points in dimension n) that is updated by the algorithms
 *  steps.
 * <p>
 * <p>
 *  The simplex update procedure ({@link Nelder_Mead_Simplex} or
 * {@link Multi_Directional_Simplex})  must be passed to the
 * {@code optimize} method.
 * </p>
 * <p>
 *  Each call to {@code optimize} will re-use the start configuration of
 *  the current simplex and move it such that its first vertex is at the
 *  provided start point of the optimization.
 *  If the {@code optimize} method is called to solve a different problem
 *  and the number of parameters change, the simplex must be re-initialized
 *  to one with the appropriate dimensions.
 * </p>
 * <p>
 *  Convergence is checked by providing the <em>worst</em> points of
 *  previous and current simplex to the convergence checker, not the best
 *  ones.
 * </p>
 * <p>
 *  This simplex optimizer implementation does not directly support constrained
 *  optimization with simple bounds; so, for such optimizations, either a more
 *  dedicated algorithm must be used like
 *  {@link CMAES_Optimizer} or {@link BOBYQA_Optimizer}, or the objective
 *  function must be wrapped in an adapter like
 *  {@link org.hipparchus.optim.nonlinear.scalar.Multivariate_FunctionMappingAdapter
 *  Multivariate_FunctionMappingAdapter} or
 *  {@link org.hipparchus.optim.nonlinear.scalar.Multivariate_FunctionPenalty_adapter
 *  Multivariate_FunctionPenalty_adapter}.
 *  <br/>
 *  The call to {@link #optimize(Optimization_data[]) optimize} will throw
 *  {@link Math_Runtime_Exception} if bounds are passed to it.
 * </p>
 *
 */
class Simplex_Optimizer extends Multivariate_Optimizer 
{
    /** Simplex update rule. */
    private Abstract_Simplex simplex;

    /**
     * @param checker Convergence checker.
     */
    public Simplex_Optimizer(Convergence_Checker<Point_valuePair> checker) 
    {
        super(checker);
    }

    /**
     * @param rel Relative threshold.
     * @param abs Absolute threshold.
     */
    public Simplex_Optimizer(double rel, double abs) 
    {
        this(new Simple_Value_Checker(rel, abs));
    }

    /**
     * {@inherit_doc}
     *
     * @param opt_data Optimization data. In addition to those documented in
     * {@link Multivariate_Optimizer#parse_optimization_data(Optimization_data[])
     * Multivariate_Optimizer}, this method will register the following data:
     * <ul>
     *  <li>{@link Abstract_Simplex}</li>
     * </ul>
     * @return {@inherit_doc}
     */
    //override
    public Point_valuePair optimize(Optimization_data... opt_data) 
    {
        // Set up base class and perform computation.
        return super.optimize(opt_data);
    }

    /** {@inherit_doc} */
    //override
    protected Point_valuePair do_optimize() 
    {
        check_parameters();

        // Indirect call to "compute_objective_value" in order to update the
        // evaluations counter.
        const Multivariate_Function eval_func
            = Multivariate_Function() 
            {
                /** {@inherit_doc} */
                //override
                public double value(std::vector<double> point) 
                {
                    return compute_objective_value(point);
                }
            };

        const bool is_minim = get_goal_type() == Goal_Type.MINIMIZE;
        const Comparator<Point_valuePair> comparator
            = Comparator<Point_valuePair>() 
            {
            /** {@inherit_doc} */
            //override
            public int compare(const Point_valuePair o1, const Point_valuePair o2) 
            {
                const double v1 = o1.get_value();
                const double v2 = o2.get_value();
                return is_minim ? Double.compare(v1, v2) : Double.compare(v2, v1);
            }
        };

        // Initialize search.
        simplex.build(get_start_point());
        simplex.evaluate(eval_func, comparator);

        Point_valuePair[] previous = NULL;
        int iteration = 0;
        const Convergence_Checker<Point_valuePair> checker = get_convergence_checker();
        while (true) 
        {
            if (get_iterations() > 0) 
            {
                bool converged = true;
                for (int i{}; i < simplex.get_size(); i++) 
                {
                    Point_valuePair prev = previous[i];
                    converged = converged &&
                        checker.converged(iteration, prev, simplex.get_point(i));
                }
                if (converged) 
                {
                    // We have found an optimum.
                    return simplex.get_point(0);
                }
            }

            // We still need to search.
            previous = simplex.get_points();
            simplex.iterate(eval_func, comparator);

            increment_iteration_count();
        }
    }

    /**
     * Scans the list of (required and optional) optimization data that
     * characterize the problem.
     *
     * @param opt_data Optimization data.
     * The following data will be looked for:
     * <ul>
     *  <li>{@link Abstract_Simplex}</li>
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
            if (data instanceof Abstract_Simplex) 
            {
                simplex = (Abstract_Simplex) data;
                // If more data must be parsed, this statement _must_ be
                // changed to "continue".
                break;
            }
        }
    }

    /**
     * @Math_Runtime_Exception if bounds were passed to the
     * {@link #optimize(Optimization_data[]) optimize} method.
     * @Null_Argument_Exception if no initial simplex was passed to the
     * {@link #optimize(Optimization_data[]) optimize} method.
     */
    private void check_parameters() 
    {
        if (simplex == NULL) 
        {
            throw Null_Argument_Exception();
        }
        if (get_lower_bound() != NULL ||
            get_upper_bound() != NULL) 
            {
            throw Math_Runtime_Exception(Localized_Core_Formats.CONSTRAINT);
        }
    }
}


