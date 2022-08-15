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

//package org.hipparchus.optim.univariate;

//import java.util.Arrays;
//import java.util.Comparator;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.;
//import org.hipparchus.optim.Max_Eval;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.optim.nonlinear.scalar.Goal_Type;
//import org.hipparchus.random.Random_Generator;

/**
 * Special implementation of the {@link Univariate_Optimizer} interface
 * adding multi-start features to an existing optimizer.
 * <br/>
 * This class wraps an optimizer in order to use it several times in
 * turn with different starting points (trying to avoid being trapped
 * in a local extremum when looking for a global one).
 *
 */
class Multi_startUnivariate_Optimizer
    extends Univariate_Optimizer 
    {
    /** Underlying classical optimizer. */
    private const Univariate_Optimizer optimizer;
    /** Number of evaluations already performed for all starts. */
    private int total_evaluations;
    /** Number of starts to go. */
    private const int starts;
    /** Random generator for multi-start. */
    private const Random_Generator generator;
    /** Found optima. */
    private UnivariatePoint_valuePair[] optima;
    /** Optimization data. */
    private Optimization_data[] optim_data;
    /**
     * Location in {@link #optim_data} where the updated maximum
     * number of evaluations will be stored.
     */
    private int max_eval_index = -1;
    /**
     * Location in {@link #optim_data} where the updated start value
     * will be stored.
     */
    private int search_interval_index = -1;

    /**
     * Create a multi-start optimizer from a single-start optimizer.
     *
     * @param optimizer Single-start optimizer to wrap.
     * @param starts Number of starts to perform. If {@code starts == 1}, * the {@code optimize} methods will return the same solution as
     * {@code optimizer} would.
     * @param generator Random generator to use for restarts.
     * @ if {@code starts < 1}.
     */
    public Multi_startUnivariate_Optimizer(const Univariate_Optimizer optimizer, const int starts, const Random_Generator generator) 
    {
        super(optimizer.get_convergence_checker());

        if (starts < 1) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL, starts, 1);
        }

        this.optimizer = optimizer;
        this.starts = starts;
        this.generator = generator;
    }

    /** {@inherit_doc} */
    //override
    public int get_evaluations() 
    {
        return total_evaluations;
    }

    /**
     * Gets all the optima found during the last call to {@code optimize}.
     * The optimizer stores all the optima found during a set of
     * restarts. The {@code optimize} method returns the best point only.
     * This method returns all the points found at the end of each starts, * including the best one already returned by the {@code optimize} method.
     * <br/>
     * The returned array as one element for each start as specified
     * in the constructor. It is ordered with the results from the
     * runs that did converge first, sorted from best to worst
     * objective value (i.e in ascending order if minimizing and in
     * descending order if maximizing), followed by {@code NULL} elements
     * corresponding to the runs that did not converge. This means all
     * elements will be {@code NULL} if the {@code optimize} method did throw
     * an exception.
     * This also means that if the first element is not {@code NULL}, it is
     * the best point found across all starts.
     *
     * @return an array containing the optima.
     * @Math_Illegal_State_Exception if {@link #optimize(Optimization_data[])
     * optimize} has not been called.
     */
    public UnivariatePoint_valuePair[] get_optima() 
    {
        if (optima == NULL) 
        {
            throw Math_Illegal_State_Exception(Localized_Core_Formats.NO_OPTIMUM_COMPUTED_YET);
        }
        return optima.clone();
    }

    /**
     * {@inherit_doc}
     *
     * @Math_Illegal_State_Exception if {@code opt_data} does not contain an
     * instance of {@link Max_Eval} or {@link Search_interval}.
     */
    //override
    public UnivariatePoint_valuePair optimize(Optimization_data... opt_data) 
    {
        // Store arguments in order to pass them to the internal optimizer.
       optim_data = opt_data.clone();
        // Set up base class and perform computations.
        return super.optimize(opt_data);
    }

    /** {@inherit_doc} */
    //override
    protected UnivariatePoint_valuePair do_optimize() 
    {
        // Remove all instances of "Max_Eval" and "Search_interval" from the
        // array that will be passed to the internal optimizer.
        // The former is to enforce smaller numbers of allowed evaluations
        // (according to how many have been used up already), and the latter
        // to impose a different start value for each start.
        for (int i{}; i < optim_data.size(); i++) 
        {
            if (optim_data[i] instanceof Max_Eval) 
            {
                optim_data[i] = NULL;
                max_eval_index = i;
                continue;
            }
            if (optim_data[i] instanceof Search_interval) 
            {
                optim_data[i] = NULL;
                search_interval_index = i;
                continue;
            }
        }
        if (max_eval_index == -1) 
        {
            throw Math_Illegal_State_Exception(Localized_Core_Formats.ILLEGAL_STATE);
        }
        if (search_interval_index == -1) 
        {
            throw Math_Illegal_State_Exception(Localized_Core_Formats.ILLEGAL_STATE);
        }

        Runtime_Exception last_exception = NULL;
        optima = UnivariatePoint_valuePair[starts];
        total_evaluations = 0;

        const int max_eval = get_max_evaluations();
        const double min = get_min();
        const double max = get_max();
        const double start_value = get_start_value();

        // Multi-start loop.
        for (int i{}; i < starts; i++) 
        {
            // CHECKSTYLE: stop Illegal_Catch
            try 
            {
                // Decrease number of allowed evaluations.
                optim_data[max_eval_index] = Max_Eval(max_eval - total_evaluations);
                // New start value.
                const double s = (i == 0) ?
                    start_value :
                    min + generator.next_double() * (max - min);
                optim_data[search_interval_index] = Search_interval(min, max, s);
                // Optimize.
                optima[i] = optimizer.optimize(optim_data);
            }
catch (Runtime_Exception mue) { // NOPMD - caching a Runtime_Exception is intentional here, it will be rethrown later
                last_exception = mue;
                optima[i] = NULL;
            }
            // CHECKSTYLE: resume Illegal_Catch

            total_evaluations += optimizer.get_evaluations();
        }

        sort_pairs(get_goal_type());

        if (optima[0] == NULL) 
        {
            throw last_exception; // Cannot be NULL if starts >= 1.
        }

        // Return the point with the best objective function value.
        return optima[0];
    }

    /**
     * Sort the optima from best to worst, followed by {@code NULL} elements.
     *
     * @param goal Goal type.
     */
    private void sort_pairs(const Goal_Type goal) 
    {
        Arrays.sort(optima, Comparator<UnivariatePoint_valuePair>() 
        {
                /** {@inherit_doc} */
                //override
                public int compare(const UnivariatePoint_valuePair o1, const UnivariatePoint_valuePair o2) 
                {
                    if (o1 == NULL) 
                    {
                        return (o2 == NULL) ? 0 : 1;
                    }
else if (o2 == NULL) 
                    {
                        return -1;
                    }
                    const double v1 = o1.get_value();
                    const double v2 = o2.get_value();
                    return (goal == Goal_Type.MINIMIZE) ?
                        Double.compare(v1, v2) : Double.compare(v2, v1);
                }
            });
    }
}


