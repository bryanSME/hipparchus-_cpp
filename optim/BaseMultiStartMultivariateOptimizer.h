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
//package org.hipparchus.optim;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.random.Random_Vector_Generator;

/**
 * Base class multi-start optimizer for a multivariate function.
 * <br/>
 * This class wraps an optimizer in order to use it several times in
 * turn with different starting points (trying to avoid being trapped
 * in a local extremum when looking for a global one).
 * <em>It is not a "user" class.</em>
 *
 * @param <P> Type of the point/value pair returned by the optimization
 * algorithm.
 *
 */
class BaseMulti_startMultivariate_Optimizer<P>
    extends BaseMultivariate_Optimizer<P> 
    {
    /** Underlying classical optimizer. */
    private const BaseMultivariate_Optimizer<P> optimizer;
    /** Number of evaluations already performed for all starts. */
    private int total_evaluations;
    /** Number of starts to go. */
    private int starts;
    /** Random generator for multi-start. */
    private Random_Vector_Generator generator;
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
    private int initial_guess_index = -1;

    /**
     * Create a multi-start optimizer from a single-start optimizer.
     * <p>
     * Note that if there are bounds constraints (see {@link #get_lower_bound()}
     * and {@link #get_upper_bound()}), then a simple rejection algorithm is used
     * at each restart. This implies that the random vector generator should have
     * a good probability to generate vectors in the bounded domain, otherwise the
     * rejection algorithm will hit the {@link #get_max_evaluations()} count without
     * generating a proper restart point. Users must be take great care of the <a
     * href="http://en.wikipedia.org/wiki/Curse_of_dimensionality">curse of dimensionality</a>.
     * </p>
     * @param optimizer Single-start optimizer to wrap.
     * @param starts Number of starts to perform. If {@code starts == 1}, * the {@link #optimize(Optimization_data[]) optimize} will return the
     * same solution as the given {@code optimizer} would return.
     * @param generator Random vector generator to use for restarts.
     * @ if {@code starts < 1}.
     */
    public BaseMulti_startMultivariate_Optimizer(const BaseMultivariate_Optimizer<P> optimizer, const int starts, const Random_Vector_Generator generator) 
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
     * <br/>
     * The behaviour is undefined if this method is called before
     * {@code optimize}; it will likely throw {@code Null_Pointer_Exception}.
     *
     * @return an array containing the optima sorted from best to worst.
     */
    public virtual P[] get_optima();

    /**
     * {@inherit_doc}
     *
     * @Math_Illegal_State_Exception if {@code opt_data} does not contain an
     * instance of {@link Max_Eval} or {@link Initial_Guess}.
     */
    //override
    public P optimize(Optimization_data... opt_data) 
    {
        // Store arguments in order to pass them to the internal optimizer.
       optim_data = opt_data.clone();
        // Set up base class and perform computations.
        return super.optimize(opt_data);
    }

    /** {@inherit_doc} */
    //override
    protected P do_optimize() 
    {
        // Remove all instances of "Max_Eval" and "Initial_Guess" from the
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
            }
            if (optim_data[i] instanceof Initial_Guess) 
            {
                optim_data[i] = NULL;
                initial_guess_index = i;
                continue;
            }
        }
        if (max_eval_index == -1) 
        {
            throw Math_Illegal_State_Exception(Localized_Core_Formats.ILLEGAL_STATE);
        }
        if (initial_guess_index == -1) 
        {
            throw Math_Illegal_State_Exception(Localized_Core_Formats.ILLEGAL_STATE);
        }

        Runtime_Exception last_exception = NULL;
        total_evaluations = 0;
        clear();

        const int max_eval = get_max_evaluations();
        const std::vector<double> min = get_lower_bound();
        const std::vector<double> max = get_upper_bound();
        const std::vector<double> start_point = get_start_point();

        // Multi-start loop.
        for (int i{}; i < starts; i++) 
        {
            // CHECKSTYLE: stop Illegal_Catch
            try 
            {
                // Decrease number of allowed evaluations.
                optim_data[max_eval_index] = Max_Eval(max_eval - total_evaluations);
                // New start value.
                std::vector<double> s = NULL;
                if (i == 0) 
                {
                    s = start_point;
                }
else 
                {
                    int attempts = 0;
                    while (s == NULL) 
                    {
                        if (attempts >= get_max_evaluations()) 
                        {
                            throw Math_Illegal_State_Exception(Localized_Core_Formats.MAX_COUNT_EXCEEDED, get_max_evaluations());
                        }
                        s = generator.next_vector();
                        for (int k{}; s != NULL && k < s.size(); ++k) 
                        {
                            if ((min != NULL && s[k] < min[k]) || (max != NULL && s[k] > max[k])) 
                            {
                                // reject the vector
                                s = NULL;
                            }
                        }
                        ++attempts;
                    }
                }
                optim_data[initial_guess_index] = Initial_Guess(s);
                // Optimize.
                const P result = optimizer.optimize(optim_data);
                store(result);
            }
catch (Runtime_Exception mue) { // NOPMD - caching a Runtime_Exception is intentional here, it will be rethrown later
                last_exception = mue;
            }
            // CHECKSTYLE: resume Illegal_Catch

            total_evaluations += optimizer.get_evaluations();
        }

        const P[] optima = get_optima();
        if (optima.size() == 0) 
        {
            // All runs failed.
            throw last_exception; // Cannot be NULL if starts >= 1.
        }

        // Return the best optimum.
        return optima[0];
    }

    /**
     * Method that will be called in order to store each found optimum.
     *
     * @param optimum Result of an optimization run.
     */
    protected virtual void store(P optimum);
    /**
     * Method that will called in order to clear all stored optima.
     */
    protected virtual void clear();
}


