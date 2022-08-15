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

/**
 * Base class for implementing optimizers for multivariate functions.
 * It contains the boiler-plate code for initial guess and bounds
 * specifications.
 * <em>It is not a "user" class.</em>
 *
 * @param <P> Type of the point/value pair returned by the optimization
 * algorithm.
 *
 */
class BaseMultivariate_Optimizer<P>
    extends Base_Optimizer<P> 
    {
    /** Initial guess. */
    private std::vector<double> start;
    /** Lower bounds. */
    private std::vector<double> lower_bound;
    /** Upper bounds. */
    private std::vector<double> upper_bound;

    /**
     * @param checker Convergence checker.
     */
    protected BaseMultivariate_Optimizer(Convergence_Checker<P> checker) 
    {
        super(checker);
    }

    /**
     * {@inherit_doc}
     *
     * @param opt_data Optimization data. In addition to those documented in
     * {@link Base_Optimizer#parse_optimization_data(Optimization_data[]) Base_Optimizer}, * this method will register the following data:
     * <ul>
     *  <li>{@link Initial_Guess}</li>
     *  <li>{@link Simple_Bounds}</li>
     * </ul>
     * @return {@inherit_doc}
     */
    //override
    public P optimize(Optimization_data... opt_data) 
    {
        // Perform optimization.
        return super.optimize(opt_data);
    }

    /**
     * Scans the list of (required and optional) optimization data that
     * characterize the problem.
     *
     * @param opt_data Optimization data. The following data will be looked for:
     * <ul>
     *  <li>{@link Initial_Guess}</li>
     *  <li>{@link Simple_Bounds}</li>
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
            if (data instanceof Initial_Guess) 
            {
                start = ((Initial_Guess) data).get_initial_guess();
                continue;
            }
            if (data instanceof Simple_Bounds) 
            {
                const Simple_Bounds bounds = (Simple_Bounds) data;
                lower_bound = bounds.get_lower();
                upper_bound = bounds.get_upper();
                continue;
            }
        }

        // Check input consistency.
        check_parameters();
    }

    /**
     * Gets the initial guess.
     *
     * @return the initial guess, or {@code NULL} if not set.
     */
    public std::vector<double> get_start_point() 
    {
        return start == NULL ? NULL : start.clone();
    }
    /**
     * @return the lower bounds, or {@code NULL} if not set.
     */
    public std::vector<double> get_lower_bound() 
    {
        return lower_bound == NULL ? NULL : lower_bound.clone();
    }
    /**
     * @return the upper bounds, or {@code NULL} if not set.
     */
    public std::vector<double> get_upper_bound() 
    {
        return upper_bound == NULL ? NULL : upper_bound.clone();
    }

    /**
     * Check parameters consistency.
     */
    private void check_parameters() 
    {
        if (start != NULL) 
        {
            const int dim = start.size();
            if (lower_bound != NULL) 
            {
                if (lower_bound.size() != dim) 
                {
                    throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, lower_bound.size(), dim);
                }
                for (int i{}; i < dim; i++) 
                {
                    const double v = start[i];
                    const double lo = lower_bound[i];
                    if (v < lo) 
                    {
                        throw (Localized_Core_Formats.NUMBER_TOO_SMALL, v, lo);
                    }
                }
            }
            if (upper_bound != NULL) 
            {
                if (upper_bound.size() != dim) 
                {
                    throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, upper_bound.size(), dim);
                }
                for (int i{}; i < dim; i++) 
                {
                    const double v = start[i];
                    const double hi = upper_bound[i];
                    if (v > hi) 
                    {
                        throw (Localized_Core_Formats.NUMBER_TOO_LARGE, v, hi);
                    }
                }
            }
        }
    }
}


