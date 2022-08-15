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

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.optim.Base_Optimizer;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.optim.nonlinear.scalar.Goal_Type;

/**
 * Base class for a univariate scalar function optimizer.
 *
 */
class Univariate_Optimizer
    extends Base_Optimizer<UnivariatePoint_valuePair> 
    {
    /** Objective function. */
    private Univariate_Function function;
    /** Type of optimization. */
    private Goal_Type goal;
    /** Initial guess. */
    private double start;
    /** Lower bound. */
    private double min;
    /** Upper bound. */
    private double max;

    /**
     * @param checker Convergence checker.
     */
    protected Univariate_Optimizer(Convergence_Checker<UnivariatePoint_valuePair> checker) 
    {
        super(checker);
    }

    /**
     * {@inherit_doc}
     *
     * @param opt_data Optimization data. In addition to those documented in
     * {@link Base_Optimizer#parse_optimization_data(Optimization_data[])
     * Base_Optimizer}, this method will register the following data:
     * <ul>
     *  <li>{@link Goal_Type}</li>
     *  <li>{@link Search_interval}</li>
     *  <li>{@link UnivariateObjective_Function}</li>
     * </ul>
     * @return {@inherit_doc}
     * @Math_Illegal_State_Exception if the maximal number of
     * evaluations is exceeded.
     */
    //override
    public UnivariatePoint_valuePair optimize(Optimization_data... opt_data)
        Math_Illegal_State_Exception 
        {
        // Perform computation.
        return super.optimize(opt_data);
    }

    /**
     * @return the optimization type.
     */
    public Goal_Type get_goal_type() 
    {
        return goal;
    }

    /**
     * Scans the list of (required and optional) optimization data that
     * characterize the problem.
     *
     * @param opt_data Optimization data.
     * The following data will be looked for:
     * <ul>
     *  <li>{@link Goal_Type}</li>
     *  <li>{@link Search_interval}</li>
     *  <li>{@link UnivariateObjective_Function}</li>
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
            if (data instanceof Search_interval) 
            {
                const Search_interval interval = (Search_interval) data;
                min = interval.get_min();
                max = interval.get_max();
                start = interval.get_start_value();
                continue;
            }
            if (data instanceof UnivariateObjective_Function) 
            {
                function = ((UnivariateObjective_Function) data).get_objective_function();
                continue;
            }
            if (data instanceof Goal_Type) 
            {
                goal = (Goal_Type) data;
                continue;
            }
        }
    }

    /**
     * @return the initial guess.
     */
    public double get_start_value() 
    {
        return start;
    }
    /**
     * @return the lower bounds.
     */
    public double get_min() 
    {
        return min;
    }
    /**
     * @return the upper bounds.
     */
    public double get_max() 
    {
        return max;
    }

    /**
     * Computes the objective function value.
     * This method <em>must</em> be called by subclasses to enforce the
     * evaluation counter limit.
     *
     * @param x Point at which the objective function must be evaluated.
     * @return the objective function value at the specified point.
     * @Math_Illegal_State_Exception if the maximal number of
     * evaluations is exceeded.
     */
    protected double compute_objective_value(double x) 
    {
        super.increment_evaluation_count();
        return function.value(x);
    }
}


