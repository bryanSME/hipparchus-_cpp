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
//package org.hipparchus.optim.nonlinear.scalar;

//import org.hipparchus.analysis.Multivariate_Function;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.optim.BaseMultivariate_Optimizer;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.optim.Point_valuePair;

/**
 * Base class for a multivariate scalar function optimizer.
 *
 */
class Multivariate_Optimizer
    extends BaseMultivariate_Optimizer<Point_valuePair> 
    {
    /** Objective function. */
    private Multivariate_Function function;
    /** Type of optimization. */
    private Goal_Type goal;

    /**
     * @param checker Convergence checker.
     */
    protected Multivariate_Optimizer(Convergence_Checker<Point_valuePair> checker) 
    {
        super(checker);
    }

    /**
     * {@inherit_doc}
     *
     * @param opt_data Optimization data. In addition to those documented in
     * {@link BaseMultivariate_Optimizer#parse_optimization_data(Optimization_data[])
     * BaseMultivariate_Optimizer}, this method will register the following data:
     * <ul>
     *  <li>{@link Objective_Function}</li>
     *  <li>{@link Goal_Type}</li>
     * </ul>
     * @return {@inherit_doc}
     * @Math_Illegal_State_Exception if the maximal number of
     * evaluations is exceeded.
     */
    //override
    public Point_valuePair optimize(Optimization_data... opt_data)
        Math_Illegal_State_Exception 
        {
        // Set up base class and perform computation.
        return super.optimize(opt_data);
    }

    /**
     * Scans the list of (required and optional) optimization data that
     * characterize the problem.
     *
     * @param opt_data Optimization data.
     * The following data will be looked for:
     * <ul>
     *  <li>{@link Objective_Function}</li>
     *  <li>{@link Goal_Type}</li>
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
            if (data instanceof Goal_Type) 
            {
                goal = (Goal_Type) data;
                continue;
            }
            if (data instanceof Objective_Function) 
            {
                function = ((Objective_Function) data).get_objective_function();
                continue;
            }
        }
    }

    /**
     * @return the optimization type.
     */
    public Goal_Type get_goal_type() 
    {
        return goal;
    }

    /**
     * Computes the objective function value.
     * This method <em>must</em> be called by subclasses to enforce the
     * evaluation counter limit.
     *
     * @param params Point at which the objective function must be evaluated.
     * @return the objective function value at the specified point.
     * @Math_Illegal_State_Exception if the maximal number of
     * evaluations is exceeded.
     */
    public double compute_objective_value(std::vector<double> params) 
    {
        super.increment_evaluation_count();
        return function.value(params);
    }
}


