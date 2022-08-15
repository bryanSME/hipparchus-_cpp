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
//package org.hipparchus.optim.linear;

//import java.util.Collection;
//import java.util.Collections;

//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.optim.Point_valuePair;
//import org.hipparchus.optim.nonlinear.scalar.Multivariate_Optimizer;

/**
 * Base class for implementing linear optimizers.
 *
 */
class Linear_Optimizer
    extends Multivariate_Optimizer 
    {
    /**
     * Linear objective function.
     */
    private Linear_Objective_Function function;
    /**
     * Linear constraints.
     */
    private Collection<Linear_Constraint> linear_constraints;
    /**
     * Whether to restrict the variables to non-negative values.
     */
    private bool non_negative;

    /**
     * Simple constructor with default settings.
     *
     */
    protected Linear_Optimizer() 
    {
        super(null); // No convergence checker.
    }

    /**
     * @return {@code true} if the variables are restricted to non-negative values.
     */
    protected bool is_restricted_to_non_negative() 
    {
        return non_negative;
    }

    /**
     * @return the optimization type.
     */
    protected Linear_Objective_Function get_function() 
    {
        return function;
    }

    /**
     * @return the optimization type.
     */
    protected Collection<Linear_Constraint> get_constraints() 
    {
        return Collections.unmodifiable_collection(linear_constraints);
    }

    /**
     * {@inherit_doc}
     *
     * @param opt_data Optimization data. In addition to those documented in
     * {@link Multivariate_Optimizer#parse_optimization_data(Optimization_data[])
     * Multivariate_Optimizer}, this method will register the following data:
     * <ul>
     *  <li>{@link Linear_Objective_Function}</li>
     *  <li>{@link Linear_ConstraintSet}</li>
     *  <li>{@link Non_Negative_Constraint}</li>
     * </ul>
     * @return {@inherit_doc}
     * @Math_Illegal_State_Exception if the maximal number of
     * iterations is exceeded.
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
     *  <li>{@link Linear_Objective_Function}</li>
     *  <li>{@link Linear_ConstraintSet}</li>
     *  <li>{@link Non_Negative_Constraint}</li>
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
            if (data instanceof Linear_Objective_Function) 
            {
                function = (Linear_Objective_Function) data;
                continue;
            }
            if (data instanceof Linear_ConstraintSet) 
            {
                linear_constraints = ((Linear_ConstraintSet) data).get_constraints();
                continue;
            }
            if  (data instanceof Non_Negative_Constraint) 
            {
                non_negative = ((Non_Negative_Constraint) data).is_restricted_to_non_negative();
                continue;
            }
        }
    }
}


