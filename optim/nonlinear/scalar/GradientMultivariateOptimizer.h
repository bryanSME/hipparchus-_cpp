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

//import org.hipparchus.analysis.MultivariateVector_function;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.optim.Point_valuePair;

/**
 * Base class for implementing optimizers for multivariate scalar
 * differentiable functions.
 * It contains boiler-plate code for dealing with gradient evaluation.
 *
 */
class GradientMultivariate_Optimizer
    extends Multivariate_Optimizer 
    {
    /**
     * Gradient of the objective function.
     */
    private MultivariateVector_function gradient;

    /**
     * @param checker Convergence checker.
     */
    protected GradientMultivariate_Optimizer(Convergence_Checker<Point_valuePair> checker) 
    {
        super(checker);
    }

    /**
     * Compute the gradient vector.
     *
     * @param params Point at which the gradient must be evaluated.
     * @return the gradient at the specified point.
     */
    protected std::vector<double> compute_objective_gradient(const std::vector<double> params) 
    {
        return gradient.value(params);
    }

    /**
     * {@inherit_doc}
     *
     * @param opt_data Optimization data. In addition to those documented in
     * {@link Multivariate_Optimizer#parse_optimization_data(Optimization_data[])
     * Multivariate_Optimizer}, this method will register the following data:
     * <ul>
     *  <li>{@link Objective_Function_Gradient}</li>
     * </ul>
     * @return {@inherit_doc}
     * @Math_Illegal_State_Exception if the maximal number of
     * evaluations (of the objective function) is exceeded.
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
     *  <li>{@link Objective_Function_Gradient}</li>
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
            if  (data instanceof Objective_Function_Gradient) 
            {
                gradient = ((Objective_Function_Gradient) data).get_objective_function_gradient();
                // If more data must be parsed, this statement _must_ be
                // changed to "continue".
                break;
            }
        }
    }
}


