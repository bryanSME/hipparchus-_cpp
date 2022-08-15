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
//import org.hipparchus.optim.Optimization_data;

/**
 * Gradient of the scalar function to be optimized.
 *
 */
class Objective_Function_Gradient : Optimization_data 
{
    /** Function to be optimized. */
    private const MultivariateVector_function gradient;

    /**
     * @param g Gradient of the function to be optimized.
     */
    public Objective_Function_Gradient(MultivariateVector_function g) 
    {
        gradient = g;
    }

    /**
     * Gets the gradient of the function to be optimized.
     *
     * @return the objective function gradient.
     */
    public MultivariateVector_function get_objective_function_gradient() 
    {
        return gradient;
    }
}


