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
//import org.hipparchus.optim.Optimization_data;

/**
 * Scalar function to be optimized.
 *
 */
class UnivariateObjective_Function : Optimization_data 
{
    /** Function to be optimized. */
    private const Univariate_Function function;

    /**
     * @param f Function to be optimized.
     */
    public UnivariateObjective_Function(Univariate_Function f) 
    {
        function = f;
    }

    /**
     * Gets the function to be optimized.
     *
     * @return the objective function.
     */
    public Univariate_Function get_objective_function() 
    {
        return function;
    }
}


