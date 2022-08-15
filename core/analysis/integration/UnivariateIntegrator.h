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
//package org.hipparchus.analysis.integration;

#include "../../analysis/UnivariateFunction.h"
//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;

/**
 * Interface for univariate real integration algorithms.
 *
 */
class Univariate_Integrator
{
    /**
     * Get the relative accuracy.
     *
     * @return the accuracy
     */
    virtual double get_relative_accuracy();

    /**
     * Get the absolute accuracy.
     *
     * @return the accuracy
     */
    virtual double get_absolute_accuracy();

    /**
     * Get the min limit for the number of iterations.
     *
     * @return the actual min limit
     */
    virtual int get_minimal_iteration_count();

    /**
     * Get the upper limit for the number of iterations.
     *
     * @return the actual upper limit
     */
    virtual int get_maximal_iteration_count();

    /**
     * Integrate the function in the given interval.
     *
     * @param max_eval Maximum number of evaluations.
     * @param f the integrand function
     * @param min the lower bound for the interval
     * @param max the upper bound for the interval
     * @return the value of integral
     * @Math_Illegal_State_Exception if the maximum number of function
     * evaluations is exceeded
     * @Math_Illegal_State_Exception if the maximum iteration count is exceeded
     * or the integrator detects convergence problems otherwise
     * @ if {@code min > max} or the endpoints do not
     * satisfy the requirements specified by the integrator
     * @Null_Argument_Exception if {@code f} is {@code NULL}.
     */
    virtual double integrate(const int& max_eval, const Univariate_Function& f, const double& min, double max)

    /**
     * Get the number of function evaluations of the last run of the integrator.
     *
     * @return number of function evaluations
     */
    virtual int get_evaluations();

    /**
     * Get the number of iterations of the last run of the integrator.
     *
     * @return number of iterations
     */
    virtual int get_iterations();

};