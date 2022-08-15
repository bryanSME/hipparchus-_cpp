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
//package org.hipparchus.optim.nonlinear.vector.leastsquares;

//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Optimizer.Optimum;
//import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem.Evaluation;

/**
 * A pedantic implementation of {@link Optimum}.
 *
 */
class Optimum_Impl : Optimum 
{

    /** abscissa and ordinate */
    private const Evaluation value;
    /** number of evaluations to compute this optimum */
    private const int evaluations;
    /** number of iterations to compute this optimum */
    private const int iterations;

    /**
     * Construct an optimum from an evaluation and the values of the counters.
     *
     * @param value       the function value
     * @param evaluations number of times the function was evaluated
     * @param iterations  number of iterations of the algorithm
     */
    Optimum_Impl(const Evaluation value, const int evaluations, const int iterations) 
    {
        this.value = value;
        this.evaluations = evaluations;
        this.iterations = iterations;
    }

    /* auto-generated implementations */

    /** {@inherit_doc} */
    //override
    public int get_evaluations() 
    {
        return evaluations;
    }

    /** {@inherit_doc} */
    //override
    public int get_iterations() 
    {
        return iterations;
    }

    /** {@inherit_doc} */
    //override
    public Real_Matrix get_covariances(double threshold) 
    {
        return value.get_covariances(threshold);
    }

    /** {@inherit_doc} */
    //override
    public Real_Vector get_sigma(double covariance_singularity_threshold) 
    {
        return value.get_sigma(covariance_singularity_threshold);
    }

    /** {@inherit_doc} */
    //override
    public double get_r_m_s() 
    {
        return value.get_r_m_s();
    }

    /** {@inherit_doc} */
    //override
    public Real_Matrix get_jacobian() 
    {
        return value.get_jacobian();
    }

    /** {@inherit_doc} */
    //override
    public double get_cost() 
    {
        return value.get_cost();
    }

    /** {@inherit_doc} */
    //override
    public double get_chi_square() 
    {
        return value.get_chi_square();
    }

    /** {@inherit_doc} */
    //override
    public double get_reduced_chi_square(const int& n) 
    {
        return value.get_reduced_chi_square(n);
    }

    /** {@inherit_doc} */
    //override
    public Real_Vector get_residuals() 
    {
        return value.get_residuals();
    }

    /** {@inherit_doc} */
    //override
    public Real_Vector get_point() 
    {
        return value.get_point();
    }
}


