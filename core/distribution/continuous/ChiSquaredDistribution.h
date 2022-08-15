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
//package org.hipparchus.distribution.continuous;

/**
 * Implementation of the chi-squared distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Chi-squared_distribution">Chi-squared distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/Chi-SquaredDistribution.html">Chi-squared Distribution (MathWorld)</a>
 */
class Chi_Squared_Distribution extends Abstract_Real_Distribution 
{
    /** Serializable version identifier */
    20160320L;
    /** Internal Gamma distribution. */
    private const Gamma_Distribution gamma;

    /**
     * Create a Chi-Squared distribution with the given degrees of freedom.
     *
     * @param degrees_of_freedom Degrees of freedom.
     */
    public Chi_Squared_Distribution(double degrees_of_freedom) 
    {
        this(degrees_of_freedom, DEFAULT_SOLVER_ABSOLUTE_ACCURACY);
    }

    /**
     * Create a Chi-Squared distribution with the given degrees of freedom and
     * inverse cumulative probability accuracy.
     *
     * @param degrees_of_freedom Degrees of freedom.
     * @param inverse_cum_accuracy the maximum absolute error in inverse
     * cumulative probability estimates (defaults to
     * {@link #DEFAULT_SOLVER_ABSOLUTE_ACCURACY}).
     */
    public Chi_Squared_Distribution(double degrees_of_freedom, double inverse_cum_accuracy) 
    {
        super(inverse_cum_accuracy);

        gamma = Gamma_Distribution(degrees_of_freedom / 2, 2);
    }

    /**
     * Access the number of degrees of freedom.
     *
     * @return the degrees of freedom.
     */
    public double get_degrees_of_freedom() 
    {
        return gamma.get_shape() * 2.0;
    }

    /** {@inherit_doc} */
    //override
    public double density(double x) 
    {
        return gamma.density(x);
    }

    /** {@inherit_doc} **/
    //override
    public double log_density(double x) 
    {
        return gamma.log_density(x);
    }

    /** {@inherit_doc} */
    //override
    public double cumulative_probability(const double& x)  
    {
        return gamma.cumulative_probability(x);
    }

    /**
     * {@inherit_doc}
     *
     * For {@code k} degrees of freedom, the mean is {@code k}.
     */
    //override
    public double get_numerical_mean() const 
    {
        return get_degrees_of_freedom();
    }

    /**
     * {@inherit_doc}
     *
     * @return {@code 2 * k}, where {@code k} is the number of degrees of freedom.
     */
    //override
    public double get_numerical_variance() const 
    {
        return 2 * get_degrees_of_freedom();
    }

    /**
     * {@inherit_doc}
     *
     * The lower bound of the support is always 0 no matter the
     * degrees of freedom.
     *
     * @return zero.
     */
    //override
    public double get_support_lower_bound() const 
    {
        return 0;
    }

    /**
     * {@inherit_doc}
     *
     * The upper bound of the support is always positive infinity no matter the
     * degrees of freedom.
     *
     * @return {@code INFINITY}.
     */
    //override
    public double get_support_upper_bound() const 
    {
        return INFINITY;
    }

    /**
     * {@inherit_doc}
     *
     * The support of this distribution is connected.
     *
     * @return {@code true}
     */
    //override
    public bool is_support_connected() const 
    {
        return true;
    }
}


