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

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Utils;

/**
 * Implementation of the uniform real distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)">
 * Uniform distribution (continuous), at Wikipedia</a>
 */
class UniformReal_Distribution extends Abstract_Real_Distribution 
{
    
    20120109L;
    /** Lower bound of this distribution (inclusive). */
    private const double lower;
    /** Upper bound of this distribution (exclusive). */
    private const double upper;

    /**
     * Create a standard uniform real distribution with lower bound (inclusive)
     * equal to zero and upper bound (exclusive) equal to one.
     */
    public UniformReal_Distribution() 
    {
        this(0, 1);
    }

    /**
     * Create a uniform real distribution using the given lower and upper
     * bounds.
     *
     * @param lower Lower bound of this distribution (inclusive).
     * @param upper Upper bound of this distribution (exclusive).
     * @ if {@code lower >= upper}.
     */
    public UniformReal_Distribution(double lower, double upper)
         
        {
        super();
        if (lower >= upper) 
        {
            throw (
                            hipparchus::exception::Localized_Core_Formats_Type::LOWER_BOUND_NOT_BELOW_UPPER_BOUND, lower, upper, false);
        }

        this.lower = lower;
        this.upper = upper;
    }

    /** {@inherit_doc} */
    //override
    public double density(double x) 
    {
        if (x < lower || x > upper) 
        {
            return 0.0;
        }
        return 1 / (upper - lower);
    }

    /** {@inherit_doc} */
    //override
    public double cumulative_probability(const double& x)  
    {
        if (x <= lower) 
        {
            return 0;
        }
        if (x >= upper) 
        {
            return 1;
        }
        return (x - lower) / (upper - lower);
    }

    /** {@inherit_doc} */
    //override
    public double inverse_cumulative_probability(const double& p)
         
        {
        Math_Utils::check_range_inclusive(p, 0, 1);
        return p * (upper - lower) + lower;
    }

    /**
     * {@inherit_doc}
     *
     * For lower bound {@code lower} and upper bound {@code upper}, the mean is
     * {@code 0.5 * (lower + upper)}.
     */
    //override
    public double get_numerical_mean() const 
    {
        return 0.5 * (lower + upper);
    }

    /**
     * {@inherit_doc}
     *
     * For lower bound {@code lower} and upper bound {@code upper}, the
     * variance is {@code (upper - lower)^2 / 12}.
     */
    //override
    public double get_numerical_variance() const 
    {
        double ul = upper - lower;
        return ul * ul / 12;
    }

    /**
     * {@inherit_doc}
     *
     * The lower bound of the support is equal to the lower bound parameter
     * of the distribution.
     *
     * @return lower bound of the support
     */
    //override
    public double get_support_lower_bound() const 
    {
        return lower;
    }

    /**
     * {@inherit_doc}
     *
     * The upper bound of the support is equal to the upper bound parameter
     * of the distribution.
     *
     * @return upper bound of the support
     */
    //override
    public double get_support_upper_bound() const 
    {
        return upper;
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


