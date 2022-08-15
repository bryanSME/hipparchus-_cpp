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
//import org.hipparchus.special.Erf;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
#include <cmath>
#include "AbstractRealDistribution.h"
/**
 * Implementation of the normal (gaussian) distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Normal_distribution">Normal distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/Normal_Distribution.html">Normal distribution (MathWorld)</a>
 */
class Normal_Distribution : Abstract_Real_Distribution 
{
private:
    /** &radic;(2) */
    static constexpr double SQRT2{ std::sqrt(2.0) };
    /** Mean of this distribution. */
    const double mean;
    /** Standard deviation of this distribution. */
    const double standard_deviation;
    /** The value of {@code log(sd) + 0.5*log(2*pi)} stored for faster computation. */
    const double log_standard_deviation_plus_half_log2_pi;

public:
    /**
     * Create a normal distribution with mean equal to zero and standard
     * deviation equal to one.
     */
    Normal_Distribution() 
    {
        Normal_Distribution(0, 1);
    }

    /**
     * Create a normal distribution using the given mean, standard deviation.
     *
     * @param mean Mean for this distribution.
     * @param sd Standard deviation for this distribution.
     * @ if {@code sd <= 0}.
     */
    Normal_Distribution(const double& mean, const double& sd)
    {
        if (sd <= 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::STANDARD_DEVIATION, sd);
        }

        mean = mean;
        standard_deviation = sd;
        log_standard_deviation_plus_half_log2_pi =
                std::log(sd) + 0.5 * std::log(2 * std::numbers::pi);
    }

    /**
     * Access the mean.
     *
     * @return the mean for this distribution.
     */
    double get_mean() const
    {
        return mean;
    }

    /**
     * Access the standard deviation.
     *
     * @return the standard deviation for this distribution.
     */
    double get_standard_deviation() const
    {
        return standard_deviation;
    }

    /** {@inherit_doc} */
    //override
    double density(const double& x)
    {
        return std::exp(log_density(x));
    }

    /** {@inherit_doc} */
    //override
    double log_density(const double& x) 
    {
        const double x0 = x - mean;
        const double x1 = x0 / standard_deviation;
        return -0.5 * x1 * x1 - log_standard_deviation_plus_half_log2_pi;
    }

    /**
     * {@inherit_doc}
     *
     * If {@code x} is more than 40 standard deviations from the mean, 0 or 1
     * is returned, as in these cases the actual value is within
     * {@code Double.MIN_VALUE} of 0 or 1.
     */
    //override
    double cumulative_probability(const double& x)  
    {
        const double dev = x - mean;
        if (std::abs(dev) > 40 * standard_deviation) 
        {
            return dev < 0
                ? 0.0
                : 1.0;
        }
        return 0.5 * Erf.erfc(-dev / (standard_deviation * SQRT2));
    }

    /** {@inherit_doc} */
    //override
    double inverse_cumulative_probability(const double& p)  
    {
        Math_Utils::check_range_inclusive(p, 0, 1);
        return mean + standard_deviation * SQRT2 * Erf.erf_inv(2 * p - 1);
    }

    /** {@inherit_doc} */
    //override
    double probability(const double& x0, const double& x1)
    {
        if (x0 > x1) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT, x0, x1, true);
        }
        const double denom = standard_deviation * SQRT2;
        const double v0 = (x0 - mean) / denom;
        const double v1 = (x1 - mean) / denom;
        return 0.5 * Erf.erf(v0, v1);
    }

    /**
     * {@inherit_doc}
     *
     * For mean parameter {@code mu}, the mean is {@code mu}.
     */
    //override
    double get_numerical_mean() const 
    {
        return get_mean();
    }

    /**
     * {@inherit_doc}
     *
     * For standard deviation parameter {@code s}, the variance is {@code s^2}.
     */
    //override
    double get_numerical_variance() const 
    {
        const double s = get_standard_deviation();
        return s * s;
    }

    /**
     * {@inherit_doc}
     *
     * The lower bound of the support is always negative infinity
     * no matter the parameters.
     *
     * @return lower bound of the support (always
     * {@code -INFINITY})
     */
    //override
    double get_support_lower_bound() const 
    {
        return -INFINITY;
    }

    /**
     * {@inherit_doc}
     *
     * The upper bound of the support is always positive infinity
     * no matter the parameters.
     *
     * @return upper bound of the support (always
     * {@code INFINITY})
     */
    //override
    double get_support_upper_bound() const 
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
    bool is_support_connected() const 
    {
        return true;
    }
};