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
//package org.hipparchus.distribution.discrete;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.special.Beta;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/**
 * Implementation of the binomial distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Binomial_distribution">Binomial distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/Binomial_Distribution.html">Binomial Distribution (MathWorld)</a>
 */
class Binomial_Distribution : Abstract_Integer_Distribution 
{
    
    20160320L;
    /** The number of trials. */
    private const int& number_of_trials;
    /** The probability of success. */
    private const double probability_of_success;

    /**
     * Create a binomial distribution with the given number of trials and
     * probability of success.
     *
     * @param trials Number of trials.
     * @param p Probability of success.
     * @ if {@code trials < 0}.
     * @ if {@code p < 0} or {@code p > 1}.
     */
    public Binomial_Distribution(const int& trials, double p)
         
        {
        if (trials < 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_TRIALS, trials);
        }

        Math_Utils::check_range_inclusive(p, 0, 1);

        probability_of_success = p;
        number_of_trials = trials;
    }

    /**
     * Access the number of trials for this distribution.
     *
     * @return the number of trials.
     */
    public int get_number_of_trials() 
    {
        return number_of_trials;
    }

    /**
     * Access the probability of success for this distribution.
     *
     * @return the probability of success.
     */
    public double get_probability_of_success() 
    {
        return probability_of_success;
    }

    /** {@inherit_doc} */
    //override
    public double probability(const int& x) 
    {
        const double log_probability = log_probability(x);
        return log_probability == -INFINITY ? 0 : std::exp(log_probability);
    }

    /** {@inherit_doc} **/
    //override
    public double log_probability(const int& x) 
    {
        if (number_of_trials == 0) 
        {
            return (x == 0) ? 0. : -INFINITY;
        }
        double ret;
        if (x < 0 || x > number_of_trials) 
        {
            ret = -INFINITY;
        }
else 
        {
            ret = Saddle_Point_Expansion.log_binomial_probability(x, number_of_trials, probability_of_success, 1.0 - probability_of_success);
        }
        return ret;
    }

    /** {@inherit_doc} */
    //override
    public double cumulative_probability(const int& x) 
    {
        double ret;
        if (x < 0) 
        {
            ret = 0.0;
        }
else if (x >= number_of_trials) 
        {
            ret = 1.0;
        }
else 
        {
            ret = 1.0 - Beta.regularized_beta(probability_of_success, x + 1.0, number_of_trials - x);
        }
        return ret;
    }

    /**
     * {@inherit_doc}
     *
     * For {@code n} trials and probability parameter {@code p}, the mean is
     * {@code n * p}.
     */
    //override
    public double get_numerical_mean() const 
    {
        return number_of_trials * probability_of_success;
    }

    /**
     * {@inherit_doc}
     *
     * For {@code n} trials and probability parameter {@code p}, the variance is
     * {@code n * p * (1 - p)}.
     */
    //override
    public double get_numerical_variance() const 
    {
        const double p = probability_of_success;
        return number_of_trials * p * (1 - p);
    }

    /**
     * {@inherit_doc}
     *
     * The lower bound of the support is always 0 except for the probability
     * parameter {@code p = 1}.
     *
     * @return lower bound of the support (0 or the number of trials)
     */
    //override
    public int get_support_lower_bound() 
    {
        return probability_of_success < 1.0 ? 0 : number_of_trials;
    }

    /**
     * {@inherit_doc}
     *
     * The upper bound of the support is the number of trials except for the
     * probability parameter {@code p = 0}.
     *
     * @return upper bound of the support (number of trials or 0)
     */
    //override
    public int get_support_upper_bound() 
    {
        return probability_of_success > 0.0 ? number_of_trials : 0;
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


