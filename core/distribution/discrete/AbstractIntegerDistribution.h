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

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.Math_Runtime_Exception;
#include <cmath>
#include <exception>
#include <limits>
#include "../../distribution/IntegerDistribution.h"
#include "../../util/MathUtils.h"

/**
 * Base class for integer-valued discrete distributions.
 * <p>
 * Default implementations are provided for some of the methods that
 * do not vary from distribution to distribution.
 */
class Abstract_Integer_Distribution : public Integer_Distribution
{
public:
    /**
     * {@inherit_doc}
     *
     * The default implementation uses the identity
     * <p>
     * {@code P(x0 < X <= x1) = P(X <= x1) - P(X <= x0)}
     */
    //override
    double probability(const int& x0, const int& x1)  
    {
        if (x1 < x0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT, x0, x1, true);
        }
        return cumulative_probability(x1) - cumulative_probability(x0);
    }

    /**
     * {@inherit_doc}
     *
     * The default implementation returns
     * <ul>
     * <li>{@link #get_support_lower_bound()} for {@code p = 0},</li>
     * <li>{@link #get_support_upper_bound()} for {@code p = 1}, and</li>
     * <li>{@link #solve_inverse_cumulative_probability(double, int, int)} for
     *     {@code 0 < p < 1}.</li>
     * </ul>
     */
    //override
    int inverse_cumulative_probability(const double& p)  
    {
        Math_Utils::check_range_inclusive(p, 0, 1);

        int lower = get_support_lower_bound();
        if (p == 0.0) 
        {
            return lower;
        }
        if (lower == std::numeric_limits<int>::min()) 
        {
            if (checked_cumulative_probability(lower) >= p) 
            {
                return lower;
            }
        }
        else 
        {
            lower -= 1; // this ensures cumulative_probability(lower) < p, which
                        // is important for the solving step
        }

        int upper = get_support_upper_bound();
        if (p == 1.0) 
        {
            return upper;
        }

        // use the one-sided Chebyshev inequality to narrow the bracket
        // cf. Abstract_Real_Distribution.inverse_cumulative_probabilitystatic_cast<double>(
        const double mu = get_numerical_mean();
        const double sigma = std::sqrt(get_numerical_variance());
        const bool chebyshev_applies =
                !(std::isinf(mu)    || std::isnan(mu)    ||
                  std::isinf(sigma) || std::isnan(sigma) ||
                  sigma == 0.0);
        if (chebyshev_applies) 
        {
            double k = std::sqrt((1.0 - p) / p);
            double tmp = mu - k * sigma;
            if (tmp > lower) 
            {
                lower = (static_cast<int>(std::ceil(tmp))) - 1;
            }
            k = 1.0 / k;
            tmp = mu + k * sigma;
            if (tmp < upper) 
            {
                upper = (static_cast<int>(std::ceil(tmp))) - 1;
            }
        }

        return solve_inverse_cumulative_probability(p, lower, upper);
    }

    /**
     * {@inherit_doc}
     * <p>
     * The default implementation simply computes the logarithm of {@code probability(x)}.
     */
    //override
    double log_probability(const int& x)
    {
        return std::log(probability(x));
    }

protected:
    /**
     * This is a utility function used by {@link
     * #inverse_cumulative_probabilitystatic_cast<double>(}. It assumes {@code 0 < p < 1} and
     * that the inverse cumulative probability lies in the bracket {@code
     * (lower, upper]}. The implementation does simple bisection to find the
     * smallest {@code p}-quantile {@code inf{x in Z | P(X<=x) >= p}}.
     *
     * @param p the cumulative probability
     * @param lower a value satisfying {@code cumulative_probability(lower) < p}
     * @param upper a value satisfying {@code p <= cumulative_probability(upper)}
     * @return the smallest {@code p}-quantile of this distribution
     */
    int solve_inverse_cumulative_probability(const double& p, int lower, int upper) 
    {
        while (lower + 1 < upper) 
        {
            int xm = (lower + upper) / 2;
            if (xm < lower || xm > upper) 
            {
                /*
                 * Overflow.
                 * There will never be an overflow in both calculation methods
                 * for xm at the same time
                 */
                xm = lower + (upper - lower) / 2;
            }

            const auto pm = checked_cumulative_probability(xm);
            if (pm >= p) 
            {
                upper = xm;
            }
            else 
            {
                lower = xm;
            }
        }
        return upper;
    }

private:

    /**
     * Computes the cumulative probability function and checks for {@code NaN}
     * values returned.
     * <p>
     * Throws {@code Math_Runtime_Exception} if the value is {@code NaN}.
     * Reany exception encountered evaluating the cumulative
     * probability function.
     * Throws {@code Math_Runtime_Exception} if the cumulative
     * probability function returns {@code NaN}.
     *
     * @param argument input value
     * @return the cumulative probability
     * @Math_Runtime_Exception if the cumulative probability is {@code NaN}
     */
    double checked_cumulative_probability(const int& argument)
    {
        double result = cumulative_probability(argument);
        if (std::isnan(result)) 
        {
            throw std::exception("Not implemented - checked_cumulative_probability");
            //throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::DISCRETE_CUMULATIVE_PROBABILITY_RETURNED_NAN, argument);
        }
        return result;
    }
};