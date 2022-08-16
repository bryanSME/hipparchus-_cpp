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
//import org.hipparchus.special.Beta;
//import org.hipparchus.special.Gamma;
//import org.hipparchus.util.FastMath;
#include "AbstractRealDistribution.h"
#include "../../special/Beta.h"
#include "../../special/Gamma.h"

/**
 * Implements the Beta distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Beta_distribution">Beta distribution</a>
 */
class Beta_Distribution : public Abstract_Real_Distribution 
{
private:
    /** First shape parameter. */
    const double my_alpha{};
    /** Second shape parameter. */
    const double my_beta{};
    /** Normalizing factor used in density computations. */
    const double my_z{};

public:
    /**
     * Build a instance.
     *
     * @param alpha First shape parameter (must be positive).
     * @param beta Second shape parameter (must be positive).
     */
    Beta_Distribution(const double& alpha, const double& beta) 
    {
        Beta_Distribution(alpha, beta, DEFAULT_SOLVER_ABSOLUTE_ACCURACY);
    }

    /**
     * Build a instance.
     *
     * @param alpha First shape parameter (must be positive).
     * @param beta Second shape parameter (must be positive).
     * @param inverse_cum_accuracy Maximum absolute error in inverse
     * cumulative probability estimates (defaults to
     * {@link #DEFAULT_SOLVER_ABSOLUTE_ACCURACY}).
     */
    Beta_Distribution(const double& alpha, const double& beta, const double& inverse_cum_accuracy) 
        : 
        my_alpha{ alpha },
        my_beta{ beta },
        my_z{ Gamma::log_gamma(alpha) +
                     Gamma::log_gamma(beta) -
                     Gamma::log_gamma(alpha + beta) }
    {
        super(inverse_cum_accuracy);
    }

    /**
     * Access the first shape parameter, {@code alpha}.
     *
     * @return the first shape parameter.
     */
    double get_alpha() const
    {
        return my_alpha;
    }

    /**
     * Access the second shape parameter, {@code beta}.
     *
     * @return the second shape parameter.
     */
    double get_beta() const
    {
        return my_beta;
    }

    /** {@inherit_doc} **/
    //override
    double log_density(const double& x) const
    {
        if (x < 0 || x > 1) 
        {
            return -INFINITY;
        }
        if (x == 0) 
        {
            if (my_alpha < 1) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::CANNOT_COMPUTE_BETA_DENSITY_AT_0_FOR_SOME_ALPHA, alpha, 1, false);
            }
            return -INFINITY;
        }
        if (x == 1) 
        {
            if (my_beta < 1) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::CANNOT_COMPUTE_BETA_DENSITY_AT_1_FOR_SOME_BETA, beta, 1, false);
            }
            return -INFINITY;
        } 
        {
            double log_x = std::log(x);
            double log1m_x = std::log1p(-x);
            return (my_alpha - 1) * log_x + (my_beta - 1) * log1m_x - my_z;
        }
    }

    /** {@inherit_doc} */
    //override
    double density(const double& x) const
    {
        const double log_density = log_density(x);
        return log_density == -INFINITY
            ? 0
            : std::exp(log_density);
    }

    /** {@inherit_doc} */
    //override
    double cumulative_probability(const double& x) const
    {
        if (x <= 0) 
        {
            return 0;
        }
        if (x >= 1) 
        {
            return 1;
        }
        return Beta::regularized_beta(x, my_alpha, my_beta);
    }

    /**
     * {@inherit_doc}
     *
     * For first shape parameter {@code alpha} and second shape parameter
     * {@code beta}, the mean is {@code alpha / (alpha + beta)}.
     */
    //override
    double get_numerical_mean() const 
    {
        const double& a = get_alpha();
        return a / (a + get_beta());
    }

    /**
     * {@inherit_doc}
     *
     * For first shape parameter {@code alpha} and second shape parameter
     * {@code beta}, the variance is
     * {@code (alpha * beta) / [(alpha + beta)^2 * (alpha + beta + 1)]}.
     */
    //override
    double get_numerical_variance() const 
    {
        const double& a = get_alpha();
        const double b = get_beta();
        const double& alphabetasum = a + b;
        return (a * b) / ((alphabetasum * alphabetasum) * (alphabetasum + 1));
    }

    /**
     * {@inherit_doc}
     *
     * The lower bound of the support is always 0 no matter the parameters.
     *
     * @return lower bound of the support (always 0)
     */
    //override
    double get_support_lower_bound() const 
    {
        return 0;
    }

    /**
     * {@inherit_doc}
     *
     * The upper bound of the support is always 1 no matter the parameters.
     *
     * @return upper bound of the support (always 1)
     */
    //override
    double get_support_upper_bound() const 
    {
        return 1;
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