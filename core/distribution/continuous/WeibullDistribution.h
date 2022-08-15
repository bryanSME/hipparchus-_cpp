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

//import org.hipparchus.exception.;
#include "AbstractRealDistribution.h"
#include "../../special/Gamma.h"
#include "../../util/MathUtils.h"
#include "../../exception/LocalizedCoreFormats.h"

/**
 * Implementation of the Weibull distribution. This implementation uses the
 * two parameter form of the distribution defined by
 * <a href="http://mathworld.wolfram.com/Weibull_Distribution.html">
 * Weibull Distribution</a>, equations (1) and (2).
 *
 * @see <a href="http://en.wikipedia.org/wiki/Weibull_distribution">Weibull distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/Weibull_Distribution.html">Weibull distribution (MathWorld)</a>
 */
class Weibull_Distribution : public Abstract_Real_Distribution 
{
private:
    /** The shape parameter. */
    const double my_shape{};
    /** The scale parameter. */
    const double my_scale{};

public:
    /**
     * Create a Weibull distribution with the given shape and scale.
     *
     * @param alpha Shape parameter.
     * @param beta Scale parameter.
     * @ if {@code alpha <= 0} or {@code beta <= 0}.
     */
    Weibull_Distribution(const double& alpha, const double& beta) : my_scale{ beta }, my_shape{ alpha }
    {
        if (alpha <= 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::SHAPE, alpha);
        }
        if (beta <= 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::SCALE, beta);
        }
    }

    /**
     * Access the shape parameter, {@code alpha}.
     *
     * @return the shape parameter, {@code alpha}.
     */
    double get_shape() const
    {
        return my_shape;
    }

    /**
     * Access the scale parameter, {@code beta}.
     *
     * @return the scale parameter, {@code beta}.
     */
    double get_scale() const
    {
        return my_scale;
    }

    /** {@inherit_doc} */
    //override
    double density(const double& x) 
    {
        if (x < 0) 
        {
            return 0;
        }

        const double xscale = x / my_scale;
        const double xscalepow = std::pow(xscale, my_shape - 1);

        /*
         * std::pow(x / scale, shape) =
         * std::pow(xscale, shape) =
         * std::pow(xscale, shape - 1) * xscale
         */
        const double xscalepowshape = xscalepow * xscale;

        return (my_shape / my_scale) * xscalepow * std::exp(-xscalepowshape);
    }

    /** {@inherit_doc} */
    //override
    double log_density(double x) 
    {
        if (x < 0) 
        {
            return -INFINITY;
        }

        const double xscale = x / my_scale;
        const double logxscalepow = std::log(xscale) * (my_shape - 1);

        /*
         * std::pow(x / scale, shape) =
         * std::pow(xscale, shape) =
         * std::pow(xscale, shape - 1) * xscale
         */
        const double xscalepowshape = std::exp(logxscalepow) * xscale;

        return std::log(my_shape / my_scale) + logxscalepow - xscalepowshape;
    }

    /** {@inherit_doc} */
    //override
    double cumulative_probability(const double& x) 
    {
        return x <= 0.0
            ? 0.0
            : 1.0 - std::exp(-std::pow(x / my_scale, my_shape));
    }

    /**
     * {@inherit_doc}
     *
     * Returns {@code 0} when {@code p == 0} and
     * {@code INFINITY} when {@code p == 1}.
     */
    //override
    double inverse_cumulative_probability(const double& p) 
    {
        Math_Utils::check_range_inclusive(p, 0, 1);

        if (p == 0) 
        {
            return 0.0;
        }
        if (p == 1) 
        {
            return INFINITY;
        }
        return my_scale * std::pow(-std::log1p(-p), 1.0 / my_shape);
    }

    /**
     * {@inherit_doc}
     *
     * The mean is {@code scale * Gamma(1 + (1 / shape))}, where {@code Gamma()}
     * is the Gamma-function.
     */
    //override
    double get_numerical_mean() const 
    {
        const double sh = get_shape();
        const double sc = get_scale();

        return sc * std::exp(Gamma::log_gamma(1 + (1 / sh)));
    }

    /**
     * {@inherit_doc}
     *
     * The variance is {@code scale^2 * Gamma(1 + (2 / shape)) - mean^2}
     * where {@code Gamma()} is the Gamma-function.
     */
    //override
    double get_numerical_variance() const 
    {
        const double sh = get_shape();
        const double sc = get_scale();
        const double mn = get_numerical_mean();

        return (sc * sc) * std::exp(Gamma::log_gamma(1 + (2 / sh))) - (mn * mn);
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