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

#include <cmath>

//import org.hipparchus.exception.;
//import org.hipparchus.special.Erf;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/**
 * This class : the <a href="http://en.wikipedia.org/wiki/L%C3%A9vy_distribution">
 * L&eacute;vy distribution</a>.
 */
class Levy_Distribution : Abstract_Real_Distribution 
{
private:

    /** Location parameter. */
    const double my_mu;

    /** Scale parameter. */
    const double my_c;  // Setting this to 1 returns a cum_prob of 1.0

    /** Half of c (for calculations). */
    const double my_half_c;

public:
    /**
     * Build a instance.
     *
     * @param mu location parameter
     * @param c scale parameter
     */
    Levy_Distribution(const double& mu, const double& c) : my_mu{ mu }, my_c{ c }, my_half_c{ 0.5 * c }
    {
        super();
    }


    /** {@inherit_doc}
    * <p>
    * From Wikipedia: The probability density function of the L&eacute;vy distribution
    * over the domain is
    * </p>
    * <pre>
    * f(x; &mu;, c) = &radic;(c / 2&pi;) * e<sup>-c / 2 (x - &mu;)</sup> / (x - &mu;)<sup>3/2</sup>
    * </pre>
    * <p>
    * For this distribution, {@code X}, this method returns {@code P(X < x)}.
    * If {@code x} is less than location parameter &mu;, {@codeNAN} is
    * returned, as in these cases the distribution is not defined.
    * </p>
    */
    //override
    double density(const double& x) 
    {
        if (x < my_mu) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }

        const double delta = x - my_mu;
        const double f     = my_half_c / delta;
        return std::sqrt(f / std::numbers::pi) * std::exp(-f) /delta;
    }

    /** {@inherit_doc}
     *
     * See documentation of {@link #densitystatic_cast<double>(} for computation details.
     */
    //override
    double log_density(double x) 
    {
        if (x < my_mu) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }

        const double delta = x - my_mu;
        const double f     = half_c / delta;
        return 0.5 * std::log(f / std::numbers::pi) - f - std::log(delta);
    }

    /** {@inherit_doc}
     * <p>
     * From Wikipedia: the cumulative distribution function is
     * </p>
     * <pre>
     * f(x; u, c) = erfc (&radic; (c / 2 (x - u )))
     * </pre>
     */
    //override
    double cumulative_probability(const double& x) 
    {
        if (x < my_mu) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return Erf.erfc(std::sqrt(half_c / (x - my_mu)));
    }

    /** {@inherit_doc} */
    //override
    double inverse_cumulative_probability(const double& p)  
    {
        Math_Utils::check_range_inclusive(p, 0, 1);
        const double t = Erf.erfc_inv(p);
        return my_mu + my_half_c / (t * t);
    }

    /** Get the scale parameter of the distribution.
     * @return scale parameter of the distribution
     */
    double get_scale() const
    {
        return my_c;
    }

    /** Get the location parameter of the distribution.
     * @return location parameter of the distribution
     */
    public double get_location() const
    {
        return my_mu;
    }

    /** {@inherit_doc} */
    //override
    public double get_numerical_mean() const 
    {
        return INFINITY;
    }

    /** {@inherit_doc} */
    //override
    public double get_numerical_variance() const 
    {
        return INFINITY;
    }

    /** {@inherit_doc} */
    //override
    public double get_support_lower_bound() const 
    {
        return my_mu;
    }

    /** {@inherit_doc} */
    //override
    public double get_support_upper_bound() const 
    {
        return INFINITY;
    }

    /** {@inherit_doc} */
    //override
    public bool is_support_connected() const 
    {
        return true;
    }

}


