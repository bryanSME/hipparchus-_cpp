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
//import org.hipparchus.special.Gamma;
//import org.hipparchus.util.FastMath;
#include <cmath>
/**
 * This class : the Nakagami distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Nakagami_distribution">Nakagami Distribution (Wikipedia)</a>
 */
class Nakagami_Distribution extends Abstract_Real_Distribution 
{

    
    20141003;

    /** The shape parameter. */
    private const double mu;
    /** The scale parameter. */
    private const double omega;

    /**
     * Build a instance.
     *
     * @param mu shape parameter
     * @param omega scale parameter (must be positive)
     * @ if {@code mu < 0.5}
     * @ if {@code omega <= 0}
     */
    public Nakagami_Distribution(double mu, double omega)
         
        {
        this(mu, omega, DEFAULT_SOLVER_ABSOLUTE_ACCURACY);
    }

    /**
     * Build a instance.
     *
     * @param mu shape parameter
     * @param omega scale parameter (must be positive)
     * @param inverse_absolute_accuracy the maximum absolute error in inverse
     * cumulative probability estimates (defaults to {@link #DEFAULT_SOLVER_ABSOLUTE_ACCURACY}).
     * @ if {@code mu < 0.5}
     * @ if {@code omega <= 0}
     */
    public Nakagami_Distribution(double mu, double omega, double inverse_absolute_accuracy)
         
        {
        super(inverse_absolute_accuracy);

        if (mu < 0.5) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, mu, 0.5);
        }
        if (omega <= 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_SCALE, omega);
        }

        this.mu = mu;
        this.omega = omega;
    }

    /**
     * Access the shape parameter, {@code mu}.
     *
     * @return the shape parameter.
     */
    public double get_shape() 
    {
        return mu;
    }

    /**
     * Access the scale parameter, {@code omega}.
     *
     * @return the scale parameter.
     */
    public double get_scale() 
    {
        return omega;
    }

    /** {@inherit_doc} */
    //override
    public double density(double x) 
    {
        if (x <= 0) 
        {
            return 0.0;
        }
        return 2.0 * std::pow(mu, mu) / (Gamma::gamma(mu) * std::pow(omega, mu)) *
                     std::pow(x, 2 * mu - 1) * std::exp(-mu * x * x / omega);
    }

    /** {@inherit_doc} */
    //override
    public double cumulative_probability(const double& x) 
    {
        return Gamma::regularized_gamma_p(mu, mu * x * x / omega);
    }

    /** {@inherit_doc} */
    //override
    public double get_numerical_mean() const 
    {
        return Gamma::gamma(mu + 0.5) / Gamma::gamma(mu) * std::sqrt(omega / mu);
    }

    /** {@inherit_doc} */
    //override
    public double get_numerical_variance() const 
    {
        double v = Gamma::gamma(mu + 0.5) / Gamma::gamma(mu);
        return omega * (1 - 1 / mu * v * v);
    }

    /** {@inherit_doc} */
    //override
    public double get_support_lower_bound() const 
    {
        return 0;
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


