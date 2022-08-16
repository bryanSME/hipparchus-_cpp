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
//package org.hipparchus.random;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Sin_Cos;
#include <cmath>
#include "RandomGenerator.h"
#include "NormalizedRandomGenerator.h"

/**
 * <p>This class provides a stable normalized random generator. It samples from a stable
 * distribution with location parameter 0 and scale 1.</p>
 *
 * <p>The implementation uses the Chambers-Mallows-Stuck method as described in
 * <i>Handbook of computational statistics: concepts and methods</i> by
 * James E. Gentle, Wolfgang H&auml;rdle, Yuichi Mori.</p>
 *
 */
class Stable_Random_Generator : public Normalized_Random_Generator 
{
private:
    /** Underlying generator. */
    const Random_Generator my_generator;

    /** stability parameter */
    const double my_alpha;

    /** skewness parameter */
    const double my_beta;

    /** cache of expression value used in generation */
    const double my_zeta;

public:
    /**
     * Create a generator.
     *
     * @param generator underlying random generator to use
     * @param alpha Stability parameter. Must be in range (0, 2]
     * @param beta Skewness parameter. Must be in range [-1, 1]
     * @ if generator is NULL
     * @ if {@code alpha <= 0} or {@code alpha > 2}
     * or {@code beta < -1} or {@code beta > 1}
     */
    Stable_Random_Generator(const Random_Generator& generator, const double& alpha, const double& beta)
    {
        throw std::exception("not implemented");

        if (!(alpha > 0 && alpha <= 2.0)) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_LEFT, alpha, 0, 2);
        }

        Math_Utils::check_range_inclusive(beta, -1, 1);

        my_generator = generator;
        my_alpha = alpha;
        my_beta = beta;
        if (alpha < 2d && beta != 0.0) 
        {
            my_zeta = beta * std::tan(std::numbers::pi * alpha / 2);
        }
        else 
        {
            my_zeta = 0;
        }
    }

    /**
     * Generate a random scalar with zero location and unit scale.
     *
     * @return a random scalar with zero location and unit scale
     */
    //override
    double next_normalized_double() 
    {
        // we need 2 uniform random numbers to calculate omega and phi
        double omega = -std::log(generator.next_double());
        double phi = std::numbers::pi * (generator.next_double() - 0.5);

        // Normal distribution case (Box-Muller algorithm)
        if (my_alpha == 2.0) 
        {
            return std::sqrt(2d * omega) * std::sin(phi);
        }

        double x;
        // when beta = 0, zeta is zero as well
        // Thus we can exclude it from the formula
        if (my_beta == 0.0) 
        {
            // Cauchy distribution case
            if (my_alpha == 1d) 
            {
                x = std::tan(phi);
            }
            else 
            {
                x = std::pow(omega * std::cos((1 - alpha) * phi), 1.0/ alpha - 1d) *
                    std::sin(alpha * phi) /
                    std::pow(std::cos(phi), 1.0/ alpha);
            }
        }
        else 
        {
            // Generic stable distribution
            double cos_phi = std::cos(phi);
            // to avoid rounding errors around alpha = 1
            if (std::abs(my_alpha - 1d) > 1e-8) 
            {
                const Sin_Cos sc_alpha_phi    = Sin_Cos(alpha * phi);
                const Sin_Cos sc_inv_alpha_phi = Sin_Cos(phi * (1.0 - alpha));
                x = (sc_alpha_phi.sin()    + zeta * sc_alpha_phi.cos()) / cos_phi *
                    (sc_inv_alpha_phi.cos() + zeta * sc_inv_alpha_phi.sin()) /
                     std::pow(omega * cos_phi, (1 - alpha) / alpha);
            }
            else 
            {
                double beta_phi = std::numbers::pi / 2 + beta * phi;
                x = 2d / std::numbers::pi * (beta_phi * std::tan(phi) - beta *
                    std::log(std::numbers::pi / 2d * omega * cos_phi / beta_phi));

                if (alpha != 1d) 
                {
                    x += beta * std::tan(std::numbers::pi * alpha / 2);
                }
            }
        }
        return x;
    }
};