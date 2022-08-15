#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
//package org.hipparchus.special.elliptic.carlson;

//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
#include <vector>

/** Duplication algorithm for Carlson R<sub>J</sub> elliptic integral.
 * @since 2.0
 */
class Rj_Real_Duplication : public Real_Duplication 
{

private:
    /** Delta product. */
    double my_delta;

    /** sₘ iteration parameter. */
    double my_s_m;

    /** Simple constructor.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param delta precomputed value of (p-x)(p-y)(p-z)
     */
    Rj_Real_Duplication(const double& x, const double& y, const double z, const double p, const double delta) 
    {
        super(x, y, z, p);
        this.delta = delta;
    }

    /** {@inherit_doc} */
    //override
    protected void initial_mean_point(const std::vector<double>& va) 
    {
        va[4] = (va[0] + va[1] + va[2] + va[3] * 2) / 5.0;
    }

    /** {@inherit_doc} */
    //override
    protected double convergence_criterion(const double r, const double max) 
    {
        return max / (std::sqrt(std::sqrt(std::sqrt(r * 0.25))));
    }

    /** {@inherit_doc} */
    //override
    protected void update(const int m, const std::vector<double>& va_m, const std::vector<double> sqrt_m, const  double four_m) 
    {
        const double dM = (sqrt_m[3] + sqrt_m[0]) * (sqrt_m[3] + sqrt_m[1]) * (sqrt_m[3] + sqrt_m[2]);
        if (m == 0) 
        {
            s_m = dM * 0.5;
        }
else 
        {
            // equation A.3 in Carlson[2000]
            const double rM = s_m * (std::sqrt(delta / (s_m * s_m * four_m) + 1.0) + 1.0);
            s_m = (dM * rM - delta / (four_m * four_m)) / ((dM + rM / four_m) * 2);
        }

        // equation 2.19 in Carlson[1995]
        const double lambda_a = sqrt_m[0] * sqrt_m[1];
        const double lambda_b = sqrt_m[0] * sqrt_m[2];
        const double lambda_c = sqrt_m[1] * sqrt_m[2];

        // equations 2.19 and 2.20 in Carlson[1995]
        va_m[0] = Math_Arrays::linear_combination(0.25, va_m[0], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // xₘ
        va_m[1] = Math_Arrays::linear_combination(0.25, va_m[1], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // yₘ
        va_m[2] = Math_Arrays::linear_combination(0.25, va_m[2], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // zₘ
        va_m[3] = Math_Arrays::linear_combination(0.25, va_m[3], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // pₘ
        va_m[4] = Math_Arrays::linear_combination(0.25, va_m[4], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // aₘ

    }

    /** {@inherit_doc} */
    //override
    protected double evaluate(const std::vector<double>& va0, const double& a_m, const  double four_m) 
    {

        // compute symmetric differences
        const double inv    = 1.0 / (a_m * four_m);
        const double big_x   = (va0[4] - va0[0]) * inv;
        const double big_y   = (va0[4] - va0[1]) * inv;
        const double big_z   = (va0[4] - va0[2]) * inv;
        const double big_p   = (big_x + big_y + big_z) * -0.5;
        const double big_p2  = big_p * big_p;

        // compute elementary symmetric functions (we already know e1 = 0 by construction)
        const double xyz    = big_x * big_y * big_z;
        const double e2     = big_x * (big_y + big_z) + big_y * big_z - big_p * big_p * 3;
        const double e3     = xyz + big_p * 2 * (e2 + big_p2 * 2);
        const double e4     = (xyz * 2 + big_p * (e2 + big_p2 * 3)) * big_p;
        const double e5     = xyz * big_p2;

        const double e2e2   = e2   * e2;
        const double e2e3   = e2   * e3;
        const double e2e4   = e2   * e4;
        const double e2e5   = e2   * e5;
        const double e3e3   = e3   * e3;
        const double e3e4   = e3   * e4;
        const double e2e2e2 = e2e2 * e2;
        const double e2e2e3 = e2e2 * e3;

        // evaluate integral using equation 19.36.1 in DLMF
        // (which add more terms than equation 2.7 in Carlson[1995])
        const double poly = ((e3e4 + e2e5) * Rd_Real_Duplication.E3_E4_P_E2_E5 +
                              e2e2e3       * Rd_Real_Duplication.E2_E2_E3 +
                              e2e4         * Rd_Real_Duplication.E2_E4 +
                              e3e3         * Rd_Real_Duplication.E3_E3 +
                              e2e2e2       * Rd_Real_Duplication.E2_E2_E2 +
                              e5           * Rd_Real_Duplication.E5 +
                              e2e3         * Rd_Real_Duplication.E2_E3 +
                              e4           * Rd_Real_Duplication.E4 +
                              e2e2         * Rd_Real_Duplication.E2_E2 +
                              e3           * Rd_Real_Duplication.E3 +
                              e2           * Rd_Real_Duplication.E2 +
                              Rd_Real_Duplication.CONSTANT) /
                             Rd_Real_Duplication.DENOMINATOR;
        const double poly_term = poly / (a_m * std::sqrt(a_m) * four_m);

        // compute a single R_C term
        const double rc_term = Rc_Real_Duplication(1.0, delta / (s_m * s_m * four_m) + 1.0).integral() * 3 / s_m;

        return poly_term + rc_term;

    }

}


