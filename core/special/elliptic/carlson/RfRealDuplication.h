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
#include <vector>
#include <cmath>
#include <algorithm>

/** Duplication algorithm for Carlson R<sub>F</sub> elliptic integral.
 * @since 2.0
 */
class Rf_Real_Duplication //extends Real_Duplication 
{
private:
    /** Max number of iterations in the AGM scale. */
    static constexpr int AGM_MAX{ 32 };

    /** Constant term in R<sub>F</sub> polynomial. */
    static constexpr double CONSTANT{ 240240 };

    /** Coefficient of E₂ in R<sub>F</sub> polynomial. */
    static constexpr double E2{ -24024 };

    /** Coefficient of E₃ in R<sub>F</sub> polynomial. */
    static constexpr double E3{ 17160 };

    /** Coefficient of E₂² in R<sub>F</sub> polynomial. */
    static constexpr double E2_E2{ 10010 };

    /** Coefficient of E₂E₃ in R<sub>F</sub> polynomial. */
    static constexpr double E2_E3{ -16380 };

    /** Coefficient of E₃² in R<sub>F</sub> polynomial. */
    static constexpr double E3_E3{ 6930 };

    /** Coefficient of E₂³ in R<sub>F</sub> polynomial. */
    static constexpr double E2_E2_E2{ -5775 };

    /** Denominator in R<sub>F</sub> polynomial. */
    static constexpr double DENOMINATOR{ 240240 };

    /** Compute Carlson complete elliptic integral R<sub>F</sub>(u, v, 0).
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @return Carlson complete elliptic integral R<sub>F</sub>(u, v, 0)
     */
    double complete_integral(const double& x, const double& y)
    {

        double x_m = std::sqrt(x);
        double yM = std::sqrt(y);

        // iterate down
        for (int i{ 1 }; i < AGM_MAX; ++i)
        {

            const double x_m1 = x_m;
            const double y_m1 = yM;

            // arithmetic mean
            x_m = (x_m1 + y_m1) * 0.5;

            // geometric mean
            yM = std::sqrt(x_m1 * y_m1);

            // convergence (by the inequality of arithmetic and geometric means, this is non-negative)
            if (std::abs(x_m - yM) <= 4 * FastMath.ulp(x_m))
            {
                // convergence has been reached
                break;
            }

        }

        return std::numbers::pi / (x_m + yM);

    }

public:
    /** Simple constructor.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     */
    Rf_Real_Duplication(const double& x, const double& y, const double& z) 
    {
        //super(x, y, z);
    }

protected:
    /** {@inherit_doc} */
    //override
    void initial_mean_point(std::vector<double>& va) 
    {
        va[3] = (va[0] + va[1] + va[2]) / 3.0;
    }

    /** {@inherit_doc} */
    //override
    double convergence_criterion(const double r, const double max) 
    {
        return max / std::sqrt(std::sqrt(std::sqrt(r * 3.0)));
    }

    /** {@inherit_doc} */
    //override
    void update(const int& m, std::vector<double>& va_m, const std::vector<double> sqrt_m, const  double four_m) 
    {

        // equation 2.3 in Carlson[1995]
        const auto lambda_a = sqrt_m[0] * sqrt_m[1];
        const auto lambda_b = sqrt_m[0] * sqrt_m[2];
        const auto lambda_c = sqrt_m[1] * sqrt_m[2];

        // equations 2.3 and 2.4 in Carlson[1995]
        va_m[0] = Math_Arrays::linear_combination(0.25, va_m[0], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // xₘ
        va_m[1] = Math_Arrays::linear_combination(0.25, va_m[1], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // yₘ
        va_m[2] = Math_Arrays::linear_combination(0.25, va_m[2], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // zₘ
        va_m[3] = Math_Arrays::linear_combination(0.25, va_m[3], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // aₘ

    }

    /** {@inherit_doc} */
    //override
    double evaluate(const std::vector<double>& va0, const double& a_m, const  double four_m) 
    {
        // compute symmetric differences
        const double inv{ 1.0 / (a_m * four_m) };
        const double big_x{ (va0[3] - va0[0]) * inv };
        const double big_y{ (va0[3] - va0[1]) * inv };
        const double big_z{ -(big_x + big_y) };

        // compute elementary symmetric functions (we already know e1 = 0 by construction)
        const double e2{ big_x * big_y - big_z * big_z };
        const double e3{ big_x * big_y * big_z };

        const double e2e2{ e2 * e2 };
        const double e2e3{ e2 * e3 };
        const double e3e3{ e3 * e3 };
        const double e2e2e2{ e2e2 * e2 };

        // evaluate integral using equation 19.36.1 in DLMF
        // (which add more terms than equation 2.7 in Carlson[1995])
        const double poly = (e2e2e2 * E2_E2_E2 +
                             e3e3   * E3_E3 +
                             e2e3   * E2_E3 +
                             e2e2   * E2_E2 +
                             e3     * E3 +
                             e2     * E2 +
                             CONSTANT) /
                        DENOMINATOR;
        return poly / std::sqrt(a_m);

    }

public:
    /** {@inherit_doc} */
    //override
    double integral() 
    {
        const auto x{ get_vi(0) };
        const auto y{ get_vi(1) };
        const auto z{ get_vi(2) };
        if (x == 0) 
        {
            return complete_integral(y, z);
        }
        if (y == 0) 
        {
            return complete_integral(x, z);
        }
        if (z == 0) 
        {
            return complete_integral(x, y);
        }
        return super.integral();
    }
};