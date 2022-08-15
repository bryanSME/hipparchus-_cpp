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

#include <cmath>
#include <vector>
#include "../../../CalculusFieldElement.hpp"


/** Duplication algorithm for Carlson R<sub>C</sub> elliptic integral.
 * @param <T> type of the field elements (really {@link std::complex<double>} or {@link Field_Complex<double>})
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Rc_Field_Duplication : public Field_Duplication<T>
{

    /** Simple constructor.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     */
    Rc_Field_Duplication(const T& x, const T& y)
    {
        super(x, y);
    }

protected:
    /** {@inherit_doc} */
    //override
    void initial_mean_point(const std::vector<T>& va)
    {
        va[2] = va[0].add(va[1].multiply(2)).divide(3.0);
    }

    /** {@inherit_doc} */
    //override
    T convergence_criterion(const T& r, const T& max)
    {
        return max.divide(std::sqrt(std::sqrt(std::sqrt(r.multiply(3.0)))));
    }

    /** {@inherit_doc} */
    //override
    void update(const int& m, const std::vector<T>& va_m, const std::vector<T>& sqrt_m, const double& four_m)
    {
        const T lambda_a = sqrt_m[0].multiply(sqrt_m[1]).multiply(2);
        const T lambda_b = va_m[1];
        va_m[0] = va_m[0].linear_combination(0.25, va_m[0], 0.25, lambda_a, 0.25, lambda_b); // xₘ
        va_m[1] = va_m[1].linear_combination(0.25, va_m[1], 0.25, lambda_a, 0.25, lambda_b); // yₘ
        va_m[2] = va_m[2].linear_combination(0.25, va_m[2], 0.25, lambda_a, 0.25, lambda_b); // aₘ
    }

    /** {@inherit_doc} */
    //override
    T evaluate(const std::vector<T>& va0, const T& a_m, const double& four_m)
    {

        // compute the single polynomial independent variable
        const T s = va0[1].subtract(va0[2]).divide(a_m.multiply(four_m));

        // evaluate integral using equation 2.13 in Carlson[1995]
        const T poly = s.multiply(Rc_Real_Duplication.S7).
            add(Rc_Real_Duplication.S6).multiply(s).
            add(Rc_Real_Duplication.S5).multiply(s).
            add(Rc_Real_Duplication.S4).multiply(s).
            add(Rc_Real_Duplication.S3).multiply(s).
            add(Rc_Real_Duplication.S2).multiply(s).
            multiply(s).
            add(Rc_Real_Duplication.S0).
            divide(Rc_Real_Duplication.DENOMINATOR);
        return poly.divide(std::sqrt(a_m));
    }
};