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

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.complex.std::complex<double>;
//import org.hipparchus.complex.Field_Complex<double>;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
#include <vector>
#include "../../../CalculusFieldElement.hpp"

/** Duplication algorithm for Carlson symmetric forms.
 * <p>
 * The algorithms are described in B. C. Carlson 1995 paper
 * "Numerical computation of real or complex elliptic integrals", with
 * improvements described in the appendix of B. C. Carlson and James FitzSimons
 * 2000 paper "Reduction theorems for elliptic integrands with the square root
 * of two quadratic factors". They are also described in
 * <a href="https://dlmf.nist.gov/19.36#i">section 19.36(i)</a>
 * of Digital Library of Mathematical Functions.
 * </p>
 * @param <T> type of the field elements (really {@link std::complex<double>} or {@link Field_Complex<double>})
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Duplication 
{
private:
    /** Max number of iterations. */
    static constexpr int M_MAX{ 16 };

    /** Symmetric variables of the integral, plus mean point. */
    const std::vector<T> my_initial_va;

    /** Convergence criterion. */
    const double my_q;

protected:
    /** Get the i<sup>th</sup> symmetric variable.
     * @param i index of the variable
     * @return i<sup>th</sup> symmetric variable
     */
    T get_vi(const int& i)
    {
        return initial_va[i];
    }

    /** Compute initial mean point.
     * <p>
     * The initial mean point is put as the last array element
     * </>
     * @param va symmetric variables of the integral (plus placeholder for initial mean point)
     */
    virtual void initial_mean_point(std::vector<T> va);

    /** Compute convergence criterion.
     * @param r relative tolerance
     * @param max max(|a0-v[i]|)
     * @return convergence criterion
     */
    virtual T convergence_criterion(T r, T max);

    /** Update reduced variables in place.
     * <ul>
     *  <li>vₘ₊₁|i] ← (vₘ[i] + λₘ) / 4</li>
     *  <li>aₘ₊₁ ← (aₘ + λₘ) / 4</li>
     * </ul>
     * @param m iteration index
     * @param va_m reduced variables and mean point (updated in place)
     * @param sqrt_m square roots of reduced variables
     * @param four_m 4<sup>m</sup>
     */
    virtual void update(const int& m, std::vector<T> va_m, std::vector<T> sqrt_m, double four_m);

    /** Evaluate integral.
     * @param va0 initial symmetric variables and mean point of the integral
     * @param a_m reduced mean point
     * @param four_m 4<sup>m</sup>
     * @return convergence criterion
     */
    virtual T evaluate(std::vector<T> va0, T a_m, double four_m);


public:
    /** Constructor.
     * @param v symmetric variables of the integral
     */
   //@Safe_Varargs
    Field_Duplication(const T... v) 
    {
        const Field<T> field = v[0].get_field();
        const int n = v.size();
        initial_va = Math_Arrays::build_array(field, n + 1);
        System.arraycopy(v, 0, initial_va, 0, n);
        initial_mean_point(initial_va);

        T max = field.get_zero();
        const T a0 = initial_va[n];
        for (const T vi : v) 
        {
            max = std::max(max, a0.subtract(vi).abs());
        }
        this.q = convergence_criterion(FastMath.ulp(field.get_one()), max).get_real();

    }


    /** Compute Carlson elliptic integral.
     * @return Carlson elliptic integral
     */
    T integral() 
    {

        // duplication iterations
        const int& n     = initial_va.size() - 1;
        const std::vector<T> va_m   = initial_va.clone();
        const std::vector<T> sqrt_m = Math_Arrays::build_array(initial_va[0].get_field(), n);
        double    four_m = 1.0;
        for (const int& m = 0; m < M_MAX; ++m) 
        {

            if (m > 0 && q < four_m * va_m[n].norm()) 
            {
                // convergence reached
                break;
            }

            // apply duplication once more
            // (we know that {Field}std::complex<double>.sqrt() returns the root with nonnegative real part)
            for (int i{}; i < n; ++i) 
            {
                sqrt_m[i] = va_m[i].sqrt();
            }
            update(m, va_m, sqrt_m, four_m);

            four_m *= 4;

        }

        return evaluate(initial_va, va_m[n], four_m);
    }
};