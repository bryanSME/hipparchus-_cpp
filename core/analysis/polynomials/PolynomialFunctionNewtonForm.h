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
//package org.hipparchus.analysis.polynomials;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.analysis.Field_Univariate_Function;
//import org.hipparchus.analysis.differentiation.Derivative;
//import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Math_Utils;
#include <vector>
#include "../../analysis/FieldUnivariateFunction.h"
#include "../../analysis/differentiation/UnivariateDifferentiableFunction.h"
//#include "../../util/MathUtils.h"
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/**
 * Implements the representation of a real polynomial function in
 * Newton Form. For reference, see <b>Elementary Numerical Analysis</b>, * ISBN 0070124477, chapter 2.
 * <p>
 * The formula of polynomial in Newton form is
 *     p(x) = a[0] + a[1](x-c[0]) + a[2](x-c[0])(x-c[1]) + ... +
 *            a[n](x-c[0])(x-c[1])...(x-c[n-1])
 * Note that the length of a[] is one more than the length of c[]</p>
 *
 */
class Polynomial_FunctionNewtonForm : public Univariate_Differentiable_Function, public Field_Univariate_Function 
{
private:
    /**
     * The coefficients of the polynomial, ordered by degree -- i.e.
     * coefficients[0] is the constant term and coefficients[n] is the
     * coefficient of x^n where n is the degree of the polynomial.
     */
    std::vector<double> my_coefficients;

    /**
     * Centers of the Newton polynomial.
     */
    const std::vector<double> my_c;

    /**
     * When all c[i] = 0, a[] becomes normal polynomial coefficients, * i.e. a[i] = coefficients[i].
     */
    const std::vector<double> my_a;

    /**
     * Whether the polynomial coefficients are available.
     */
    bool my_coefficients_computed;

protected:
    /**
     * Calculate the normal polynomial coefficients given the Newton form.
     * It also uses nested multiplication but takes O(N^2) time.
     */
    void compute_coefficients()
    {
        const int n = degree();

        my_coefficients = std::vector<double>(n + 1);
        for (int i{}; i <= n; i++)
        {
            my_coefficients[i] = 0.0;
        }

        my_coefficients[0] = my_a[n];
        for (int i = n - 1; i >= 0; i--)
        {
            for (int j = n - i; j > 0; j--)
            {
                my_coefficients[j] = my_coefficients[j - 1] - my_c[i] * my_coefficients[j];
            }
            my_coefficients[0] = my_a[i] - my_c[i] * my_coefficients[0];
        }

        my_coefficients_computed = true;
    }

    /**
     * Verifies that the input arrays are valid.
     * <p>
     * The centers must be distinct for interpolation purposes, but not
     * for general use. Thus it is not verified here.</p>
     *
     * @param a the coefficients in Newton form formula
     * @param c the centers
     * @Null_Argument_Exception if any argument is {@code NULL}.
     * @ if any array has zero length.
     * @ if the size difference between
     * {@code a} and {@code c} is not equal to 1.
     * @see org.hipparchus.analysis.interpolation.Divided_Difference_Interpolator#compute_divided_difference(std::vector<double>, * std::vector<double>)
     */
    static void verify_input_array(const std::vector<double>& a, const std::vector<double>& c)
    {
        //Math_Utils::check_not_null(a);
        //Math_Utils::check_not_null(c);
        if (a.empty() || c.empty())
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
        }
        if (a.size() != c.size() + 1)
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::ARRAY_SIZES_SHOULD_HAVE_DIFFERENCE_1, a.size(), c.size());
        }
    }

public:
    /**
     * Construct a Newton polynomial with the given a[] and c[]. The order of
     * centers are important in that if c[] shuffle, then values of a[] would
     * completely change, not just a permutation of old a[].
     * <p>
     * The constructor makes copy of the input arrays and assigns them.</p>
     *
     * @param a Coefficients in Newton form formula.
     * @param c Centers.
     * @Null_Argument_Exception if any argument is {@code NULL}.
     * @ if any array has zero length.
     * @ if the size difference between
     * {@code a} and {@code c} is not equal to 1.
     */
    Polynomial_FunctionNewtonForm(const std::vector<double>& a, const std::vector<double> c) : my_a{ std::vector<double>(a) }, my_c{ std::vector<double>(c) }
    {
        verify_input_array(a, c);
        System.arraycopy(a, 0, my_a, 0, a.size());
        System.arraycopy(c, 0, my_c, 0, c.size());
        my_coefficients_computed = false;
    }

    /**
     * Calculate the function value at the given point.
     *
     * @param z Point at which the function value is to be computed.
     * @return the function value.
     */
    //override
    double value(const double& z) 
    {
       return evaluate(my_a, my_c, z);
    }

    /**
     * {@inherit_doc}
     */
    //override
    <T extends Derivative<T>> T value(const T t) 
    {
        verify_input_array(a, c);

        const int n = c.size();
        T value = t.get_field().get_zero().add(a[n]);
        for (int i = n - 1; i >= 0; i--) 
        {
            value = t.subtract(c[i]).multiply(value).add(a[i]);
        }

        return value;

    }

    /**
     * {@inherit_doc}
     */
    //override
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    T value(const T t) 
    {
        verify_input_array(a, c);

        const int n = c.size();
        T value = t.get_field().get_zero().add(a[n]);
        for (int i = n - 1; i >= 0; i--) 
        {
            value = t.subtract(c[i]).multiply(value).add(a[i]);
        }

        return value;

    }

    /**
     * Returns the degree of the polynomial.
     *
     * @return the degree of the polynomial
     */
    int degree() const 
    {
        return c.size();
    }

    /**
     * Returns a copy of coefficients in Newton form formula.
     * <p>
     * Changes made to the returned copy will not affect the polynomial.</p>
     *
     * @return a fresh copy of coefficients in Newton form formula
     */
    std::vector<double> get_newton_coefficients() 
    {
        auto out = std::vector<double>(a.size());
        System.arraycopy(a, 0, out, 0, a.size());
        return out;
    }

    /**
     * Returns a copy of the centers array.
     * <p>
     * Changes made to the returned copy will not affect the polynomial.</p>
     *
     * @return a fresh copy of the centers array.
     */
    std::vector<double> get_centers() 
    {
        auto out = std::vector<double>(c.size());
        System.arraycopy(c, 0, out, 0, c.size());
        return out;
    }

    /**
     * Returns a copy of the coefficients array.
     * <p>
     * Changes made to the returned copy will not affect the polynomial.</p>
     *
     * @return a fresh copy of the coefficients array.
     */
    std::vector<double> get_coefficients() 
    {
        if (!coefficients_computed) 
        {
            compute_coefficients();
        }
        auto out = std::vector<double>(coefficients.size());
        System.arraycopy(coefficients, 0, out, 0, coefficients.size());
        return out;
    }

    /**
     * Evaluate the Newton polynomial using nested multiplication. It is
     * also called <a href="http://mathworld.wolfram.com/HornersRule.html">
     * Horner's Rule</a> and takes O(N) time.
     *
     * @param a Coefficients in Newton form formula.
     * @param c Centers.
     * @param z Point at which the function value is to be computed.
     * @return the function value.
     * @Null_Argument_Exception if any argument is {@code NULL}.
     * @ if any array has zero length.
     * @ if the size difference between
     * {@code a} and {@code c} is not equal to 1.
     */
    static double evaluate(const std::vector<double>& a, const std::vector<double>& c, const double& z)
    {
        verify_input_array(a, c);

        const auto n = c.size();
        auto value = a[n];
        for (int i{ n - 1 }; i >= 0; i--) 
        {
            value = a[i] + (z - c[i]) * value;
        }

        return value;
    }
};