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

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
#include "../UnivariateFunction.h"
#include "../../util/MathArrays.h"
#include <vector>
/**
 * Implements the representation of a real polynomial function in
 * <a href="http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html">
 * Lagrange Form</a>. For reference, see <b>Introduction to Numerical
 * Analysis</b>, ISBN 038795452X, chapter 2.
 * <p>
 * The approximated function should be smooth enough for Lagrange polynomial
 * to work well. Otherwise, consider using splines instead.</p>
 *
 */
class Polynomial_Function_Lagrange_Form : Univariate_Function 
{
private:
    /**
     * The coefficients of the polynomial, ordered by degree -- i.e.
     * coefficients[0] is the constant term and coefficients[n] is the
     * coefficient of x^n where n is the degree of the polynomial.
     */
    std::vector<double> my_coefficients;
    /**
     * Interpolating points (abscissas).
     */
    const std::vector<double> my_x;
    /**
     * Function values at interpolating points.
     */
    const std::vector<double> my_y;
    /**
     * Whether the polynomial coefficients are available.
     */
    bool my_coefficients_computed;

    /**
     * Evaluate the Lagrange polynomial using
     * <a href="http://mathworld.wolfram.com/NevillesAlgorithm.html">
     * Neville's Algorithm</a>. It takes O(n^2) time.
     *
     * @param x Interpolating points array.
     * @param y Interpolating values array.
     * @param z Point at which the function value is to be computed.
     * @return the function value.
     * @ if {@code x} and {@code y} have
     * different lengths.
     * @org.hipparchus.exception.
     * if {@code x} is not sorted in strictly increasing order.
     * @ if the size of {@code x} is less
     * than 2.
     */
    static double evaluate_internal(const std::vector<double>& x, const std::vector<double>& y, const double& z)
    {
        int nearest{};
        const int n = x.size();
        auto c = std::vector<double>(n);
        auto d = std::vector<double>(n);
        double min_dist{ INFINITY };
        for (int i{}; i < n; i++)
        {
            // initialize the difference arrays
            c[i] = y[i];
            d[i] = y[i];
            // find out the abscissa closest to z
            const double dist = std::abs(z - x[i]);
            if (dist < min_dist)
            {
                nearest = i;
                min_dist = dist;
            }
        }

        // initial approximation to the function value at z
        double value = y[nearest];

        for (int i{ 1 }; i < n; i++)
        {
            for (int j{}; j < n - i; j++)
            {
                const double tc = x[j] - z;
                const double td = x[i + j] - z;
                const double divider = x[j] - x[i + j];
                // update the difference arrays
                const double w = (c[j + 1] - d[j]) / divider;
                c[j] = tc * w;
                d[j] = td * w;
            }
            // sum up the difference terms to get the const value
            if (nearest < 0.5 * (n - i + 1))
            {
                value += c[nearest];    // fork down
            }
            else
            {
                nearest--;
                value += d[nearest];    // fork up
            }
        }

        return value;
    }

public:
    /**
     * Construct a Lagrange polynomial with the given abscissas and function
     * values. The order of interpolating points are not important.
     * <p>
     * The constructor makes copy of the input arrays and assigns them.</p>
     *
     * @param x interpolating points
     * @param y function values at interpolating points
     * @ if the array lengths are different.
     * @ if the number of points is less than 2.
     * @
     * if two abscissae have the same value.
     */
    Polynomial_Function_Lagrange_Form(const std::vector<double>& x, const std::vector<double>& y)
        : my_x{ std::vector<double>(x.size()) }, my_y{std::vector<double>(y.size())}
    {
        System.arraycopy(x, 0, my_x, 0, x.size());
        System.arraycopy(y, 0, my_y, 0, y.size());
        coefficients_computed = false;

        if (!verify_interpolation_array(x, y, false)) 
        {
            Math_Arrays::sort_in_place(my_x, my_y);
            // Second check in case some abscissa is duplicated.
            verify_interpolation_array(my_x, my_y, true);
        }
    }

    /**
     * Calculate the function value at the given point.
     *
     * @param z Point at which the function value is to be computed.
     * @return the function value.
     * @ if {@code x} and {@code y} have
     * different lengths.
     * @org.hipparchus.exception.
     * if {@code x} is not sorted in strictly increasing order.
     * @ if the size of {@code x} is less
     * than 2.
     */
    //override
    double value(const double& z) 
    {
        return evaluate_internal(my_x, my_y, z);
    }

    /**
     * Returns the degree of the polynomial.
     *
     * @return the degree of the polynomial
     */
    int degree() const 
    {
        return my_x.size() - 1;
    }

    /**
     * Returns a copy of the interpolating points array.
     * <p>
     * Changes made to the returned copy will not affect the polynomial.</p>
     *
     * @return a fresh copy of the interpolating points array
     */
    std::vector<double> get_interpolating_points() const
    {
        auto out = std::vector<double>(my_x.size());
        System.arraycopy(my_x, 0, out, 0, my_x.size());
        return out;
    }

    /**
     * Returns a copy of the interpolating values array.
     * <p>
     * Changes made to the returned copy will not affect the polynomial.</p>
     *
     * @return a fresh copy of the interpolating values array
     */
    std::vector<double> get_interpolating_values() 
    {
        auto out = std::vector<double>(my_y.size());
        System.arraycopy(my_y, 0, out, 0, y.size());
        return out;
    }

    /**
     * Returns a copy of the coefficients array.
     * <p>
     * Changes made to the returned copy will not affect the polynomial.</p>
     * <p>
     * Note that coefficients computation can be ill-conditioned. Use with caution
     * and only when it is necessary.</p>
     *
     * @return a fresh copy of the coefficients array
     */
    std::vector<double> get_coefficients() 
    {
        if (!my_coefficients_computed) 
        {
            compute_coefficients();
        }
        auto out = std::vector<double>(my_coefficients.size());
        System.arraycopy(my_coefficients, 0, out, 0, my_coefficients.size());
        return out;
    }

    /**
     * Evaluate the Lagrange polynomial using
     * <a href="http://mathworld.wolfram.com/NevillesAlgorithm.html">
     * Neville's Algorithm</a>. It takes O(n^2) time.
     *
     * @param x Interpolating points array.
     * @param y Interpolating values array.
     * @param z Point at which the function value is to be computed.
     * @return the function value.
     * @ if {@code x} and {@code y} have
     * different lengths.
     * @
     * if {@code x} is not sorted in strictly increasing order.
     * @ if the size of {@code x} is less
     * than 2.
     */
    static double evaluate(const std::vector<double>& x, const std::vector<double>& y, const double& z)
    {
        if (verify_interpolation_array(x, y, false)) 
        {
            return evaluate_internal(x, y, z);
        }

        // Array is not sorted.
        auto x_new = std::vector<double>(x.size());
        auto y_new = std::vector<double>(y.size());
        System.arraycopy(x, 0, x_new, 0, x.size());
        System.arraycopy(y, 0, y_new, 0, y.size());

        Math_Arrays::sort_in_place(x_new, y_new);
        // Second check in case some abscissa is duplicated.
        verify_interpolation_array(x_new, y_new, true);
        return evaluate_internal(x_new, y_new, z);
    }


    /**
     * Check that the interpolation arrays are valid.
     * The arrays features checked by this method are that both arrays have the
     * same length and this length is at least 2.
     *
     * @param x Interpolating points array.
     * @param y Interpolating values array.
     * @param abort Whether to throw an exception if {@code x} is not sorted.
     * @ if the array lengths are different.
     * @ if the number of points is less than 2.
     * @org.hipparchus.exception.
     * if {@code x} is not sorted in strictly increasing order and {@code abort}
     * is {@code true}.
     * @return {@code false} if the {@code x} is not sorted in increasing order, * {@code true} otherwise.
     * @see #evaluate(std::vector<double>, std::vector<double>, double)
     * @see #compute_coefficients()
     */
    static bool verify_interpolation_array(const std::vector<double>& x, const std::vector<double>& y, bool abort)
    {
        Math_Arrays::check_equal_length(x, y);
        if (x.size() < 2) 
        {
            throw std::exception("hipparchus::exception::Localized_Core_Formats_Type::WRONG_NUMBER_OF_POINTS, 2, x.size(), true);");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::WRONG_NUMBER_OF_POINTS, 2, x.size(), true);
        }

        return Math_Arrays::check_order(x, Math_Arrays::Order_Direction::INCREASING, true, abort);
    }

protected:
        /**
         * Calculate the coefficients of Lagrange polynomial from the
         * interpolation data. It takes O(n^2) time.
         * Note that this computation can be ill-conditioned: Use with caution
         * and only when it is necessary.
         */
        void compute_coefficients()
        {
            const int n = degree() + 1;
            my_coefficients = std::vector<double>(n);
            for (int i{}; i < n; i++)
            {
                my_coefficients[i] = 0;
            }

            // c[] are the coefficients of P(x) = (x-x[0])(x-x[1])...(x-x[n-1])
            auto c = std::vector<double>(n + 1);
            c[0] = 1.0;
            for (int i{}; i < n; i++)
            {
                for (int j = i; j > 0; j--)
                {
                    c[j] = c[j - 1] - c[j] * my_x[i];
                }
                c[0] *= -my_x[i];
                c[i + 1] = 1;
            }

            auto tc = std::vector<double>(n);
            for (int i{}; i < n; i++)
            {
                // d = (x[i]-x[0])...(x[i]-x[i-1])(x[i]-x[i+1])...(x[i]-x[n-1])
                double d{ 1 };
                for (int j{}; j < n; j++)
                {
                    if (i != j)
                    {
                        d *= my_x[i] - my_x[j];
                    }
                }
                const double t = my_y[i] / d;
                // Lagrange polynomial is the sum of n terms, each of which is a
                // polynomial of degree n-1. tc[] are the coefficients of the i-th
                // numerator Pi(x) = (x-x[0])...(x-x[i-1])(x-x[i+1])...(x-x[n-1]).
                tc[n - 1] = c[n];     // actually c[n] = 1
                my_coefficients[n - 1] += t * tc[n - 1];
                for (int j{ n - 2 }; j >= 0; j--)
                {
                    tc[j] = c[j + 1] + tc[j + 1] * my_x[i];
                    my_coefficients[j] += t * tc[j];
                }
            }

            my_coefficients_computed = true;
        }
};