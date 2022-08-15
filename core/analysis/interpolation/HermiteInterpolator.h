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
//package org.hipparchus.analysis.interpolation;

//import java.util.Array_list;
//import java.util.Arrays;
//import java.util.List;
#include <cmath>

//import org.hipparchus.analysis.differentiation.Derivative;
//import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Vector_Function;
//import org.hipparchus.analysis.polynomials.Polynomial_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Combinatorics_Utils;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include "../../analysis/differentiation/UnivariateDifferentiableVectorFunction.h"


/** Polynomial interpolator using both sample values and sample derivatives.
 * <p>
 * The interpolation polynomials match all sample points, including both values
 * and provided derivatives. There is one polynomial for each component of
 * the values vector. All polynomials have the same degree. The degree of the
 * polynomials depends on the number of points and number of derivatives at each
 * point. For example the interpolation polynomials for n sample points without
 * any derivatives all have degree n-1. The interpolation polynomials for n
 * sample points with the two extreme points having value and first derivative
 * and the remaining points having value only all have degree n+1. The
 * interpolation polynomial for n sample points with value, first and second
 * derivative for all points all have degree 3n-1.
 * </p>
 *
 */
class Hermite_Interpolator : Univariate_Differentiable_Vector_Function 
{
private:
    /** Sample abscissae. */
    const std::vector<double> my_abscissae;

    /** Top diagonal of the divided differences array. */
    const std::vector<std::vector<double>> my_top_diagonal;

    /** Bottom diagonal of the divided differences array. */
    const std::vector<std::vector<double>> my_bottom_diagonal;

    /** Check interpolation can be performed.
     * @exception  if interpolation cannot be performed
     * because sample is empty
     */
    void check_interpolation()
    {
        if (my_abscissae.empty())
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_INTERPOLATION_SAMPLE);
        }
    }

    /** Create a polynomial from its coefficients.
     * @param c polynomials coefficients
     * @return polynomial
     */
    Polynomial_Function polynomial(double ... c)
    {
        return Polynomial_Function(c);
    }

public:
    /** Create an empty interpolator.
     */
    Hermite_Interpolator() = default;

    /** Add a sample point.
     * <p>
     * This method must be called once for each sample point. It is allowed to
     * mix some calls with values only with calls with values and first
     * derivatives.
     * </p>
     * <p>
     * The point abscissae for all calls <em>must</em> be different.
     * </p>
     * @param x abscissa of the sample point
     * @param value value and derivatives of the sample point
     * (if only one row is passed, it is the value, if two rows are
     * passed the first one is the value and the second the derivative
     * and so on)
     * @exception  if the abscissa difference between added point
     * and a previous point is zero (i.e. the two points are at same abscissa)
     * @exception Math_Runtime_Exception if the number of derivatives is larger
     * than 20, which prevents computation of a factorial
     */
    void add_sample_point(const double& x, const std::vector<double>& value)
    {
        for (int i{}; i < value.size(); ++i) 
        {
            auto y = value[i];
            if (i > 1) 
            {
                double inv = 1.0 / Combinatorics_Utils.factorial(i);
                for (int j{}; j < y.size(); ++j) 
                {
                    y[j] *= inv;
                }
            }

            // update the bottom diagonal of the divided differences array
            const int n = my_abscissae.size();
            my_bottom_diagonal.add(n - i, y);
            std::vector<double> bottom0 = y;
            for (int j = i; j < n; ++j) 
            {
                auto bottom1 = my_bottom_diagonal.at(n - (j + 1));
                const double inv = 1.0 / (x - my_abscissae.at(n - (j + 1)));
                if (std::isinfinite(inv)) 
                {
                    throw (hipparchus::exception::Localized_Core_Formats_Type::DUPLICATED_ABSCISSA_DIVISION_BY_ZERO, x);
                }
                for (int k{}; k < y.size(); ++k) 
                {
                    bottom1[k] = inv * (bottom0[k] - bottom1[k]);
                }
                bottom0 = bottom1;
            }

            // update the top diagonal of the divided differences array
            my_top_diagonal.add(bottom0);

            // update the abscissae array
            my_abscissae.add(x);

        }

    }

    /** Compute the interpolation polynomials.
     * @return interpolation polynomials array
     * @exception  if sample is empty
     */
    std::vector<Polynomial_Function> get_polynomials()
    {

        // safety check
        check_interpolation();

        // iteration initialization
        const auto zero = polynomial(0);
        auto polynomials = Polynomial_Function[my_top_diagonal.at(0).size()];
        for (int i{}; i < polynomials.size(); ++i) 
        {
            polynomials[i] = zero;
        }
        auto coeff = polynomial(1);

        // build the polynomials by iterating on the top diagonal of the divided differences array
        for (int i{}; i < my_top_diagonal.size(); ++i) 
        {
            auto tdi = my_top_diagonal.at(i);
            for (int k{}; k < polynomials.size(); ++k) 
            {
                polynomials[k] = polynomials[k].add(coeff.multiply(polynomial(tdi[k])));
            }
            coeff = coeff.multiply(polynomial(-my_abscissae.at(i), 1.0));
        }

        return polynomials;

    }

    /** Interpolate value at a specified abscissa.
     * <p>
     * Calling this method is equivalent to call the {@link Polynomial_Function#valuestatic_cast<double>(
     * value} methods of all polynomials returned by {@link #get_polynomials() get_polynomials}, * except it does not build the intermediate polynomials, so this method is faster and
     * numerically more stable.
     * </p>
     * @param x interpolation abscissa
     * @return interpolated value
     * @exception  if sample is empty
     */
    //override
    std::vector<double> value(double x)  
    {

        // safety check
        check_interpolation();

        auto value = std::vector<double>(my_top_diagonal.at(0).size());
        double value_coeff = 1;
        for (int i{}; i < my_top_diagonal.size(); ++i) 
        {
            auto divided_difference = my_top_diagonal.at(i);
            for (int k{}; k < value.size(); ++k) 
            {
                value[k] += divided_difference[k] * value_coeff;
            }
            const double delta_x = x - my_abscissae.at(i);
            value_coeff *= delta_x;
        }

        return value;

    }

    /** {@inherit_doc}. */
    //override
    //<T extends Derivative<T>> std::vector<T> value(T x)
    //{

    //    // safety check
    //    check_interpolation();

    //    const std::vector<T> value = Math_Arrays::build_array(x.get_field(), top_diagonal.get(0).size());
    //    Arrays.fill(value, x.get_field().get_zero());
    //    T value_coeff = x.get_field().get_one();
    //    for (int i{}; i < top_diagonal.size(); ++i) 
    //    {
    //        std::vector<double> divided_difference = top_diagonal.get(i);
    //        for (int k{}; k < value.size(); ++k) 
    //        {
    //            value[k] = value[k].add(value_coeff.multiply(divided_difference[k]));
    //        }
    //        const T delta_x = x.subtract(abscissae.get(i));
    //        value_coeff = value_coeff.multiply(delta_x);
    //    }

    //    return value;

    //}

    /** Interpolate value and first derivatives at a specified abscissa.
     * @param x interpolation abscissa
     * @param order maximum derivation order
     * @return interpolated value and derivatives (value in row 0, * 1<sup>st</sup> derivative in row 1, ... n<sup>th</sup> derivative in row n)
     * @exception  if sample is empty
     * @Null_Argument_Exception if x is NULL
     */
    std::vector<std::vector<double>> derivatives(const double& x, const int& order)
    {
        // safety check
        //Math_Utils::check_not_null(x);
        if (my_abscissae.empty()) 
        {
            throw std::exception("HermiteInterpolator derivatives exception not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_INTERPOLATION_SAMPLE);
        }

        auto tj = std::vector<double>(order + 1);
        tj[0] = 0;
        for (int i{}; i < order; ++i) 
        {
            tj[i + 1] = tj[i] + 1;
        }

        auto derivatives = std::vector<std::vector<double>>(order + 1, std::vector<double>(my_top_diagonal.at(0).size()));
        auto value_coeff = std::vector<double>(order + 1);
        value_coeff[0] = 1.0;
        for (int i{}; i < my_top_diagonal.size(); ++i) 
        {
            auto divided_difference = my_top_diagonal.at(i);
            const double delta_x = x - my_abscissae.at(i);
            for(int j = order; j >= 0; --j)
            {
                for (int k{}; k < derivatives[j].size(); ++k) 
                {
                    derivatives[j][k] += divided_difference[k] * value_coeff[j];
                }
                value_coeff[j] *= delta_x;
                if (j > 0) 
                {
                    value_coeff[j] += tj[j] * value_coeff[j - 1];
                }
            }
        }

        return derivatives;
    }
};