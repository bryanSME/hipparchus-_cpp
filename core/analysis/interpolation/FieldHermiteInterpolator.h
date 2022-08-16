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

#include <vector>
//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;

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
 * @param <T> Type of the field elements.
 *
 */
class Field_Hermite_Interpolator<T extends Field_Element<T>> 
{
private:
    /** Sample abscissae. */
    const std::vector<T> my_abscissae{};

    /** Top diagonal of the divided differences array. */
    const std::vector<std::vector<T>> my_top_diagonal{};

    /** Bottom diagonal of the divided differences array. */
    const std::vector<std::vector<T>> my_bottom_diagonal{};

public:
    /** Create an empty interpolator.
     */
    Field_Hermite_Interpolator() = default;

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
     * @ if derivative structures are inconsistent
     * @ if x is NULL
     */
    //@Safe_Varargs
    const void add_sample_point(const T& x, const std::vector<T> ... value)
    {

        //Math_Utils::check_not_null(x);
        T factorial = x.get_field().get_one();
        for (int i{}; i < value.size(); ++i) 
        {

            const std::vector<T> y = value[i].clone();
            if (i > 1) 
            {
                factorial = factorial.multiply(i);
                const T inv = factorial.reciprocal();
                for (int j{}; j < y.size(); ++j) 
                {
                    y[j] = y[j].multiply(inv);
                }
            }

            // update the bottom diagonal of the divided differences array
            const int n = abscissae.size();
            bottom_diagonal.add(n - i, y);
            std::vector<T> bottom0 = y;
            for (int j{ i }; j < n; ++j)
            {
                const std::vector<T> bottom1 = bottom_diagonal.get(n - (j + 1));
                if (x.equals(abscissae.get(n - (j + 1)))) 
                {
                    throw std::exception("not implmented");
                    //throw (hipparchus::exception::Localized_Core_Formats_Type::DUPLICATED_ABSCISSA_DIVISION_BY_ZERO, x);
                }
                const T inv = x.subtract(abscissae.get(n - (j + 1))).reciprocal();
                for (int k{}; k < y.size(); ++k) 
                {
                    bottom1[k] = inv.multiply(bottom0[k].subtract(bottom1[k]));
                }
                bottom0 = bottom1;
            }

            // update the top diagonal of the divided differences array
            my_top_diagonal.add(bottom0.clone());

            // update the abscissae array
            abscissae.add(x);

        }

    }

    /** Interpolate value at a specified abscissa.
     * @param x interpolation abscissa
     * @return interpolated value
     * @exception  if sample is empty
     * @ if x is NULL
     */
    std::vector<T> value(T x)  
    {
        // safety check
        //Math_Utils::check_not_null(x);
        if (abscissae.is_empty()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_INTERPOLATION_SAMPLE);
        }

        const std::vector<T> value = Math_Arrays::build_array(x.get_field(), top_diagonal.get(0).size());
        T value_coeff = x.get_field().get_one();
        for (int i{}; i < top_diagonal.size(); ++i) 
        {
            std::vector<T> divided_difference = top_diagonal.get(i);
            for (int k{}; k < value.size(); ++k) 
            {
                value[k] = value[k].add(divided_difference[k].multiply(value_coeff));
            }
            const T delta_x = x.subtract(abscissae.get(i));
            value_coeff = value_coeff.multiply(delta_x);
        }

        return value;

    }

    /** Interpolate value and first derivatives at a specified abscissa.
     * @param x interpolation abscissa
     * @param order maximum derivation order
     * @return interpolated value and derivatives (value in row 0, * 1<sup>st</sup> derivative in row 1, ... n<sup>th</sup> derivative in row n)
     * @exception  if sample is empty
     * @ if x is NULL
     */
    std::vector<std::vector<T>> derivatives(T x, int order)  
    {
        // safety check
        //Math_Utils::check_not_null(x);
        if (abscissae.is_empty()) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_INTERPOLATION_SAMPLE);
        }

        const T zero = x.get_field().get_zero();
        const T one  = x.get_field().get_one();
        const std::vector<T> tj = Math_Arrays::build_array(x.get_field(), order + 1);
        tj[0] = zero;
        for (int i{}; i < order; ++i) 
        {
            tj[i + 1] = tj[i].add(one);
        }

        const std::vector<std::vector<T>> derivatives =
                Math_Arrays::build_array(x.get_field(), order + 1, top_diagonal.get(0).size());
        const std::vector<T> value_coeff = Math_Arrays::build_array(x.get_field(), order + 1);
        value_coeff[0] = x.get_field().get_one();
        for (int i{}; i < top_diagonal.size(); ++i) 
        {
            std::vector<T> divided_difference = top_diagonal.get(i);
            const T delta_x = x.subtract(abscissae.get(i));
            for (int j = order; j >= 0; --j) 
            {
                for (int k{}; k < derivatives[j].size(); ++k) 
                {
                    derivatives[j][k] =
                            derivatives[j][k].add(divided_difference[k].multiply(value_coeff[j]));
                }
                value_coeff[j] = value_coeff[j].multiply(delta_x);
                if (j > 0) 
                {
                    value_coeff[j] = value_coeff[j].add(tj[j].multiply(value_coeff[j - 1]));
                }
            }
        }

        return derivatives;
    }

};