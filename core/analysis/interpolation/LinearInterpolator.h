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

//import java.lang.reflect.Array;
#include<vector>

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.analysis.polynomials.Field_Polynomial_Function;
//import org.hipparchus.analysis.polynomials.Field_Polynomial_Spline_Function;
//import org.hipparchus.analysis.polynomials.Polynomial_Function;
//import org.hipparchus.analysis.polynomials.Polynomial_Spline_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/**
 * Implements a linear function for interpolation of real univariate functions.
 *
 */
class Linear_Interpolator : public Univariate_Interpolator, public Field_Univariate_Interpolator 
{

    /**
     * Computes a linear interpolating function for the data set.
     *
     * @param x the arguments for the interpolation points
     * @param y the values for the interpolation points
     * @return a function which interpolates the data set
     * @ if {@code x} and {@code y}
     * have different sizes.
     * @ if {@code x} is not sorted in
     * strict increasing order.
     * @ if the size of {@code x} is smaller
     * than 2.
     */
    //override
    public Polynomial_Spline_Function interpolate(const std::vector<double>& x, const std::vector<double>& y)
    {
        //Math_Utils::check_not_null(x);
        //Math_Utils::check_not_null(y);
        Math_Arrays::check_equal_length(x, y);

        if (x.size() < 2) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, x.size(), 2, true);
        }

        // Number of intervals.  The number of data points is n + 1.
        int n = x.size() - 1;

        Math_Arrays::check_order(x);

        // Slope of the lines between the datapoints.
        const double m[] = std::vector<double>(n];
        for (int i{}; i < n; i++) 
        {
            m[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        }

        const Polynomial_Function polynomials[] = Polynomial_Function[n];
        const double coefficients[] = std::vector<double>(2);
        for (int i{}; i < n; i++) 
        {
            coefficients[0] = y[i];
            coefficients[1] = m[i];
            polynomials[i] = Polynomial_Function(coefficients);
        }

        return Polynomial_Spline_Function(x, polynomials);
    }

    /**
     * Computes a linear interpolating function for the data set.
     *
     * @param x the arguments for the interpolation points
     * @param y the values for the interpolation points
     * @param <T> the type of the field elements
     * @return a function which interpolates the data set
     * @ if {@code x} and {@code y}
     * have different sizes.
     * @ if {@code x} is not sorted in
     * strict increasing order.
     * @ if the size of {@code x} is smaller
     * than 2.
     * @since 1.5
     */
    //override
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public  Field_Polynomial_Spline_Function<T> interpolate(const T x[], const T y[])
    {
        //Math_Utils::check_not_null(x);
        //Math_Utils::check_not_null(y);
        Math_Arrays::check_equal_length(x, y);

        if (x.size() < 2) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, x.size(), 2, true);
        }

        // Number of intervals.  The number of data points is n + 1.
        int n = x.size() - 1;

        Math_Arrays::check_order(x);

        // Slope of the lines between the datapoints.
        const T m[] = Math_Arrays::build_array(x[0].get_field(), n);
        for (int i{}; i < n; i++) 
        {
            m[i] = y[i + 1].subtract(y[i]).divide(x[i + 1].subtract(x[i]));
        }

        //@Suppress_Warnings("unchecked")
        const Field_Polynomial_Function<T> polynomials[] =
                        (Field_Polynomial_Function<T>[]) Array.new_instance(Field_Polynomial_Function.class, n);
        const T coefficients[] = Math_Arrays::build_array(x[0].get_field(), 2);
        for (int i{}; i < n; i++) 
        {
            coefficients[0] = y[i];
            coefficients[1] = m[i];
            polynomials[i] = Field_Polynomial_Function<>(coefficients);
        }

        return Field_Polynomial_Spline_Function<>(x, polynomials);
    }

}


