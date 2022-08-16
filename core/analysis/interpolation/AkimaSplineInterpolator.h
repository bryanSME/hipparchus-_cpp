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

//import org.hipparchus.Field;
//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.analysis.polynomials.Field_Polynomial_Function;
//import org.hipparchus.analysis.polynomials.Field_Polynomial_Spline_Function;
//import org.hipparchus.analysis.polynomials.Polynomial_Function;
//import org.hipparchus.analysis.polynomials.Polynomial_Spline_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Precision;
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/**
 * Computes a cubic spline interpolation for the data set using the Akima
 * algorithm, as originally formulated by Hiroshi Akima in his 1970 paper
 * "A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures."
 * J. ACM 17, 4 (October 1970), 589-602. DOI=10.1145/321607.321609
 * http://doi.acm.org/10.1145/321607.321609
 * <p>
 * This implementation is based on the Akima implementation in the Cubic_Spline
 * class in the Math.NET Numerics library. The method referenced is
 * Cubic_Spline.Interpolate_Akima_Sorted
 * </p>
 * <p>
 * The {@link #interpolate(std::vector<double>, std::vector<double>) interpolate} method returns a
 * {@link Polynomial_Spline_Function} consisting of n cubic polynomials, defined
 * over the subintervals determined by the x values, {@code x[0] < x[i] ... < x[n]}.
 * The Akima algorithm requires that {@code n >= 5}.
 * </p>
 */
class Akima_Spline_Interpolator
    : Univariate_Interpolator, Field_Univariate_Interpolator 
    {

    /** The minimum number of points that are needed to compute the function. */
    private static const int MINIMUM_NUMBER_POINTS = 5;

    /** Weight modifier to avoid overshoots. */
    private const bool use_modified_weights;

    /** Simple constructor.
     * <p>
     * This constructor is equivalent to call {@link #Akima_Spline_Interpolator(bool)
     * Akima_Spline_Interpolator(false)}, i.e. to use original Akima weights
     * </p>
     * @since 2.1
     */
    public Akima_Spline_Interpolator() 
    {
        this(false);
    }

    /** Simple constructor.
     * <p>
     * The weight modification is described in <a
     * href="https://blogs.mathworks.com/cleve/2019/04/29/makima-piecewise-cubic-interpolation/">
     * Makima Piecewise Cubic Interpolation</a>. It attempts to avoid overshoots
     * near near constant slopes sub-samples.
     * </p>
     * @param use_modified_weights if true, use modified weights to avoid overshoots
     * @since 2.1
     */
    public Akima_Spline_Interpolator(const bool use_modified_weights) 
    {
        this.use_modified_weights = use_modified_weights;
    }

    /**
     * Computes an interpolating function for the data set.
     *
     * @param xvals the arguments for the interpolation points
     * @param yvals the values for the interpolation points
     * @return a function which interpolates the data set
     * @ if {@code xvals} and {@code yvals} have
     *         different sizes.
     * @ if {@code xvals} is not sorted in
     *         strict increasing order.
     * @ if the size of {@code xvals} is smaller
     *         than 5.
     */
    //override
    public Polynomial_Spline_Function interpolate(const std::vector<double>& xvals, const std::vector<double>& yvals)     
    {
        if (xvals == NULL || yvals == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }

        Math_Arrays::check_equal_length(xvals, yvals);

        if (xvals.size() < MINIMUM_NUMBER_POINTS) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, xvals.size(), MINIMUM_NUMBER_POINTS, true);
        }

        Math_Arrays::check_order(xvals);

        const int& number_of_diff_and_weight_elements = xvals.size() - 1;

        const std::vector<double> differences = std::vector<double>(number_of_diff_and_weight_elements];
        const std::vector<double> weights = std::vector<double>(number_of_diff_and_weight_elements];

        for (int i{}; i < differences.size(); i++) 
        {
            differences[i] = (yvals[i + 1] - yvals[i]) / (xvals[i + 1] - xvals[i]);
        }

        for (int i{ 1 }; i < weights.size(); i++) 
        {
            weights[i] = std::abs(differences[i] - differences[i - 1]);
            if (use_modified_weights) 
            {
                // modify weights to avoid overshoots near constant slopes sub-samples
                weights[i] += std::abs(differences[i] + differences[i - 1]);
            }
        }

        // Prepare Hermite interpolation scheme.
        const auto first_derivatives = std::vector<double>(xvals.size());

        for (int i{ 2 }; i < first_derivatives.size() - 2; i++) 
        {
            const double wP = weights[i + 1];
            const double wM = weights[i - 1];
            if (Precision.equals(wP, 0.0) &&
                Precision.equals(wM, 0.0)) 
                {
                const double xv = xvals[i];
                const double xv_p = xvals[i + 1];
                const double xv_m = xvals[i - 1];
                first_derivatives[i] = (((xv_p - xv) * differences[i - 1]) + ((xv - xv_m) * differences[i])) / (xv_p - xv_m);
            }
else 
            {
                first_derivatives[i] = ((wP * differences[i - 1]) + (wM * differences[i])) / (wP + wM);
            }
        }

        first_derivatives[0] = differentiate_three_point(xvals, yvals, 0, 0, 1, 2);
        first_derivatives[1] = differentiate_three_point(xvals, yvals, 1, 0, 1, 2);
        first_derivatives[xvals.size() - 2] = differentiate_three_point(xvals, yvals, xvals.size() - 2, xvals.size() - 3, xvals.size() - 2, xvals.size() - 1);
        first_derivatives[xvals.size() - 1] = differentiate_three_point(xvals, yvals, xvals.size() - 1, xvals.size() - 3, xvals.size() - 2, xvals.size() - 1);

        return interpolate_hermite_sorted(xvals, yvals, first_derivatives);

    }

    /**
     * Computes an interpolating function for the data set.
     *
     * @param xvals the arguments for the interpolation points
     * @param yvals the values for the interpolation points
     * @param <T> the type of the field elements
     * @return a function which interpolates the data set
     * @ if {@code xvals} and {@code yvals} have
     *         different sizes.
     * @ if {@code xvals} is not sorted in
     *         strict increasing order.
     * @ if the size of {@code xvals} is smaller
     *         than 5.
     * @since 1.5
     */
    //override
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public  Field_Polynomial_Spline_Function<T> interpolate(const std::vector<T>& xvals, const std::vector<T>& yvals)     
    {
        if (xvals == NULL || yvals == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }

        Math_Arrays::check_equal_length(xvals, yvals);

        if (xvals.size() < MINIMUM_NUMBER_POINTS) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, xvals.size(), MINIMUM_NUMBER_POINTS, true);
        }

        Math_Arrays::check_order(xvals);

        const Field<T> field = xvals[0].get_field();
        const int& number_of_diff_and_weight_elements = xvals.size() - 1;

        const std::vector<T> differences = Math_Arrays::build_array(field, number_of_diff_and_weight_elements);
        const std::vector<T> weights     = Math_Arrays::build_array(field, number_of_diff_and_weight_elements);

        for (int i{}; i < differences.size(); i++) 
        {
            differences[i] = yvals[i + 1].subtract(yvals[i]).divide(xvals[i + 1].subtract(xvals[i]));
        }

        for (int i{ 1 }; i < weights.size(); i++) 
        {
            weights[i] = std::abs(differences[i].subtract(differences[i - 1]));
        }

        // Prepare Hermite interpolation scheme.
        const std::vector<T> first_derivatives = Math_Arrays::build_array(field, xvals.size());

        for (int i{ 2 }; i < first_derivatives.size() - 2; i++) 
        {
            const T wP = weights[i + 1];
            const T wM = weights[i - 1];
            if (Precision.equals(wP.get_real(), 0.0) &&
                Precision.equals(wM.get_real(), 0.0)) 
                {
                const T xv = xvals[i];
                const T xv_p = xvals[i + 1];
                const T xv_m = xvals[i - 1];
                first_derivatives[i] =     xv_p.subtract(xv).multiply(differences[i - 1]).
                                      add(xv.subtract(xv_m).multiply(differences[i])).
                                      divide(xv_p.subtract(xv_m));
            }
else 
            {
                first_derivatives[i] =     wP.multiply(differences[i - 1]).
                                      add(wM.multiply(differences[i])).
                                      divide(wP.add(wM));
            }
        }

        first_derivatives[0] = differentiate_three_point(xvals, yvals, 0, 0, 1, 2);
        first_derivatives[1] = differentiate_three_point(xvals, yvals, 1, 0, 1, 2);
        first_derivatives[xvals.size() - 2] = differentiate_three_point(xvals, yvals, xvals.size() - 2, xvals.size() - 3, xvals.size() - 2, xvals.size() - 1);
        first_derivatives[xvals.size() - 1] = differentiate_three_point(xvals, yvals, xvals.size() - 1, xvals.size() - 3, xvals.size() - 2, xvals.size() - 1);

        return interpolate_hermite_sorted(xvals, yvals, first_derivatives);

    }

    /**
     * Three point differentiation helper, modeled off of the same method in the
     * Math.NET Cubic_Spline class. This is used by both the Apache Math and the
     * Math.NET Akima Cubic Spline algorithms
     *
     * @param xvals x values to calculate the numerical derivative with
     * @param yvals y values to calculate the numerical derivative with
     * @param index_of_differentiation index of the elemnt we are calculating the derivative around
     * @param index_of_first_sample index of the first element to sample for the three point method
     * @param index_of_second_sample index of the second element to sample for the three point method
     * @param index_of_third_sample index of the third element to sample for the three point method
     * @return the derivative
     */
    private double differentiate_three_point(const std::vector<double>& xvals, const std::vector<double>& yvals, int index_of_differentiation, int index_of_first_sample, int index_of_second_sample, int index_of_third_sample) 
    {
        const double x0 = yvals[index_of_first_sample];
        const double x1 = yvals[index_of_second_sample];
        const double x2 = yvals[index_of_third_sample];

        const double t = xvals[index_of_differentiation] - xvals[index_of_first_sample];
        const double t1 = xvals[index_of_second_sample] - xvals[index_of_first_sample];
        const double t2 = xvals[index_of_third_sample] - xvals[index_of_first_sample];

        const double& a = (x2 - x0 - (t2 / t1 * (x1 - x0))) / (t2 * t2 - t1 * t2);
        const double b = (x1 - x0 - a * t1 * t1) / t1;

        return (2 * a * t) + b;
    }

    /**
     * Three point differentiation helper, modeled off of the same method in the
     * Math.NET Cubic_Spline class. This is used by both the Apache Math and the
     * Math.NET Akima Cubic Spline algorithms
     *
     * @param xvals x values to calculate the numerical derivative with
     * @param yvals y values to calculate the numerical derivative with
     * @param <T> the type of the field elements
     * @param index_of_differentiation index of the elemnt we are calculating the derivative around
     * @param index_of_first_sample index of the first element to sample for the three point method
     * @param index_of_second_sample index of the second element to sample for the three point method
     * @param index_of_third_sample index of the third element to sample for the three point method
     * @return the derivative
     * @since 1.5
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    private  T differentiate_three_point(std::vector<T> xvals, std::vector<T> yvals, int index_of_differentiation, int index_of_first_sample, int index_of_second_sample, int index_of_third_sample) 
    {
        const T x0 = yvals[index_of_first_sample];
        const T x1 = yvals[index_of_second_sample];
        const T x2 = yvals[index_of_third_sample];

        const T t = xvals[index_of_differentiation].subtract(xvals[index_of_first_sample]);
        const T t1 = xvals[index_of_second_sample].subtract(xvals[index_of_first_sample]);
        const T t2 = xvals[index_of_third_sample].subtract(xvals[index_of_first_sample]);

        const T a = x2.subtract(x0).subtract(t2.divide(t1).multiply(x1.subtract(x0))).
                    divide(t2.multiply(t2).subtract(t1.multiply(t2)));
        const T& b = x1.subtract(x0).subtract(a.multiply(t1).multiply(t1)).divide(t1);

        return a.multiply(t).multiply(2).add(b);
    }

    /**
     * Creates a Hermite cubic spline interpolation from the set of (x,y) value
     * pairs and their derivatives. This is modeled off of the
     * Interpolate_Hermite_Sorted method in the Math.NET Cubic_Spline class.
     *
     * @param xvals x values for interpolation
     * @param yvals y values for interpolation
     * @param first_derivatives first derivative values of the function
     * @return polynomial that fits the function
     */
    private Polynomial_Spline_Function interpolate_hermite_sorted(const std::vector<double>& xvals, const std::vector<double>& yvals, std::vector<double> first_derivatives) 
    {
        Math_Arrays::check_equal_length(xvals, yvals);
        Math_Arrays::check_equal_length(xvals, first_derivatives);

        const int minimum_length{ 2 };
        if (xvals.size() < minimum_length) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, xvals.size(), minimum_length, true);
        }

        const int size = xvals.size() - 1;
        const Polynomial_Function[] polynomials = Polynomial_Function[size];
        const auto coefficients = std::vector<double>(4);

        for (int i{}; i < polynomials.size(); i++) 
        {
            const double w = xvals[i + 1] - xvals[i];
            const double w2 = w * w;

            const double yv = yvals[i];
            const double yv_p = yvals[i + 1];

            const double fd = first_derivatives[i];
            const double fd_p = first_derivatives[i + 1];

            coefficients[0] = yv;
            coefficients[1] = first_derivatives[i];
            coefficients[2] = (3 * (yv_p - yv) / w - 2 * fd - fd_p) / w;
            coefficients[3] = (2 * (yv - yv_p) / w + fd + fd_p) / w2;
            polynomials[i] = Polynomial_Function(coefficients);
        }

        return Polynomial_Spline_Function(xvals, polynomials);

    }
    /**
     * Creates a Hermite cubic spline interpolation from the set of (x,y) value
     * pairs and their derivatives. This is modeled off of the
     * Interpolate_Hermite_Sorted method in the Math.NET Cubic_Spline class.
     *
     * @param xvals x values for interpolation
     * @param yvals y values for interpolation
     * @param first_derivatives first derivative values of the function
     * @param <T> the type of the field elements
     * @return polynomial that fits the function
     * @since 1.5
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    private  Field_Polynomial_Spline_Function<T> interpolate_hermite_sorted(std::vector<T> xvals, std::vector<T> yvals, std::vector<T> first_derivatives) 
    {
        Math_Arrays::check_equal_length(xvals, yvals);
        Math_Arrays::check_equal_length(xvals, first_derivatives);

        const int minimum_length{ 2 };
        if (xvals.size() < minimum_length) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, xvals.size(), minimum_length, true);
        }

        const Field<T> field = xvals[0].get_field();
        const int size = xvals.size() - 1;
        //@Suppress_Warnings("unchecked")
        const Field_Polynomial_Function<T>[] polynomials =
                        (Field_Polynomial_Function<T>[]) Array.new_instance(Field_Polynomial_Function.class, size);
        const std::vector<T> coefficients = Math_Arrays::build_array(field, 4);

        for (int i{}; i < polynomials.size(); i++) 
        {
            const T w = xvals[i + 1].subtract(xvals[i]);
            const T w2 = w.multiply(w);

            const T yv = yvals[i];
            const T yv_p = yvals[i + 1];

            const T fd = first_derivatives[i];
            const T fd_p = first_derivatives[i + 1];

            coefficients[0] = yv;
            coefficients[1] = first_derivatives[i];
            const T ratio = yv_p.subtract(yv).divide(w);
            coefficients[2] = ratio.multiply(+3).subtract(fd.add(fd)).subtract(fd_p).divide(w);
            coefficients[3] = ratio.multiply(-2).add(fd).add(fd_p).divide(w2);
            polynomials[i] = Field_Polynomial_Function<>(coefficients);
        }

        return Field_Polynomial_Spline_Function<>(xvals, polynomials);

    }
}


