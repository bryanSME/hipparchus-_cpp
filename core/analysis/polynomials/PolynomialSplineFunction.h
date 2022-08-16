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

//import java.util.Arrays;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.analysis.Field_Univariate_Function;
//import org.hipparchus.analysis.differentiation.Derivative;
//import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <vector>
#include "PolynomialFunction.h"
#include "../../util/MathUtils.h"
#include "../../util/MathArrays.h"
#include "../differentiation/UnivariateDifferentiableFunction.h"
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/**
 * Represents a polynomial spline function.
 * <p>
 * A <strong>polynomial spline function</strong> consists of a set of
 * <i>interpolating polynomials</i> and an ascending array of domain
 * <i>knot points</i>, determining the intervals over which the spline function
 * is defined by the constituent polynomials.  The polynomials are assumed to
 * have been computed to match the values of another function at the knot
 * points.  The value consistency constraints are not currently enforced by
 * <code>Polynomial_Spline_Function</code> itself, but are assumed to hold among
 * the polynomials and knot points passed to the constructor.</p>
 * <p>
 * N.B.:  The polynomials in the <code>polynomials</code> property must be
 * centered on the knot points to compute the spline function values.
 * See below.</p>
 * <p>
 * The domain of the polynomial spline function is
 * <code>[smallest knot, largest knot]</code>.  Attempts to evaluate the
 * function at values outside of this range generate Illegal_Argument_Exceptions.
 * </p>
 * <p>
 * The value of the polynomial spline function for an argument <code>x</code>
 * is computed as follows:
 * <ol>
 * <li>The knot array is searched to find the segment to which <code>x</code>
 * belongs.  If <code>x</code> is less than the smallest knot point or greater
 * than the largest one, an <code>Illegal_Argument_Exception</code>
 * is thrown.</li>
 * <li> Let <code>j</code> be the index of the largest knot point that is less
 * than or equal to <code>x</code>.  The value returned is
 * {@code polynomials[j](x - knot[j])}</li></ol>
 *
 */
class Polynomial_Spline_Function : Univariate_Differentiable_Function, Field_Univariate_Function 
{
private:
    /**
     * Spline segment interval delimiters (knots).
     * Size is n + 1 for n segments.
     */
    std::vector<double> my_knots;
    /**
     * The polynomial functions that make up the spline.  The first element
     * determines the value of the spline over the first subinterval, the
     * second over the second, etc.   Spline function values are determined by
     * evaluating these functions at {@code (x - knot[i])} where i is the
     * knot segment to which x belongs.
     */
    std::vector<Polynomial_Function> my_polynomials;
    /**
     * Number of spline segments. It is equal to the number of polynomials and
     * to the number of partition points - 1.
     */
    int my_n;


    /**
     * Construct a polynomial spline function with the given segment delimiters
     * and interpolating polynomials.
     * The constructor copies both arrays and assigns the copies to the knots
     * and polynomials properties, respectively.
     *
     * @param knots Spline segment interval delimiters.
     * @param polynomials Polynomial functions that make up the spline.
     * @Null_Argument_Exception if either of the input arrays is {@code NULL}.
     * @ if knots has length less than 2.
     * @ if {@code polynomials.size() != knots.size() - 1}.
     * @ if the {@code knots} array is not strictly increasing.
     *
     */
    Polynomial_Spline_Function(const std::vector<double>& knots, const std::vector<Polynomial_Function>& polynomials)
    {
        if (knots.size() < 2) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_ENOUGH_POINTS_IN_SPLINE_PARTITION, 2, knots.size(), false);
        }
        Math_Utils::check_dimension(polynomials.size(), knots.size() - 1);
        Math_Arrays::check_order(knots);

        my_n = knots.size() -1;
        my_knots = std::vector<double>(n + 1);
        System.arraycopy(knots, 0, my_knots, 0, n + 1);
        my_polynomials = std::vector<Polynomial_Function>(n);
        System.arraycopy(polynomials, 0, my_polynomials, 0, n);
    }

    /**
     * Compute the value for the function.
     * See {@link Polynomial_Spline_Function} for details on the algorithm for
     * computing the value of the function.
     *
     * @param v Point for which the function value should be computed.
     * @return the value.
     * @ if {@code v} is outside of the domain of the
     * spline function (smaller than the smallest knot point or larger than the
     * largest knot point).
     */
    //override
    double value(double v) 
    {
        Math_Utils::check_range_inclusive(v, my_knots[0], my_knots[n]);
        int i = Arrays.binary_search(my_knots, v);
        if (i < 0) 
        {
            i = -i - 2;
        }
        // This will handle the case where v is the last knot value
        // There are only n-1 polynomials, so if v is the last knot
        // then we will use the last polynomial to calculate the value.
        if ( i >= my_polynomials.size() ) 
        {
            i--;
        }
        return my_polynomials[i].value(v - my_knots[i]);
    }

    /**
     * Get the derivative of the polynomial spline function.
     *
     * @return the derivative function.
     */
    Polynomial_Spline_Function polynomial_spline_derivative() 
    {
        auto derivative_polynomials = std::vector<Polynomial_Function>(my_n);
        for (int i{}; i < my_n; i++) 
        {
            derivative_polynomials[i] = my_polynomials[i].polynomial_derivative();
        }
        return Polynomial_Spline_Function(my_knots, derivative_polynomials);
    }


    /** {@inherit_doc}
     */
    //override
    <T extends Derivative<T>> T value(const T t) 
    {
        const double t0 = t.get_real();
        Math_Utils::check_range_inclusive(t0, knots[0], knots[n]);
        int i = Arrays.binary_search(knots, t0);
        if (i < 0) 
        {
            i = -i - 2;
        }
        // This will handle the case where t is the last knot value
        // There are only n-1 polynomials, so if t is the last knot
        // then we will use the last polynomial to calculate the value.
        if ( i >= polynomials.size() ) 
        {
            i--;
        }
        return polynomials[i].value(t.subtract(knots[i]));
    }

    /**
     * {@inherit_doc}
     */
    //override
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    T value(const T t) 
    {
        const double t0 = t.get_real();
        Math_Utils::check_range_inclusive(t0, knots[0], knots[n]);
        int i = Arrays.binary_search(knots, t0);
        if (i < 0) 
        {
            i = -i - 2;
        }
        // This will handle the case where t is the last knot value
        // There are only n-1 polynomials, so if t is the last knot
        // then we will use the last polynomial to calculate the value.
        if ( i >= polynomials.size() ) 
        {
            i--;
        }
        return polynomials[i].value(t.subtract(knots[i]));
    }

    /**
     * Get the number of spline segments.
     * It is also the number of polynomials and the number of knot points - 1.
     *
     * @return the number of spline segments.
     */
    int get_n() const
    {
        return my_n;
    }

    /**
     * Get a copy of the interpolating polynomials array.
     * It returns a fresh copy of the array. Changes made to the copy will
     * not affect the polynomials property.
     *
     * @return the interpolating polynomials.
     */
    std::vector<Polynomial_Function> get_polynomials() 
    {
        Polynomial_Function p[] = Polynomial_Function[n];
        System.arraycopy(polynomials, 0, p, 0, n);
        return p;
    }

    /**
     * Get an array copy of the knot points.
     * It returns a fresh copy of the array. Changes made to the copy
     * will not affect the knots property.
     *
     * @return the knot points.
     */
    std::vector<double> get_knots() 
    {
        auto out = std::vector<double>(n + 1);
        System.arraycopy(knots, 0, out, 0, n + 1);
        return out;
    }

    /**
     * Indicates whether a point is within the interpolation range.
     *
     * @param x Point.
     * @return {@code true} if {@code x} is a valid point.
     */
    bool is_valid_point(double x) const
    {
        return !(x < knots[0] || x > knots[n]);
    }
}


