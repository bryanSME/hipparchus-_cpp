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
#include <vector>
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

//import org.hipparchus.Field;
//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.analysis.polynomials.Field_Polynomial_Function;
//import org.hipparchus.analysis.polynomials.Field_Polynomial_Spline_Function;
//import org.hipparchus.analysis.polynomials.Polynomial_Function;
//import org.hipparchus.analysis.polynomials.Polynomial_Spline_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include "UnivariateInterpolator.h"
#include "FieldUnivariateInterpolator.h"
#include "../../util/MathArrays.h"
#include "../polynomials/PolynomialSplineFunction.h"

/**
 * Computes a natural (also known as "free", "unclamped") cubic spline interpolation for the data set.
 * <p>
 * The {@link #interpolate(std::vector<double>, std::vector<double>)} method returns a {@link Polynomial_Spline_Function}
 * consisting of n cubic polynomials, defined over the subintervals determined by the x values, * {@code x[0] < x[i] ... < x[n].}  The x values are referred to as "knot points."</p>
 * <p>
 * The value of the Polynomial_Spline_Function at a point x that is greater than or equal to the smallest
 * knot point and strictly less than the largest knot point is computed by finding the subinterval to which
 * x belongs and computing the value of the corresponding polynomial at <code>x - x[i] </code> where
 * <code>i</code> is the index of the subinterval.  See {@link Polynomial_Spline_Function} for more details.
 * </p>
 * <p>
 * The interpolating polynomials satisfy: <ol>
 * <li>The value of the Polynomial_Spline_Function at each of the input x values equals the
 *  corresponding y value.</li>
 * <li>Adjacent polynomials are equal through two derivatives at the knot points (i.e., adjacent polynomials
 *  "match up" at the knot points, as do their first and second derivatives).</li>
 * </ol></p>
 * <p>
 * The cubic spline interpolation algorithm implemented is as described in R.L. Burden, J.D. Faires, * <u>Numerical Analysis</u>, 4th Ed., 1989, PWS-Kent, ISBN 0-53491-585-X, pp 126-131.
 * </p>
 *
 */
class Spline_Interpolator : public Univariate_Interpolator, public Field_Univariate_Interpolator
{
	/**
	 * Computes an interpolating function for the data set.
	 * @param x the arguments for the interpolation points
	 * @param y the values for the interpolation points
	 * @return a function which interpolates the data set
	 * @ if {@code x} and {@code y}
	 * have different sizes.
	 * @ if {@code x} is not sorted in
	 * strict increasing order.
	 * @ if the size of {@code x} is smaller
	 * than 3.
	 */
	 //override
	public Polynomial_Spline_Function interpolate(const std::vector<double>& x, const std::vector<double>& y)
	{
		//Math_Utils::check_not_null(x);
		//Math_Utils::check_not_null(y);
		Math_Arrays::check_equal_length(x, y);
		if (x.size() < 3)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, x.size(), 3, true);
		}

		// Number of intervals.  The number of data points is n + 1.
		const int n = x.size() - 1;

		Math_Arrays::check_order(x);

		// Differences between knot points
		const double h[] = std::vector<double>(n];
		for (int i{}; i < n; i++)
		{
			h[i] = x[i + 1] - x[i];
		}

		const double mu[] = std::vector<double>(n];
		const double z[] = std::vector<double>(n + 1];
		mu[0] = 0;
		z[0] = 0;
		double g;
		for (int i{ 1 }; i < n; i++)
		{
			g = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
			mu[i] = h[i] / g;
			z[i] = (3 * (y[i + 1] * h[i - 1] - y[i] * (x[i + 1] - x[i - 1]) + y[i - 1] * h[i]) /
				(h[i - 1] * h[i]) - h[i - 1] * z[i - 1]) / g;
		}

		// cubic spline coefficients --  b is linear, c quadratic, d is cubic (original y's are constants)
		const std::vector<double>& b = std::vector<double>(n];
		const double c[] = std::vector<double>(n + 1];
		const double d[] = std::vector<double>(n];

		z[n] = 0;
		c[n] = 0;

		for (int j = n - 1; j >= 0; j--)
		{
			c[j] = z[j] - mu[j] * c[j + 1];
			b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
			d[j] = (c[j + 1] - c[j]) / (3d * h[j]);
		}

		const Polynomial_Function polynomials[] = Polynomial_Function[n];
		auto coefficients = std::vector<double>(4);
		for (int i{}; i < n; i++)
		{
			coefficients[0] = y[i];
			coefficients[1] = b[i];
			coefficients[2] = c[i];
			coefficients[3] = d[i];
			polynomials[i] = Polynomial_Function(coefficients);
		}

		return Polynomial_Spline_Function(x, polynomials);
	}

	/**
	 * Computes an interpolating function for the data set.
	 * @param x the arguments for the interpolation points
	 * @param y the values for the interpolation points
	 * @param <T> the type of the field elements
	 * @return a function which interpolates the data set
	 * @ if {@code x} and {@code y}
	 * have different sizes.
	 * @ if {@code x} is not sorted in
	 * strict increasing order.
	 * @ if the size of {@code x} is smaller
	 * than 3.
	 * @since 1.5
	 */
	 //override
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public  Field_Polynomial_Spline_Function<T> interpolate(T x[], T y[])
	{
		//Math_Utils::check_not_null(x);
		//Math_Utils::check_not_null(y);
		Math_Arrays::check_equal_length(x, y);
		if (x.size() < 3)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, x.size(), 3, true);
		}

		// Number of intervals.  The number of data points is n + 1.
		const int n = x.size() - 1;

		Math_Arrays::check_order(x);

		// Differences between knot points
		const Field<T> field = x[0].get_field();
		const std::vector<T> h = Math_Arrays::build_array(field, n);
		for (int i{}; i < n; i++)
		{
			h[i] = x[i + 1].subtract(x[i]);
		}

		std::vector<T> mu = Math_Arrays::build_array(field, n);
		std::vector<T> z = Math_Arrays::build_array(field, n + 1);
		mu[0] = field.get_zero();
		z[0] = field.get_zero();
		for (int i{ 1 }; i < n; i++)
		{
			const T g = x[i + 1].subtract(x[i - 1]).multiply(2).subtract(h[i - 1].multiply(mu[i - 1]));
			mu[i] = h[i].divide(g);
			z[i] = y[i + 1].multiply(h[i - 1]).
				subtract(y[i].multiply(x[i + 1].subtract(x[i - 1]))).
				add(y[i - 1].multiply(h[i])).
				multiply(3).
				divide(h[i - 1].multiply(h[i])).
				subtract(h[i - 1].multiply(z[i - 1])).
				divide(g);
		}

		// cubic spline coefficients --  b is linear, c quadratic, d is cubic (original y's are constants)
		const std::vector<T> b = Math_Arrays::build_array(field, n);
		std::vector<T> c = Math_Arrays::build_array(field, n + 1);
		const std::vector<T> d = Math_Arrays::build_array(field, n);

		z[n] = field.get_zero();
		c[n] = field.get_zero();

		for (int j = n - 1; j >= 0; j--)
		{
			c[j] = z[j].subtract(mu[j].multiply(c[j + 1]));
			b[j] = y[j + 1].subtract(y[j]).divide(h[j]).
				subtract(h[j].multiply(c[j + 1].add(c[j]).add(c[j])).divide(3));
			d[j] = c[j + 1].subtract(c[j]).divide(h[j].multiply(3));
		}

		//@Suppress_Warnings("unchecked")
		const Field_Polynomial_Function<T> polynomials[] =
			(Field_Polynomial_Function<T>[]) Array.new_instance(Field_Polynomial_Function.class, n);
		std::vector<T> coefficients = Math_Arrays::build_array(field, 4);
		for (int i{}; i < n; i++)
		{
			coefficients[0] = y[i];
			coefficients[1] = b[i];
			coefficients[2] = c[i];
			coefficients[3] = d[i];
			polynomials[i] = Field_Polynomial_Function<>(coefficients);
		}

		return Field_Polynomial_Spline_Function<>(x, polynomials);
	}
};