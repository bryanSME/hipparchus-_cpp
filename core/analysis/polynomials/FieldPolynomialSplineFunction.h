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
  * @param <T> the type of the field elements
  * @since 1.5
  */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Polynomial_Spline_Function : Calculus_Field_Univariate_Function<T>
{
private:
	/**
	 * Spline segment interval delimiters (knots).
	 * Size is n + 1 for n segments.
	 */
	const std::vector<T> my_knots;

	/**
	 * The polynomial functions that make up the spline.  The first element
	 * determines the value of the spline over the first subinterval, the
	 * second over the second, etc.   Spline function values are determined by
	 * evaluating these functions at {@code (x - knot[i])} where i is the
	 * knot segment to which x belongs.
	 */
	const std::vector<Field_Polynomial_Function<T>> polynomials;

	/**
	 * Number of spline segments. It is equal to the number of polynomials and
	 * to the number of partition points - 1.
	 */
	const int my_n;

public:
	/**
	 * Construct a polynomial spline function with the given segment delimiters
	 * and interpolating polynomials.
	 * The constructor copies both arrays and assigns the copies to the knots
	 * and polynomials properties, respectively.
	 *
	 * @param knots Spline segment interval delimiters.
	 * @param polynomials Polynomial functions that make up the spline.
	 * @ if either of the input arrays is {@code NULL}.
	 * @ if knots has length less than 2.
	 * @ if {@code polynomials.size() != knots.size() - 1}.
	 * @ if the {@code knots} array is not strictly increasing.
	 *
	 */
	 //@Suppress_Warnings("unchecked")
	Field_Polynomial_Spline_Function(const std::vector<T>& knots, const Field_Polynomial_Function<T> polynomials[])
	{
		throw std::exception("not implmented");
		//if (knots == NULL || polynomials == NULL)
		//{
		//    throw std::exception("not implemented");
		//    //throw ();
		//}
		//if (knots.size() < 2)
		//{
		//    throw std::exception("not implemented");
		//    //throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_ENOUGH_POINTS_IN_SPLINE_PARTITION, 2, knots.size(), false);
		//}
		//Math_Utils::check_dimension(polynomials.size(), knots.size() - 1);
		//Math_Arrays::check_order(knots);

		//my_n = knots.size() -1;
		//my_knots = knots;
		//my_polynomials = (std::vector<Field_Polynomial_Function<T>>) Array.new_instance(Field_Polynomial_Function.class, n);
		//System.arraycopy(polynomials, 0, this.polynomials, 0, n);
	}

	/** Get the {@link Field} to which the instance belongs.
	 * @return {@link Field} to which the instance belongs
	 */
	Field<T> get_field()
	{
		return knots[0].get_field();
	}

	/**
	 * Compute the value for the function.
	 * See {@link Field_Polynomial_Spline_Function} for details on the algorithm for
	 * computing the value of the function.
	 *
	 * @param v Point for which the function value should be computed.
	 * @return the value.
	 * @ if {@code v} is outside of the domain of the
	 * spline function (smaller than the smallest knot point or larger than the
	 * largest knot point).
	 */
	T value(const double& v)
	{
		return value(get_field().get_zero().add(v));
	}

	/**
	 * Compute the value for the function.
	 * See {@link Field_Polynomial_Spline_Function} for details on the algorithm for
	 * computing the value of the function.
	 *
	 * @param v Point for which the function value should be computed.
	 * @return the value.
	 * @ if {@code v} is outside of the domain of the
	 * spline function (smaller than the smallest knot point or larger than the
	 * largest knot point).
	 */
	 //override
	T value(const T& v)
	{
		Math_Utils::check_range_inclusive(v.get_real(), knots[0].get_real(), knots[n].get_real());
		int i = Arrays.binary_search(knots, v);
		if (i < 0)
		{
			i = -i - 2;
		}
		// This will handle the case where v is the last knot value
		// There are only n-1 polynomials, so if v is the last knot
		// then we will use the last polynomial to calculate the value.
		if (i >= polynomials.size())
		{
			i--;
		}
		return polynomials[i].value(v.subtract(knots[i]));
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
	std::vector<Field_Polynomial_Function<T>> get_polynomials() const
	{
		return polynomials;
	}

	/**
	 * Get an array copy of the knot points.
	 * It returns a fresh copy of the array. Changes made to the copy
	 * will not affect the knots property.
	 *
	 * @return the knot points.
	 */
	std::vector<T> get_knots() const
	{
		return my_knots;
	}

	/**
	 * Indicates whether a point is within the interpolation range.
	 *
	 * @param x Point.
	 * @return {@code true} if {@code x} is a valid point.
	 */
	bool is_valid_point(const T& x) const
	{
		return !(x.get_real() < knots[0].get_real() || x.get_real() > knots[n].get_real());
	}
	/**
	 * Get the derivative of the polynomial spline function.
	 *
	 * @return the derivative function.
	 */
	 //@Suppress_Warnings("unchecked")
	Field_Polynomial_Spline_Function<T> polynomial_spline_derivative()
	{
		Field_Polynomial_Function<T> derivative_polynomials[] = (Field_Polynomial_Function<T>[]) Array.new_instance(Field_Polynomial_Function.class, n);
		for (int i{}; i < n; i++)
		{
			derivative_polynomials[i] = polynomials[i].polynomial_derivative();
		}
		return Field_Polynomial_Spline_Function<>(knots, derivative_polynomials);
	}
};