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
 //package org.hipparchus.analysis.polynomials;

 //import org.hipparchus.Calculus_Field_Element;
 //import org.hipparchus.Field;
 //import org.hipparchus.analysis.Calculus_Field_Univariate_Function;
 //import org.hipparchus.exception.Localized_Core_Formats;
 //import org.hipparchus.exception.;
 //import org.hipparchus.exception.;
 //import org.hipparchus.util.FastMath;
 //import org.hipparchus.util.Math_Arrays;
 //import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/**
 * Immutable representation of a real polynomial function with real coefficients.
 * <p>
 * <a href="http://mathworld.wolfram.com/Horners_method.html">Horner's Method</a>
 * is used to evaluate the function.</p>
 * @param <T> the type of the field elements
 * @since 1.5
 *
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Polynomial_Function : Calculus_Field_Univariate_Function<T>
{
	/**
	 * The coefficients of the polynomial, ordered by degree -- i.e., * coefficients[0] is the constant term and coefficients[n] is the
	 * coefficient of x^n where n is the degree of the polynomial.
	 */
	private const T coefficients[];

	/**
	 * Construct a polynomial with the given coefficients.  The first element
	 * of the coefficients array is the constant term.  Higher degree
	 * coefficients follow in sequence.  The degree of the resulting polynomial
	 * is the index of the last non-null element of the array, or 0 if all elements
	 * are NULL.
	 * <p>
	 * The constructor makes a copy of the input array and assigns the copy to
	 * the coefficients property.</p>
	 *
	 * @param c Polynomial coefficients.
	 * @ if {@code c} is {@code NULL}.
	 * @ if {@code c} is empty.
	 */
	public Field_Polynomial_Function(const T c[])
	{
		super();
		th_Utils.check_not_null(c);
		int n = c.size();
		if (n == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
		}
		while ((n > 1) && (c[n - 1].get_real() == 0))
		{
			--n;
		}
		this.coefficients = Math_Arrays::build_array(c[0].get_field(), n);
		System.arraycopy(c, 0, this.coefficients, 0, n);
	}

	/**
	 * Compute the value of the function for the given argument.
	 * <p>
	 *  The value returned is </p><p>
	 *  {@code coefficients[n] * x^n + ... + coefficients[1] * x  + coefficients[0]}
	 * </p>
	 *
	 * @param x Argument for which the function value should be computed.
	 * @return the value of the polynomial at the given point.
	 *
	 * @see org.hipparchus.analysis.Univariate_Function#valuestatic_cast<double>(
	 */
	public T value(double x)
	{
		return evaluate(coefficients, get_field().get_zero().add(x));
	}

	/**
	 * Compute the value of the function for the given argument.
	 * <p>
	 *  The value returned is </p><p>
	 *  {@code coefficients[n] * x^n + ... + coefficients[1] * x  + coefficients[0]}
	 * </p>
	 *
	 * @param x Argument for which the function value should be computed.
	 * @return the value of the polynomial at the given point.
	 *
	 * @see org.hipparchus.analysis.Univariate_Function#valuestatic_cast<double>(
	 */
	 //override
	public T value(T x)
	{
		return evaluate(coefficients, x);
	}

	/** Get the {@link Field} to which the instance belongs.
	 * @return {@link Field} to which the instance belongs
	 */
	public Field<T> get_field()
	{
		return coefficients[0].get_field();
	}

	/**
	 * Returns the degree of the polynomial.
	 *
	 * @return the degree of the polynomial.
	 */
	public int degree() const
	{
		return coefficients.size() - 1;
	}

	/**
	 * Returns a copy of the coefficients array.
	 * <p>
	 * Changes made to the returned copy will not affect the coefficients of
	 * the polynomial.</p>
	 *
	 * @return a fresh copy of the coefficients array.
	 */
	public std::vector<T> get_coefficients()
	{
		return coefficients.clone();
	}

	/**
	 * Uses Horner's Method to evaluate the polynomial with the given coefficients at
	 * the argument.
	 *
	 * @param coefficients Coefficients of the polynomial to evaluate.
	 * @param argument Input value.
	 * @param <T> the type of the field elements
	 * @return the value of the polynomial.
	 * @ if {@code coefficients} is empty.
	 * @ if {@code coefficients} is {@code NULL}.
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	protected static  T evaluate(std::vector<T> coefficients, T argument)
	{
		//Math_Utils::check_not_null(coefficients);
		int n = coefficients.size();
		if (n == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
		}
		T result = coefficients[n - 1];
		for (int j = n - 2; j >= 0; j--)
		{
			result = argument.multiply(result).add(coefficients[j]);
		}
		return result;
	}

	/**
	 * Add a polynomial to the instance.
	 *
	 * @param p Polynomial to add.
	 * @return a polynomial which is the sum of the instance and {@code p}.
	 */
	public Field_Polynomial_Function<T> add(const Field_Polynomial_Function<T> p)
	{
		// identify the lowest degree polynomial
		const int low_length = std::min(coefficients.size(), p.coefficients.size());
		const int high_length = std::max(coefficients.size(), p.coefficients.size());

		// build the coefficients array
		std::vector<T> new_coefficients = Math_Arrays::build_array(get_field(), high_length);
		for (int i{}; i < low_length; ++i)
		{
			new_coefficients[i] = coefficients[i].add(p.coefficients[i]);
		}
		System.arraycopy((coefficients.size() < p.coefficients.size()) ?
			p.coefficients : coefficients, low_length, new_coefficients, low_length, high_length - low_length);

		return Field_Polynomial_Function<>(new_coefficients);
	}

	/**
	 * Subtract a polynomial from the instance.
	 *
	 * @param p Polynomial to subtract.
	 * @return a polynomial which is the instance minus {@code p}.
	 */
	public Field_Polynomial_Function<T> subtract(const Field_Polynomial_Function<T> p)
	{
		// identify the lowest degree polynomial
		int low_length = std::min(coefficients.size(), p.coefficients.size());
		int high_length = std::max(coefficients.size(), p.coefficients.size());

		// build the coefficients array
		std::vector<T> new_coefficients = Math_Arrays::build_array(get_field(), high_length);
		for (int i{}; i < low_length; ++i)
		{
			new_coefficients[i] = coefficients[i].subtract(p.coefficients[i]);
		}
		if (coefficients.size() < p.coefficients.size())
		{
			for (int i = low_length; i < high_length; ++i)
			{
				new_coefficients[i] = p.coefficients[i].negate();
			}
		}
		else
		{
			System.arraycopy(coefficients, low_length, new_coefficients, low_length, high_length - low_length);
		}

		return Field_Polynomial_Function<>(new_coefficients);
	}

	/**
	 * Negate the instance.
	 *
	 * @return a polynomial with all coefficients negated
	 */
	public Field_Polynomial_Function<T> negate()
	{
		const std::vector<T> new_coefficients = Math_Arrays::build_array(get_field(), coefficients.size());
		for (int i{}; i < coefficients.size(); ++i)
		{
			new_coefficients[i] = coefficients[i].negate();
		}
		return Field_Polynomial_Function<>(new_coefficients);
	}

	/**
	 * Multiply the instance by a polynomial.
	 *
	 * @param p Polynomial to multiply by.
	 * @return a polynomial equal to this times {@code p}
	 */
	public Field_Polynomial_Function<T> multiply(const Field_Polynomial_Function<T> p)
	{
		const Field<T> field = get_field();
		const std::vector<T> new_coefficients = Math_Arrays::build_array(field, coefficients.size() + p.coefficients.size() - 1);

		for (int i{}; i < new_coefficients.size(); ++i)
		{
			new_coefficients[i] = field.get_zero();
			for (int j = std::max(0, i + 1 - p.coefficients.size());
				j < std::min(coefficients.size(), i + 1);
				++j)
			{
				new_coefficients[i] = new_coefficients[i].add(coefficients[j].multiply(p.coefficients[i - j]));
			}
		}

		return Field_Polynomial_Function<>(new_coefficients);
	}

	/**
	 * Returns the coefficients of the derivative of the polynomial with the given coefficients.
	 *
	 * @param coefficients Coefficients of the polynomial to differentiate.
	 * @param <T> the type of the field elements
	 * @return the coefficients of the derivative or {@code NULL} if coefficients has length 1.
	 * @ if {@code coefficients} is empty.
	 * @ if {@code coefficients} is {@code NULL}.
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	protected static  std::vector<T> differentiate(std::vector<T> coefficients)
	{
		//Math_Utils::check_not_null(coefficients);
		int n = coefficients.size();
		if (n == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
		}
		const Field<T> field = coefficients[0].get_field();
		const std::vector<T> result = Math_Arrays::build_array(field, std::max(1, n - 1));
		if (n == 1)
		{
			result[0] = field.get_zero();
		}
		else
		{
			for (int i = n - 1; i > 0; i--)
			{
				result[i - 1] = coefficients[i].multiply(i);
			}
		}
		return result;
	}

	/**
	 * Returns an anti-derivative of this polynomial, with 0 constant term.
	 *
	 * @return a polynomial whose derivative has the same coefficients as this polynomial
	 */
	public Field_Polynomial_Function<T> anti_derivative()
	{
		const Field<T> field = get_field();
		const int d = degree();
		const std::vector<T> anti = Math_Arrays::build_array(field, d + 2);
		anti[0] = field.get_zero();
		for (int i{ 1 }; i <= d + 1; i++)
		{
			anti[i] = coefficients[i - 1].multiply(1.0 / i);
		}
		return Field_Polynomial_Function<>(anti);
	}

	/**
	 * Returns the definite integral of this polymomial over the given interval.
	 * <p>
	 * [lower, upper] must describe a finite interval (neither can be infinite
	 * and lower must be less than or equal to upper).
	 *
	 * @param lower lower bound for the integration
	 * @param upper upper bound for the integration
	 * @return the integral of this polymomial over the given interval
	 * @ if the bounds do not describe a finite interval
	 */
	public T integrate(const double lower, const double upper)
	{
		const T zero = get_field().get_zero();
		return integrate(zero.add(lower), zero.add(upper));
	}

	/**
	 * Returns the definite integral of this polymomial over the given interval.
	 * <p>
	 * [lower, upper] must describe a finite interval (neither can be infinite
	 * and lower must be less than or equal to upper).
	 *
	 * @param lower lower bound for the integration
	 * @param upper upper bound for the integration
	 * @return the integral of this polymomial over the given interval
	 * @ if the bounds do not describe a finite interval
	 */
	public T integrate(const T lower, const T upper)
	{
		if (std::isinf(lower.get_real()) || std::isinf(upper.get_real()))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INFINITE_BOUND);
		}
		if (lower.get_real() > upper.get_real())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::LOWER_BOUND_NOT_BELOW_UPPER_BOUND);
		}
		const Field_Polynomial_Function<T> anti = anti_derivative();
		return anti.value(upper).subtract(anti.value(lower));
	}

	/**
	 * Returns the derivative as a {@link Field_Polynomial_Function}.
	 *
	 * @return the derivative polynomial.
	 */
	public Field_Polynomial_Function<T> polynomial_derivative()
	{
		return Field_Polynomial_Function<>(differentiate(coefficients));
	}
}
