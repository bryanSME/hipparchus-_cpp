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
 //package org.hipparchus.analysis.integration.gauss;

 //import org.hipparchus.Calculus_Field_Element;
 //import org.hipparchus.Field;
 //import org.hipparchus.exception.;
 //import org.hipparchus.util.Math_Arrays;
 //import org.hipparchus.util.Pair;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/**
 * Factory that creates a
 * <a href="http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature">
 * Gauss-type quadrature rule using Hermite polynomials</a>
 * of the first kind.
 * Such a quadrature rule allows the calculation of improper integrals
 * of a function
 * <p>
 *  \\(f(x) e^{-x^2}\\)
 * </p><p>
 * Recurrence relation and weights computation follow
 * <a href="http://en.wikipedia.org/wiki/Abramowitz_and_Stegun">
 * Abramowitz and Stegun, 1964</a>.
 * </p><p>
 * The coefficients of the standard Hermite polynomials grow very rapidly.
 * In order to avoid overflows, each Hermite polynomial is normalized with
 * respect to the underlying scalar product.
 * @param <T> Type of the number used to represent the points and weights of
 * the quadrature rules.
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldHermiteRule_Factory : public FieldAbstractRule_Factory<T>
{
public:
	/** Simple constructor
	 * @param field field to which rule coefficients belong
	 */
	FieldHermiteRule_Factory(const Field<T> field)
	{
		super(field);
	}

protected:
	/** {@inherit_doc} */
	//override
	Pair<std::vector<T>, std::vector<T>> compute_rule(const int& number_of_points)
	{
		const Field<T> field = get_field();
		const T        sqrtPi = field.get_zero().get_pi().sqrt();

		if (number_of_points == 1)
		{
			// Break recursion.
			const std::vector<T> points = Math_Arrays::build_array(field, number_of_points);
			const std::vector<T> weights = Math_Arrays::build_array(field, number_of_points);
			points[0] = field.get_zero();
			weights[0] = sqrtPi;
			return Pair<>(points, weights);
		}

		// find nodes as roots of Hermite polynomial
		const std::vector<T> points = find_roots(number_of_points, Hermite<>(field, number_of_points)::ratio);
		enforce_symmetry(points);

		// compute weights
		const std::vector<T> weights = Math_Arrays::build_array(field, number_of_points);
		const Hermite<T> hm1 = Hermite<>(field, number_of_points - 1);
		for (int i{}; i < number_of_points; i++)
		{
			const T y = hm1.hNhNm1(points[i])[0];
			weights[i] = sqrtPi.divide(y.multiply(y).multiply(number_of_points));
		}

		return Pair<>(points, weights);
	}

	/** Hermite polynomial, normalized to avoid overflow.
	 * <p>
	 * The regular Hermite polynomials and associated weights are given by:
	 *   <pre>
	 *     H\xe2\x82\x80(x)   = 1
	 *     H\xe2\x82\x81(x)   = 2 x
	 *     H\xe2\x82\x99\xe2\x82\x8a\xe2\x82\x81(x) = 2x H\xe2\x82\x99(x) - 2n H\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x), and H'\xe2\x82\x99(x) = 2n H\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x)
	 *     w\xe2\x82\x99(x\xe1\xb5\xa2) = [2\xe2\x81\xbf\xe2\x81\xbb\xc2\xb9 n! \xe2\x88\x9a\xcf\x80]/[n H\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x\xe1\xb5\xa2)]\xc2\xb2
	 *   </pre>
	 * </p>
	 * <p>
	 * In order to avoid overflow with normalize the polynomials h\xe2\x82\x99(x) = H\xe2\x82\x99(x) / \xe2\x88\x9a[2\xe2\x81\xbf n!]
	 * so the recurrence relations and weights become:
	 *   <pre>
	 *     h\xe2\x82\x80(x)   = 1
	 *     h\xe2\x82\x81(x)   = \xe2\x88\x9a2 x
	 *     h\xe2\x82\x99\xe2\x82\x8a\xe2\x82\x81(x) = [\xe2\x88\x9a2 x h\xe2\x82\x99(x) - \xe2\x88\x9an h\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x)]/\xe2\x88\x9a(n+1), and h'\xe2\x82\x99(x) = 2n h\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x)
	 *     u\xe2\x82\x99(x\xe1\xb5\xa2) = \xe2\x88\x9a\xcf\x80/[n N\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x\xe1\xb5\xa2)\xc2\xb2]
	 *   </pre>
	 * </p>
	 * @param <T> Type of the field elements.
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static class Hermite
	{
	private:
		/** \xe2\x88\x9a2. */
		const T my_sqrt2;

		/** Degree. */
		const int my_degree;

		/** Compute N\xe2\x82\x99(x) and N\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x).
		 * @param x point at which polynomials are evaluated
		 * @return array containing N\xe2\x82\x99(x) at index 0 and N\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x) at index 1
		 */
		std::vector<T> hNhNm1(const T& x)
		{
			std::vector<T> h = Math_Arrays::build_array(x.get_field(), 2);
			h[0] = sqrt2.multiply(x);
			h[1] = x.get_field().get_one();
			T sqrt_n = x.get_field().get_one();
			for (const int n{ 1 }; n < degree; n++)
			{
				// apply recurrence relation h\xe2\x82\x99\xe2\x82\x8a\xe2\x82\x81(x) = [\xe2\x88\x9a2 x h\xe2\x82\x99(x) - \xe2\x88\x9an h\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x)]/\xe2\x88\x9a(n+1)
				const T sqrt_np = x.get_field().get_zero().new_instance(n + 1).sqrt();
				const T hp = (h[0].multiply(x).multiply(sqrt2).subtract(h[1].multiply(sqrt_n))).divide(sqrt_np);
				h[1] = h[0];
				h[0] = hp;
				sqrt_n = sqrt_np;
			}
			return h;
		}

	public:
		/** Simple constructor.
		 * @param field field to which rule coefficients belong
		 * @param degree polynomial degree
		 */
		Hermite(const Field<T>& field, const int& degree)
			:
			my_sqrt2{ field.get_zero().new_instance(2).sqrt() },
			my_degree{ degree }
		{};

		/** Compute ratio H(x)/H'(x).
		 * @param x point at which ratio must be computed
		 * @return ratio H(x)/H'(x)
		 */
		T ratio(const T& x)
		{
			std::vector<T> h = hNhNm1(x);
			return h[0].divide(h[1].multiply(2 * degree));
		}
	}
};