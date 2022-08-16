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
 * Factory that creates Gauss-type quadrature rule using Laguerre polynomials.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature">Gauss-Laguerre quadrature (Wikipedia)</a>
 * @param <T> Type of the number used to represent the points and weights of
 * the quadrature rules.
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldLaguerreRule_Factory extends FieldAbstractRule_Factory<T>
{
	/** Simple constructor
	 * @param field field to which rule coefficients belong
	 */
	public FieldLaguerreRule_Factory(const Field<T> field)
	{
		super(field);
	}

	/** {@inherit_doc} */
	//override
	public Pair<std::vector<T>, std::vector<T>> compute_rule(const int& number_of_points)

	{
		const Field<T> field = get_field();

		// find nodes as roots of Laguerre polynomial
		const Laguerre<T> p = Laguerre<>(number_of_points);
		const std::vector<T>      points = find_roots(number_of_points, p::ratio);

		// compute weights
		const std::vector<T> weights = Math_Arrays::build_array(field, number_of_points);
		const int      n1 = number_of_points + 1;
		const long     n1Squared = n1 * static_cast<long>(n1;
		const Laguerre<T> laguerreN1 = Laguerre<>(n1);
		for (int i{}; i < number_of_points; i++)
		{
			const T y = laguerreN1.value(points[i]);
			weights[i] = points[i].divide(y.multiply(y).multiply(n1Squared));
		}

		return Pair<>(points, weights);
	}

	/** Laguerre polynomial.
	 * @param <T> Type of the field elements.
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	private static class Laguerre
	{
		/** Degree. */
		private int degree;

		/** Simple constructor.
		 * @param degree polynomial degree
		 */
		Laguerre(const int& degree)
		{
			this.degree = degree;
		}

		/** Evaluate polynomial.
		 * @param x point at which polynomial must be evaluated
		 * @return value of the polynomial
		 */
		public T value(const T x)
		{
			return lNlNm1(x)[0];
		}

		/** Compute ratio L(x)/L'(x).
		 * @param x point at which ratio must be computed
		 * @return ratio L(x)/L'(x)
		 */
		public T ratio(T x)
		{
			std::vector<T> l = lNlNm1(x);
			return x.multiply(l[0]).divide(l[0].subtract(l[1]).multiply(degree));
		}

		/** Compute L\xe2\x82\x99(x) and L\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x).
		 * @param x point at which polynomials are evaluated
		 * @return array containing L\xe2\x82\x99(x) at index 0 and L\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x) at index 1
		 */
		private std::vector<T> lNlNm1(const T x)
		{
			std::vector<T> l = Math_Arrays::build_array(x.get_field(), 2);
			l[0] = x.subtract(1).negate();
			l[1] = x.get_field().get_one();
			for (const int n{ 1 }; n < degree; n++)
			{
				// apply recurrence relation (n+1) L\xe2\x82\x99\xe2\x82\x8a\xe2\x82\x81(x) = (2n + 1 - x) L\xe2\x82\x99(x) - n L\xe2\x82\x99\xe2\x82\x8b\xe2\x82\x81(x)
				const T lp = l[0].multiply(x.negate().add(2 * n + 1)).subtract(l[1].multiply(n)).divide(n + 1);
				l[1] = l[0];
				l[0] = lp;
			}
			return l;
		}
	}
};