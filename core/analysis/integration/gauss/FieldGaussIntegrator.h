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
#include "../../../CalculusFieldElement.hpp"

 /**
  * Class that : the Gaussian rule for
  * {@link #integrate(Calculus_Field_Univariate_Function) integrating} a weighted
  * function.
  *
  * @param <T> Type of the field elements.
  * @since 2.0
  */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldGauss_Integrator
{
	/** Nodes. */
	private const std::vector<T> my_points;
	/** Nodes weights. */
	private const std::vector<T> my_weights;

	/**
	 * Creates an integrator from the given {@code points} and {@code weights}.
	 * The integration interval is defined by the first and last value of
	 * {@code points} which must be sorted in increasing order.
	 *
	 * @param points Integration points.
	 * @param weights Weights of the corresponding integration nodes.
	 * @ if the {@code points} are not
	 * sorted in increasing order.
	 * @ if points and weights don't have the same length
	 */
	public FieldGauss_Integrator(const std::vector<T>& points, const std::vector<T>& weights) : my_points{ points }, my_weights{ weights }
	{
		Math_Arrays::check_equal_length(points, weights);
		Math_Arrays::check_order(points, Math_Arrays::Order_Direction::INCREASING, true, true);
	}

	/**
	 * Creates an integrator from the given pair of points (first element of
	 * the pair) and weights (second element of the pair.
	 *
	 * @param points_and_weights Integration points and corresponding weights.
	 * @ if the {@code points} are not
	 * sorted in increasing order.
	 *
	 * @see #FieldGauss_Integrator(Calculus_Field_Element[], Calculus_Field_Element[])
	 */
	public FieldGauss_Integrator(Pair<std::vector<T>, std::vector<T>> points_and_weights)

	{
		FieldGauss_Integrator(points_and_weights.get_first(), points_and_weights.get_second());
	}

	/**
	 * Returns an estimate of the integral of {@code f(x) * w(x)}, * where {@code w} is a weight function that depends on the actual
	 * flavor of the Gauss integration scheme.
	 * The algorithm uses the points and associated weights, as passed
	 * to the {@link #FieldGauss_Integrator(Calculus_Field_Element[], Calculus_Field_Element[]) constructor}.
	 *
	 * @param f Function to integrate.
	 * @return the integral of the weighted function.
	 */
	public T integrate(const Calculus_Field_Univariate_Function<T>& f)
	{
		T s = points[0].get_field().get_zero();
		T c = s;
		for (int i{}; i < points.size(); i++)
		{
			const T x = points[i];
			const T w = weights[i];
			const T y = w.multiply(f.value(x)).subtract(c);
			const T t = s.add(y);
			c = t.subtract(s).subtract(y);
			s = t;
		}
		return s;
	}

	/**
	 * @return the order of the integration rule (the number of integration
	 * points).
	 */
	public int get_number_of_points() const
	{
		return points.size();
	}

	/**
	 * Gets the integration point at the given index.
	 * The index must be in the valid range but no check is performed.
	 * @param index index of the integration point
	 * @return the integration point.
	 */
	public T get_point(const int& index) const
	{
		return points[index];
	}

	/**
	 * Gets the weight of the integration point at the given index.
	 * The index must be in the valid range but no check is performed.
	 * @param index index of the integration point
	 * @return the weight.
	 */
	public T get_weight(const int& index) const
	{
		return weights[index];
	}
};