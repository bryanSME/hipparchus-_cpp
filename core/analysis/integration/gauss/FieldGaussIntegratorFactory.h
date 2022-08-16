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
  //package org.hipparchus.analysis.integration.gauss;

  //import org.hipparchus.Calculus_Field_Element;
  //import org.hipparchus.Field;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Pair;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/**
 * Class that provides different ways to compute the nodes and weights to be
 * used by the {@link Gauss_Integrator Gaussian integration rule}.
 * @param <T> Type of the field elements.
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldGauss_Integrator_factory
{
	/** Generator of Gauss-Legendre integrators. */
	private const FieldRule_Factory<T> legendre;
	/** Generator of Gauss-Hermite integrators. */
	private const FieldRule_Factory<T> hermite;
	/** Generator of Gauss-Laguerre integrators. */
	private const FieldRule_Factory<T> laguerre;

	/** Simple constructor.
	 * @param field field to which function argument and value belong
	 */
	public FieldGauss_Integrator_factory(const Field<T> field)
	{
		legendre = FieldLegendreRule_Factory<>(field);
		hermite = FieldHermiteRule_Factory<>(field);
		laguerre = FieldLaguerreRule_Factory<>(field);
	}

	/**
	 * Creates a Gauss-Laguerre integrator of the given order.
	 * The call to the
	 * {@link Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method will perform an integration on the interval
	 * \\([0, +\\infty)\\): the computed value is the improper integral of
	 * \\(e^{-x} f(x)\\)
	 * where \\(f(x)\\) is the function passed to the
	 * {@link Symmetric_Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method.
	 *
	 * @param number_of_points Order of the integration rule.
	 * @return a Gauss-Legendre integrator.
	 */
	public FieldGauss_Integrator<T> laguerre(const int& number_of_points)
	{
		return FieldGauss_Integrator<>(laguerre.get_rule(number_of_points));
	}

	/**
	 * Creates a Gauss-Legendre integrator of the given order.
	 * The call to the
	 * {@link FieldGauss_Integrator#integrate(org.hipparchus.analysis.Calculus_Field_Univariate_Function)
	 * integrate} method will perform an integration on the natural interval
	 * {@code [-1 , 1]}.
	 *
	 * @param number_of_points Order of the integration rule.
	 * @return a Gauss-Legendre integrator.
	 */
	public FieldGauss_Integrator<T> legendre(const int& number_of_points)
	{
		return FieldGauss_Integrator<>(legendre.get_rule(number_of_points));
	}

	/**
	 * Creates a Gauss-Legendre integrator of the given order.
	 * The call to the
	 * {@link FieldGauss_Integrator#integrate(org.hipparchus.analysis.Calculus_Field_Univariate_Function)
	 * integrate} method will perform an integration on the given interval.
	 *
	 * @param number_of_points Order of the integration rule.
	 * @param lower_bound Lower bound of the integration interval.
	 * @param upper_bound Upper bound of the integration interval.
	 * @return a Gauss-Legendre integrator.
	 * @ if number of points is not positive
	 */
	public FieldGauss_Integrator<T> legendre(const int& number_of_points, T lower_bound, T upper_bound)

	{
		return FieldGauss_Integrator<>(transform(legendre.get_rule(number_of_points), lower_bound, upper_bound));
	}

	/**
	 * Creates a Gauss-Hermite integrator of the given order.
	 * The call to the
	 * {@link Symmetric_Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method will perform a weighted integration on the interval
	 * \\([-\\infty, +\\infty]\\): the computed value is the improper integral of
	 * \\(e^{-x^2}f(x)\\)
	 * where \\(f(x)\\) is the function passed to the
	 * {@link Symmetric_Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method.
	 *
	 * @param number_of_points Order of the integration rule.
	 * @return a Gauss-Hermite integrator.
	 */
	public SymmetricFieldGauss_Integrator<T> hermite(const int& number_of_points)
	{
		return SymmetricFieldGauss_Integrator<>(hermite.get_rule(number_of_points));
	}

	/**
	 * Performs a change of variable so that the integration can be performed
	 * on an arbitrary interval {@code [a, b]}.
	 * It is assumed that the natural interval is {@code [-1, 1]}.
	 *
	 * @param rule Original points and weights.
	 * @param a Lower bound of the integration interval.
	 * @param b Lower bound of the integration interval.
	 * @return the points and weights adapted to the interval.
	 */
	private Pair<std::vector<T>, std::vector<T>> transform(Pair<std::vector<T>, std::vector<T>> rule, T a, T b)
	{
		const std::vector<T> points = rule.get_first();
		const std::vector<T> weights = rule.get_second();

		// Scaling
		const T scale = b.subtract(a).multiply(0.5);
		const T shift = a.add(scale);

		for (int i{}; i < points.size(); i++)
		{
			points[i] = points[i].multiply(scale).add(shift);
			weights[i] = weights[i].multiply(scale);
		}

		return Pair<>(points, weights);
	}
};