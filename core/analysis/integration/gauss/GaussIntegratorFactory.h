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

  //import org.hipparchus.dfp.DFP_Field;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Pair;

  /**
   * Class that provides different ways to compute the nodes and weights to be
   * used by the {@link Gauss_Integrator Gaussian integration rule}.
   */
class Gauss_Integrator_factory
{
	/** Number of digits for Legendre high precision. */
	public static const int DEFAULT_DECIMAL_DIGITS = 40;

	/** Generator of Gauss-Legendre integrators. */
	private const Rule_Factory legendre;
	/** Generator of Gauss-Legendre integrators. */
	private const Rule_Factory legendre_high_precision;
	/** Generator of Gauss-Hermite integrators. */
	private const Rule_Factory hermite;
	/** Generator of Gauss-Laguerre integrators. */
	private const Rule_Factory laguerre;

	/** Simple constructor.
	 */
	public Gauss_Integrator_factory()
	{
		this(DEFAULT_DECIMAL_DIGITS);
	}

	/** Simple constructor.
	 * @param decimal_digits minimum number of decimal digits for {@link #legendre_high_precisionstatic_cast<int>(}
	 */
	public Gauss_Integrator_factory(const int decimal_digits)
	{
		legendre = LegendreRule_Factory();
		legendre_high_precision = ConvertingRule_Factory<>(new FieldLegendreRule_Factory<>(new DFP_Field(decimal_digits)));
		hermite = HermiteRule_Factory();
		laguerre = LaguerreRule_Factory();
	}

	/**
	 * Creates a Gauss-Laguerre integrator of the given order.
	 * The call to the
	 * {@link Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method will perform an integration on the interval
	 * \([0, +\infty)\): the computed value is the improper integral of
	 * \(e^{-x} f(x)\)
	 * where \(f(x)\) is the function passed to the
	 * {@link Symmetric_Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method.
	 *
	 * @param number_of_points Order of the integration rule.
	 * @return a Gauss-Legendre integrator.
	 */
	public Gauss_Integrator laguerre(const int& number_of_points)
	{
		return Gauss_Integrator(laguerre.get_rule(number_of_points));
	}

	/**
	 * Creates a Gauss-Legendre integrator of the given order.
	 * The call to the
	 * {@link Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method will perform an integration on the natural interval
	 * {@code [-1 , 1]}.
	 *
	 * @param number_of_points Order of the integration rule.
	 * @return a Gauss-Legendre integrator.
	 */
	public Gauss_Integrator legendre(const int& number_of_points)
	{
		return Gauss_Integrator(legendre.get_rule(number_of_points));
	}

	/**
	 * Creates a Gauss-Legendre integrator of the given order.
	 * The call to the
	 * {@link Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method will perform an integration on the given interval.
	 *
	 * @param number_of_points Order of the integration rule.
	 * @param lower_bound Lower bound of the integration interval.
	 * @param upper_bound Upper bound of the integration interval.
	 * @return a Gauss-Legendre integrator.
	 * @ if number of points is not positive
	 */
	public Gauss_Integrator legendre(const int& number_of_points, double lower_bound, double upper_bound)

	{
		return Gauss_Integrator(transform(legendre.get_rule(number_of_points), lower_bound, upper_bound));
	}

	/**
	 * Creates a Gauss-Legendre integrator of the given order.
	 * The call to the
	 * {@link Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method will perform an integration on the natural interval
	 * {@code [-1 , 1]}.
	 *
	 * @param number_of_points Order of the integration rule.
	 * @return a Gauss-Legendre integrator.
	 * @ if number of points is not positive
	 */
	public Gauss_Integrator legendre_high_precision(const int& number_of_points)

	{
		return Gauss_Integrator(legendre_high_precision.get_rule(number_of_points));
	}

	/**
	 * Creates an integrator of the given order, and whose call to the
	 * {@link Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method will perform an integration on the given interval.
	 *
	 * @param number_of_points Order of the integration rule.
	 * @param lower_bound Lower bound of the integration interval.
	 * @param upper_bound Upper bound of the integration interval.
	 * @return a Gauss-Legendre integrator.
	 * @ if number of points is not positive
	 */
	public Gauss_Integrator legendre_high_precision(const int& number_of_points, double lower_bound, double upper_bound)

	{
		return Gauss_Integrator(transform(legendre_high_precision.get_rule(number_of_points), lower_bound, upper_bound));
	}

	/**
	 * Creates a Gauss-Hermite integrator of the given order.
	 * The call to the
	 * {@link Symmetric_Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method will perform a weighted integration on the interval
	 * \([-\infty, +\infty]\): the computed value is the improper integral of
	 * \(e^{-x^2}f(x)\)
	 * where \(f(x)\) is the function passed to the
	 * {@link Symmetric_Gauss_Integrator#integrate(org.hipparchus.analysis.Univariate_Function)
	 * integrate} method.
	 *
	 * @param number_of_points Order of the integration rule.
	 * @return a Gauss-Hermite integrator.
	 */
	public Symmetric_Gauss_Integrator hermite(const int& number_of_points)
	{
		return Symmetric_Gauss_Integrator(hermite.get_rule(number_of_points));
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
	private Pair<std::vector<double>, std::vector<double>> transform(Pair<std::vector<double>, std::vector<double>> rule, const double& a, double b)
	{
		const std::vector<double> points = rule.get_first();
		const std::vector<double> weights = rule.get_second();

		// Scaling
		const double scale = (b - a) / 2;
		const double shift = a + scale;

		for (int i{}; i < points.size(); i++)
		{
			points[i] = points[i] * scale + shift;
			weights[i] *= scale;
		}

		return Pair<>(points, weights);
	}
}
