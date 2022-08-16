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

  //package org.hipparchus.distribution.continuous;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.special.Beta;
  //import org.hipparchus.util.FastMath;
#include "AbstractRealDistribution.h"

  /**
   * Implementation of the F-distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/F-distribution">F-distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/F-Distribution.html">F-distribution (MathWorld)</a>
   */
class F_Distribution : Abstract_Real_Distribution
{
private:
	/** The numerator degrees of freedom. */
	const double my_numerator_degrees_of_freedom;
	/** The numerator degrees of freedom. */
	const double my_denominator_degrees_of_freedom;
	/** Cached numerical variance */
	const double my_numerical_variance;

	/**
	 * Calculates the numerical variance.
	 *
	 * @return the variance of this distribution
	 */
	double calculate_numerical_variance()
	{
		const double denominator_d_f = get_denominator_degrees_of_freedom();

		if (denominator_d_f > 4)
		{
			const double numerator_d_f = get_numerator_degrees_of_freedom();
			const double denom_d_f_minus_two = denominator_d_f - 2;

			return (2 * (denominator_d_f * denominator_d_f) * (numerator_d_f + denominator_d_f - 2)) /
				((numerator_d_f * (denom_d_f_minus_two * denom_d_f_minus_two) * (denominator_d_f - 4)));
		}

		return std::numeric_limits<double>::quiet_NaN();
	}

public:
	/**
	 * Creates an F distribution using the given degrees of freedom.
	 *
	 * @param numerator_degrees_of_freedom Numerator degrees of freedom.
	 * @param denominator_degrees_of_freedom Denominator degrees of freedom.
	 * @ if
	 * {@code numerator_degrees_of_freedom <= 0} or
	 * {@code denominator_degrees_of_freedom <= 0}.
	 */
	F_Distribution(const double& numerator_degrees_of_freedom, const double& denominator_degrees_of_freedom)
	{
		F_Distribution(numerator_degrees_of_freedom, denominator_degrees_of_freedom, DEFAULT_SOLVER_ABSOLUTE_ACCURACY);
	}

	/**
	 * Creates an F distribution.
	 *
	 * @param numerator_degrees_of_freedom Numerator degrees of freedom.
	 * @param denominator_degrees_of_freedom Denominator degrees of freedom.
	 * @param inverse_cum_accuracy the maximum absolute error in inverse
	 * cumulative probability estimates.
	 * @ if {@code numerator_degrees_of_freedom <= 0} or
	 * {@code denominator_degrees_of_freedom <= 0}.
	 */
	F_Distribution(const double& numerator_degrees_of_freedom, const double& denominator_degrees_of_freedom, const double& inverse_cum_accuracy)
		:
		my_numerator_degrees_of_freedom{ numerator_degrees_of_freedom },
		my_denominator_degrees_of_freedom{ denominator_degrees_of_freedom },
		my_numerical_variance{ calculate_numerical_variance() }
	{
		super(inverse_cum_accuracy);

		if (numerator_degrees_of_freedom <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DEGREES_OF_FREEDOM, numerator_degrees_of_freedom);
		}
		if (denominator_degrees_of_freedom <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DEGREES_OF_FREEDOM, denominator_degrees_of_freedom);
		}

	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	double density(const double& x)
	{
		return std::exp(log_density(x));
	}

	/** {@inherit_doc} **/
	//override
	double log_density(const double& x)
	{
		const double nhalf = numerator_degrees_of_freedom / 2;
		const double mhalf = denominator_degrees_of_freedom / 2;
		const double logx = std::log(x);
		const double logn = std::log(numerator_degrees_of_freedom);
		const double logm = std::log(denominator_degrees_of_freedom);
		const double lognxm = std::log(numerator_degrees_of_freedom * x + denominator_degrees_of_freedom);
		return nhalf * logn + nhalf * logx - logx +
			mhalf * logm - nhalf * lognxm - mhalf * lognxm -
			Beta::log_beta(nhalf, mhalf);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The implementation of this method is based on
	 * <ul>
	 *  <li>
	 *   <a href="http://mathworld.wolfram.com/F-Distribution.html">
	 *   F-Distribution</a>, equation (4).
	 *  </li>
	 * </ul>
	 */
	 //override
	double cumulative_probability(const double& x)
	{
		if (x <= 0)
		{
			return 0;
		}
		double n = my_numerator_degrees_of_freedom;
		double m = my_denominator_degrees_of_freedom;

		return Beta::regularized_beta((n * x) / (m + n * x), 0.5 * n, 0.5 * m);
	}

	/**
	 * Access the numerator degrees of freedom.
	 *
	 * @return the numerator degrees of freedom.
	 */
	double get_numerator_degrees_of_freedom() const
	{
		return my_numerator_degrees_of_freedom;
	}

	/**
	 * Access the denominator degrees of freedom.
	 *
	 * @return the denominator degrees of freedom.
	 */
	double get_denominator_degrees_of_freedom() const
	{
		return my_denominator_degrees_of_freedom;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For denominator degrees of freedom parameter {@code b}, the mean is
	 * <ul>
	 *  <li>if {@code b > 2} then {@code b / (b - 2)},</li>
	 *  <li>else undefined ({@codeNAN}).
	 * </ul>
	 */
	 //override
	double get_numerical_mean() const
	{
		const double denominator_d_f = get_denominator_degrees_of_freedom();

		if (denominator_d_f > 2)
		{
			return denominator_d_f / (denominator_d_f - 2);
		}

		return std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * {@inherit_doc}
	 *
	 * For numerator degrees of freedom parameter {@code a} and denominator
	 * degrees of freedom parameter {@code b}, the variance is
	 * <ul>
	 *  <li>
	 *    if {@code b > 4} then
	 *    {@code [2 * b^2 * (a + b - 2)] / [a * (b - 2)^2 * (b - 4)]}, *  </li>
	 *  <li>else undefined ({@codeNAN}).
	 * </ul>
	 */
	 //override
	double get_numerical_variance() const
	{
		return my_numerical_variance;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always 0 no matter the parameters.
	 *
	 * @return lower bound of the support (always 0)
	 */
	 //override
	double get_support_lower_bound() const
	{
		return 0;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is always positive infinity
	 * no matter the parameters.
	 *
	 * @return upper bound of the support (always INFINITY)
	 */
	 //override
	double get_support_upper_bound() const
	{
		return INFINITY;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The support of this distribution is connected.
	 *
	 * @return {@code true}
	 */
	 //override
	bool is_support_connected() const
	{
		return true;
	}
};