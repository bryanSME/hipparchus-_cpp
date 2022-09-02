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
  //import org.hipparchus.special.Gamma;
  //import org.hipparchus.util.FastMath;

  /**
   * Implementation of Student's t-distribution.
   */
class T_Distribution : public Abstract_Real_Distribution
{
public:
	/** The degrees of freedom. */
	double my_degrees_of_freedom;
	/** Static computation factor based on degrees_of_freedom. */
	double my_factor;

	/**
	 * Create a t distribution using the given degrees of freedom.
	 *
	 * @param degrees_of_freedom Degrees of freedom.
	 * @ if {@code degrees_of_freedom <= 0}
	 */
	T_Distribution(const double& degrees_of_freedom)
	{
		T_Distribution(degrees_of_freedom, DEFAULT_SOLVER_ABSOLUTE_ACCURACY);
	}

	/**
	 * Create a t distribution using the given degrees of freedom and the
	 * specified inverse cumulative probability absolute accuracy.
	 *
	 * @param degrees_of_freedom Degrees of freedom.
	 * @param inverse_cum_accuracy the maximum absolute error in inverse
	 * cumulative probability estimates
	 * (defaults to {@link #DEFAULT_SOLVER_ABSOLUTE_ACCURACY}).
	 * @ if {@code degrees_of_freedom <= 0}
	 */
	T_Distribution(const double& degrees_of_freedom, const double& inverse_cum_accuracy)
	{
		Abstract_Real_Distribution(inverse_cum_accuracy);

		if (degrees_of_freedom <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DEGREES_OF_FREEDOM, degrees_of_freedom);
		}
		my_degrees_of_freedom = degrees_of_freedom;

		const double n = degrees_of_freedom;
		const double n_plus_1_over2 = (n + 1) / 2;
		my_factor = Gamma::log_gamma(n_plus_1_over2) - 0.5 * (std::log(std::numbers::pi) + std::log(n)) - Gamma::log_gamma(n / 2);
	}

	/**
	 * Access the degrees of freedom.
	 *
	 * @return the degrees of freedom.
	 */
	double get_degrees_of_freedom() const
	{
		return my_degrees_of_freedom;
	}

	/** {@inherit_doc} */
	//override
	double density(const double& x)
	{
		return std::exp(log_density(x));
	}

	/** {@inherit_doc} */
	//override
	double log_density(const double& x)
	{
		const double n = my_degrees_of_freedom;
		const double n_plus_1_over2 = (n + 1) / 2;
		return my_factor - n_plus_1_over2 * std::log(1 + x * x / n);
	}

	/** {@inherit_doc} */
	//override
	double cumulative_probability(const double& x)
	{
		if (x == 0)
		{
			return 0.5;
		}
		double t = Beta::regularized_beta(my_degrees_of_freedom / (my_degrees_of_freedom + (x * x)), 0.5 * my_degrees_of_freedom, 0.5);
		if (x < 0.0)
		{
			return 0.5 * t;
		}
		return 1.0 - 0.5 * t;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For degrees of freedom parameter {@code df}, the mean is
	 * <ul>
	 *  <li>if {@code df > 1} then {@code 0},</li>
	 * <li>else undefined ({@codeNAN}).</li>
	 * </ul>
	 */
	 //override
	double get_numerical_mean() const
	{
		const double df = get_degrees_of_freedom();

		if (df > 1)
		{
			return 0;
		}

		return std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * {@inherit_doc}
	 *
	 * For degrees of freedom parameter {@code df}, the variance is
	 * <ul>
	 *  <li>if {@code df > 2} then {@code df / (df - 2)},</li>
	 *  <li>if {@code 1 < df <= 2} then positive infinity
	 *  ({@code INFINITY}),</li>
	 *  <li>else undefined ({@codeNAN}).</li>
	 * </ul>
	 */
	 //override
	double get_numerical_variance() const
	{
		const double df = get_degrees_of_freedom();

		if (df > 2)
		{
			return df / (df - 2);
		}

		if (df > 1 && df <= 2)
		{
			return INFINITY;
		}

		return std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always negative infinity no matter the
	 * parameters.
	 *
	 * @return lower bound of the support (always
	 * {@code -INFINITY})
	 */
	 //override
	double get_support_lower_bound() const
	{
		return -INFINITY;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is always positive infinity no matter the
	 * parameters.
	 *
	 * @return upper bound of the support (always
	 * {@code INFINITY})
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