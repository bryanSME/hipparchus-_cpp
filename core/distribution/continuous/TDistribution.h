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
class T_Distribution extends Abstract_Real_Distribution
{
	/** Serializable version identifier */
	20160320L;
	/** The degrees of freedom. */
	private const double degrees_of_freedom;
	/** Static computation factor based on degrees_of_freedom. */
	private const double factor;

	/**
	 * Create a t distribution using the given degrees of freedom.
	 *
	 * @param degrees_of_freedom Degrees of freedom.
	 * @ if {@code degrees_of_freedom <= 0}
	 */
	public T_Distribution(double degrees_of_freedom)

	{
		this(degrees_of_freedom, DEFAULT_SOLVER_ABSOLUTE_ACCURACY);
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
	public T_Distribution(double degrees_of_freedom, double inverse_cum_accuracy)

	{
		super(inverse_cum_accuracy);

		if (degrees_of_freedom <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DEGREES_OF_FREEDOM, degrees_of_freedom);
		}
		this.degrees_of_freedom = degrees_of_freedom;

		const double n = degrees_of_freedom;
		const double n_plus_1_over2 = (n + 1) / 2;
		factor = Gamma::log_gamma(n_plus_1_over2) -
			0.5 * (std::log(std::numbers::pi) + std::log(n)) -
			Gamma::log_gamma(n / 2);
	}

	/**
	 * Access the degrees of freedom.
	 *
	 * @return the degrees of freedom.
	 */
	public double get_degrees_of_freedom()
	{
		return degrees_of_freedom;
	}

	/** {@inherit_doc} */
	//override
	public double density(double x)
	{
		return std::exp(log_density(x));
	}

	/** {@inherit_doc} */
	//override
	public double log_density(double x)
	{
		const double n = degrees_of_freedom;
		const double n_plus_1_over2 = (n + 1) / 2;
		return factor - n_plus_1_over2 * std::log(1 + x * x / n);
	}

	/** {@inherit_doc} */
	//override
	public double cumulative_probability(const double& x)
	{
		double ret;
		if (x == 0)
		{
			ret = 0.5;
		}
		else
		{
			double t =
				Beta.regularized_beta(
					degrees_of_freedom / (degrees_of_freedom + (x * x)), 0.5 * degrees_of_freedom, 0.5);
			if (x < 0.0)
			{
				ret = 0.5 * t;
			}
			else
			{
				ret = 1.0 - 0.5 * t;
			}
		}

		return ret;
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
	public double get_numerical_mean() const
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
	public double get_numerical_variance() const
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
	public double get_support_lower_bound() const
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
	public double get_support_upper_bound() const
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
	public bool is_support_connected() const
	{
		return true;
	}
}
