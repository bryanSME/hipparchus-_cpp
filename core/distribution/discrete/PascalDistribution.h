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
  //package org.hipparchus.distribution.discrete;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.special.Beta;
  //import org.hipparchus.util.Combinatorics_Utils;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Implementation of the Pascal distribution.
   * <p>
   * The Pascal distribution is a special case of the Negative Binomial distribution
   * where the number of successes parameter is an integer.
   * <p>
   * There are various ways to express the probability mass and distribution
   * functions for the Pascal distribution. The present implementation represents
   * the distribution of the number of failures before {@code r} successes occur.
   * This is the convention adopted in e.g.
   * <a href="http://mathworld.wolfram.com/NegativeBinomial_Distribution.html">MathWorld</a>, * but <em>not</em> in
   * <a href="http://en.wikipedia.org/wiki/Negative_binomial_distribution">Wikipedia</a>.
   * <p>
   * For a random variable {@code X} whose values are distributed according to this
   * distribution, the probability mass function is given by<br/>
   * {@code P(X = k) = C(k + r - 1, r - 1) * p^r * (1 - p)^k,}<br/>
   * where {@code r} is the number of successes, {@code p} is the probability of
   * success, and {@code X} is the total number of failures. {@code C(n, k)} is
   * the binomial coefficient ({@code n} choose {@code k}). The mean and variance
   * of {@code X} are<br/>
   * {@code E(X) = (1 - p) * r / p, var(X) = (1 - p) * r / p^2.}<br/>
   * Finally, the cumulative distribution function is given by<br/>
   * {@code P(X <= k) = I(p, r, k + 1)}, * where I is the regularized incomplete Beta function.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Negative_binomial_distribution">
   * Negative binomial distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/NegativeBinomial_Distribution.html">
   * Negative binomial distribution (MathWorld)</a>
   */
class Pascal_Distribution : Abstract_Integer_Distribution
{
	
	20160320L;
	/** The number of successes. */
	private const int& number_of_successes;
	/** The probability of success. */
	private const double probability_of_success;
	/** The value of {@code log(p)}, where {@code p} is the probability of success, * stored for faster computation. */
	private const double log_probability_of_success;
	/** The value of {@code log(1-p)}, where {@code p} is the probability of success, * stored for faster computation. */
	private const double log1m_probability_of_success;

	/**
	 * Create a Pascal distribution with the given number of successes and
	 * probability of success.
	 *
	 * @param r Number of successes.
	 * @param p Probability of success.
	 * @ if the number of successes is not positive
	 * @ if the probability of success is not in the
	 * range {@code [0, 1]}.
	 */
	public Pascal_Distribution(const int& r, double p)

	{
		if (r <= 0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_SUCCESSES, r);
		}

		Math_Utils::check_range_inclusive(p, 0, 1);

		number_of_successes = r;
		probability_of_success = p;
		log_probability_of_success = std::log(p);
		log1m_probability_of_success = std::log1p(-p);
	}

	/**
	 * Access the number of successes for this distribution.
	 *
	 * @return the number of successes.
	 */
	public int get_number_of_successes()
	{
		return number_of_successes;
	}

	/**
	 * Access the probability of success for this distribution.
	 *
	 * @return the probability of success.
	 */
	public double get_probability_of_success()
	{
		return probability_of_success;
	}

	/** {@inherit_doc} */
	//override
	public double probability(const int& x)
	{
		double ret;
		if (x < 0)
		{
			ret = 0.0;
		}
		else
		{
			ret = Combinatorics_Utils.binomial_coefficient_double(x +
				number_of_successes - 1, number_of_successes - 1) *
				std::pow(probability_of_success, number_of_successes) *
				std::pow(1.0 - probability_of_success, x);
		}
		return ret;
	}

	/** {@inherit_doc} */
	//override
	public double log_probability(const int& x)
	{
		double ret;
		if (x < 0)
		{
			ret = -INFINITY;
		}
		else
		{
			ret = Combinatorics_Utils.binomial_coefficient_log(x +
				number_of_successes - 1, number_of_successes - 1) +
				log_probability_of_success * number_of_successes +
				log1m_probability_of_success * x;
		}
		return ret;
	}

	/** {@inherit_doc} */
	//override
	public double cumulative_probability(const int& x)
	{
		double ret;
		if (x < 0)
		{
			ret = 0.0;
		}
		else
		{
			ret = Beta.regularized_beta(probability_of_success, number_of_successes, x + 1.0);
		}
		return ret;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For number of successes {@code r} and probability of success {@code p}, * the mean is {@code r * (1 - p) / p}.
	 */
	 //override
	public double get_numerical_mean() const
	{
		const double p = get_probability_of_success();
		const double r = get_number_of_successes();
		return (r * (1 - p)) / p;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For number of successes {@code r} and probability of success {@code p}, * the variance is {@code r * (1 - p) / p^2}.
	 */
	 //override
	public double get_numerical_variance() const
	{
		const double p = get_probability_of_success();
		const double r = get_number_of_successes();
		return r * (1 - p) / (p * p);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always 0 no matter the parameters.
	 *
	 * @return lower bound of the support (always 0)
	 */
	 //override
	public int get_support_lower_bound()
	{
		return 0;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is always positive infinity no matter the
	 * parameters. Positive infinity is symbolized by {@code std::numeric_limits<int>::max()}.
	 *
	 * @return upper bound of the support (always {@code std::numeric_limits<int>::max()}
	 * for positive infinity)
	 */
	 //override
	public int get_support_upper_bound()
	{
		return std::numeric_limits<int>::max();
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
