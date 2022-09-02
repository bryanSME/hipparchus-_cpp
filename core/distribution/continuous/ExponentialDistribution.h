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

#include "../../util/MathUtils.h"

  /**
   * Implementation of the exponential distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Exponential_distribution">Exponential distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/Exponential_Distribution.html">Exponential distribution (MathWorld)</a>
   */
class Exponential_Distribution : public Abstract_Real_Distribution
{
private:
	/** The mean of this distribution. */
	const double my_mean;
	/** The logarithm of the mean, stored to reduce computing time. **/
	const double my_log_mean;

public:
	/**
	 * Create an exponential distribution with the given mean.
	 *
	 * @param mean Mean of this distribution.
	 * @ if {@code mean <= 0}.
	 */
	Exponential_Distribution(const double& mean)
		: 
		my_mean{ mean },
		my_log_mean{ std::log(mean) }
	{
		if (mean <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::MEAN, mean);
		}
	}

	/**
	 * Access the mean.
	 *
	 * @return the mean.
	 */
	double get_mean() const
	{
		return my_mean;
	}

	/** {@inherit_doc} */
	//override
	double density(const double& x)
	{
		const double log_density = Abstract_Real_Distribution::log_density(x);
		return std::isinf(log_density)
			? 0
			: std::exp(log_density);
	}

	/** {@inherit_doc} **/
	//override
	double log_density(const double& x)
	{
		if (x < 0)
		{
			return -INFINITY;
		}
		return -x / my_mean - my_log_mean;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The implementation of this method is based on:
	 * <ul>
	 * <li>
	 * <a href="http://mathworld.wolfram.com/Exponential_Distribution.html">
	 * Exponential Distribution</a>, equation (1).</li>
	 * </ul>
	 */
	 //override
	double cumulative_probability(const double& x) const
	{
		return (x <= 0.0)
			? 0.0
			: 1.0 - std::exp(-x / my_mean);
	}

	/**
	 * {@inherit_doc}
	 *
	 * Returns {@code 0} when {@code p= = 0} and
	 * {@code INFINITY} when {@code p == 1}.
	 */
	 //override
	double inverse_cumulative_probability(const double& p)
	{
		Math_Utils::check_range_inclusive(p, 0, 1);

		return (p == 1.0)
			? INFINITY
			: -my_mean * std::log(1.0 - p);
	}

	/**
	 * {@inherit_doc}
	 *
	 * For mean parameter {@code k}, the mean is {@code k}.
	 */
	 //override
	double get_numerical_mean() const
	{
		return get_mean();
	}

	/**
	 * {@inherit_doc}
	 *
	 * For mean parameter {@code k}, the variance is {@code k^2}.
	 */
	 //override
	double get_numerical_variance() const
	{
		const double m = get_mean();
		return m * m;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always 0 no matter the mean parameter.
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
	 * no matter the mean parameter.
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