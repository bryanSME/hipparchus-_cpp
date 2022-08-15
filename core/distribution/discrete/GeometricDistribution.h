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
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
#include <cmath>
#include "../../util/MathUtils.h"
#include "AbstractIntegerDistribution.h"

  /**
   * Implementation of the geometric distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Geometric_distribution">Geometric distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/Geometric_Distribution.html">Geometric Distribution (MathWorld)</a>
   */
class Geometric_Distribution : Abstract_Integer_Distribution
{
private:
	/** The probability of success. */
	const double my_probability_of_success;
	/** {@code log(p)} where p is the probability of success. */
	const double my_log_probability_of_success;
	/** {@code log(1 - p)} where p is the probability of success. */
	const double my_log1m_probability_of_success;

public:
	/**
	 * Create a geometric distribution with the given probability of success.
	 *
	 * @param p probability of success.
	 * @ if {@code p <= 0} or {@code p > 1}.
	 */
	Geometric_Distribution(const double& p)
		:
		my_probability_of_success{ p },
		my_log_probability_of_success{ std::log(p) },
		my_log1m_probability_of_success{ std::log1p(-p) }
	{
		if (p <= 0 || p > 1)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_LEFT, p, 0, 1);
		}
	}

	/**
	 * Access the probability of success for this distribution.
	 *
	 * @return the probability of success.
	 */
	double get_probability_of_success() const
	{
		return my_probability_of_success;
	}

	/** {@inherit_doc} */
	//override
	double probability(const int& x) const
	{
		return x < 0
			? 0
			: std::exp(my_log1m_probability_of_success * x) * my_probability_of_success;
	}

	/** {@inherit_doc} */
	//override
	double log_probability(const int& x) const
	{
		return x < 0
			? -INFINITY
			: x * my_log1m_probability_of_success + my_log_probability_of_success;
	}

	/** {@inherit_doc} */
	//override
	double cumulative_probability(const int& x) const
	{
		return x < 0
			? 0
			: -std::expm1(my_log1m_probability_of_success * (x + 1));
	}

	/**
	 * {@inherit_doc}
	 *
	 * For probability parameter {@code p}, the mean is {@code (1 - p) / p}.
	 */
	 //override
	double get_numerical_mean() const
	{
		return (1 - my_probability_of_success) / my_probability_of_success;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For probability parameter {@code p}, the variance is
	 * {@code (1 - p) / (p * p)}.
	 */
	 //override
	double get_numerical_variance() const
	{
		return (1 - my_probability_of_success) / (my_probability_of_success * my_probability_of_success);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always 0.
	 *
	 * @return lower bound of the support (always 0)
	 */
	 //override
	int get_support_lower_bound() const
	{
		return 0;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is infinite (which we approximate as
	 * {@code std::numeric_limits<int>::max()}).
	 *
	 * @return upper bound of the support (always std::numeric_limits<int>::max())
	 */
	 //override
	int get_support_upper_bound() const
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
	bool is_support_connected() const
	{
		return true;
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	int inverse_cumulative_probability(const double& p) const
	{
		Math_Utils::check_range_inclusive(p, 0, 1);

		if (p == 1)
		{
			return std::numeric_limits<int>::max();
		}
		if (p == 0)
		{
			return 0;
		}
		return std::max(0, static_cast<int>(std::ceil(std::log1p(-p) / my_log1m_probability_of_success - 1)));
	}
};