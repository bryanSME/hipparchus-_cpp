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

  /**
   * Implementation of the uniform integer distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Uniform_distribution_(discrete)">
   * Uniform distribution (discrete), at Wikipedia</a>
   */
class UniformInteger_Distribution : Abstract_Integer_Distribution
{
	/** Lower bound (inclusive) of this distribution. */
	private const int lower;
	/** Upper bound (inclusive) of this distribution. */
	private const int upper;

	/**
	 * Creates a uniform integer distribution using the given lower and
	 * upper bounds (both inclusive).
	 *
	 * @param lower Lower bound (inclusive) of this distribution.
	 * @param upper Upper bound (inclusive) of this distribution.
	 * @ if {@code lower >= upper}.
	 */
	public UniformInteger_Distribution(const int& lower, int upper)

	{
		if (lower > upper)
		{
			throw (
				hipparchus::exception::Localized_Core_Formats_Type::LOWER_BOUND_NOT_BELOW_UPPER_BOUND, lower, upper, true);
		}
		this.lower = lower;
		this.upper = upper;
	}

	/** {@inherit_doc} */
	//override
	public double probability(const int& x)
	{
		if (x < lower || x > upper)
		{
			return 0;
		}
		return 1.0 / (upper - lower + 1);
	}

	/** {@inherit_doc} */
	//override
	public double cumulative_probability(const int& x)
	{
		if (x < lower)
		{
			return 0;
		}
		if (x > upper)
		{
			return 1;
		}
		return (x - lower + 1.0) / (upper - lower + 1.0);
	}

	/**
	 * {@inherit_doc}
	 *
	 * For lower bound {@code lower} and upper bound {@code upper}, the mean is
	 * {@code 0.5 * (lower + upper)}.
	 */
	 //override
	public double get_numerical_mean() const
	{
		return 0.5 * (lower + upper);
	}

	/**
	 * {@inherit_doc}
	 *
	 * For lower bound {@code lower} and upper bound {@code upper}, and
	 * {@code n = upper - lower + 1}, the variance is {@code (n^2 - 1) / 12}.
	 */
	 //override
	public double get_numerical_variance() const
	{
		double n = upper - lower + 1;
		return (n * n - 1) / 12.0;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is equal to the lower bound parameter
	 * of the distribution.
	 *
	 * @return lower bound of the support
	 */
	 //override
	public int get_support_lower_bound()
	{
		return lower;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is equal to the upper bound parameter
	 * of the distribution.
	 *
	 * @return upper bound of the support
	 */
	 //override
	public int get_support_upper_bound()
	{
		return upper;
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
