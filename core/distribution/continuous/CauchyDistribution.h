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
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Implementation of the Cauchy distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Cauchy_distribution">Cauchy distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/Cauchy_Distribution.html">Cauchy Distribution (MathWorld)</a>
   */
class Cauchy_Distribution : public Abstract_Real_Distribution
{
private:
	/** The median of this distribution. */
	const double my_median;
	/** The scale of this distribution. */
	const double my_scale;

public:
	/**
	 * Creates a Cauchy distribution with the median equal to zero and scale
	 * equal to one.
	 */
	Cauchy_Distribution()
	{
		Cauchy_Distribution(0, 1);
	}

	/**
	 * Creates a Cauchy distribution.
	 *
	 * @param median Median for this distribution
	 * @param scale Scale parameter for this distribution
	 * @ if {@code scale <= 0}
	 */
	Cauchy_Distribution(const double& median, const double& scale)
		:
		my_scale{ scale },
		my_median{ median }
	{
		if (scale <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SCALE, scale);
		}
	}

	/** {@inherit_doc} */
	//override
	double cumulative_probability(const double& x)
	{
		return 0.5 + (std::atan((x - my_median) / my_scale) / std::numbers::pi);
	}

	/**
	 * Access the median.
	 *
	 * @return the median for this distribution.
	 */
	double get_median() const
	{
		return my_median;
	}

	/**
	 * Access the scale parameter.
	 *
	 * @return the scale parameter for this distribution.
	 */
	double get_scale() const
	{
		return my_scale;
	}

	/** {@inherit_doc} */
	//override
	double density(const double& x)
	{
		const double dev = x - my_median;
		return (1 / std::numbers::pi) * (my_scale / (dev * dev + my_scale * my_scale));
	}

	/**
	 * {@inherit_doc}
	 *
	 * Returns {@code -INFINITY} when {@code p == 0}
	 * and {@code INFINITY} when {@code p == 1}.
	 */
	 //override
	double inverse_cumulative_probability(const double& p)
	{
		Math_Utils::check_range_inclusive(p, 0, 1);

		if (p == 0)
		{
			return -INFINITY;
		}
		if (p == 1)
		{
			return INFINITY;
		}
		return my_median + my_scale * std::tan(std::numbers::pi * (p - .5));
	}

	/**
	 * {@inherit_doc}
	 *
	 * The mean is always undefined no matter the parameters.
	 *
	 * @return mean (alwaysNAN)
	 */
	 //override
	double get_numerical_mean() const
	{
		return std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * {@inherit_doc}
	 *
	 * The variance is always undefined no matter the parameters.
	 *
	 * @return variance (alwaysNAN)
	 */
	 //override
	double get_numerical_variance() const
	{
		return std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always negative infinity no matter
	 * the parameters.
	 *
	 * @return lower bound of the support (always -INFINITY)
	 */
	 //override
	double get_support_lower_bound() const
	{
		return -INFINITY;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is always positive infinity no matter
	 * the parameters.
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