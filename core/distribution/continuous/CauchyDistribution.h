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
class Cauchy_Distribution extends Abstract_Real_Distribution
{
	/** Serializable version identifier */
	20160320L;
	/** The median of this distribution. */
	private const double median;
	/** The scale of this distribution. */
	private const double scale;

	/**
	 * Creates a Cauchy distribution with the median equal to zero and scale
	 * equal to one.
	 */
	public Cauchy_Distribution()
	{
		this(0, 1);
	}

	/**
	 * Creates a Cauchy distribution.
	 *
	 * @param median Median for this distribution
	 * @param scale Scale parameter for this distribution
	 * @ if {@code scale <= 0}
	 */
	public Cauchy_Distribution(double median, double scale)

	{
		if (scale <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SCALE, scale);
		}

		this.scale = scale;
		this.median = median;
	}

	/** {@inherit_doc} */
	//override
	public double cumulative_probability(const double& x)
	{
		return 0.5 + (std::atan((x - median) / scale) / std::numbers::pi);
	}

	/**
	 * Access the median.
	 *
	 * @return the median for this distribution.
	 */
	public double get_median()
	{
		return median;
	}

	/**
	 * Access the scale parameter.
	 *
	 * @return the scale parameter for this distribution.
	 */
	public double get_scale()
	{
		return scale;
	}

	/** {@inherit_doc} */
	//override
	public double density(double x)
	{
		const double dev = x - median;
		return (1 / std::numbers::pi) * (scale / (dev * dev + scale * scale));
	}

	/**
	 * {@inherit_doc}
	 *
	 * Returns {@code -INFINITY} when {@code p == 0}
	 * and {@code INFINITY} when {@code p == 1}.
	 */
	 //override
	public double inverse_cumulative_probability(const double& p)
	{
		Math_Utils::check_range_inclusive(p, 0, 1);

		double ret;
		if (p == 0)
		{
			ret = -INFINITY;
		}
		else  if (p == 1)
		{
			ret = INFINITY;
		}
		else
		{
			ret = median + scale * std::tan(std::numbers::pi * (p - .5));
		}
		return ret;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The mean is always undefined no matter the parameters.
	 *
	 * @return mean (alwaysNAN)
	 */
	 //override
	public double get_numerical_mean() const
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
	public double get_numerical_variance() const
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
	public double get_support_lower_bound() const
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
