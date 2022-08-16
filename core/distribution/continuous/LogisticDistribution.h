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
   * This class : the Logistic distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Logistic_distribution">Logistic Distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/Logistic_Distribution.html">Logistic Distribution (Mathworld)</a>
   */
class Logistic_Distribution extends Abstract_Real_Distribution
{
	20141003L;

	/** The location parameter. */
	private const double mu;
	/** The scale parameter. */
	private const double s;

	/**
	 * Build a instance.
	 *
	 * @param mu location parameter
	 * @param s scale parameter (must be positive)
	 * @ if {@code beta <= 0}
	 */
	public Logistic_Distribution(double mu, double s)

	{
		if (s <= 0.0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_SCALE, s);
		}

		this.mu = mu;
		this.s = s;
	}

	/**
	 * Access the location parameter, {@code mu}.
	 *
	 * @return the location parameter.
	 */
	public double get_location()
	{
		return mu;
	}

	/**
	 * Access the scale parameter, {@code s}.
	 *
	 * @return the scale parameter.
	 */
	public double get_scale()
	{
		return s;
	}

	/** {@inherit_doc} */
	//override
	public double density(double x)
	{
		double z = (x - mu) / s;
		double v = std::exp(-z);
		return 1 / s * v / ((1.0 + v) * (1.0 + v));
	}

	/** {@inherit_doc} */
	//override
	public double cumulative_probability(const double& x)
	{
		double z = 1 / s * (x - mu);
		return 1.0 / (1.0 + std::exp(-z));
	}

	/** {@inherit_doc} */
	//override
	public double inverse_cumulative_probability(const double& p)
	{
		Math_Utils::check_range_inclusive(p, 0, 1);

		if (p == 0)
		{
			return 0.0;
		}
		else if (p == 1)
		{
			return INFINITY;
		}
		return s * Math.log(p / (1.0 - p)) + mu;
	}

	/** {@inherit_doc} */
	//override
	public double get_numerical_mean() const
	{
		return mu;
	}

	/** {@inherit_doc} */
	//override
	public double get_numerical_variance() const
	{
		return (Math_Utils::PI_SQUARED / 3.0) * (1.0 / (s * s));
	}

	/** {@inherit_doc} */
	//override
	public double get_support_lower_bound() const
	{
		return -INFINITY;
	}

	/** {@inherit_doc} */
	//override
	public double get_support_upper_bound() const
	{
		return INFINITY;
	}

	/** {@inherit_doc} */
	//override
	public bool is_support_connected() const
	{
		return true;
	}
}
