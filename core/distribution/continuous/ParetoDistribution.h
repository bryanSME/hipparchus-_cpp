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

  /**
   * Implementation of the Pareto distribution.
   * <p>
   * <strong>Parameters:</strong>
   * The probability distribution function of {@code X} is given by (for {@code x >= k}):
   * <pre>
   *  α * k^α / x^(α + 1)
   * </pre>
   * <p>
   * <ul>
   * <li>{@code k} is the <em>scale</em> parameter: this is the minimum possible value of {@code X},</li>
   * <li>{@code α} is the <em>shape</em> parameter: this is the Pareto index</li>
   * </ul>
   *
   * @see <a href="http://en.wikipedia.org/wiki/Pareto_distribution">
   * Pareto distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/Pareto_Distribution.html">
   * Pareto distribution (MathWorld)</a>
   */
class Pareto_Distribution : Abstract_Real_Distribution
{
private:

	/** The scale parameter of this distribution. */
	const double my_scale{};
	/** The shape parameter of this distribution. */
	const double my_shape{};

public:
	/**
	 * Create a Pareto distribution with a scale of {@code 1} and a shape of {@code 1}.
	 */
	Pareto_Distribution()
	{
		Pareto_Distribution(1, 1);
	}

	/**
	 * Create a Pareto distribution using the specified scale and shape.
	 *
	 * @param scale the scale parameter of this distribution
	 * @param shape the shape parameter of this distribution
	 * @ if {@code scale <= 0} or {@code shape <= 0}.
	 */
	Pareto_Distribution(const double& scale, const double& shape)
	{
		Pareto_Distribution(scale, shape, DEFAULT_SOLVER_ABSOLUTE_ACCURACY);
	}

	/**
	 * Creates a Pareto distribution.
	 *
	 * @param scale Scale parameter of this distribution.
	 * @param shape Shape parameter of this distribution.
	 * @param inverse_cum_accuracy Inverse cumulative probability accuracy.
	 * @ if {@code scale <= 0} or {@code shape <= 0}.
	 */
	Pareto_Distribution(const double& scale, const double& shape, const double& inverse_cum_accuracy)
		:
		my_scale{ scale },
		my_shape{ shape }
	{
		super(inverse_cum_accuracy);

		if (scale <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SCALE, scale);
		}

		if (shape <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SHAPE, shape);
		}
	}

	/**
	 * Returns the scale parameter of this distribution.
	 *
	 * @return the scale parameter
	 */
	double get_scale() const
	{
		return my_scale;
	}

	/**
	 * Returns the shape parameter of this distribution.
	 *
	 * @return the shape parameter
	 */
	double get_shape() const
	{
		return my_shape;
	}

	/**
	 * {@inherit_doc}
	 * <p>
	 * For scale {@code k}, and shape {@code α} of this distribution, the PDF
	 * is given by
	 * <ul>
	 * <li>{@code 0} if {@code x < k},</li>
	 * <li>{@code α * k^α / x^(α + 1)} otherwise.</li>
	 * </ul>
	 */
	 //override
	double density(const double& x) const
	{
		if (x < my_scale)
		{
			return 0;
		}
		return std::pow(my_scale, my_shape) / std::pow(x, my_shape + 1) * my_shape;
	}

	/** {@inherit_doc}
	 *
	 * See documentation of {@link #densitystatic_cast<double>(} for computation details.
	 */
	 //override
	double log_density(const double& x) const
	{
		if (x < my_scale)
		{
			return -INFINITY;
		}
		return std::log(my_scale) * my_shape - std::log(x) * (my_shape + 1) + std::log(my_shape);
	}

	/**
	 * {@inherit_doc}
	 * <p>
	 * For scale {@code k}, and shape {@code α} of this distribution, the CDF is given by
	 * <ul>
	 * <li>{@code 0} if {@code x < k},</li>
	 * <li>{@code 1 - (k / x)^α} otherwise.</li>
	 * </ul>
	 */
	 //override
	double cumulative_probability(const double& x) const
	{
		if (x <= my_scale)
		{
			return 0;
		}
		return 1 - std::pow(my_scale / x, my_shape);
	}

	/**
	 * {@inherit_doc}
	 * <p>
	 * For scale {@code k} and shape {@code α}, the mean is given by
	 * <ul>
	 * <li>{@code ∞} if {@code α <= 1},</li>
	 * <li>{@code α * k / (α - 1)} otherwise.</li>
	 * </ul>
	 */
	 //override
	double get_numerical_mean() const
	{
		if (my_shape <= 1)
		{
			return INFINITY;
		}
		return my_shape * my_scale / (my_shape - 1);
	}

	/**
	 * {@inherit_doc}
	 * <p>
	 * For scale {@code k} and shape {@code α}, the variance is given by
	 * <ul>
	 * <li>{@code ∞} if {@code 1 < α <= 2},</li>
	 * <li>{@code k^2 * α / ((α - 1)^2 * (α - 2))} otherwise.</li>
	 * </ul>
	 */
	 //override
	double get_numerical_variance() const
	{
		if (my_shape <= 2)
		{
			return INFINITY;
		}
		double s = my_shape - 1;
		return my_scale * my_scale * my_shape / (s * s) / (my_shape - 2);
	}

	/**
	 * {@inherit_doc}
	 * <p>
	 * The lower bound of the support is equal to the scale parameter {@code k}.
	 *
	 * @return lower bound of the support
	 */
	 //override
	double get_support_lower_bound() const
	{
		return my_scale;
	}

	/**
	 * {@inherit_doc}
	 * <p>
	 * The upper bound of the support is always positive infinity no matter the parameters.
	 *
	 * @return upper bound of the support (always {@code INFINITY})
	 */
	 //override
	double get_support_upper_bound() const
	{
		return INFINITY;
	}

	/**
	 * {@inherit_doc}
	 * <p>
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