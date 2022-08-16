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
  //import org.hipparchus.special.Erf;
  //import org.hipparchus.util.FastMath;
#include <cmath>
/**
 * Implementation of the log-normal (gaussian) distribution.
 * <p>
 * <strong>Parameters:</strong>
 * {@code X} is log-normally distributed if its natural logarithm {@code log(X)}
 * is normally distributed. The probability distribution function of {@code X}
 * is given by (for {@code x > 0})
 * <p>
 * {@code exp(-0.5 * ((ln(x) - m) / s)^2) / (s * sqrt(2 * pi) * x)}
 * <ul>
 * <li>{@code m} is the <em>location</em> parameter: this is the mean of the
 * normally distributed natural logarithm of this distribution,</li>
 * <li>{@code s} is the <em>shape</em> parameter: this is the standard
 * deviation of the normally distributed natural logarithm of this
 * distribution.
 * </ul>
 *
 * @see <a href="http://en.wikipedia.org/wiki/Log-normal_distribution">
 * Log-normal distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/Log_Normal_Distribution.html">
 * Log Normal distribution (MathWorld)</a>
 */
class Log_Normal_Distribution : Abstract_Real_Distribution
{
	/** &radic;(2 &pi;) */
	private static const double SQRT2PI = std::sqrt(2 * std::numbers::pi);

	/** &radic;(2) */
	private static const double SQRT2 = std::sqrt(2.0);

	/** The location parameter of this distribution (named m in MathWorld and µ in Wikipedia). */
	private const double location;

	/** The shape parameter of this distribution. */
	private const double shape;
	/** The value of {@code log(shape) + 0.5 * log(2*PI)} stored for faster computation. */
	private const double log_shape_plus_half_log2_pi;

	/**
	 * Create a log-normal distribution, where the mean and standard deviation
	 * of the {@link Normal_Distribution normally distributed} natural
	 * logarithm of the log-normal distribution are equal to zero and one
	 * respectively. In other words, the location of the returned distribution is
	 * {@code 0}, while its shape is {@code 1}.
	 */
	public Log_Normal_Distribution()
	{
		this(0, 1);
	}

	/**
	 * Create a log-normal distribution using the specified location and shape.
	 *
	 * @param location the location parameter of this distribution
	 * @param shape the shape parameter of this distribution
	 * @ if {@code shape <= 0}.
	 */
	public Log_Normal_Distribution(double location, double shape)

	{
		this(location, shape, DEFAULT_SOLVER_ABSOLUTE_ACCURACY);
	}

	/**
	 * Creates a log-normal distribution.
	 *
	 * @param location Location parameter of this distribution.
	 * @param shape Shape parameter of this distribution.
	 * @param inverse_cum_accuracy Inverse cumulative probability accuracy.
	 * @ if {@code shape <= 0}.
	 */
	public Log_Normal_Distribution(const double& location, const double& shape, const double& inverse_cum_accuracy)
	{
		super(inverse_cum_accuracy);

		if (shape <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SHAPE, shape);
		}

		this.location = location;
		this.shape = shape;
		this.log_shape_plus_half_log2_pi = std::log(shape) + 0.5 * std::log(2 * std::numbers::pi);
	}

	/**
	 * Returns the location parameter of this distribution.
	 *
	 * @return the location parameter
	 * @since 1.4
	 */
	public double get_location()
	{
		return location;
	}

	/**
	 * Returns the shape parameter of this distribution.
	 *
	 * @return the shape parameter
	 */
	public double get_shape()
	{
		return shape;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For location {@code m}, and shape {@code s} of this distribution, the PDF
	 * is given by
	 * <ul>
	 * <li>{@code 0} if {@code x <= 0},</li>
	 * <li>{@code exp(-0.5 * ((ln(x) - m) / s)^2) / (s * sqrt(2 * pi) * x)}
	 * otherwise.</li>
	 * </ul>
	 */
	 //override
	public double density(double x)
	{
		if (x <= 0)
		{
			return 0;
		}
		const double x0 = std::log(x) - location;
		const double x1 = x0 / shape;
		return std::exp(-0.5 * x1 * x1) / (shape * SQRT2PI * x);
	}

	/** {@inherit_doc}
	 *
	 * See documentation of {@link #densitystatic_cast<double>(} for computation details.
	 */
	 //override
	public double log_density(double x)
	{
		if (x <= 0)
		{
			return -INFINITY;
		}
		const double log_x = std::log(x);
		const double x0 = log_x - location;
		const double x1 = x0 / shape;
		return -0.5 * x1 * x1 - (log_shape_plus_half_log2_pi + log_x);
	}

	/**
	 * {@inherit_doc}
	 *
	 * For location {@code m}, and shape {@code s} of this distribution, the CDF
	 * is given by
	 * <ul>
	 * <li>{@code 0} if {@code x <= 0},</li>
	 * <li>{@code 0} if {@code ln(x) - m < 0} and {@code m - ln(x) > 40 * s}, as
	 * in these cases the actual value is within {@code Double.MIN_VALUE} of 0, * <li>{@code 1} if {@code ln(x) - m >= 0} and {@code ln(x) - m > 40 * s}, * as in these cases the actual value is within {@code Double.MIN_VALUE} of 1,</li>
	 * <li>{@code 0.5 + 0.5 * erf((ln(x) - m) / (s * sqrt(2))} otherwise.</li>
	 * </ul>
	 */
	 //override
	public double cumulative_probability(const double& x)
	{
		if (x <= 0)
		{
			return 0;
		}
		const double dev = std::log(x) - location;
		if (std::abs(dev) > 40 * shape)
		{
			return dev < 0 ? 0.0 : 1.0;
		}
		return 0.5 + 0.5 * Erf.erf(dev / (shape * SQRT2));
	}

	/** {@inherit_doc} */
	//override
	public double probability(double x0, double x1)

	{
		if (x0 > x1)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT, x0, x1, true);
		}
		if (x0 <= 0 || x1 <= 0)
		{
			return super.probability(x0, x1);
		}
		const double denom = shape * SQRT2;
		const double v0 = (std::log(x0) - location) / denom;
		const double v1 = (std::log(x1) - location) / denom;
		return 0.5 * Erf.erf(v0, v1);
	}

	/**
	 * {@inherit_doc}
	 *
	 * For location {@code m} and shape {@code s}, the mean is
	 * {@code exp(m + s^2 / 2)}.
	 */
	 //override
	public double get_numerical_mean() const
	{
		double s = shape;
		return std::exp(location + (s * s / 2));
	}

	/**
	 * {@inherit_doc}
	 *
	 * For location {@code m} and shape {@code s}, the variance is
	 * {@code (exp(s^2) - 1) * exp(2 * m + s^2)}.
	 */
	 //override
	public double get_numerical_variance() const
	{
		const double s = shape;
		const double ss = s * s;
		return (std::expm1(ss)) * std::exp(2 * location + ss);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always 0 no matter the parameters.
	 *
	 * @return lower bound of the support (always 0)
	 */
	 //override
	public double get_support_lower_bound() const
	{
		return 0;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is always positive infinity
	 * no matter the parameters.
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
