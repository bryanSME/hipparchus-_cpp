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
#include <cmath>
#include "../continuous/NormalDistribution.h"
#include "../discrete/SaddlePointExpansion.h"

  //import org.hipparchus.distribution.continuous.Normal_Distribution;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.special.Gamma;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Implementation of the Poisson distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Poisson_distribution">Poisson distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/Poisson_Distribution.html">Poisson distribution (MathWorld)</a>
   */
class Poisson_Distribution : Abstract_Integer_Distribution
{
public:
	/** Default maximum number of iterations for cumulative probability calculations. */
	static constexpr int DEFAULT_MAX_ITERATIONS{ 10000000 };
	/** Default convergence criterion. */
	static constexpr double DEFAULT_EPSILON{ 1e-12 };

private:
	/** Distribution used to compute normal approximation. */
	Normal_Distribution my_normal;
	/** Mean of the distribution. */
	double my_mean;

	/**
	 * Maximum number of iterations for cumulative probability. Cumulative
	 * probabilities are estimated using either Lanczos series approximation
	 * of {@link Gamma#regularized_gamma_p(double, double, double, int)}
	 * or continued fraction approximation of
	 * {@link Gamma#regularized_gamma_q(double, double, double, int)}.
	 */
	int my_max_iterations;

	/** Convergence criterion for cumulative probability. */
	double my_epsilon;

public:
	/**
	 * Creates a Poisson distribution with specified mean.
	 *
	 * @param p the Poisson mean
	 * @ if {@code p <= 0}.
	 */
	Poisson_Distribution(double p)
	{
		Poisson_Distribution(p, DEFAULT_EPSILON, DEFAULT_MAX_ITERATIONS);
	}

	/**
	 * Creates a Poisson distribution with specified mean, convergence
	 * criterion and maximum number of iterations.
	 *
	 * @param p Poisson mean.
	 * @param epsilon Convergence criterion for cumulative probabilities.
	 * @param max_iterations the maximum number of iterations for cumulative
	 * probabilities.
	 * @ if {@code p <= 0}.
	 */
	Poisson_Distribution(const double& p, const double& epsilon, const int& max_iterations)
	{
		if (p <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::MEAN, p);
		}
		my_mean = p;
		my_epsilon = epsilon;
		my_max_iterations = max_iterations;

		// Use the same RNG instance as the parent class.
		my_normal = Normal_Distribution(p, std::sqrt(p));
	}

	/**
	 * Creates a Poisson distribution with the specified mean and
	 * convergence criterion.
	 *
	 * @param p Poisson mean.
	 * @param epsilon Convergence criterion for cumulative probabilities.
	 * @ if {@code p <= 0}.
	 */
	Poisson_Distribution(const double& p, const double& epsilon)
	{
		Poisson_Distribution(p, epsilon, DEFAULT_MAX_ITERATIONS);
	}

	/**
	 * Creates a Poisson distribution with the specified mean and maximum
	 * number of iterations.
	 *
	 * @param p Poisson mean.
	 * @param max_iterations Maximum number of iterations for cumulative probabilities.
	 */
	Poisson_Distribution(const double& p, const int& max_iterations)
	{
		Poisson_Distribution(p, DEFAULT_EPSILON, max_iterations);
	}

	/**
	 * Get the mean for the distribution.
	 *
	 * @return the mean for the distribution.
	 */
	double get_mean() const
	{
		return my_mean;
	}

	/** {@inherit_doc} */
	//override
	double probability(const int& x)
	{
		const double log_probability = log_probability(x);
		return log_probability == -INFINITY ? 0 : std::exp(log_probability);
	}

	/** {@inherit_doc} */
	//override
	double log_probability(const int& x)
	{
		if (x < 0 || x == std::numeric_limits<int>::max())
		{
			return -INFINITY;
		}
		if (x == 0)
		{
			return -my_mean;
		}
		return -Saddle_Point_Expansion::get_stirling_error(x) -
			Saddle_Point_Expansion::get_deviance_part(x, my_mean) -
			0.5 * std::log(Math_Utils::TWO_PI) - 0.5 * std::log(x);
	}

	/** {@inherit_doc} */
	//override
	double cumulative_probability(const int& x)
	{
		if (x < 0)
		{
			return 0;
		}
		if (x == std::numeric_limits<int>::max())
		{
			return 1;
		}
		return Gamma::regularized_gamma_q(static_cast<double>(x + 1, my_mean, my_epsilon, my_max_iterations);
	}

	/**
	 * Calculates the Poisson distribution function using a normal
	 * approximation. The {@code N(mean, sqrt(mean))} distribution is used
	 * to approximate the Poisson distribution. The computation uses
	 * "half-correction" (evaluating the normal distribution function at
	 * {@code x + 0.5}).
	 *
	 * @param x Upper bound, inclusive.
	 * @return the distribution function value calculated using a normal
	 * approximation.
	 */
	double normal_approximate_probability(const int& x)
	{
		// calculate the probability using half-correction
		return my_normal.cumulative_probability(x + 0.5);
	}

	/**
	 * {@inherit_doc}
	 *
	 * For mean parameter {@code p}, the mean is {@code p}.
	 */
	 //override
	double get_numerical_mean() const
	{
		return get_mean();
	}

	/**
	 * {@inherit_doc}
	 *
	 * For mean parameter {@code p}, the variance is {@code p}.
	 */
	 //override
	double get_numerical_variance() const
	{
		return get_mean();
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always 0 no matter the mean parameter.
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
	 * The upper bound of the support is positive infinity, * regardless of the parameter values. There is no integer infinity, * so this method returns {@code std::numeric_limits<int>::max()}.
	 *
	 * @return upper bound of the support (always {@code std::numeric_limits<int>::max()} for
	 * positive infinity)
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
};