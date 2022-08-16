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

  /**
   * Implementation of the hypergeometric distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Hypergeometric_distribution">Hypergeometric distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/Hypergeometric_Distribution.html">Hypergeometric distribution (MathWorld)</a>
   */
class Hypergeometric_Distribution : Abstract_Integer_Distribution
{
	/** The number of successes in the population. */
	private const int& number_of_successes;
	/** The population size. */
	private const int population_size;
	/** The sample size. */
	private const int sample_size;
	/** Cached numerical variance */
	private const double numerical_variance;

	/**
	 * Construct a hypergeometric distribution with the specified population
	 * size, number of successes in the population, and sample size.
	 *
	 * @param population_size Population size.
	 * @param number_of_successes Number of successes in the population.
	 * @param sample_size Sample size.
	 * @ if {@code number_of_successes < 0}.
	 * @ if {@code population_size <= 0}.
	 * @ if {@code number_of_successes > population_size}, * or {@code sample_size > population_size}.
	 */
	public Hypergeometric_Distribution(const int& population_size, int number_of_successes, int sample_size)

	{
		if (population_size <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::POPULATION_SIZE, population_size);
		}
		if (number_of_successes < 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_SUCCESSES, number_of_successes);
		}
		if (sample_size < 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_SAMPLES, sample_size);
		}

		if (number_of_successes > population_size)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_SUCCESS_LARGER_THAN_POPULATION_SIZE, number_of_successes, population_size, true);
		}
		if (sample_size > population_size)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SAMPLE_SIZE_LARGER_THAN_POPULATION_SIZE, sample_size, population_size, true);
		}

		this.number_of_successes = number_of_successes;
		this.population_size = population_size;
		this.sample_size = sample_size;
		this.numerical_variance = calculate_numerical_variance();
	}

	/** {@inherit_doc} */
	//override
	public double cumulative_probability(const int& x)
	{
		double ret;

		std::vector<int> domain = get_domain(population_size, number_of_successes, sample_size);
		if (x < domain[0])
		{
			ret = 0.0;
		}
		else if (x >= domain[1])
		{
			ret = 1.0;
		}
		else
		{
			ret = inner_cumulative_probability(domain[0], x, 1);
		}

		return ret;
	}

	/**
	 * Return the domain for the given hypergeometric distribution parameters.
	 *
	 * @param n Population size.
	 * @param m Number of successes in the population.
	 * @param k Sample size.
	 * @return a two element array containing the lower and upper bounds of the
	 * hypergeometric distribution.
	 */
	private std::vector<int> get_domain(const int& n, int m, const int& k)
	{
		return std::vector<int> { get_lower_domain(n, m, k), get_upper_domain(m, k) };
	}

	/**
	 * Return the lowest domain value for the given hypergeometric distribution
	 * parameters.
	 *
	 * @param n Population size.
	 * @param m Number of successes in the population.
	 * @param k Sample size.
	 * @return the lowest domain value of the hypergeometric distribution.
	 */
	private int get_lower_domain(const int& n, int m, const int& k)
	{
		return std::max(0, m - (n - k));
	}

	/**
	 * Access the number of successes.
	 *
	 * @return the number of successes.
	 */
	public int get_number_of_successes()
	{
		return number_of_successes;
	}

	/**
	 * Access the population size.
	 *
	 * @return the population size.
	 */
	public int get_population_size()
	{
		return population_size;
	}

	/**
	 * Access the sample size.
	 *
	 * @return the sample size.
	 */
	public int get_sample_size()
	{
		return sample_size;
	}

	/**
	 * Return the highest domain value for the given hypergeometric distribution
	 * parameters.
	 *
	 * @param m Number of successes in the population.
	 * @param k Sample size.
	 * @return the highest domain value of the hypergeometric distribution.
	 */
	private int get_upper_domain(const int& m, const int& k)
	{
		return std::min(k, m);
	}

	/** {@inherit_doc} */
	//override
	public double probability(const int& x)
	{
		const double log_probability = log_probability(x);
		return log_probability == -INFINITY ? 0 : std::exp(log_probability);
	}

	/** {@inherit_doc} */
	//override
	public double log_probability(const int& x)
	{
		double ret;

		std::vector<int> domain = get_domain(population_size, number_of_successes, sample_size);
		if (x < domain[0] || x > domain[1])
		{
			ret = -INFINITY;
		}
		else
		{
			double p = static_cast<double>(sample_size / static_cast<double>(population_size;
			double q = static_cast<double>((population_size - sample_size) / static_cast<double>(population_size;
			double p1 = Saddle_Point_Expansion.log_binomial_probability(x, number_of_successes, p, q);
			double p2 =
				Saddle_Point_Expansion.log_binomial_probability(sample_size - x, population_size - number_of_successes, p, q);
			double p3 =
				Saddle_Point_Expansion.log_binomial_probability(sample_size, population_size, p, q);
			ret = p1 + p2 - p3;
		}

		return ret;
	}

	/**
	 * For this distribution, {@code X}, this method returns {@code P(X >= x)}.
	 *
	 * @param x Value at which the CDF is evaluated.
	 * @return the upper tail CDF for this distribution.
	 */
	public double upper_cumulative_probability(const int& x)
	{
		double ret;

		const std::vector<int> domain = get_domain(population_size, number_of_successes, sample_size);
		if (x <= domain[0])
		{
			ret = 1.0;
		}
		else if (x > domain[1])
		{
			ret = 0.0;
		}
		else
		{
			ret = inner_cumulative_probability(domain[1], x, -1);
		}

		return ret;
	}

	/**
	 * For this distribution, {@code X}, this method returns
	 * {@code P(x0 <= X <= x1)}.
	 * This probability is computed by summing the point probabilities for the
	 * values {@code x0, x0 + 1, x0 + 2, ..., x1}, in the order directed by
	 * {@code dx}.
	 *
	 * @param x0 Inclusive lower bound.
	 * @param x1 Inclusive upper bound.
	 * @param dx Direction of summation (1 indicates summing from x0 to x1, and
	 * 0 indicates summing from x1 to x0).
	 * @return {@code P(x0 <= X <= x1)}.
	 */
	private double inner_cumulative_probability(const int& x0, int x1, int dx)
	{
		double ret = probability(x0);
		while (x0 != x1)
		{
			x0 += dx;
			ret += probability(x0);
		}
		return ret;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For population size {@code N}, number of successes {@code m}, and sample
	 * size {@code n}, the mean is {@code n * m / N}.
	 */
	 //override
	public double get_numerical_mean() const
	{
		return get_sample_size() * (get_number_of_successes() / static_cast<double>(get_population_size());
	}

	/**
	 * {@inherit_doc}
	 *
	 * For population size {@code N}, number of successes {@code m}, and sample
	 * size {@code n}, the variance is
	 * {@code [n * m * (N - n) * (N - m)] / [N^2 * (N - 1)]}.
	 */
	 //override
	public double get_numerical_variance() const
	{
		return numerical_variance;
	}

	/**
	 * Calculate the numerical variance.
	 *
	 * @return the variance of this distribution
	 */
	private double calculate_numerical_variance()
	{
		const double N = get_population_size();
		const double m = get_number_of_successes();
		const double n = get_sample_size();
		return (n * m * (N - n) * (N - m)) / (N * N * (N - 1));
	}

	/**
	 * {@inherit_doc}
	 *
	 * For population size {@code N}, number of successes {@code m}, and sample
	 * size {@code n}, the lower bound of the support is
	 * {@code max(0, n + m - N)}.
	 *
	 * @return lower bound of the support
	 */
	 //override
	public int get_support_lower_bound()
	{
		return std::max(0, get_sample_size() + get_number_of_successes() - get_population_size());
	}

	/**
	 * {@inherit_doc}
	 *
	 * For number of successes {@code m} and sample size {@code n}, the upper
	 * bound of the support is {@code min(m, n)}.
	 *
	 * @return upper bound of the support
	 */
	 //override
	public int get_support_upper_bound()
	{
		return std::min(get_number_of_successes(), get_sample_size());
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
