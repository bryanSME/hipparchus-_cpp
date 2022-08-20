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

#include <exception>
#include <cmath>

  /**
   * Implementation of the binomial distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Binomial_distribution">Binomial distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/Binomial_Distribution.html">Binomial Distribution (MathWorld)</a>
   */
class Binomial_Distribution : public Abstract_Integer_Distribution
{
private:
	/** The number of trials. */
	const int my_number_of_trials;
	/** The probability of success. */
	const double my_probability_of_success;

public:
	/**
	 * Create a binomial distribution with the given number of trials and
	 * probability of success.
	 *
	 * @param trials Number of trials.
	 * @param p Probability of success.
	 * @ if {@code trials < 0}.
	 * @ if {@code p < 0} or {@code p > 1}.
	 */
	Binomial_Distribution(const int& trials, const double& p)
		:
		my_probability_of_success{ p },
		my_number_of_trials{ trials }
	{
		if (trials < 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_TRIALS, trials);
		}

		Math_Utils::check_range_inclusive(p, 0, 1);
	}

	/**
	 * Access the number of trials for this distribution.
	 *
	 * @return the number of trials.
	 */
	int get_number_of_trials() const
	{
		return my_number_of_trials;
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
	double probability(const int& x)
	{
		const double _log_probability = log_probability(x);
		return _log_probability == -INFINITY 
			? 0 
			: std::exp(_log_probability);
	}

	/** {@inherit_doc} **/
	//override
	double log_probability(const int& x)
	{
		if (my_number_of_trials == 0)
		{
			return (x == 0) 
				? 0. 
				: -INFINITY;
		}
		if (x < 0 || x > my_number_of_trials)
		{
			return -INFINITY;
		}
		return Saddle_Point_Expansion::log_binomial_probability(x, my_number_of_trials, my_probability_of_success, 1.0 - my_probability_of_success);
	}

	/** {@inherit_doc} */
	//override
	double cumulative_probability(const int& x)
	{
		if (x < 0)
		{
			return 0.0;
		}
		if (x >= my_number_of_trials)
		{
			return 1.0;
		}
		return 1.0 - Beta::regularized_beta(my_probability_of_success, x + 1.0, my_number_of_trials - x);
	}

	/**
	 * {@inherit_doc}
	 *
	 * For {@code n} trials and probability parameter {@code p}, the mean is
	 * {@code n * p}.
	 */
	 //override
	double get_numerical_mean() const
	{
		return my_number_of_trials * my_probability_of_success;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For {@code n} trials and probability parameter {@code p}, the variance is
	 * {@code n * p * (1 - p)}.
	 */
	 //override
	double get_numerical_variance() const
	{
		const double p = my_probability_of_success;
		return my_number_of_trials * p * (1 - p);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always 0 except for the probability
	 * parameter {@code p = 1}.
	 *
	 * @return lower bound of the support (0 or the number of trials)
	 */
	 //override
	int get_support_lower_bound()
	{
		return my_probability_of_success < 1.0
			? 0
			: my_number_of_trials;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is the number of trials except for the
	 * probability parameter {@code p = 0}.
	 *
	 * @return upper bound of the support (number of trials or 0)
	 */
	 //override
	int get_support_upper_bound()
	{
		return my_probability_of_success > 0.0
			? my_number_of_trials
			: 0;
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