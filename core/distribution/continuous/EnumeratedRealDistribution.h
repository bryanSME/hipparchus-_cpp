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

  //import java.util.Array_list;
  //import java.util.Hash_Map;
  //import java.util.List;
  //import java.util.Map;
  //import java.util.Map.Entry;

  //import org.hipparchus.distribution.Enumerated_Distribution;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.std::pair;
#include "AbstractRealDistribution.h"
#include "../EnumeratedDistribution.hpp"
#include "../../util/MathArrays.h"
#include "../../util/MathUtils.h"

/**
 * Implementation of a real-valued {@link Enumerated_Distribution}.
 * <p>
 * Values with zero-probability are allowed but they do not extend the
 * support.
 * <p>
 * Duplicate values are allowed. Probabilities of duplicate values are
 * combined when computing cumulative probabilities and statistics.
 */
class Enumerated_Real_Distribution : Abstract_Real_Distribution
{
private:
	/**
	 * {@link Enumerated_Distribution} (using the {@link double} wrapper)
	 * used to generate the pmf.
	 */
	Enumerated_Distribution<double> my_inner_distribution;

	/**
	 * Create the list of std::pairs representing the distribution from singletons and probabilities.
	 *
	 * @param singletons values
	 * @param probabilities probabilities
	 * @return list of value/probability pairs
	 * @ if probabilities contains negative, infinite or NaN values or only 0's
	 */
	static std::vector<std::pair<double, double>> create_distribution(const std::vector<double>& singletons, const std::vector<double>& probabilities)
	{
		Math_Arrays::check_equal_length(singletons, probabilities);
		auto samples = std::vector<double, double>(singletons.size());

		const std::vector<double> normalized_probabilities = Enumerated_Distribution::check_and_normalize(probabilities);
		for (int i{}; i < singletons.size(); i++)
		{
			samples.push_back(std::pair<double, double>(singletons[i], normalized_probabilities[i]));
		}
		return samples;
	}

public:
	/**
	 * Create a discrete real-valued distribution from the input data.  Values are assigned
	 * mass based on their frequency.  For example, [0,1,1,2] as input creates a distribution
	 * with values 0, 1 and 2 having probability masses 0.25, 0.5 and 0.25 respectively, *
	 * @param data input dataset
	 */
	Enumerated_Real_Distribution(const std::vector<double>& data)
	{
		super();
		auto data_map = std::unordered_map<double, int>(data.size());
		for (double value : data)
		{
			int count = data_map.at(value);
			if (count == NULL)
			{
				count = 0;
			}
			data_map.try_emplace(value, ++count);
		}
		const int mass_points = data_map.size();
		const double denom = data.size();
		auto values = std::vector<double>(mass_points);
		auto probabilities = std::vector<double>(mass_points);
		int index{};
		for (const auto& [key, value] : data_map)
		{
			values[index] = key;
			probabilities[index] = value / denom;
			index++;
		}
		my_inner_distribution = Enumerated_Distribution<double>(create_distribution(values, probabilities));
	}

	/**
	 * Create a discrete real-valued distribution using the given probability mass function
	 * enumeration.
	 *
	 * @param singletons array of random variable values.
	 * @param probabilities array of probabilities.
	 * @ if
	 * {@code singletons.size() != probabilities.size()}
	 * @ if any of the probabilities are negative.
	 * @ if any of the probabilities are NaN.
	 * @ if any of the probabilities are infinite.
	 */
	Enumerated_Real_Distribution(const std::vector<double> singletons, const std::vector<double> probabilities)
	{
		super();
		my_inner_distribution = Enumerated_Distribution<double>(create_distribution(singletons, probabilities));
	}

	/**
	 * For a random variable {@code X} whose values are distributed according to
	 * this distribution, this method returns {@code P(X = x)}. In other words, * this method represents the probability mass function (PMF) for the
	 * distribution.
	 * <p>
	 * Note that if {@code x1} and {@code x2} satisfy {@code x1.equals(x2)}, * or both are NULL, then {@code probability(x1) = probability(x2)}.
	 *
	 * @param x the point at which the PMF is evaluated
	 * @return the value of the probability mass function at {@code x}
	 */
	double probability(const double& x)
	{
		return my_inner_distribution.probability(x);
	}

	/**
	 * For a random variable {@code X} whose values are distributed according to
	 * this distribution, this method returns {@code P(X = x)}. In other words, * this method represents the probability mass function (PMF) for the
	 * distribution.
	 *
	 * @param x the point at which the PMF is evaluated
	 * @return the value of the probability mass function at point {@code x}
	 */
	 //override
	double density(const double& x)
	{
		return probability(x);
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	double cumulative_probability(const double& x)
	{
		double probability{};

		for (const auto& [key, value] : my_inner_distribution.get_pmf())
		{
			if (key <= x)
			{
				probability += value;
			}
		}

		return probability;
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	double inverse_cumulative_probability(const double& p)
	{
		Math_Utils::check_range_inclusive(p, 0, 1);

		double probability{};
		double x = get_support_lower_bound();
		for (const auto& [key, value] : my_inner_distribution.get_pmf())
		{
			if (value == 0.0)
			{
				continue;
			}

			probability += value;
			x = key;

			if (probability >= p)
			{
				break;
			}
		}

		return x;
	}

	/**
	 * {@inherit_doc}
	 *
	 * @return {@code sum(singletons[i] * probabilities[i])}
	 */
	 //override
	double get_numerical_mean() const
	{
		double mean{};

		for (const auto& [key, value] : my_inner_distribution.get_pmf())
		{
			mean += value * key;
		}

		return mean;
	}

	/**
	 * {@inherit_doc}
	 *
	 * @return {@code sum((singletons[i] - mean) ^ 2 * probabilities[i])}
	 */
	 //override
	double get_numerical_variance() const
	{
		double mean{};
		double mean_of_squares{};

		for (const auto& [key, value] : my_inner_distribution.get_pmf())
		{
			mean += value * key;
			mean_of_squares += value * key * key;
		}

		return mean_of_squares - mean * mean;
	}

	/**
	 * {@inherit_doc}
	 *
	 * Returns the lowest value with non-zero probability.
	 *
	 * @return the lowest value with non-zero probability.
	 */
	 //override
	double get_support_lower_bound() const
	{
		double min = INFINITY;
		for (const auto& sample : my_inner_distribution.get_pmf())
		{
			if (sample.get_key() < min && sample.get_value() > 0)
			{
				min = sample.get_key();
			}
		}
		return min;
	}

	/**
	 * {@inherit_doc}
	 *
	 * Returns the highest value with non-zero probability.
	 *
	 * @return the highest value with non-zero probability.
	 */
	 //override
	double get_support_upper_bound() const
	{
		double max = -INFINITY;
		for (const auto& [key, value] : my_inner_distribution.get_pmf())
		{
			if (key > max && value > 0)
			{
				max = key;
			}
		}

		return max;
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
	 * Return the probability mass function as a list of (value, probability) pairs.
	 *
	 * @return the probability mass function.
	 */
	std::vector<std::pair<double, double>> get_pmf()
	{
		return my_inner_distribution.get_pmf();
	}
};