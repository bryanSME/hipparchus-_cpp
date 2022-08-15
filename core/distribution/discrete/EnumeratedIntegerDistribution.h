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

  //import java.util.Array_list;
  //import java.util.Hash_Map;
  //import java.util.List;
  //import java.util.Map;
  //import java.util.Map.Entry;

  //import org.hipparchus.distribution.Enumerated_Distribution;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Pair;
#include <vector>
#include "AbstractIntegerDistribution.h"
#include <unordered_map>
#include "../EnumeratedDistribution.hpp"

  /**
   * Implementation of an integer-valued {@link Enumerated_Distribution}.
   * <p>
   * Values with zero-probability are allowed but they do not extend the
   * support.
   * <p>
   * Duplicate values are allowed. Probabilities of duplicate values are combined
   * when computing cumulative probabilities and statistics.
   */
class Enumerated_Integer_Distribution : public Abstract_Integer_Distribution
{
private:
	/**
	 * {@link Enumerated_Distribution} instance (using the {@link Integer} wrapper)
	 * used to generate the pmf.
	 */
	Enumerated_Distribution<int> my_inner_distribution;

	/**
	* Create the list of Pairs representing the distribution from singletons and probabilities.
	 *
	 * @param singletons values
	 * @param probabilities probabilities
	 * @return list of value/probability pairs
	 * @ if probabilities contains negative, infinite or NaN values or only 0's
	 */
	static std::vector<std::pair<int, double>> create_distribution(const std::vector<int>& singletons, const std::vector<double>& probabilities)
	{
		Math_Utils::check_dimension(singletons.size(), probabilities.size());
		auto samples = std::vector<std::pair<int, double>>(singletons.size());

		const std::vector<double> normalized_probabilities = Enumerated_Distribution<int>::check_and_normalize(probabilities);
		for (int i{}; i < singletons.size(); i++)
		{
			samples.push_back(std::pair<int, double>(singletons[i], normalized_probabilities[i]));
		}
		return samples;
	}

public:
	/**
	 * Create a discrete distribution using the given probability mass function
	 * definition.
	 *
	 * @param singletons array of random variable values.
	 * @param probabilities array of probabilities.
	 * @ if
	 * {@code singletons.size() != probabilities.size()}
	 * @ if probabilities contains negative, infinite or NaN values or only 0's
	 */
	Enumerated_Integer_Distribution(const std::vector<int>& singletons, const std::vector<double>& probabilities)
	{
		inner_distribution = Enumerated_Distribution<int>(create_distribution(singletons, probabilities));
	}

	/**
	 * Create a discrete integer-valued distribution from the input data.  Values are assigned
	 * mass based on their frequency.  For example, [0,1,1,2] as input creates a distribution
	 * with values 0, 1 and 2 having probability masses 0.25, 0.5 and 0.25 respectively, *
	 * @param data input dataset
	 */
	Enumerated_Integer_Distribution(const std::vector<int>& data)
	{
		auto data_map = std::unordered_map<int, int>();
		for (const int& value : data)
		{
			auto count = data_map.at(value);
			if (count == NULL)
			{
				count = 0;
			}
			data_map.try_emplace(value, ++count);
		}
		const int mass_points = data_map.size();
		const double denom = data.size();
		auto values = std::vector<int>(mass_points);
		auto probabilities = std::vector<double>(mass_points);
		int index{};
		for (const auto& [first, second] : data_map)
		{
			values[index] = first;
			probabilities[index] = second / denom;
			index++;
		}
		my_inner_distribution = Enumerated_Distribution<int>(create_distribution(values, probabilities));
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	double probability(const int& x)
	{
		return my_inner_distribution.probability(x);
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	double cumulative_probability(const int& x)
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
	int get_support_lower_bound()
	{
		auto min = std::numeric_limits<int>::max();
		for (const auto& [key, value] : my_inner_distribution.get_pmf())
		{
			if (key < min && value > 0)
			{
				min = key;
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
	int get_support_upper_bound()
	{
		auto max = std::numeric_limits<int>::min();
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
	std::vector<std::pair<int, double>> get_pmf()
	{
		return my_inner_distribution.get_pmf();
	}
};