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
  //package org.hipparchus.distribution.multivariate;

  //import java.util.Array_list;
  //import java.util.List;

  //import org.hipparchus.distribution.Multivariate_Real_Distribution;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.random.Random_Generator;
  //import org.hipparchus.random.Well19937c;
  //import org.hipparchus.util.Pair;
#include <type_traits>
#include <vector>
#include "../../random/Well19937c.h"

/**
 * Class for representing <a href="http://en.wikipedia.org/wiki/Mixture_model">
 * mixture model</a> distributions.
 *
 * @param <T> Type of the mixture components.
 */
template<typename T, typename std::enable_if<std::is_base_of<Multivariate_Real_Distribution, T>::value>::type* = nullptr>
class Mixture_Multivariate_Real_Distribution : public AbstractMultivariate_Real_Distribution
{
private:
	/** Normalized weight of each mixture component. */
	const std::vector<double> my_weight;
	/** Mixture components. */
	const std::vector<T> my_distribution;

public:
	/**
	 * Creates a mixture model from a list of distributions and their
	 * associated weights.
	 * <p>
	 * <b>Note:</b> this constructor will implicitly create an instance of
	 * {@link Well19937c} as random generator to be used for sampling only (see
	 * {@link #sample()} and {@link #samplestatic_cast<int>(}). In case no sampling is
	 * needed for the created distribution, it is advised to pass {@code NULL}
	 * as random generator via the appropriate constructors to avoid the
	 * additional initialisation overhead.
	 *
	 * @param components List of (weight, distribution) pairs from which to sample.
	 */
	Mixture_Multivariate_Real_Distribution(const std::vector<std::pair<double, T>>& components)
	{
		Mixture_Multivariate_Real_Distribution(Well19937c(), components);
	}

	/**
	 * Creates a mixture model from a list of distributions and their
	 * associated weights.
	 *
	 * @param rng Random number generator.
	 * @param components Distributions from which to sample.
	 * @ if any of the weights is negative.
	 * @ if not all components have the same
	 * number of variables.
	 */
	Mixture_Multivariate_Real_Distribution(const Random_Generator& rng, const std::vector<std::pair<double, T>>& components)
	{
		super(rng, components.get(0).get_second().get_dimension());

		const int num_comp = components.size();
		const int dim = get_dimension();
		double weight_sum{};
		for (int i{}; i < num_comp; i++)
		{
			const Pair<Double, T> comp = components.get(i);
			if (comp.get_second().get_dimension() != dim)
			{
				throw std::exception("not implemented");
				// throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, comp.get_second().get_dimension(), dim);
			}
			if (comp.get_first() < 0)
			{
				throw std::exception("not implemented");
				// throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, comp.get_first(), 0);
			}
			weight_sum += comp.get_first();
		}

		// Check for overflow.
		if (std::isinf(weight_sum))
		{
			throw std::exception("not implemented");
			//throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW);
		}

		// Store each distribution and its normalized weight.
		distribution = std::vector < std::pair<double, T>();
		weight = std::vector<double>(num_comp);
		for (int i{}; i < num_comp; i++)
		{
			const std::pair<double, T> comp = components.get(i);
			weight[i] = comp.get_first() / weight_sum;
			distribution.add(comp.get_second());
		}
	}

	/** {@inherit_doc} */
	//override
	double density(const std::vector<double>& values)
	{
		double p{};
		for (int i{}; i < weight.size(); i++)
		{
			p += weight[i] * distribution.get(i).density(values);
		}
		return p;
	}

	/** {@inherit_doc} */
	//override
	public std::vector<double> sample()
	{
		// Sampled values.
		std::vector<double> vals;

		// Determine which component to sample from.
		const double random_value = random.next_double();
		double sum{};

		for (int i{}; i < weight.size(); i++)
		{
			sum += weight[i];
			if (random_value <= sum)
			{
				// pick model i
				vals = distribution.get(i).sample();
				break;
			}
		}

		if (vals.empty())
		{
			// This should never happen, but it ensures we won't return a NULL in
			// case the loop above has some floating point inequality problem on
			// the const iteration.
			return distribution.get(weight.size() - 1).sample();
		}

		return vals;
	}

	/** {@inherit_doc} */
	//override
	void reseed_random_generator(const long& seed)
	{
		// Seed needs to be propagated to underlying components
		// in order to maintain consistency between runs.
		super.reseed_random_generator(seed);

		for (int i{}; i < distribution.size(); i++)
		{
			// Make each component's seed different in order to avoid
			// using the same sequence of random numbers.
			distribution.get(i).reseed_random_generator(i + 1 + seed);
		}
	}

	/**
	 * Gets the distributions that make up the mixture model.
	 *
	 * @return the component distributions and associated weights.
	 */
	std::vector<std::pair<double, T>> get_components()
	{
		auto list = std::vector<std::pair<double, T>>(weight.size());

		for (int i{}; i < weight.size(); i++)
		{
			list.add(std::pair<double, T>(weight[i], distribution.get(i)));
		}

		return list;
	}
};