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

  //import org.hipparchus.exception.;
  //import org.hipparchus.random.Random_Generator;
  //import org.hipparchus.util.Pair;
#include "MixtureMultivariateRealDistribution.hpp"

/**
 * Multivariate normal mixture distribution.
 * This class is mainly syntactic sugar.
 *
 * @see Mixture_Multivariate_Real_Distribution
 */
class Mixture_Multivariate_Normal_Distribution : public Mixture_Multivariate_Real_Distribution<Multivariate_Normal_Distribution>
{
public:
	/**
	 * Creates a multivariate normal mixture distribution.
	 * <p>
	 * <b>Note:</b> this constructor will implicitly create an instance of
	 * {@link org.hipparchus.random.Well19937c Well19937c} as random
	 * generator to be used for sampling only (see {@link #sample()} and
	 * {@link #samplestatic_cast<int>(}). In case no sampling is needed for the created
	 * distribution, it is advised to pass {@code NULL} as random generator via
	 * the appropriate constructors to avoid the additional initialisation
	 * overhead.
	 *
	 * @param weights Weights of each component.
	 * @param means Mean vector for each component.
	 * @param covariances Covariance matrix for each component.
	 */
	Mixture_Multivariate_Normal_Distribution(std::vector<double> weights, std::vector<std::vector<double>> means, const std::vector<std::vector<double>>& covariances)
	{
		super(create_components(weights, means, covariances));
	}

	/**
	 * Creates a mixture model from a list of distributions and their
	 * associated weights.
	 * <p>
	 * <b>Note:</b> this constructor will implicitly create an instance of
	 * {@link org.hipparchus.random.Well19937c Well19937c} as random
	 * generator to be used for sampling only (see {@link #sample()} and
	 * {@link #samplestatic_cast<int>(}). In case no sampling is needed for the created
	 * distribution, it is advised to pass {@code NULL} as random generator via
	 * the appropriate constructors to avoid the additional initialisation
	 * overhead.
	 *
	 * @param components List of (weight, distribution) pairs from which to sample.
	 */
	Mixture_Multivariate_Normal_Distribution(std::vector<std::pair<double, Multivariate_Normal_Distribution>> components)
	{
		super(components);
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
	Mixture_Multivariate_Normal_Distribution(Random_Generator rng, std::vector<std::pair<double, Multivariate_Normal_Distribution>> components)
	{
		super(rng, components);
	}

private:
	/**
	 * @param weights Weights of each component.
	 * @param means Mean vector for each component.
	 * @param covariances Covariance matrix for each component.
	 * @return the list of components.
	 */
	static std::vector<std::pair<double, Multivariate_Normal_Distribution>> create_components(const std::vector<double>& weights, const std::vector<std::vector<double>>& means, const std::vector<std::vector<double>>& covariances)
	{
		auto mvns = std::vector<std::pair<double, Multivariate_Normal_Distribution>>(weights.size());

		for (int i{}; i < weights.size(); i++)
		{
			const auto dist = Multivariate_Normal_Distribution(means[i], covariances[i]);

			mvns.push_back(std::pair<double, Multivariate_Normal_Distribution>(weights[i], dist));
		}

		return mvns;
	}
};