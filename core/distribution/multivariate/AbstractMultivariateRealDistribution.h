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

  //import org.hipparchus.distribution.Multivariate_Real_Distribution;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.random.Random_Generator;
#include  <vector>
#include "../MultivariateRealDistribution.h"
#include "../../random/RandomGenerator.h"

/**
 * Base class for multivariate probability distributions.
 */
class AbstractMultivariate_Real_Distribution : public Multivariate_Real_Distribution
{
private:
	/** The number of dimensions or columns in the multivariate distribution. */
	int my_dimension;

protected:
	/** RNG instance used to generate samples from the distribution. */
	Random_Generator my_random;

	/**
	 * @param rng Random number generator.
	 * @param n Number of dimensions.
	 */
	AbstractMultivariate_Real_Distribution(const Random_Generator& rng, const int& n)
		:
		my_random{ rng },
		my_dimension{ n }
	{};

public:
	/** {@inherit_doc} */
	//override
	void reseed_random_generator(const long& seed)
	{
		my_random.set_seed(seed);
	}

	/** {@inherit_doc} */
	//override
	int get_dimension() const
	{
		return my_dimension;
	}

	/** {@inherit_doc} */
	//override
	virtual std::vector<double> sample();

	/** {@inherit_doc} */
	//override
	std::vector<std::vector<double>> sample(const int& sample_size)
	{
		if (sample_size <= 0)
		{
			throw std::exception("Not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_SAMPLES, sample_size);
		}
		auto out = std::vector<std::vector<double>>(sample_size, std::vector<double>(my_dimension);
		for (int i{}; i < sample_size; i++)
		{
			out[i] = sample();
		}
		return out;
	}
};