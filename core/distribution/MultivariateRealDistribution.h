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

#include <vector>

  /**
   * Base interface for multivariate continuous distributions.
   * <p>
   * This is based largely on the Real_Distribution interface, but cumulative
   * distribution functions are not required because they are often quite
   * difficult to compute for multivariate distributions.
   */
class Multivariate_Real_Distribution
{
	/**
	 * Returns the probability density function (PDF) of this distribution
	 * evaluated at the specified point {@code x}. In general, the PDF is the
	 * derivative of the cumulative distribution function. If the derivative
	 * does not exist at {@code x}, then an appropriate replacement should be
	 * returned, e.g. {@code INFINITY}, {@codeNAN}, or
	 * the limit inferior or limit superior of the difference quotient.
	 *
	 * @param x Point at which the PDF is evaluated.
	 * @return the value of the probability density function at point {@code x}.
	 */
	virtual double density(std::vector<double> x) = 0;

	/**
	 * Reseeds the random generator used to generate samples.
	 *
	 * @param seed Seed with which to initialize the random number generator.
	 */
	virtual void reseed_random_generator(long seed) = 0;

	/**
	 * Gets the number of random variables of the distribution.
	 * It is the size of the array returned by the {@link #sample() sample}
	 * method.
	 *
	 * @return the number of variables.
	 */
	virtual int get_dimension() = 0;

	/**
	 * Generates a random value vector sampled from this distribution.
	 *
	 * @return a random value vector.
	 */
	virtual std::vector<double> sample() = 0;

	/**
	 * Generates a list of a random value vectors from the distribution.
	 *
	 * @param sample_size the number of random vectors to generate.
	 * @return an array representing the random samples.
	 * @org.hipparchus.exception.
	 * if {@code sample_size} is not positive.
	 *
	 * @see #sample()
	 */
	virtual std::vector<std::vector<double>> sample(const int& sample_size) = 0;
};