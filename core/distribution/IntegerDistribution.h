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

  /**
   * Interface for discrete distributions.
   */
class Integer_Distribution
{
public:
	/**
	 * For a random variable {@code X} whose values are distributed according to
	 * this distribution, this method returns {@code log(P(X = x))}, where
	 * {@code log} is the natural logarithm. In other words, this method
	 * represents the logarithm of the probability mass function (PMF) for the
	 * distribution. Note that due to the floating point precision and
	 * under/overflow issues, this method will for some distributions be more
	 * precise and faster than computing the logarithm of
	 * {@link #probabilitystatic_cast<int>(}.
	 *
	 * @param x the point at which the PMF is evaluated
	 * @return the logarithm of the value of the probability mass function at {@code x}
	 */
	double log_probability(const int& x);

	/**
	 * For a random variable {@code X} whose values are distributed according
	 * to this distribution, this method returns {@code P(X = x)}. In other
	 * words, this method represents the probability mass function (PMF)
	 * for the distribution.
	 *
	 * @param x the point at which the PMF is evaluated
	 * @return the value of the probability mass function at {@code x}
	 */
	double probability(const int& x);

	/**
	 * For a random variable {@code X} whose values are distributed according
	 * to this distribution, this method returns {@code P(x0 < X <= x1)}.
	 *
	 * @param x0 the exclusive lower bound
	 * @param x1 the inclusive upper bound
	 * @return the probability that a random variable with this distribution
	 * will take a value between {@code x0} and {@code x1}, * excluding the lower and including the upper endpoint
	 * @ if {@code x0 > x1}
	 */
	double probability(const int& x0, int x1);

	/**
	 * For a random variable {@code X} whose values are distributed according
	 * to this distribution, this method returns {@code P(X <= x)}.  In other
	 * words, this method represents the (cumulative) distribution function
	 * (CDF) for this distribution.
	 *
	 * @param x the point at which the CDF is evaluated
	 * @return the probability that a random variable with this
	 * distribution takes a value less than or equal to {@code x}
	 */
	double cumulative_probability(const int& x);

	/**
	 * Computes the quantile function of this distribution.
	 * For a random variable {@code X} distributed according to this distribution, * the returned value is
	 * <ul>
	 * <li><code>inf{x in Z | P(X&lt;=x) &gt;= p}</code> for {@code 0 < p <= 1},</li>
	 * <li><code>inf{x in Z | P(X&lt;=x) &gt; 0}</code> for {@code p = 0}.</li>
	 * </ul>
	 * If the result exceeds the range of the data type {@code int}, * then {@code std::numeric_limits<int>::min()} or {@code std::numeric_limits<int>::max()} is returned.
	 *
	 * @param p the cumulative probability
	 * @return the smallest {@code p}-quantile of this distribution
	 * (largest 0-quantile for {@code p = 0})
	 * @ if {@code p < 0} or {@code p > 1}
	 */
	int inverse_cumulative_probability(const double& p);

	/**
	 * Use this method to get the numerical value of the mean of this
	 * distribution.
	 *
	 * @return the mean or {@codeNAN} if it is not defined
	 */
	double get_numerical_mean();

	/**
	 * Use this method to get the numerical value of the variance of this
	 * distribution.
	 *
	 * @return the variance (possibly {@code INFINITY} or
	 * {@codeNAN} if it is not defined)
	 */
	double get_numerical_variance();

	/**
	 * Access the lower bound of the support. This method must return the same
	 * value as {@code inverse_cumulative_probability(0)}. In other words, this
	 * method must return
	 * <p><code>inf {x in Z | P(X &lt;= x) &gt; 0}</code>.</p>
	 *
	 * @return lower bound of the support ({@code std::numeric_limits<int>::min()}
	 * for negative infinity)
	 */
	int get_support_lower_bound();

	/**
	 * Access the upper bound of the support. This method must return the same
	 * value as {@code inverse_cumulative_probability(1)}. In other words, this
	 * method must return
	 * <p><code>inf {x in R | P(X &lt;= x) = 1}</code>.</p>
	 *
	 * @return upper bound of the support ({@code std::numeric_limits<int>::max()}
	 * for positive infinity)
	 */
	int get_support_upper_bound();

	/**
	 * Use this method to get information about whether the support is
	 * connected, i.e. whether all integers between the lower and upper bound of
	 * the support are included in the support.
	 *
	 * @return whether the support is connected or not
	 */
	bool is_support_connected();
};