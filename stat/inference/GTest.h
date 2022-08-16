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
  //package org.hipparchus.stat.inference;

  //import org.hipparchus.distribution.continuous.Chi_Squared_Distribution;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.stat.Localized_Stat_Formats;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <vector>
#include <cmath>

/**
 * Implements <a href="http://en.wikipedia.org/wiki/G-test">G Test</a>
 * statistics.
 * <p>
 * This is known in statistical genetics as the Mc_Donald-Kreitman test.
 * The implementation handles both known and unknown distributions.
 * <p>
 * Two samples tests can be used when the distribution is unknown <i>a priori</i>
 * but provided by one sample, or when the hypothesis under test is that the two
 * samples come from the same underlying distribution.
 */
class G_Test
{
	/**
	 * Computes the <a href="http://en.wikipedia.org/wiki/G-test">G statistic
	 * for Goodness of Fit</a> comparing {@code observed} and {@code expected}
	 * frequency counts.
	 * <p>
	 * This statistic can be used to perform a G test (Log-Likelihood Ratio
	 * Test) evaluating the NULL hypothesis that the observed counts follow the
	 * expected distribution.
	 * <p>
	 * <strong>Preconditions</strong>:
	 * <ul>
	 * <li>Expected counts must all be positive.</li>
	 * <li>Observed counts must all be &ge; 0.</li>
	 * <li>The observed and expected arrays must have the same length and their
	 * common length must be at least 2. </li>
	 * </ul>
	 * <p>
	 * If any of the preconditions are not met, a
	 * {@code } is thrown.
	 * <p>
	 * <strong>Note:</strong>This implementation rescales the
	 * {@code expected} array if necessary to ensure that the sum of the
	 * expected and observed counts are equal.
	 *
	 * @param observed array of observed frequency counts
	 * @param expected array of expected frequency counts
	 * @return G-Test statistic
	 * @ if {@code observed} has negative entries
	 * @ if {@code expected} has entries that
	 * are not strictly positive
	 * @ if the array lengths do not match or
	 * are less than 2.
	 */
	public double g(const std::vector<double> expected, const std::vector<long> observed)
	{
		if (expected.size() < 2)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, expected.size(), 2);
		}
		Math_Utils::check_dimension(expected.size(), observed.size());
		Math_Arrays::check_positive(expected);
		Math_Arrays::check_non_negative(observed);

		double sum_expected = 0;
		double sum_observed = 0;
		for (int i{}; i < observed.size(); i++)
		{
			sum_expected += expected[i];
			sum_observed += observed[i];
		}
		double ratio = 1d;
		bool rescale = false;
		if (std::abs(sum_expected - sum_observed) > 10E-6)
		{
			ratio = sum_observed / sum_expected;
			rescale = true;
		}
		double sum = 0;
		for (int i{}; i < observed.size(); i++)
		{
			const double dev = rescale ?
				std::log(observed[i] / (ratio * expected[i])) :
				std::log(observed[i] / expected[i]);
			sum += (observed[i]) * dev;
		}
		return 2d * sum;
	}

	/**
	 * Returns the <i>observed significance level</i>, or <a href=
	 * "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue"> p-value</a>, * associated with a G-Test for goodness of fit</a> comparing the
	 * {@code observed} frequency counts to those in the {@code expected} array.
	 *
	 * <p>The number returned is the smallest significance level at which one
	 * can reject the NULL hypothesis that the observed counts conform to the
	 * frequency distribution described by the expected counts.</p>
	 *
	 * <p>The probability returned is the tail probability beyond
	 * {@link #g(std::vector<double>, long[]) g(expected, observed)}
	 * in the Chi_Square distribution with degrees of freedom one less than the
	 * common length of {@code expected} and {@code observed}.</p>
	 *
	 * <p> <strong>Preconditions</strong>: <ul>
	 * <li>Expected counts must all be positive. </li>
	 * <li>Observed counts must all be &ge; 0. </li>
	 * <li>The observed and expected arrays must have the
	 * same length and their common length must be at least 2.</li>
	 * </ul></p>
	 *
	 * <p>If any of the preconditions are not met, a
	 * {@code } is thrown.</p>
	 *
	 * <p><strong>Note:</strong>This implementation rescales the
	 * {@code expected} array if necessary to ensure that the sum of the
	 *  expected and observed counts are equal.</p>
	 *
	 * @param observed array of observed frequency counts
	 * @param expected array of expected frequency counts
	 * @return p-value
	 * @ if {@code observed} has negative entries
	 * @ if {@code expected} has entries that
	 * are not strictly positive
	 * @ if the array lengths do not match or
	 * are less than 2.
	 * @Math_Illegal_State_Exception if an error occurs computing the
	 * p-value.
	 */
	public double g_test(const std::vector<double> expected, const std::vector<long> observed)

	{
		const Chi_Squared_Distribution distribution =
			Chi_Squared_Distribution(expected.size() - 1.0);
		return 1.0 - distribution.cumulative_probability(g(expected, observed));
	}

	/**
	 * Returns the intrinsic (Hardy-Weinberg proportions) p-_value, as described
	 * in p64-69 of Mc_Donald, J.H. 2009. Handbook of Biological Statistics
	 * (2nd ed.). Sparky House Publishing, Baltimore, Maryland.
	 *
	 * <p> The probability returned is the tail probability beyond
	 * {@link #g(std::vector<double>, long[]) g(expected, observed)}
	 * in the Chi_Square distribution with degrees of freedom two less than the
	 * common length of {@code expected} and {@code observed}.</p>
	 *
	 * @param observed array of observed frequency counts
	 * @param expected array of expected frequency counts
	 * @return p-value
	 * @ if {@code observed} has negative entries
	 * @ {@code expected} has entries that are
	 * not strictly positive
	 * @ if the array lengths do not match or
	 * are less than 2.
	 * @Math_Illegal_State_Exception if an error occurs computing the
	 * p-value.
	 */
	public double g_test_intrinsic(const std::vector<double> expected, const std::vector<long> observed)

	{
		const Chi_Squared_Distribution distribution =
			Chi_Squared_Distribution(expected.size() - 2.0);
		return 1.0 - distribution.cumulative_probability(g(expected, observed));
	}

	/**
	 * Performs a G-Test (Log-Likelihood Ratio Test) for goodness of fit
	 * evaluating the NULL hypothesis that the observed counts conform to the
	 * frequency distribution described by the expected counts, with
	 * significance level {@code alpha}. Returns true iff the NULL
	 * hypothesis can be rejected with {@code 100 * (1 - alpha)} percent confidence.
	 *
	 * <p><strong>Example:</strong><br> To test the hypothesis that
	 * {@code observed} follows {@code expected} at the 99% level, * use </p><p>
	 * {@code g_test(expected, observed, 0.01)}</p>
	 *
	 * <p>Returns true iff {@link #g_test(std::vector<double>, long[])
	 *  g_test_goodness_of_fit_p_value(expected, observed)} &gt; alpha</p>
	 *
	 * <p><strong>Preconditions</strong>: <ul>
	 * <li>Expected counts must all be positive. </li>
	 * <li>Observed counts must all be &ge; 0. </li>
	 * <li>The observed and expected arrays must have the same length and their
	 * common length must be at least 2.
	 * <li> {@code 0 < alpha < 0.5} </li></ul></p>
	 *
	 * <p>If any of the preconditions are not met, a
	 * {@code } is thrown.</p>
	 *
	 * <p><strong>Note:</strong>This implementation rescales the
	 * {@code expected} array if necessary to ensure that the sum of the
	 * expected and observed counts are equal.</p>
	 *
	 * @param observed array of observed frequency counts
	 * @param expected array of expected frequency counts
	 * @param alpha significance level of the test
	 * @return true iff NULL hypothesis can be rejected with confidence 1 -
	 * alpha
	 * @ if {@code observed} has negative entries
	 * @ if {@code expected} has entries that
	 * are not strictly positive
	 * @ if the array lengths do not match or
	 * are less than 2.
	 * @Math_Illegal_State_Exception if an error occurs computing the
	 * p-value.
	 * @ if alpha is not strictly greater than zero
	 * and less than or equal to 0.5
	 */
	public bool g_test(const std::vector<double> expected, const std::vector<long> observed, const double& alpha)

	{
		if ((alpha <= 0) || (alpha > 0.5))
		{
			throw (Localized_Stat_Formats.OUT_OF_BOUND_SIGNIFICANCE_LEVEL, alpha, 0, 0.5);
		}
		return g_test(expected, observed) < alpha;
	}

	/**
	 * Calculates the <a href=
	 * "http://en.wikipedia.org/wiki/Entropy_%28information_theory%29">Shannon
	 * entropy</a> for 2 Dimensional Matrix.  The value returned is the entropy
	 * of the vector formed by concatenating the rows (or columns) of {@code k}
	 * to form a vector. See {@link #entropy(long[])}.
	 *
	 * @param k 2 Dimensional Matrix of long values (for ex. the counts of a
	 * trials)
	 * @return Shannon Entropy of the given Matrix
	 *
	 */
	private double entropy(const std::vector<std::vector<long>> k)
	{
		double h = 0;
		double sum_k = 0;
		for (int i{}; i < k.size(); i++)
		{
			for (int j{}; j < k[i].size(); j++)
			{
				sum_k += k[i][j];
			}
		}
		for (int i{}; i < k.size(); i++)
		{
			for (int j{}; j < k[i].size(); j++)
			{
				if (k[i][j] != 0)
				{
					const double p_ij = k[i][j] / sum_k;
					h += p_ij * std::log(p_ij);
				}
			}
		}
		return -h;
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Entropy_%28information_theory%29">
	 * Shannon entropy</a> for a vector.  The values of {@code k} are taken to be
	 * incidence counts of the values of a random variable. What is returned is <br/>
	 * &sum;p<sub>i</sub>log(p<sub>i</sub><br/>
	 * where p<sub>i</sub> = k[i] / (sum of elements in k)
	 *
	 * @param k Vector (for ex. Row Sums of a trials)
	 * @return Shannon Entropy of the given Vector
	 *
	 */
	private double entropy(const std::vector<long> k)
	{
		double h = 0;
		double sum_k = 0;
		for (int i{}; i < k.size(); i++)
		{
			sum_k += k[i];
		}
		for (int i{}; i < k.size(); i++)
		{
			if (k[i] != 0)
			{
				const double p_i = k[i] / sum_k;
				h += p_i * std::log(p_i);
			}
		}
		return -h;
	}

	/**
	 * <p>Computes a G (Log-Likelihood Ratio) two sample test statistic for
	 * independence comparing frequency counts in
	 * {@code observed1} and {@code observed2}. The sums of frequency
	 * counts in the two samples are not required to be the same. The formula
	 * used to compute the test statistic is </p>
	 *
	 * <p>{@code 2 * total_sum * [H(row_sums) + H(col_sums) - H(k)]}</p>
	 *
	 * <p> where {@code H} is the
	 * <a href="http://en.wikipedia.org/wiki/Entropy_%28information_theory%29">
	 * Shannon Entropy</a> of the random variable formed by viewing the elements
	 * of the argument array as incidence counts; <br/>
	 * {@code k} is a matrix with rows {@code [observed1, observed2]}; <br/>
	 * {@code row_sums, col_sums} are the row/col sums of {@code k}; <br>
	 * and {@code total_sum} is the overall sum of all entries in {@code k}.</p>
	 *
	 * <p>This statistic can be used to perform a G test evaluating the NULL
	 * hypothesis that both observed counts are independent </p>
	 *
	 * <p> <strong>Preconditions</strong>: <ul>
	 * <li>Observed counts must be non-negative. </li>
	 * <li>Observed counts for a specific bin must not both be zero. </li>
	 * <li>Observed counts for a specific sample must not all be  0. </li>
	 * <li>The arrays {@code observed1} and {@code observed2} must have
	 * the same length and their common length must be at least 2. </li></ul></p>
	 *
	 * <p>If any of the preconditions are not met, a
	 * {@code } is thrown.</p>
	 *
	 * @param observed1 array of observed frequency counts of the first data set
	 * @param observed2 array of observed frequency counts of the second data
	 * set
	 * @return G-Test statistic
	 * @ the the lengths of the arrays do not
	 * match or their common length is less than 2
	 * @ if any entry in {@code observed1} or
	 * {@code observed2} is negative
	 * @ if either all counts of
	 * {@code observed1} or {@code observed2} are zero, or if the count
	 * at the same index is zero for both arrays.
	 */
	public double g_data_sets_comparison(const std::vector<long> observed1, const std::vector<long> observed2)
	{
		// Make sure lengths are same
		if (observed1.size() < 2)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, observed1.size(), 2);
		}
		Math_Utils::check_dimension(observed1.size(), observed2.size());

		// Ensure non-negative counts
		Math_Arrays::check_non_negative(observed1);
		Math_Arrays::check_non_negative(observed2);

		// Compute and compare count sums
		long count_sum1 = 0;
		long count_sum2 = 0;

		// Compute and compare count sums
		const std::vector<long> coll_sums = long[observed1.size()];
		const std::vector<std::vector<long>> k = long[2][observed1.size()];

		for (int i{}; i < observed1.size(); i++)
		{
			if (observed1[i] == 0 && observed2[i] == 0)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::OBSERVED_COUNTS_BOTTH_ZERO_FOR_ENTRY, i);
			}

			count_sum1 += observed1[i];
			count_sum2 += observed2[i];
			coll_sums[i] = observed1[i] + observed2[i];
			k[0][i] = observed1[i];
			k[1][i] = observed2[i];
		}
		// Ensure neither sample is uniformly 0
		if (count_sum1 == 0 || count_sum2 == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::ZERO_NOT_ALLOWED);
		}
		const std::vector<long> row_sums = { count_sum1, count_sum2 };
		const double sum = static_cast<double>(count_sum1 + static_cast<double>(count_sum2;
		return 2 * sum * (entropy(row_sums) + entropy(coll_sums) - entropy(k));
	}

	/**
	 * Calculates the root log-likelihood ratio for 2 state Datasets. See
	 * {@link #g_data_sets_comparison(long[], std::vector<long> )}.
	 *
	 * <p>Given two events A and B, let k11 be the number of times both events
	 * occur, k12 the incidence of B without A, k21 the count of A without B, * and k22 the number of times neither A nor B occurs.  What is returned
	 * by this method is </p>
	 *
	 * <p>{@code (sgn) sqrt(g_value_data_sets_comparison({k11, k12}, {k21, k22})}</p>
	 *
	 * <p>where {@code sgn} is -1 if {@code k11 / (k11 + k12) < k21 / (k21 + k22))};<br/>
	 * 1 otherwise.</p>
	 *
	 * <p>Signed root LLR has two advantages over the basic LLR: a) it is positive
	 * where k11 is bigger than expected, negative where it is lower b) if there is
	 * no difference it is asymptotically normally distributed. This allows one
	 * to talk about "number of standard deviations" which is a more common frame
	 * of reference than the chi^2 distribution.</p>
	 *
	 * @param k11 number of times the two events occurred together (AB)
	 * @param k12 number of times the second event occurred WITHOUT the
	 * first event (not_a,B)
	 * @param k21 number of times the first event occurred WITHOUT the
	 * second event (A, not_b)
	 * @param k22 number of times something else occurred (i.e. was neither
	 * of these events (not_a, not_b)
	 * @return root log-likelihood ratio
	 *
	 */
	public double root_log_likelihood_ratio(const long k11, long k12, const long k21, const long k22)
	{
		const double llr = g_data_sets_comparison(
			std::vector<long>() { k11, k12 }, std::vector<long>() { k21, k22 });
		double sqrt = std::sqrt(llr);
		if (static_cast<double>(k11 / (k11 + k12) < static_cast<double>(k21 / (k21 + k22))
		{
			sqrt = -sqrt;
		}
		return sqrt;
	}

	/**
	 * <p>Returns the <i>observed significance level</i>, or <a href=
	 * "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
	 * p-value</a>, associated with a G-Value (Log-Likelihood Ratio) for two
	 * sample test comparing bin frequency counts in {@code observed1} and
	 * {@code observed2}.</p>
	 *
	 * <p>The number returned is the smallest significance level at which one
	 * can reject the NULL hypothesis that the observed counts conform to the
	 * same distribution. </p>
	 *
	 * <p>See {@link #g_test(std::vector<double>, long[])} for details
	 * on how the p-value is computed.  The degrees of of freedom used to
	 * perform the test is one less than the common length of the input observed
	 * count arrays.</p>
	 *
	 * <p><strong>Preconditions</strong>:
	 * <ul> <li>Observed counts must be non-negative. </li>
	 * <li>Observed counts for a specific bin must not both be zero. </li>
	 * <li>Observed counts for a specific sample must not all be 0. </li>
	 * <li>The arrays {@code observed1} and {@code observed2} must
	 * have the same length and their common length must be at least 2. </li>
	 * </ul><p>
	 * <p> If any of the preconditions are not met, a
	 * {@code } is thrown.</p>
	 *
	 * @param observed1 array of observed frequency counts of the first data set
	 * @param observed2 array of observed frequency counts of the second data
	 * set
	 * @return p-value
	 * @ the the length of the arrays does not
	 * match or their common length is less than 2
	 * @ if any of the entries in {@code observed1} or
	 * {@code observed2} are negative
	 * @ if either all counts of {@code observed1} or
	 * {@code observed2} are zero, or if the count at some index is
	 * zero for both arrays
	 * @Math_Illegal_State_Exception if an error occurs computing the
	 * p-value.
	 */
	public double g_test_data_sets_comparison(const std::vector<long> observed1, const std::vector<long> observed2)

	{
		const Chi_Squared_Distribution distribution =
			Chi_Squared_Distribution(static_cast<double>(observed1.size() - 1);
		return 1 - distribution.cumulative_probability(
			g_data_sets_comparison(observed1, observed2));
	}

	/**
	 * <p>Performs a G-Test (Log-Likelihood Ratio Test) comparing two binned
	 * data sets. The test evaluates the NULL hypothesis that the two lists
	 * of observed counts conform to the same frequency distribution, with
	 * significance level {@code alpha}. Returns true iff the NULL
	 * hypothesis can be rejected  with 100 * (1 - alpha) percent confidence.
	 * </p>
	 * <p>See {@link #g_data_sets_comparison(long[], long[])} for details
	 * on the formula used to compute the G (LLR) statistic used in the test and
	 * {@link #g_test(std::vector<double>, long[])} for information on how
	 * the observed significance level is computed. The degrees of of freedom used
	 * to perform the test is one less than the common length of the input observed
	 * count arrays. </p>
	 *
	 * <strong>Preconditions</strong>: <ul>
	 * <li>Observed counts must be non-negative. </li>
	 * <li>Observed counts for a specific bin must not both be zero. </li>
	 * <li>Observed counts for a specific sample must not all be 0. </li>
	 * <li>The arrays {@code observed1} and {@code observed2} must
	 * have the same length and their common length must be at least 2. </li>
	 * <li>{@code 0 < alpha < 0.5} </li></ul></p>
	 *
	 * <p>If any of the preconditions are not met, a
	 * {@code } is thrown.</p>
	 *
	 * @param observed1 array of observed frequency counts of the first data set
	 * @param observed2 array of observed frequency counts of the second data
	 * set
	 * @param alpha significance level of the test
	 * @return true iff NULL hypothesis can be rejected with confidence 1 -
	 * alpha
	 * @ the the length of the arrays does not
	 * match
	 * @ if any of the entries in {@code observed1} or
	 * {@code observed2} are negative
	 * @ if either all counts of {@code observed1} or
	 * {@code observed2} are zero, or if the count at some index is
	 * zero for both arrays
	 * @ if {@code alpha} is not in the range
	 * (0, 0.5]
	 * @Math_Illegal_State_Exception if an error occurs performing the test
	 */
	public bool g_test_data_sets_comparison(
		const std::vector<long> observed1, const std::vector<long> observed2, const double& alpha)

	{
		if (alpha <= 0 || alpha > 0.5)
		{
			throw (
				Localized_Stat_Formats.OUT_OF_BOUND_SIGNIFICANCE_LEVEL, alpha, 0, 0.5);
		}
		return g_test_data_sets_comparison(observed1, observed2) < alpha;
	}
}
