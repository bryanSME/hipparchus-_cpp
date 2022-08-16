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

  //import java.util.Array_list;
  //import java.util.List;

  //import org.hipparchus.distribution.continuous.Normal_Distribution;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.exception.;
  //import org.hipparchus.stat.ranking.NaN_Strategy;
  //import org.hipparchus.stat.ranking.Natural_Ranking;
  //import org.hipparchus.stat.ranking.Ties_Strategy;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
#include "../ranking/NaturalRanking.h"
#include "../../core/distribution/continuous/NormalDistribution.h"
#include <vector>
#include "../../core/util/MathArrays.h"

/**
 * An implementation of the Wilcoxon signed-rank test.
 *
 * This implementation currently handles only paired (equal length) samples
 * and discards tied pairs from the analysis. The latter behavior differs from
 * the R implementation of wilcox.test and corresponds to the "wilcox"
 * zero_method configurable in scipy.stats.wilcoxon.
 */
class Wilcoxon_Signed_Rank_Test
{
private:
	/** Ranking algorithm. */
	const Natural_Ranking my_natural_ranking;

	/**
	 * Ensures that the provided arrays fulfills the assumptions. Also computes
	 * and returns the number of tied pairs (i.e., zero differences).
	 *
	 * @param x first sample
	 * @param y second sample
	 * @return the number of indices where x[i] == y[i]
	 * @ if {@code x} or {@code y} are {@code NULL}.
	 * @ if {@code x} or {@code y} are
	 *         zero-length
	 * @ if {@code x} and {@code y} do not
	 *         have the same length.
	 * @ if all pairs are tied (i.e., if no
	 *         data remains when tied pairs have been removed.
	 */
	int ensure_data_conformance(const std::vector<double>& x, const std::vector<double>& y)
	{
		if (x.empty() || y.empty())
		{
			throw std::exception("not implemented");
			// throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
		}
		Math_Arrays::check_equal_length(y, x);
		int n_ties = 0;
		for (int i{}; i < x.size(); i++)
		{
			if (x[i] == y[i])
			{
				n_ties++;
			}
		}
		if (x.size() - n_ties == 0)
		{
			throw std::exception("not implemented");
			// throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_DATA);
		}
		return n_ties;
	}

	/**
	 * Calculates y[i] - x[i] for all i, discarding ties.
	 *
	 * @param x first sample
	 * @param y second sample
	 * @return z = y - x (minus tied values)
	 */
	std::vector<double> calculate_differences(const std::vector<double>& x, const std::vector<double>& y)
	{
		std::vector<double> differences{};
		for (int i{}; i < x.size(); ++i)
		{
			if (y[i] != x[i])
			{
				differences.push_back(y[i] - x[i]);
			}
		}
		const int n_diff = differences.size();
		auto z = std::vector<double>(n_diff);
		for (int i{}; i < n_diff; i++)
		{
			z[i] = differences.at(i);
		}
		return z;
	}

	/**
	 * Calculates |z[i]| for all i
	 *
	 * @param z sample
	 * @return |z|
	 * @ if {@code z} is {@code NULL}
	 * @ if {@code z} is zero-length.
	 */
	std::vector<double> calculate_absolute_differences(const std::vector<double>& z)
	{
		if (z.empty())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
		}

		auto z_abs = std::vector<double>(z.size());

		for (int i{}; i < z.size(); ++i)
		{
			z_abs[i] = std::abs(z[i]);
		}

		return z_abs;
	}

	/**
	 * Calculates the p-value associated with a Wilcoxon signed rank statistic
	 * by enumerating all possible rank sums and counting the number that exceed
	 * the given value.
	 *
	 * @param stat Wilcoxon signed rank statistic value
	 * @param n number of subjects (corresponding to x.size())
	 * @return two-sided exact p-value
	 */
	double calculate_exact_p_value(const double& stat, const int& n)
	{
		const int m = 1 << n;
		int larger_rank_sums{};
		for (int i{}; i < m; ++i)
		{
			int rank_sum{};

			// Generate all possible rank sums
			for (int j{}; j < n; ++j)
			{
				// (i >> j) & 1 extract i's j-th bit from the right
				if (((i >> j) & 1) == 1)
				{
					rank_sum += j + 1;
				}
			}

			if (rank_sum >= stat)
			{
				++larger_rank_sums;
			}
		}

		/*
		 * larger_rank_sums / m gives the one-sided p-value, so it's multiplied
		 * with 2 to get the two-sided p-value
		 */
		return 2 * (static_cast<double>(larger_rank_sums)) / m;
	}

	/**
	 * Computes an estimate of the (2-sided) p-value using the normal
	 * approximation. Includes a continuity correction in computing the
	 * correction factor.
	 *
	 * @param stat Wilcoxon rank sum statistic
	 * @param n number of subjects (corresponding to x.size() minus any tied ranks)
	 * @return two-sided asymptotic p-value
	 */
	double calculate_asymptotic_p_value(const double& stat, const int& n)
	{
		const double ES = n * (n + 1) / 4.0;

		/*
		 * Same as (but saves computations): const double Var_W = (static_cast<double>( (N *
		 * (N + 1) * (2*N + 1))) / 24;
		 */
		const double Var_S = ES * ((2 * n + 1) / 6.0);

		double z = stat - ES;
		const double t = FastMath.signum(z);
		z = (z - t * 0.5) / std::sqrt(Var_S);

		// want 2-sided tail probability, so make sure z < 0
		if (z > 0)
		{
			z = -z;
		}
		const Normal_Distribution standard_normal = Normal_Distribution(0, 1);
		return 2 * standard_normal.cumulative_probability(z);
	}

public:
	/**
	 * Create a test instance where NaN's are left in place and ties get the
	 * average of applicable ranks.
	 */
	Wilcoxon_Signed_Rank_Test()
		:
		my_natural_ranking{ Natural_Ranking(NaN_Strategy::FIXED, Ties_Strategy::AVERAGE) }
	{
	}

	/**
	 * Create a test instance using the given strategies for NaN's and ties.
	 *
	 * @param nan_strategy specifies the strategy that should be used for
	 *       NAN's
	 * @param ties_strategy specifies the strategy that should be used for ties
	 */
	Wilcoxon_Signed_Rank_Test(const NaN_Strategy& nan_strategy, const Ties_Strategy& ties_strategy)
		:
		my_natural_ranking{ Natural_Ranking(nan_strategy, ties_strategy) }
	{
	}

	/**
	 * Computes the
	 * <a href="http://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test">
	 * Wilcoxon signed ranked statistic</a> comparing means for two related
	 * samples or repeated measurements on a single sample.
	 * <p>
	 * This statistic can be used to perform a Wilcoxon signed ranked test
	 * evaluating the NULL hypothesis that the two related samples or repeated
	 * measurements on a single sample have equal mean.
	 * </p>
	 * <p>
	 * Let X<sub>i</sub> denote the i'th individual of the first sample and
	 * Y<sub>i</sub> the related i'th individual in the second sample. Let
	 * Z<sub>i</sub> = Y<sub>i</sub> - X<sub>i</sub>.
	 * </p>
	 * <p>
	 * <strong>Preconditions</strong>:
	 * <ul>
	 * <li>The differences Z<sub>i</sub> must be independent.</li>
	 * <li>Each Z<sub>i</sub> comes from a continuous population (they must be
	 * identical) and is symmetric about a common median.</li>
	 * <li>The values that X<sub>i</sub> and Y<sub>i</sub> represent are
	 * ordered, so the comparisons greater than, less than, and equal to are
	 * meaningful.</li>
	 * </ul>
	 * </p>
	 *
	 * @param x the first sample
	 * @param y the second sample
	 * @return wilcoxon_signed_rank statistic (the larger of W+ and W-)
	 * @ if {@code x} or {@code y} are {@code NULL}.
	 * @ if {@code x} or {@code y} are
	 *         zero-length.
	 * @ if {@code x} and {@code y} do not
	 *         have the same length.
	 */
	double wilcoxon_signed_rank(const std::vector<double>& x, const std::vector<double>& y)
	{
		ensure_data_conformance(x, y);

		auto z = calculate_differences(x, y);
		const auto z_abs = calculate_absolute_differences(z);

		const std::vector<double> ranks = my_natural_ranking.rank(z_abs);

		double Wplus{};

		for (int i{}; i < z.size(); ++i)
		{
			if (z[i] > 0)
			{
				Wplus += ranks[i];
			}
		}

		const int n{ z.size() };
		const double Wminus = ((n * (n + 1)) / 2.0) - Wplus;

		return std::max(Wplus, Wminus);
	}

	/**
	 * Returns the <i>observed significance level</i>, or
	 * <a href= "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
	 * p-value</a>, associated with a
	 * <a href="http://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test">
	 * Wilcoxon signed ranked statistic</a> comparing mean for two related
	 * samples or repeated measurements on a single sample.
	 * <p>
	 * Let X<sub>i</sub> denote the i'th individual of the first sample and
	 * Y<sub>i</sub> the related i'th individual in the second sample. Let
	 * Z<sub>i</sub> = Y<sub>i</sub> - X<sub>i</sub>.
	 * </p>
	 * <p>
	 * <strong>Preconditions</strong>:
	 * <ul>
	 * <li>The differences Z<sub>i</sub> must be independent.</li>
	 * <li>Each Z<sub>i</sub> comes from a continuous population (they must be
	 * identical) and is symmetric about a common median.</li>
	 * <li>The values that X<sub>i</sub> and Y<sub>i</sub> represent are
	 * ordered, so the comparisons greater than, less than, and equal to are
	 * meaningful.</li>
	 * </ul>
	 * <strong>Implementation notes</strong>:
	 * <ul>
	 * <li>Tied pairs are discarded from the data.</li>
	 * <li>When {@code exact_p_value} is false, the normal approximation is used
	 * to estimate the p-value including a continuity correction factor.
	 * {@code wilcoxon_signed_rank_test(x, y, true)} should give the same results
	 * as {@code  wilcox.test(x, y, alternative = "two.sided", mu = 0, *     paired = TRUE, exact = FALSE, correct = TRUE)} in R (as long as
	 * there are no tied pairs in the data).</li>
	 * </ul>
	 * </p>
	 *
	 * @param x the first sample
	 * @param y the second sample
	 * @param exact_p_value if the exact p-value is wanted (only works for
	 *        x.size() &lt;= 30, if true and x.size() &gt; 30,  is thrown)
	 * @return p-value
	 * @ if {@code x} or {@code y} are {@code NULL}.
	 * @ if {@code x} or {@code y} are
	 *         zero-length or for all i, x[i] == y[i]
	 * @ if {@code x} and {@code y} do not
	 *         have the same length.
	 * @ if {@code exact_p_value} is
	 *         {@code true} and {@code x.size()} &gt; 30
	 * @Math_Illegal_State_Exception if the p-value can not be computed due
	 *         to a convergence error
	 * @Math_Illegal_State_Exception if the maximum number of iterations is
	 *         exceeded
	 */
	double wilcoxon_signed_rank_test(const std::vector<double>& x, const std::vector<double>& y, const bool exact_p_value)
	{
		const auto n_ties = ensure_data_conformance(x, y);

		const int n = x.size() - n_ties;
		const double stat = wilcoxon_signed_rank(x, y);

		if (exact_p_value && n > 30)
		{
			throw std::exception("not implemented");
			// throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, n, 30);
		}

		return exact_p_value
			? calculate_exact_p_value(stat, n)
			: calculate_asymptotic_p_value(stat, n);
	}
};