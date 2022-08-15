//#pragma once
///*
// * Licensed to the Apache Software Foundation (ASF) under one or more
// * contributor license agreements.  See the NOTICE file distributed with
// * this work for additional information regarding copyright ownership.
// * The ASF licenses this file to You under the Apache License, Version 2.0
// * (the "License"); you may not use this file except in compliance with
// * the License.  You may obtain a copy of the License at
// *
// *      http://www.apache.org/licenses/LICENSE-2.0
// *
// * Unless required by applicable law or agreed to in writing, software
// * distributed under the License is distributed on an "AS IS" BASIS, * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// * See the License for the specific language governing permissions and
// * limitations under the License.
// */
//
///*
// * This is not the original file distributed by the Apache Software Foundation
// * It has been modified by the Hipparchus project
// */
////package org.hipparchus.stat.inference;
//
////import java.util.Map;
////import java.util.Tree_Map;
////import java.util.stream.Long_Stream;
//
////import org.hipparchus.distribution.continuous.Normal_Distribution;
////import org.hipparchus.exception.Localized_Core_Formats;
////import org.hipparchus.exception.;
////import org.hipparchus.exception.Math_Illegal_State_Exception;
////import org.hipparchus.exception.Null_Argument_Exception;
////import org.hipparchus.stat.Localized_Stat_Formats;
////import org.hipparchus.stat.ranking.NaN_Strategy;
////import org.hipparchus.stat.ranking.Natural_Ranking;
////import org.hipparchus.stat.ranking.Ties_Strategy;
////import org.hipparchus.util.FastMath;
////import org.hipparchus.util.Precision;
//#include "../ranking/NaturalRanking.h"
//#include "../../core/distribution/continuous/NormalDistribution.h"
//
///**
// * An implementation of the Mann-Whitney U test.
// * <p>
// * The definitions and computing formulas used in this implementation follow
// * those in the article, * <a href="http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U"> Mann-Whitney U
// * Test</a>
// * <p>
// * In general, results correspond to (and have been tested against) the R
// * wilcox.test function, with {@code exact} meaning the same thing in both APIs
// * and {@code CORRECT} uniformly true in this implementation. For example, * wilcox.test(x, y, alternative = "two.sided", mu = 0, paired = FALSE, exact = FALSE
// * correct = TRUE) will return the same p-value as mann_whitney_u_test(x, y, * false). The minimum of the W value returned by R for wilcox.test(x, y...) and
// * wilcox.test(y, x...) should equal mann_whitney_u(x, y...).
// */
//class Mann_Whitney_U_Test 
//{
//private:
//    /**
//     * If the combined dataset contains no more values than this, test defaults to
//     * exact test.
//     */
//    static const int SMALL_SAMPLE_SIZE{ 50 };
//
//    /** Ranking algorithm. */
//    const Natural_Ranking my_natural_ranking;
//
//    /** Normal distribution */
//    const Normal_Distribution my_standard_normal;
//
//    /**
//     * Create a test instance using where NaN's are left in place and ties get
//     * the average of applicable ranks.
//     */
//    public Mann_Whitney_U_Test() 
//    {
//        natural_ranking = Natural_Ranking(NaN_Strategy.FIXED, Ties_Strategy.AVERAGE);
//        standard_normal = Normal_Distribution(0, 1);
//    }
//
//    /**
//     * Create a test instance using the given strategies for NaN's and ties.
//     *
//     * @param nan_strategy specifies the strategy that should be used for
//     *       NAN's
//     * @param ties_strategy specifies the strategy that should be used for ties
//     */
//    public Mann_Whitney_U_Test(const NaN_Strategy nan_strategy, const Ties_Strategy ties_strategy) 
//    {
//        natural_ranking = Natural_Ranking(nan_strategy, ties_strategy);
//        standard_normal = Normal_Distribution(0, 1);
//    }
//
//    /**
//     * Computes the
//     * <a href="http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U">
//     * Mann-Whitney U statistic</a> comparing means for two independent samples
//     * possibly of different lengths.
//     * <p>
//     * This statistic can be used to perform a Mann-Whitney U test evaluating
//     * the NULL hypothesis that the two independent samples have equal mean.
//     * <p>
//     * Let X<sub>i</sub> denote the i'th individual of the first sample and
//     * Y<sub>j</sub> the j'th individual in the second sample. Note that the
//     * samples can have different lengths.
//     * <p>
//     * <strong>Preconditions</strong>:
//     * <ul>
//     * <li>All observations in the two samples are independent.</li>
//     * <li>The observations are at least ordinal (continuous are also
//     * ordinal).</li>
//     * </ul>
//     *
//     * @param x the first sample
//     * @param y the second sample
//     * @return Mann-Whitney U statistic (minimum of U<sup>x</sup> and
//     *         U<sup>y</sup>)
//     * @Null_Argument_Exception if {@code x} or {@code y} are {@code NULL}.
//     * @ if {@code x} or {@code y} are
//     *         zero-length.
//     */
//    public double mann_whitney_u(const std::vector<double> x, const std::vector<double> y)
//        , Null_Argument_Exception 
//        {
//
//        ensure_data_conformance(x, y);
//
//        const std::vector<double> z = concatenate_samples(x, y);
//        const std::vector<double> ranks = natural_ranking.rank(z);
//
//        double sum_rank_x = 0;
//
//        /*
//         * The ranks for x is in the first x.size() entries in ranks because x
//         * is in the first x.size() entries in z
//         */
//        for (int i{}; i < x.size(); ++i) 
//        {
//            sum_rank_x += ranks[i];
//        }
//
//        /*
//         * U1 = R1 - (n1 * (n1 + 1)) / 2 where R1 is sum of ranks for sample 1, * e.g. x, n1 is the number of observations in sample 1.
//         */
//        const double U1 = sum_rank_x - (static_cast<long>( x.size() * (x.size() + 1)) / 2;
//
//        /*
//         * U1 + U2 = n1 * n2
//         */
//        const double U2 = static_cast<long>( x.size() * y.size() - U1;
//
//        return std::min(U1, U2);
//    }
//
//    /**
//     * Concatenate the samples into one array.
//     *
//     * @param x first sample
//     * @param y second sample
//     * @return concatenated array
//     */
//    private std::vector<double> concatenate_samples(const std::vector<double> x, const std::vector<double> y) 
//    {
//        const std::vector<double> z = std::vector<double>(x.size() + y.size()];
//
//        System.arraycopy(x, 0, z, 0, x.size());
//        System.arraycopy(y, 0, z, x.size(), y.size());
//
//        return z;
//    }
//
//    /**
//     * Returns the asymptotic <i>observed significance level</i>, or
//     * <a href="http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
//     * p-value</a>, associated with a <a href=
//     * "http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U">Mann-Whitney U
//     * Test</a> comparing means for two independent samples.
//     * <p>
//     * Let X<sub>i</sub> denote the i'th individual of the first sample and
//     * Y<sub>j</sub> the j'th individual in the second sample.
//     * <p>
//     * <strong>Preconditions</strong>:
//     * <ul>
//     * <li>All observations in the two samples are independent.</li>
//     * <li>The observations are at least ordinal.</li>
//     * </ul>
//     * <p>
//     * If there are no ties in the data and both samples are small (less than or
//     * equal to 50 values in the combined dataset), an exact test is performed;
//     * otherwise the test uses the normal approximation (with continuity
//     * correction).
//     * <p>
//     * If the combined dataset contains ties, the variance used in the normal
//     * approximation is bias-adjusted using the formula in the reference above.
//     *
//     * @param x the first sample
//     * @param y the second sample
//     * @return approximate 2-sized p-value
//     * @Null_Argument_Exception if {@code x} or {@code y} are {@code NULL}.
//     * @ if {@code x} or {@code y} are
//     *         zero-length
//     */
//    public double mann_whitney_u_test(const std::vector<double> x, const std::vector<double> y)
//        , Null_Argument_Exception 
//        {
//        ensure_data_conformance(x, y);
//
//        // If samples are both small and there are no ties, perform exact test
//        if (x.size() + y.size() <= SMALL_SAMPLE_SIZE &&
//            ties_map(x, y).is_empty()) 
//            {
//            return mann_whitney_u_test(x, y, true);
//        }
//else { // Normal approximation
//            return mann_whitney_u_test(x, y, false);
//        }
//    }
//
//    /**
//     * Returns the asymptotic <i>observed significance level</i>, or
//     * <a href="http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
//     * p-value</a>, associated with a <a href=
//     * "http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U">Mann-Whitney U
//     * Test</a> comparing means for two independent samples.
//     * <p>
//     * Let X<sub>i</sub> denote the i'th individual of the first sample and
//     * Y<sub>j</sub> the j'th individual in the second sample.
//     * <p>
//     * <strong>Preconditions</strong>:
//     * <ul>
//     * <li>All observations in the two samples are independent.</li>
//     * <li>The observations are at least ordinal.</li>
//     * </ul>
//     * <p>
//     * If {@code exact} is {@code true}, the p-value reported is exact, computed
//     * using the exact distribution of the U statistic. The computation in this
//     * case requires storage on the order of the product of the two sample
//     * sizes, so this should not be used for large samples.
//     * <p>
//     * If {@code exact} is {@code false}, the normal approximation is used to
//     * estimate the p-value.
//     * <p>
//     * If the combined dataset contains ties and {@code exact} is {@code true}, *  is thrown. If {@code exact} is {@code false}
//     * and the ties are present, the variance used to compute the approximate
//     * p-value in the normal approximation is bias-adjusted using the formula in
//     * the reference above.
//     *
//     * @param x the first sample
//     * @param y the second sample
//     * @param exact true means compute the p-value exactly, false means use the
//     *        normal approximation
//     * @return approximate 2-sided p-value
//     * @Null_Argument_Exception if {@code x} or {@code y} are {@code NULL}.
//     * @ if {@code x} or {@code y} are
//     *         zero-length or if {@code exact} is {@code true} and ties are
//     *         present in the data
//     */
//    public double mann_whitney_u_test(const std::vector<double> x, const std::vector<double> y, const bool exact)
//        , Null_Argument_Exception 
//        {
//        ensure_data_conformance(x, y);
//        const Map<Double, Integer> ties_map = ties_map(x, y);
//        const double u = mann_whitney_u(x, y);
//        if (exact) 
//        {
//            if (!ties_map.is_empty()) 
//            {
//                throw (Localized_Stat_Formats.TIES_ARE_NOT_ALLOWED);
//            }
//            return exact_p(x.size(), y.size(), u);
//        }
//
//        return approximate_p(u, x.size(), y.size(), var_u(x.size(), y.size(), ties_map));
//    }
//
//    /**
//     * Ensures that the provided arrays fulfills the assumptions.
//     *
//     * @param x first sample
//     * @param y second sample
//     * @Null_Argument_Exception if {@code x} or {@code y} are {@code NULL}.
//     * @ if {@code x} or {@code y} are
//     *         zero-length.
//     */
//    private void ensure_data_conformance(const std::vector<double> x, const std::vector<double> y)
//        , Null_Argument_Exception 
//        {
//
//        if (x == NULL || y == NULL) 
//        {
//            throw Null_Argument_Exception();
//        }
//        if (x.size() == 0 || y.size() == 0) 
//        {
//            throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
//        }
//    }
//
//    /**
//     * Estimates the 2-sided p-value associated with a Mann-Whitney U statistic
//     * value using the normal approximation.
//     * <p>
//     * The variance passed in is assumed to be corrected for ties. Continuity
//     * correction is applied to the normal approximation.
//     *
//     * @param u Mann-Whitney U statistic
//     * @param n1 number of subjects in first sample
//     * @param n2 number of subjects in second sample
//     * @param var_u variance of U (corrected for ties if these exist)
//     * @return two-sided asymptotic p-value
//     * @Math_Illegal_State_Exception if the p-value can not be computed due
//     *         to a convergence error
//     * @Math_Illegal_State_Exception if the maximum number of iterations is
//     *         exceeded
//     */
//    private double approximate_p(const double u, const int& n1, const int& n2, const double var_u)
//        Math_Illegal_State_Exception 
//        {
//
//        const double mu = static_cast<long>( n1 * n2 / 2.0;
//
//        // If u == mu, return 1
//        if (Precision.equals(mu, u)) 
//        {
//            return 1;
//        }
//
//        // Force z <= 0 so we get tail probability. Also apply continuity
//        // correction
//        const double z = -Math.abs((u - mu) + 0.5) / std::sqrt(var_u);
//
//        return 2 * standard_normal.cumulative_probability(z);
//    }
//
//    /**
//     * Calculates the (2-sided) p-value associated with a Mann-Whitney U
//     * statistic.
//     * <p>
//     * To compute the p-value, the probability densities for each value of U up
//     * to and including u are summed and the resulting tail probability is
//     * multiplied by 2.
//     * <p>
//     * The result of this computation is only valid when the combined n + m
//     * sample has no tied values.
//     * <p>
//     * This method should not be used for large values of n or m as it maintains
//     * work arrays of size n*m.
//     *
//     * @param u Mann-Whitney U statistic value
//     * @param n first sample size
//     * @param m second sample size
//     * @return two-sided exact p-value
//     */
//    private double exact_p(const int& n, const int m, const double u) 
//    {
//        const double nm = m * n;
//        if (u > nm) { // Quick exit if u is out of range
//            return 1;
//        }
//        // Need to convert u to a mean deviation, so cumulative probability is
//        // tail probability
//        const double crit = u < nm / 2 ? u : nm / 2 - u;
//
//        double cum = 0;
//        for (const int& ct = 0; ct <= crit; ct++) 
//        {
//            cum += u_density(n, m, ct);
//        }
//        return 2 * cum;
//    }
//
//    /**
//     * Computes the probability density function for the Mann-Whitney U
//     * statistic.
//     * <p>
//     * This method should not be used for large values of n or m as it maintains
//     * work arrays of size n*m.
//     *
//     * @param n first sample size
//     * @param m second sample size
//     * @param u U-statistic value
//     * @return the probability that a U statistic derived from random samples of
//     *         size n and m (containing no ties) equals u
//     */
//    private double u_density(const int& n, const int m, double u) 
//    {
//        if (u < 0 || u > m * n) 
//        {
//            return 0;
//        }
//        const std::vector<long> freq = u_frequencies(n, m);
//        return freq[static_cast<int>( std::round(u + 1)] /
//               static_cast<double>( Long_Stream.of(freq).sum();
//    }
//
//    /**
//     * Computes frequency counts for values of the Mann-Whitney U statistc. If
//     * freq[] is the returned array, freq[u + 1] counts the frequency of U = u
//     * among all possible n-m orderings. Therefore, P(u = U) = freq[u + 1] / sum
//     * where sum is the sum of the values in the returned array.
//     * <p>
//     * Implements the algorithm presented in "Algorithm AS 62: A Generator for
//     * the Sampling Distribution of the Mann-Whitney U Statistic", L. C. Dinneen
//     * and B. C. Blakesley Journal of the Royal Statistical Society. Series C
//     * (Applied Statistics) Vol. 22, No. 2 (1973), pp. 269-273.
//     *
//     * @param n first sample size
//     * @param m second sample size
//     * @return array of U statistic value frequencies
//     */
//    private std::vector<long> u_frequencies(const int& n, const int m) 
//    {
//        const int max = std::max(m, n);
//        if (max > 100) 
//        {
//            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, max, 100);
//        }
//        const int min = std::min(m, n);
//        const std::vector<long> out = long[n * m + 2];
//        const std::vector<long> work = long[n * m + 2];
//        for (int i{ 1 }; i < out.size(); i++) 
//        {
//            out[i] = (i <= (max + 1)) ? 1 : 0;
//        }
//        work[1] = 0;
//        int in = max;
//        for (int i{ 2 }; i <= min; i++) 
//        {
//            work[i] = 0;
//            in = in + max;
//            int n1 = in + 2;
//            long l = 1 + in / 2;
//            int k = i;
//            for (int j{ 1 }; j <= l; j++) 
//            {
//                k++;
//                n1 = n1 - 1;
//                const long sum = out[j] + work[j];
//                out[j] = sum;
//                work[k] = sum - out[n1];
//                out[n1] = sum;
//            }
//        }
//        return out;
//    }
//
//    /**
//     * Computes the variance for a U-statistic associated with samples of
//     * sizes{@code n} and {@code m} and ties described by {@code ties_map}. If
//     * {@code ties_map} is non-empty, the multiplicity counts in its values set
//     * are used to adjust the variance.
//     *
//     * @param n first sample size
//     * @param m second sample size
//     * @param ties_map map of <value, multiplicity>
//     * @return ties-adjusted variance
//     */
//    private double var_u(const int& n, const int m, Map<Double, Integer> ties_map) 
//    {
//        const double nm = static_cast<long>( n * m;
//        if (ties_map.is_empty()) 
//        {
//            return nm * (n + m + 1) / 12.0;
//        }
//        const long t_sum = ties_map.entry_set().stream()
//            .map_to_long(e -> e.get_value() * e.get_value() * e.get_value() -
//                            e.get_value())
//            .sum();
//        const double total_n = n + m;
//        return (nm / 12) * (total_n + 1 - t_sum / (total_n * (total_n - 1)));
//
//    }
//
//    /**
//     * Creates a map whose keys are values occurring more than once in the
//     * combined dataset formed from x and y. Map entry values are the number of
//     * occurrences. The returned map is empty iff there are no ties in the data.
//     *
//     * @param x first dataset
//     * @param y second dataset
//     * @return map of <value, number of times it occurs> for values occurring
//     *         more than once or an empty map if there are no ties (the returned
//     *         map is <em>not</em> thread-safe, which is OK in the context of the callers)
//     */
//    private Map<Double, Integer> ties_map(const std::vector<double> x, const std::vector<double> y) 
//    {
//        const Map<Double, Integer> ties_map = Tree_Map<>(); // NOPMD - no concurrent access in the callers context
//        for (int i{}; i < x.size(); i++) 
//        {
//            ties_map.merge(x[i], 1, Integer::sum);
//        }
//        for (int i{}; i < y.size(); i++) 
//        {
//            ties_map.merge(y[i], 1, Integer::sum);
//        }
//        ties_map.entry_set().remove_if(e -> e.get_value() == 1);
//        return ties_map;
//    }
//}
//
//
