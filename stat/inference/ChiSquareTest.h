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
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.stat.Localized_Stat_Formats;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <vector>
#include "../../core/util/MathArrays.h"
#include "../../core/util/MathArrays.h"

/**
 * Implements Chi-Square test statistics.
 * <p>
 * This implementation handles both known and unknown distributions.
 * <p>
 * Two samples tests can be used when the distribution is unknown <i>a priori</i>
 * but provided by one sample, or when the hypothesis under test is that the two
 * samples come from the same underlying distribution.
 */
class Chi_Square_Test 
{
public:
    /**
     * Computes the <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm">
     * Chi-Square statistic</a> comparing <code>observed</code> and <code>expected</code>
     * frequency counts.
     * <p>
     * This statistic can be used to perform a Chi-Square test evaluating the NULL
     * hypothesis that the observed counts follow the expected distribution.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Expected counts must all be positive.</li>
     * <li>Observed counts must all be &ge; 0.</li>
     * <li>The observed and expected arrays must have the same length and
     * their common length must be at least 2.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>Illegal_Argument_Exception</code> is thrown.
     * <p>
     * <strong>Note: </strong>This implementation rescales the
     * <code>expected</code> array if necessary to ensure that the sum of the
     * expected and observed counts are equal.
     *
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @return chi_square test statistic
     * @ if <code>observed</code> has negative entries
     * @ if <code>expected</code> has entries that are
     * not strictly positive
     * @ if the arrays length is less than 2
     */
    double chi_square(const std::vector<double> expected, const std::vector<long>& observed)
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
        double ratio = 1.0;
        bool rescale = false;
        if (std::abs(sum_expected - sum_observed) > 10E-6) 
        {
            ratio = sum_observed / sum_expected;
            rescale = true;
        }
        double sum_sq = 0.0;
        for (int i{}; i < observed.size(); i++) 
        {
            if (rescale) 
            {
                const double dev = observed[i] - ratio * expected[i];
                sum_sq += dev * dev / (ratio * expected[i]);
            }
            else 
            {
                const double dev = observed[i] - expected[i];
                sum_sq += dev * dev / expected[i];
            }
        }
        return sum_sq;
    }

    /**
     * Returns the <i>observed significance level</i>, or <a href=
     * "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
     * p-value</a>, associated with a
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm">
     * Chi-square goodness of fit test</a> comparing the <code>observed</code>
     * frequency counts to those in the <code>expected</code> array.
     * <p>
     * The number returned is the smallest significance level at which one can reject
     * the NULL hypothesis that the observed counts conform to the frequency distribution
     * described by the expected counts.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Expected counts must all be positive.</li>
     * <li>Observed counts must all be &ge; 0.</li>
     * <li>The observed and expected arrays must have the same length and
     * their common length must be at least 2.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>Illegal_Argument_Exception</code> is thrown.
     * <p>
     * <strong>Note: </strong>This implementation rescales the
     * <code>expected</code> array if necessary to ensure that the sum of the
     * expected and observed counts are equal.
     *
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @return p-value
     * @ if <code>observed</code> has negative entries
     * @ if <code>expected</code> has entries that are
     * not strictly positive
     * @ if the arrays length is less than 2
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    double chi_square_test(const std::vector<double> expected, const std::vector<long> observed)
    {
        const Chi_Squared_Distribution distribution = Chi_Squared_Distribution(expected.size() - 1.0);
        return 1.0 - distribution.cumulative_probability(chi_square(expected, observed));
    }

    /**
     * Performs a <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm">
     * Chi-square goodness of fit test</a> evaluating the NULL hypothesis that the
     * observed counts conform to the frequency distribution described by the expected
     * counts, with significance level <code>alpha</code>.  Returns true iff the NULL
     * hypothesis can be rejected with 100 * (1 - alpha) percent confidence.
     * <p>
     * <strong>Example:</strong><br>
     * To test the hypothesis that <code>observed</code> follows
     * <code>expected</code> at the 99% level, use
     * <code>chi_square_test(expected, observed, 0.01)</code>
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Expected counts must all be positive.</li>
     * <li>Observed counts must all be &ge; 0.</li>
     * <li>The observed and expected arrays must have the same length and
     * their common length must be at least 2.</li>
     * <li><code> 0 &lt; alpha &lt; 0.5</code></li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>Illegal_Argument_Exception</code> is thrown.
     * <p>
     * <strong>Note: </strong>This implementation rescales the
     * <code>expected</code> array if necessary to ensure that the sum of the
     * expected and observed counts are equal.
     *
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @param alpha significance level of the test
     * @return true iff NULL hypothesis can be rejected with confidence
     * 1 - alpha
     * @ if <code>observed</code> has negative entries
     * @ if <code>expected</code> has entries that are
     * not strictly positive
     * @ if the arrays length is less than 2
     * @ if <code>alpha</code> is not in the range (0, 0.5]
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    bool chi_square_test(const std::vector<double> expected, const std::vector<long> observed, const double& alpha)
    {

        if ((alpha <= 0) || (alpha > 0.5)) 
        {
            throw (Localized_Stat_Formats.OUT_OF_BOUND_SIGNIFICANCE_LEVEL, alpha, 0, 0.5);
        }
        return chi_square_test(expected, observed) < alpha;

    }

    /**
     * Computes the Chi-Square statistic associated with a
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section4/prc45.htm">
     * chi-square test of independence</a> based on the input <code>counts</code>
     * array, viewed as a two-way table.
     * <p>
     * The rows of the 2-way table are
     * <code>count[0], ... , count[count.size() - 1] </code>
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>All counts must be &ge; 0.</li>
     * <li>The count array must be rectangular (i.e. all count[i] subarrays
     * must have the same length).</li>
     * <li>The 2-way table represented by <code>counts</code> must have at
     * least 2 columns and at least 2 rows.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>Illegal_Argument_Exception</code> is thrown.
     *
     * @param counts array representation of 2-way table
     * @return chi_square test statistic
     * @Null_Argument_Exception if the array is NULL
     * @ if the array is not rectangular
     * @ if {@code counts} has negative entries
     */
    double chi_square(const std::vector<std::vector<long>> counts)
    {

        check_array(counts);
        int n_rows = counts.size();
        int n_cols = counts[0].size();

        // compute row, column and total sums
        std::vector<double> row_sum = std::vector<double>(n_rows];
        std::vector<double> col_sum = std::vector<double>(n_cols];
        double total = 0.0;
        for (int row{}; row < n_rows; row++) 
        {
            for (int col{};  col < n_cols; col++) 
            {
                row_sum[row] += counts[row][col];
                col_sum[col] += counts[row][col];
                total += counts[row][col];
            }
        }

        // compute expected counts and chi-square
        double sum_sq = 0.0;
        for (int row{}; row < n_rows; row++) 
        {
            for (int col{};  col < n_cols; col++) 
            {
                const double expected = (row_sum[row] * col_sum[col]) / total;
                sum_sq += ((counts[row][col] - expected) *
                        (counts[row][col] - expected)) / expected;
            }
        }
        return sum_sq;
    }

    /**
     * Returns the <i>observed significance level</i>, or <a href=
     * "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
     * p-value</a>, associated with a
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section4/prc45.htm">
     * chi-square test of independence</a> based on the input <code>counts</code>
     * array, viewed as a two-way table.
     * <p>
     * The rows of the 2-way table are
     * <code>count[0], ... , count[count.size() - 1] </code>
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>All counts must be &ge; 0.</li>
     * <li>The count array must be rectangular (i.e. all count[i] subarrays must have
     * the same length).</li>
     * <li>The 2-way table represented by <code>counts</code> must have at least 2
     * columns and at least 2 rows.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>Illegal_Argument_Exception</code> is thrown.
     *
     * @param counts array representation of 2-way table
     * @return p-value
     * @Null_Argument_Exception if the array is NULL
     * @ if the array is not rectangular
     * @ if {@code counts} has negative entries
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public double chi_square_test(const std::vector<std::vector<long>> counts)
    {
        check_array(counts);
        double df = (static_cast<double>( counts.size() -1) * (static_cast<double>( counts[0].size() - 1);
        const Chi_Squared_Distribution distribution = Chi_Squared_Distribution(df);
        return 1 - distribution.cumulative_probability(chi_square(counts));
    }

    /**
     * Performs a <a href="http://www.itl.nist.gov/div898/handbook/prc/section4/prc45.htm">
     * chi-square test of independence</a> evaluating the NULL hypothesis that the
     * classifications represented by the counts in the columns of the input 2-way table
     * are independent of the rows, with significance level <code>alpha</code>.
     * Returns true iff the NULL hypothesis can be rejected with 100 * (1 - alpha) percent
     * confidence.
     * <p>
     * The rows of the 2-way table are
     * <code>count[0], ... , count[count.size() - 1] </code>
     * <p>
     * <strong>Example:</strong><br>
     * To test the NULL hypothesis that the counts in
     * <code>count[0], ... , count[count.size() - 1] </code>
     * all correspond to the same underlying probability distribution at the 99% level, * use <code>chi_square_test(counts, 0.01)</code>.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>All counts must be &ge; 0.</li>
     * <li>The count array must be rectangular (i.e. all count[i] subarrays must have the
     * same length).</li>
     * <li>The 2-way table represented by <code>counts</code> must have at least 2 columns and
     * at least 2 rows.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>Illegal_Argument_Exception</code> is thrown.
     *
     * @param counts array representation of 2-way table
     * @param alpha significance level of the test
     * @return true iff NULL hypothesis can be rejected with confidence
     * 1 - alpha
     * @Null_Argument_Exception if the array is NULL
     * @ if the array is not rectangular
     * @ if {@code counts} has any negative entries
     * @ if <code>alpha</code> is not in the range (0, 0.5]
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    bool chi_square_test(const std::vector<std::vector<long>>& counts, const double& alpha)
    {
        if ((alpha <= 0) || (alpha > 0.5)) 
        {
            throw (Localized_Stat_Formats.OUT_OF_BOUND_SIGNIFICANCE_LEVEL, alpha, 0, 0.5);
        }
        return chi_square_test(counts) < alpha;
    }

    /**
     * Computes a
     * <a href="http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/chi2samp.htm">
     * Chi-Square two sample test statistic</a> comparing bin frequency counts
     * in <code>observed1</code> and <code>observed2</code>.
     * <p>
     * The sums of frequency counts in the two samples are not required to be the
     * same. The formula used to compute the test statistic is
     * <p>
     * <code>
     * &sum;[(K * observed1[i] - observed2[i]/K)<sup>2</sup> / (observed1[i] + observed2[i])]
     * </code>
     * <p>
     * where
     * <p>
     * <code>K = &sqrt;[&sum;(observed2 / &sum;(observed1)]</code>
     * <p>
     * This statistic can be used to perform a Chi-Square test evaluating the
     * NULL hypothesis that both observed counts follow the same distribution.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Observed counts must be non-negative.</li>
     * <li>Observed counts for a specific bin must not both be zero.</li>
     * <li>Observed counts for a specific sample must not all be 0.</li>
     * <li>The arrays <code>observed1</code> and <code>observed2</code> must have
     * the same length and their common length must be at least 2.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>Illegal_Argument_Exception</code> is thrown.
     *
     * @param observed1 array of observed frequency counts of the first data set
     * @param observed2 array of observed frequency counts of the second data set
     * @return chi_square test statistic
     * @ the the length of the arrays does not match
     * @ if any entries in <code>observed1</code> or
     * <code>observed2</code> are negative
     * @ if either all counts of <code>observed1</code> or
     * <code>observed2</code> are zero, or if the count at some index is zero
     * for both arrays
     */
    double chi_square_data_sets_comparison(std::vector<long>& observed1, std::vector<long>& observed2)
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
        for (int i{}; i < observed1.size(); i++) 
        {
            count_sum1 += observed1[i];
            count_sum2 += observed2[i];
        }
        // Ensure neither sample is uniformly 0
        if (count_sum1 == 0 || count_sum2 == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::ZERO_NOT_ALLOWED);
        }
        // Compare and compute weight only if different
        double weight = 0.0;
        bool unequal_counts = count_sum1 != count_sum2;
        if (unequal_counts) 
        {
            weight = std::sqrt(static_cast<double>( count_sum1 / static_cast<double>( count_sum2);
        }
        // Compute Chi_Square statistic
        double sum_sq = 0.0;
        for (int i{}; i < observed1.size(); i++) 
        {
            if (observed1[i] == 0 && observed2[i] == 0) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::OBSERVED_COUNTS_BOTTH_ZERO_FOR_ENTRY, i);
            }
            else 
            {
                const double obs1 = observed1[i];
                const double obs2 = observed2[i];
                const double dev;
                if (unequal_counts) { // apply weights
                    dev = obs1/weight - obs2 * weight;
                }
                else 
                {
                    dev = obs1 - obs2;
                }
                sum_sq += (dev * dev) / (obs1 + obs2);
            }
        }
        return sum_sq;
    }

    /**
     * Returns the <i>observed significance level</i>, or <a href=
     * "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
     * p-value</a>, associated with a Chi-Square two sample test comparing
     * bin frequency counts in <code>observed1</code> and
     * <code>observed2</code>.
     * <p>
     * The number returned is the smallest significance level at which one
     * can reject the NULL hypothesis that the observed counts conform to the
     * same distribution.
     * <p>
     * See {@link #chi_square_data_sets_comparison(long[], long[])} for details
     * on the formula used to compute the test statistic. The degrees of
     * of freedom used to perform the test is one less than the common length
     * of the input observed count arrays.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Observed counts must be non-negative.</li>
     * <li>Observed counts for a specific bin must not both be zero.</li>
     * <li>Observed counts for a specific sample must not all be 0.</li>
     * <li>The arrays <code>observed1</code> and <code>observed2</code> must
     * have the same length and their common length must be at least 2.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>Illegal_Argument_Exception</code> is thrown.
     *
     * @param observed1 array of observed frequency counts of the first data set
     * @param observed2 array of observed frequency counts of the second data set
     * @return p-value
     * @ the the length of the arrays does not match
     * @ if any entries in <code>observed1</code> or
     * <code>observed2</code> are negative
     * @ if either all counts of <code>observed1</code> or
     * <code>observed2</code> are zero, or if the count at the same index is zero
     * for both arrays
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    double chi_square_test_data_sets_comparison(std::vector<long> observed1, std::vector<long> observed2)
    {
        const Chi_Squared_Distribution distribution =
                Chi_Squared_Distribution(static_cast<double>( observed1.size() - 1);
        return 1 - distribution.cumulative_probability(
                chi_square_data_sets_comparison(observed1, observed2));
    }

    /**
     * Performs a Chi-Square two sample test comparing two binned data
     * sets. The test evaluates the NULL hypothesis that the two lists of
     * observed counts conform to the same frequency distribution, with
     * significance level <code>alpha</code>.  Returns true iff the NULL
     * hypothesis can be rejected with 100 * (1 - alpha) percent confidence.
     * <p>
     * See {@link #chi_square_data_sets_comparison(long[], long[])} for
     * details on the formula used to compute the Chisquare statistic used
     * in the test. The degrees of of freedom used to perform the test is
     * one less than the common length of the input observed count arrays.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Observed counts must be non-negative.</li>
     * <li>Observed counts for a specific bin must not both be zero.</li>
     * <li>Observed counts for a specific sample must not all be 0.</li>
     * <li>The arrays <code>observed1</code> and <code>observed2</code> must
     * have the same length and their common length must be at least 2.</li>
     * <li><code> 0 &lt; alpha &lt; 0.5</code></li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>Illegal_Argument_Exception</code> is thrown.
     *
     * @param observed1 array of observed frequency counts of the first data set
     * @param observed2 array of observed frequency counts of the second data set
     * @param alpha significance level of the test
     * @return true iff NULL hypothesis can be rejected with confidence
     * 1 - alpha
     * @ the the length of the arrays does not match
     * @ if any entries in <code>observed1</code> or
     * <code>observed2</code> are negative
     * @ if either all counts of <code>observed1</code> or
     * <code>observed2</code> are zero, or if the count at the same index is zero
     * for both arrays
     * @ if <code>alpha</code> is not in the range (0, 0.5]
     * @Math_Illegal_State_Exception if an error occurs performing the test
     */
    bool chi_square_test_data_sets_comparison(const std::vector<long> observed1, const std::vector<long> observed2, const double& alpha)
    {

        if (alpha <= 0 ||
            alpha > 0.5) 
            {
            throw (Localized_Stat_Formats.OUT_OF_BOUND_SIGNIFICANCE_LEVEL, alpha, 0, 0.5);
        }
        return chi_square_test_data_sets_comparison(observed1, observed2) < alpha;

    }

private:
    /**
     * Checks to make sure that the input std::vector<std::vector<long>> array is rectangular, * has at least 2 rows and 2 columns, and has all non-negative entries.
     *
     * @param in input 2-way table to check
     * @Null_Argument_Exception if the array is NULL
     * @ if the array is not valid
     * @ if the array contains any negative entries
     */
    void check_array(const std::vector<std::vector<long>> in)
    {

        if (in.size() < 2) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, in.size(), 2);
        }

        if (in[0].size() < 2) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, in[0].size(), 2);
        }

        Math_Arrays::check_rectangular(in);
        Math_Arrays::check_non_negative(in);
    }

};