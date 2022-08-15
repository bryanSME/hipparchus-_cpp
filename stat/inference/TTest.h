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

//import org.hipparchus.distribution.continuous.T_Distribution;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.stat.Localized_Stat_Formats;
//import org.hipparchus.stat.Stat_Utils;
//import org.hipparchus.stat.descriptive.Statistical_Summary;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/**
 * An implementation for Student's t-tests.
 * <p>
 * Tests can be:
 * <ul>
 * <li>One-sample or two-sample</li>
 * <li>One-sided or two-sided</li>
 * <li>Paired or unpaired (for two-sample tests)</li>
 * <li>Homoscedastic (equal variance assumption) or heteroscedastic
 * (for two sample tests)</li>
 * <li>Fixed significance level (bool-valued) or returning p-values.</li>
 * </ul>
 * <p>
 * Test statistics are available for all tests.  Methods including "Test" in
 * in their names perform tests, all other methods return t-statistics.  Among
 * the "Test" methods, <code>double-</code>valued methods return p-values;
 * <code>bool-</code>valued methods perform fixed significance level tests.
 * Significance levels are always specified as numbers between 0 and 0.5
 * (e.g. tests at the 95% level  use <code>alpha=0.05</code>).
 * <p>
 * Input to tests can be either <code>std::vector<double></code> arrays or
 * {@link Statistical_Summary} instances.
 * <p>
 * Uses Hipparchus {@link org.hipparchus.distribution.continuous.T_Distribution}
 * implementation to estimate exact p-values.
 */
class T_Test 
{

    /**
     * Computes a paired, 2-sample t-statistic based on the data in the input
     * arrays.  The t-statistic returned is equivalent to what would be returned by
     * computing the one-sample t-statistic {@link #t(double, std::vector<double>)}, with
     * <code>mu = 0</code> and the sample array consisting of the (signed)
     * differences between corresponding entries in <code>sample1</code> and
     * <code>sample2.</code>
     * <p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The input arrays must have the same length and their common length
     * must be at least 2.
     * </li></ul></p>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return t statistic
     * @Null_Argument_Exception if the arrays are <code>null</code>
     * @ if the arrays are empty
     * @ if the length of the arrays is not equal
     * @ if the length of the arrays is &lt; 2
     */
    public double paired_t(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception 
        {

        check_sample_data(sample1);
        check_sample_data(sample2);
        double mean_difference = Stat_Utils.mean_difference(sample1, sample2);
        return t(mean_difference, 0, Stat_Utils.variance_difference(sample1, sample2, mean_difference), sample1.size());
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i> p-value</i>, associated with a paired, two-sample, two-tailed t-test
     * based on the data in the input arrays.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the NULL hypothesis that the mean of the paired
     * differences is 0 in favor of the two-sided alternative that the mean paired
     * difference is not equal to 0. For a one-sided test, divide the returned
     * value by 2.</p>
     * <p>
     * This test is equivalent to a one-sample t-test computed using
     * {@link #t_test(double, std::vector<double>)} with <code>mu = 0</code> and the sample
     * array consisting of the signed differences between corresponding elements of
     * <code>sample1</code> and <code>sample2.</code></p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the p-value depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The input array lengths must be the same and their common length must
     * be at least 2.
     * </li></ul></p>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return p-value for t-test
     * @Null_Argument_Exception if the arrays are <code>null</code>
     * @ if the arrays are empty
     * @ if the length of the arrays is not equal
     * @ if the length of the arrays is &lt; 2
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public double paired_t_test(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        double mean_difference = Stat_Utils.mean_difference(sample1, sample2);
        return t_test(mean_difference, 0, Stat_Utils.variance_difference(sample1, sample2, mean_difference), sample1.size());
    }

    /**
     * Performs a paired t-test evaluating the NULL hypothesis that the
     * mean of the paired differences between <code>sample1</code> and
     * <code>sample2</code> is 0 in favor of the two-sided alternative that the
     * mean paired difference is not equal to 0, with significance level
     * <code>alpha</code>.
     * <p>
     * Returns <code>true</code> iff the NULL hypothesis can be rejected with
     * confidence <code>1 - alpha</code>.  To perform a 1-sided test, use
     * <code>alpha * 2</code></p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The input array lengths must be the same and their common length
     * must be at least 2.
     * </li>
     * <li> <code> 0 &lt; alpha &lt; 0.5 </code>
     * </li></ul></p>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @param alpha significance level of the test
     * @return true if the NULL hypothesis can be rejected with
     * confidence 1 - alpha
     * @Null_Argument_Exception if the arrays are <code>null</code>
     * @ if the arrays are empty
     * @ if the length of the arrays is not equal
     * @ if the length of the arrays is &lt; 2
     * @ if <code>alpha</code> is not in the range (0, 0.5]
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public bool paired_t_test(const std::vector<double> sample1, const std::vector<double> sample2, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_significance_level(alpha);
        return paired_t_test(sample1, sample2) < alpha;

    }

    /**
     * Computes a <a href="http://www.itl.nist.gov/div898/handbook/prc/section2/prc22.htm#formula">
     * t statistic </a> given observed values and a comparison constant.
     * <p>
     * This statistic can be used to perform a one sample t-test for the mean.
     * </p><p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The observed array length must be at least 2.
     * </li></ul></p>
     *
     * @param mu comparison constant
     * @param observed array of values
     * @return t statistic
     * @Null_Argument_Exception if <code>observed</code> is <code>null</code>
     * @ if the length of <code>observed</code> is &lt; 2
     */
    public double t(const double& mu, const std::vector<double> observed)
        , Null_Argument_Exception 
        {

        check_sample_data(observed);
        // No try-catch or advertised exception because args have just been checked
        return t(Stat_Utils.mean(observed), mu, Stat_Utils.variance(observed), observed.size());
    }

    /**
     * Computes a <a href="http://www.itl.nist.gov/div898/handbook/prc/section2/prc22.htm#formula">
     * t statistic </a> to use in comparing the mean of the dataset described by
     * <code>sample_stats</code> to <code>mu</code>.
     * <p>
     * This statistic can be used to perform a one sample t-test for the mean.
     * </p><p>
     * <strong>Preconditions</strong>: <ul>
     * <li><code>observed.get_n() &ge; 2</code>.
     * </li></ul></p>
     *
     * @param mu comparison constant
     * @param sample_stats Descriptive_Statistics holding sample summary statitstics
     * @return t statistic
     * @Null_Argument_Exception if <code>sample_stats</code> is <code>null</code>
     * @ if the number of samples is &lt; 2
     */
    public double t(const double& mu, const Statistical_Summary sample_stats)
        , Null_Argument_Exception 
        {

        check_sample_data(sample_stats);
        return t(sample_stats.get_mean(), mu, sample_stats.get_variance(), sample_stats.get_n());
    }

    /**
     * Computes a 2-sample t statistic,  under the hypothesis of equal
     * subpopulation variances.  To compute a t-statistic without the
     * equal variances hypothesis, use {@link #t(std::vector<double>, std::vector<double>)}.
     * <p>
     * This statistic can be used to perform a (homoscedastic) two-sample
     * t-test to compare sample means.</p>
     * <p>
     * The t-statistic is</p>
     * <p>
     * &nbsp;&nbsp;<code>  t = (m1 - m2) / (sqrt(1/n1 +1/n2) sqrt(var))</code>
     * </p><p>
     * where <strong><code>n1</code></strong> is the size of first sample;
     * <strong><code> n2</code></strong> is the size of second sample;
     * <strong><code> m1</code></strong> is the mean of first sample;
     * <strong><code> m2</code></strong> is the mean of second sample</li>
     * </ul>
     * and <strong><code>var</code></strong> is the pooled variance estimate:
     * </p><p>
     * <code>var = sqrt(((n1 - 1)var1 + (n2 - 1)var2) / ((n1-1) + (n2-1)))</code>
     * </p><p>
     * with <strong><code>var1</code></strong> the variance of the first sample and
     * <strong><code>var2</code></strong> the variance of the second sample.
     * </p><p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li></ul></p>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return t statistic
     * @Null_Argument_Exception if the arrays are <code>null</code>
     * @ if the length of the arrays is &lt; 2
     */
    public double homoscedastic_t(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception 
        {

        check_sample_data(sample1);
        check_sample_data(sample2);
        // No try-catch or advertised exception because args have just been checked
        return homoscedastic_t(Stat_Utils.mean(sample1), Stat_Utils.mean(sample2), Stat_Utils.variance(sample1), Stat_Utils.variance(sample2), sample1.size(), sample2.size());
    }

    /**
     * Computes a 2-sample t statistic, without the hypothesis of equal
     * subpopulation variances.  To compute a t-statistic assuming equal
     * variances, use {@link #homoscedastic_t(std::vector<double>, std::vector<double>)}.
     * <p>
     * This statistic can be used to perform a two-sample t-test to compare
     * sample means.</p>
     * <p>
     * The t-statistic is</p>
     * <p>
     * &nbsp;&nbsp; <code>  t = (m1 - m2) / sqrt(var1/n1 + var2/n2)</code>
     * </p><p>
     *  where <strong><code>n1</code></strong> is the size of the first sample
     * <strong><code> n2</code></strong> is the size of the second sample;
     * <strong><code> m1</code></strong> is the mean of the first sample;
     * <strong><code> m2</code></strong> is the mean of the second sample;
     * <strong><code> var1</code></strong> is the variance of the first sample;
     * <strong><code> var2</code></strong> is the variance of the second sample;
     * </p><p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li></ul></p>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return t statistic
     * @Null_Argument_Exception if the arrays are <code>null</code>
     * @ if the length of the arrays is &lt; 2
     */
    public double t(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception 
        {

        check_sample_data(sample1);
        check_sample_data(sample2);
        // No try-catch or advertised exception because args have just been checked
        return t(Stat_Utils.mean(sample1), Stat_Utils.mean(sample2), Stat_Utils.variance(sample1), Stat_Utils.variance(sample2), sample1.size(), sample2.size());
    }

    /**
     * Computes a 2-sample t statistic </a>, comparing the means of the datasets
     * described by two {@link Statistical_Summary} instances, without the
     * assumption of equal subpopulation variances.  Use
     * {@link #homoscedastic_t(Statistical_Summary, Statistical_Summary)} to
     * compute a t-statistic under the equal variances assumption.
     * <p>
     * This statistic can be used to perform a two-sample t-test to compare
     * sample means.</p>
     * <p>
      * The returned  t-statistic is</p>
     * <p>
     * &nbsp;&nbsp; <code>  t = (m1 - m2) / sqrt(var1/n1 + var2/n2)</code>
     * </p><p>
     * where <strong><code>n1</code></strong> is the size of the first sample;
     * <strong><code> n2</code></strong> is the size of the second sample;
     * <strong><code> m1</code></strong> is the mean of the first sample;
     * <strong><code> m2</code></strong> is the mean of the second sample
     * <strong><code> var1</code></strong> is the variance of the first sample;
     * <strong><code> var2</code></strong> is the variance of the second sample
     * </p><p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The datasets described by the two Univariates must each contain
     * at least 2 observations.
     * </li></ul></p>
     *
     * @param sample_stats1 Statistical_Summary describing data from the first sample
     * @param sample_stats2 Statistical_Summary describing data from the second sample
     * @return t statistic
     * @Null_Argument_Exception if the sample statistics are <code>null</code>
     * @ if the number of samples is &lt; 2
     */
    public double t(const Statistical_Summary sample_stats1, const Statistical_Summary sample_stats2)
        , Null_Argument_Exception 
        {

        check_sample_data(sample_stats1);
        check_sample_data(sample_stats2);
        return t(sample_stats1.get_mean(), sample_stats2.get_mean(), sample_stats1.get_variance(), sample_stats2.get_variance(), sample_stats1.get_n(), sample_stats2.get_n());
    }

    /**
     * Computes a 2-sample t statistic, comparing the means of the datasets
     * described by two {@link Statistical_Summary} instances, under the
     * assumption of equal subpopulation variances.  To compute a t-statistic
     * without the equal variances assumption, use
     * {@link #t(Statistical_Summary, Statistical_Summary)}.
     * <p>
     * This statistic can be used to perform a (homoscedastic) two-sample
     * t-test to compare sample means.</p>
     * <p>
     * The t-statistic returned is</p>
     * <p>
     * &nbsp;&nbsp;<code>  t = (m1 - m2) / (sqrt(1/n1 +1/n2) sqrt(var))</code>
     * </p><p>
     * where <strong><code>n1</code></strong> is the size of first sample;
     * <strong><code> n2</code></strong> is the size of second sample;
     * <strong><code> m1</code></strong> is the mean of first sample;
     * <strong><code> m2</code></strong> is the mean of second sample
     * and <strong><code>var</code></strong> is the pooled variance estimate:
     * </p><p>
     * <code>var = sqrt(((n1 - 1)var1 + (n2 - 1)var2) / ((n1-1) + (n2-1)))</code>
     * </p><p>
     * with <strong><code>var1</code></strong> the variance of the first sample and
     * <strong><code>var2</code></strong> the variance of the second sample.
     * </p><p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The datasets described by the two Univariates must each contain
     * at least 2 observations.
     * </li></ul></p>
     *
     * @param sample_stats1 Statistical_Summary describing data from the first sample
     * @param sample_stats2 Statistical_Summary describing data from the second sample
     * @return t statistic
     * @Null_Argument_Exception if the sample statistics are <code>null</code>
     * @ if the number of samples is &lt; 2
     */
    public double homoscedastic_t(const Statistical_Summary sample_stats1, const Statistical_Summary sample_stats2)
        , Null_Argument_Exception 
        {

        check_sample_data(sample_stats1);
        check_sample_data(sample_stats2);
        return homoscedastic_t(sample_stats1.get_mean(), sample_stats2.get_mean(), sample_stats1.get_variance(), sample_stats2.get_variance(), sample_stats1.get_n(), sample_stats2.get_n());
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a one-sample, two-tailed t-test
     * comparing the mean of the input array with the constant <code>mu</code>.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the NULL hypothesis that the mean equals
     * <code>mu</code> in favor of the two-sided alternative that the mean
     * is different from <code>mu</code>. For a one-sided test, divide the
     * returned value by 2.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">here</a>
     * </p><p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The observed array length must be at least 2.
     * </li></ul></p>
     *
     * @param mu constant value to compare sample mean against
     * @param sample array of sample data values
     * @return p-value
     * @Null_Argument_Exception if the sample array is <code>null</code>
     * @ if the length of the array is &lt; 2
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public double t_test(const double& mu, const std::vector<double> sample)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_sample_data(sample);
        // No try-catch or advertised exception because args have just been checked
        return t_test(Stat_Utils.mean(sample), mu, Stat_Utils.variance(sample), sample.size());
    }

    /**
     * Performs a <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">
     * two-sided t-test</a> evaluating the NULL hypothesis that the mean of the population from
     * which <code>sample</code> is drawn equals <code>mu</code>.
     * <p>
     * Returns <code>true</code> iff the NULL hypothesis can be
     * rejected with confidence <code>1 - alpha</code>.  To
     * perform a 1-sided test, use <code>alpha * 2</code></p>
     * <p>
     * <strong>Examples:</strong><br><ol>
     * <li>To test the (2-sided) hypothesis <code>sample mean = mu </code> at
     * the 95% level, use <br><code>t_test(mu, sample, 0.05) </code>
     * </li>
     * <li>To test the (one-sided) hypothesis <code> sample mean &lt; mu </code>
     * at the 99% level, first verify that the measured sample mean is less
     * than <code>mu</code> and then use
     * <br><code>t_test(mu, sample, 0.02) </code>
     * </li></ol></p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the one-sample
     * parametric t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/sg_glos.html#one-sample">here</a>
     * </p><p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The observed array length must be at least 2.
     * </li></ul></p>
     *
     * @param mu constant value to compare sample mean against
     * @param sample array of sample data values
     * @param alpha significance level of the test
     * @return p-value
     * @Null_Argument_Exception if the sample array is <code>null</code>
     * @ if the length of the array is &lt; 2
     * @ if <code>alpha</code> is not in the range (0, 0.5]
     * @Math_Illegal_State_Exception if an error computing the p-value
     */
    public bool t_test(const double& mu, const std::vector<double> sample, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_significance_level(alpha);
        return t_test(mu, sample) < alpha;
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a one-sample, two-tailed t-test
     * comparing the mean of the dataset described by <code>sample_stats</code>
     * with the constant <code>mu</code>.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the NULL hypothesis that the mean equals
     * <code>mu</code> in favor of the two-sided alternative that the mean
     * is different from <code>mu</code>. For a one-sided test, divide the
     * returned value by 2.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The sample must contain at least 2 observations.
     * </li></ul></p>
     *
     * @param mu constant value to compare sample mean against
     * @param sample_stats Statistical_Summary describing sample data
     * @return p-value
     * @Null_Argument_Exception if <code>sample_stats</code> is <code>null</code>
     * @ if the number of samples is &lt; 2
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public double t_test(const double& mu, const Statistical_Summary sample_stats)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_sample_data(sample_stats);
        return t_test(sample_stats.get_mean(), mu, sample_stats.get_variance(), sample_stats.get_n());
    }

    /**
     * Performs a <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">
     * two-sided t-test</a> evaluating the NULL hypothesis that the mean of the
     * population from which the dataset described by <code>stats</code> is
     * drawn equals <code>mu</code>.
     * <p>
     * Returns <code>true</code> iff the NULL hypothesis can be rejected with
     * confidence <code>1 - alpha</code>.  To  perform a 1-sided test, use
     * <code>alpha * 2.</code></p>
     * <p>
     * <strong>Examples:</strong><br><ol>
     * <li>To test the (2-sided) hypothesis <code>sample mean = mu </code> at
     * the 95% level, use <br><code>t_test(mu, sample_stats, 0.05) </code>
     * </li>
     * <li>To test the (one-sided) hypothesis <code> sample mean &lt; mu </code>
     * at the 99% level, first verify that the measured sample mean is less
     * than <code>mu</code> and then use
     * <br><code>t_test(mu, sample_stats, 0.02) </code>
     * </li></ol></p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the one-sample
     * parametric t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/sg_glos.html#one-sample">here</a>
     * </p><p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The sample must include at least 2 observations.
     * </li></ul></p>
     *
     * @param mu constant value to compare sample mean against
     * @param sample_stats Statistical_Summary describing sample data values
     * @param alpha significance level of the test
     * @return p-value
     * @Null_Argument_Exception if <code>sample_stats</code> is <code>null</code>
     * @ if the number of samples is &lt; 2
     * @ if <code>alpha</code> is not in the range (0, 0.5]
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public bool t_test(const double& mu, const Statistical_Summary sample_stats, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_significance_level(alpha);
        return t_test(mu, sample_stats) < alpha;
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a two-sample, two-tailed t-test
     * comparing the means of the input arrays.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the NULL hypothesis that the two means are
     * equal in favor of the two-sided alternative that they are different.
     * For a one-sided test, divide the returned value by 2.</p>
     * <p>
     * The test does not assume that the underlying popuation variances are
     * equal  and it uses approximated degrees of freedom computed from the
     * sample data to compute the p-value.  The t-statistic used is as defined in
     * {@link #t(std::vector<double>, std::vector<double>)} and the Welch-_Satterthwaite approximation
     * to the degrees of freedom is used, * as described
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section3/prc31.htm">
     * here.</a>  To perform the test under the assumption of equal subpopulation
     * variances, use {@link #homoscedastic_t_test(std::vector<double>, std::vector<double>)}.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the p-value depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li></ul></p>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return p-value for t-test
     * @Null_Argument_Exception if the arrays are <code>null</code>
     * @ if the length of the arrays is &lt; 2
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public double t_test(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_sample_data(sample1);
        check_sample_data(sample2);
        // No try-catch or advertised exception because args have just been checked
        return t_test(Stat_Utils.mean(sample1), Stat_Utils.mean(sample2), Stat_Utils.variance(sample1), Stat_Utils.variance(sample2), sample1.size(), sample2.size());
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a two-sample, two-tailed t-test
     * comparing the means of the input arrays, under the assumption that
     * the two samples are drawn from subpopulations with equal variances.
     * To perform the test without the equal variances assumption, use
     * {@link #t_test(std::vector<double>, std::vector<double>)}.</p>
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the NULL hypothesis that the two means are
     * equal in favor of the two-sided alternative that they are different.
     * For a one-sided test, divide the returned value by 2.</p>
     * <p>
     * A pooled variance estimate is used to compute the t-statistic.  See
     * {@link #homoscedastic_t(std::vector<double>, std::vector<double>)}. The sum of the sample sizes
     * minus 2 is used as the degrees of freedom.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the p-value depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li></ul></p>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return p-value for t-test
     * @Null_Argument_Exception if the arrays are <code>null</code>
     * @ if the length of the arrays is &lt; 2
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public double homoscedastic_t_test(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_sample_data(sample1);
        check_sample_data(sample2);
        // No try-catch or advertised exception because args have just been checked
        return homoscedastic_t_test(Stat_Utils.mean(sample1), Stat_Utils.mean(sample2), Stat_Utils.variance(sample1), Stat_Utils.variance(sample2), sample1.size(), sample2.size());
    }

    /**
     * Performs a
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">
     * two-sided t-test</a> evaluating the NULL hypothesis that <code>sample1</code>
     * and <code>sample2</code> are drawn from populations with the same mean, * with significance level <code>alpha</code>.  This test does not assume
     * that the subpopulation variances are equal.  To perform the test assuming
     * equal variances, use
     * {@link #homoscedastic_t_test(std::vector<double>, std::vector<double>, double)}.
     * <p>
     * Returns <code>true</code> iff the NULL hypothesis that the means are
     * equal can be rejected with confidence <code>1 - alpha</code>.  To
     * perform a 1-sided test, use <code>alpha * 2</code></p>
     * <p>
     * See {@link #t(std::vector<double>, std::vector<double>)} for the formula used to compute the
     * t-statistic.  Degrees of freedom are approximated using the
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section3/prc31.htm">
     * Welch-_Satterthwaite approximation.</a></p>
     * <p>
     * <strong>Examples:</strong><br><ol>
     * <li>To test the (2-sided) hypothesis <code>mean 1 = mean 2 </code> at
     * the 95% level,  use
     * <br><code>t_test(sample1, sample2, 0.05). </code>
     * </li>
     * <li>To test the (one-sided) hypothesis <code> mean 1 &lt; mean 2 </code>, * at the 99% level, first verify that the measured  mean of <code>sample 1</code>
     * is less than the mean of <code>sample 2</code> and then use
     * <br><code>t_test(sample1, sample2, 0.02) </code>
     * </li></ol></p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li>
     * <li> <code> 0 &lt; alpha &lt; 0.5 </code>
     * </li></ul></p>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @param alpha significance level of the test
     * @return true if the NULL hypothesis can be rejected with
     * confidence 1 - alpha
     * @Null_Argument_Exception if the arrays are <code>null</code>
     * @ if the length of the arrays is &lt; 2
     * @ if <code>alpha</code> is not in the range (0, 0.5]
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public bool t_test(const std::vector<double> sample1, const std::vector<double> sample2, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_significance_level(alpha);
        return t_test(sample1, sample2) < alpha;
    }

    /**
     * Performs a
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">
     * two-sided t-test</a> evaluating the NULL hypothesis that <code>sample1</code>
     * and <code>sample2</code> are drawn from populations with the same mean, * with significance level <code>alpha</code>,  assuming that the
     * subpopulation variances are equal.  Use
     * {@link #t_test(std::vector<double>, std::vector<double>, double)} to perform the test without
     * the assumption of equal variances.
     * <p>
     * Returns <code>true</code> iff the NULL hypothesis that the means are
     * equal can be rejected with confidence <code>1 - alpha</code>.  To
     * perform a 1-sided test, use <code>alpha * 2.</code>  To perform the test
     * without the assumption of equal subpopulation variances, use
     * {@link #t_test(std::vector<double>, std::vector<double>, double)}.</p>
     * <p>
     * A pooled variance estimate is used to compute the t-statistic. See
     * {@link #t(std::vector<double>, std::vector<double>)} for the formula. The sum of the sample
     * sizes minus 2 is used as the degrees of freedom.</p>
     * <p>
     * <strong>Examples:</strong><br><ol>
     * <li>To test the (2-sided) hypothesis <code>mean 1 = mean 2 </code> at
     * the 95% level, use <br><code>t_test(sample1, sample2, 0.05). </code>
     * </li>
     * <li>To test the (one-sided) hypothesis <code> mean 1 &lt; mean 2, </code>
     * at the 99% level, first verify that the measured mean of
     * <code>sample 1</code> is less than the mean of <code>sample 2</code>
     * and then use
     * <br><code>t_test(sample1, sample2, 0.02) </code>
     * </li></ol></p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li>
     * <li> <code> 0 &lt; alpha &lt; 0.5 </code>
     * </li></ul></p>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @param alpha significance level of the test
     * @return true if the NULL hypothesis can be rejected with
     * confidence 1 - alpha
     * @Null_Argument_Exception if the arrays are <code>null</code>
     * @ if the length of the arrays is &lt; 2
     * @ if <code>alpha</code> is not in the range (0, 0.5]
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public bool homoscedastic_t_test(const std::vector<double> sample1, const std::vector<double> sample2, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_significance_level(alpha);
        return homoscedastic_t_test(sample1, sample2) < alpha;
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a two-sample, two-tailed t-test
     * comparing the means of the datasets described by two Statistical_Summary
     * instances.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the NULL hypothesis that the two means are
     * equal in favor of the two-sided alternative that they are different.
     * For a one-sided test, divide the returned value by 2.</p>
     * <p>
     * The test does not assume that the underlying population variances are
     * equal  and it uses approximated degrees of freedom computed from the
     * sample data to compute the p-value.   To perform the test assuming
     * equal variances, use
     * {@link #homoscedastic_t_test(Statistical_Summary, Statistical_Summary)}.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the p-value depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The datasets described by the two Univariates must each contain
     * at least 2 observations.
     * </li></ul></p>
     *
     * @param sample_stats1  Statistical_Summary describing data from the first sample
     * @param sample_stats2  Statistical_Summary describing data from the second sample
     * @return p-value for t-test
     * @Null_Argument_Exception if the sample statistics are <code>null</code>
     * @ if the number of samples is &lt; 2
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public double t_test(const Statistical_Summary sample_stats1, const Statistical_Summary sample_stats2)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_sample_data(sample_stats1);
        check_sample_data(sample_stats2);
        return t_test(sample_stats1.get_mean(), sample_stats2.get_mean(), sample_stats1.get_variance(), sample_stats2.get_variance(), sample_stats1.get_n(), sample_stats2.get_n());
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a two-sample, two-tailed t-test
     * comparing the means of the datasets described by two Statistical_Summary
     * instances, under the hypothesis of equal subpopulation variances. To
     * perform a test without the equal variances assumption, use
     * {@link #t_test(Statistical_Summary, Statistical_Summary)}.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the NULL hypothesis that the two means are
     * equal in favor of the two-sided alternative that they are different.
     * For a one-sided test, divide the returned value by 2.</p>
     * <p>
     * See {@link #homoscedastic_t(std::vector<double>, std::vector<double>)} for the formula used to
     * compute the t-statistic. The sum of the  sample sizes minus 2 is used as
     * the degrees of freedom.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the p-value depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">here</a>
     * </p><p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The datasets described by the two Univariates must each contain
     * at least 2 observations.
     * </li></ul></p>
     *
     * @param sample_stats1  Statistical_Summary describing data from the first sample
     * @param sample_stats2  Statistical_Summary describing data from the second sample
     * @return p-value for t-test
     * @Null_Argument_Exception if the sample statistics are <code>null</code>
     * @ if the number of samples is &lt; 2
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public double homoscedastic_t_test(const Statistical_Summary sample_stats1, const Statistical_Summary sample_stats2)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_sample_data(sample_stats1);
        check_sample_data(sample_stats2);
        return homoscedastic_t_test(sample_stats1.get_mean(), sample_stats2.get_mean(), sample_stats1.get_variance(), sample_stats2.get_variance(), sample_stats1.get_n(), sample_stats2.get_n());
    }

    /**
     * Performs a
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">
     * two-sided t-test</a> evaluating the NULL hypothesis that
     * <code>sample_stats1</code> and <code>sample_stats2</code> describe
     * datasets drawn from populations with the same mean, with significance
     * level <code>alpha</code>.   This test does not assume that the
     * subpopulation variances are equal.  To perform the test under the equal
     * variances assumption, use
     * {@link #homoscedastic_t_test(Statistical_Summary, Statistical_Summary)}.
     * <p>
     * Returns <code>true</code> iff the NULL hypothesis that the means are
     * equal can be rejected with confidence <code>1 - alpha</code>.  To
     * perform a 1-sided test, use <code>alpha * 2</code></p>
     * <p>
     * See {@link #t(std::vector<double>, std::vector<double>)} for the formula used to compute the
     * t-statistic.  Degrees of freedom are approximated using the
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section3/prc31.htm">
     * Welch-_Satterthwaite approximation.</a></p>
     * <p>
     * <strong>Examples:</strong><br><ol>
     * <li>To test the (2-sided) hypothesis <code>mean 1 = mean 2 </code> at
     * the 95%, use
     * <br><code>t_test(sample_stats1, sample_stats2, 0.05) </code>
     * </li>
     * <li>To test the (one-sided) hypothesis <code> mean 1 &lt; mean 2 </code>
     * at the 99% level,  first verify that the measured mean of
     * <code>sample 1</code> is less than  the mean of <code>sample 2</code>
     * and then use
     * <br><code>t_test(sample_stats1, sample_stats2, 0.02) </code>
     * </li></ol></p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>: <ul>
     * <li>The datasets described by the two Univariates must each contain
     * at least 2 observations.
     * </li>
     * <li> <code> 0 &lt; alpha &lt; 0.5 </code>
     * </li></ul></p>
     *
     * @param sample_stats1 Statistical_Summary describing sample data values
     * @param sample_stats2 Statistical_Summary describing sample data values
     * @param alpha significance level of the test
     * @return true if the NULL hypothesis can be rejected with
     * confidence 1 - alpha
     * @Null_Argument_Exception if the sample statistics are <code>null</code>
     * @ if the number of samples is &lt; 2
     * @ if <code>alpha</code> is not in the range (0, 0.5]
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     */
    public bool t_test(const Statistical_Summary sample_stats1, const Statistical_Summary sample_stats2, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {

        check_significance_level(alpha);
        return t_test(sample_stats1, sample_stats2) < alpha;
    }

    //----------------------------------------------- Protected methods

    /**
     * Computes approximate degrees of freedom for 2-sample t-test.
     *
     * @param v1 first sample variance
     * @param v2 second sample variance
     * @param n1 first sample n
     * @param n2 second sample n
     * @return approximate degrees of freedom
     */
    protected double df(double v1, double v2, double n1, double n2) 
    {
        return (((v1 / n1) + (v2 / n2)) * ((v1 / n1) + (v2 / n2))) /
        ((v1 * v1) / (n1 * n1 * (n1 - 1d)) + (v2 * v2) /
                (n2 * n2 * (n2 - 1d)));
    }

    /**
     * Computes t test statistic for 1-sample t-test.
     *
     * @param m sample mean
     * @param mu constant to test against
     * @param v sample variance
     * @param n sample n
     * @return t test statistic
     */
    protected double t(const double m, const double& mu, const double v, const double n) 
    {
        return (m - mu) / std::sqrt(v / n);
    }

    /**
     * Computes t test statistic for 2-sample t-test.
     * <p>
     * Does not assume that subpopulation variances are equal.</p>
     *
     * @param m1 first sample mean
     * @param m2 second sample mean
     * @param v1 first sample variance
     * @param v2 second sample variance
     * @param n1 first sample n
     * @param n2 second sample n
     * @return t test statistic
     */
    protected double t(const double m1, const double m2, const double v1, const double v2, const double n1, const double n2)  
    {
        return (m1 - m2) / std::sqrt((v1 / n1) + (v2 / n2));
    }

    /**
     * Computes t test statistic for 2-sample t-test under the hypothesis
     * of equal subpopulation variances.
     *
     * @param m1 first sample mean
     * @param m2 second sample mean
     * @param v1 first sample variance
     * @param v2 second sample variance
     * @param n1 first sample n
     * @param n2 second sample n
     * @return t test statistic
     */
    protected double homoscedastic_t(const double m1, const double m2, const double v1, const double v2, const double n1, const double n2)  
    {
        const double pooled_variance = ((n1  - 1) * v1 + (n2 -1) * v2 ) / (n1 + n2 - 2);
        return (m1 - m2) / std::sqrt(pooled_variance * (1.0/ n1 + 1.0/ n2));
    }

    /**
     * Computes p-value for 2-sided, 1-sample t-test.
     *
     * @param m sample mean
     * @param mu constant to test against
     * @param v sample variance
     * @param n sample n
     * @return p-value
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     * @ if n is not greater than 1
     */
    protected double t_test(const double m, const double& mu, const double v, const double n)
        , Math_Illegal_State_Exception 
        {

        const double t = std::abs(t(m, mu, v, n));
        const T_Distribution distribution = T_Distribution( n - 1);
        return 2.0 * distribution.cumulative_probability(-t);

    }

    /**
     * Computes p-value for 2-sided, 2-sample t-test.
     * <p>
     * Does not assume subpopulation variances are equal. Degrees of freedom
     * are estimated from the data.</p>
     *
     * @param m1 first sample mean
     * @param m2 second sample mean
     * @param v1 first sample variance
     * @param v2 second sample variance
     * @param n1 first sample n
     * @param n2 second sample n
     * @return p-value
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     * @ if the estimated degrees of freedom is not
     * strictly positive
     */
    protected double t_test(const double m1, const double m2, const double v1, const double v2, const double n1, const double n2)
        , Math_Illegal_State_Exception 
        {

        const double t = std::abs(t(m1, m2, v1, v2, n1, n2));
        const double degrees_of_freedom = df(v1, v2, n1, n2);
        const T_Distribution distribution = T_Distribution(degrees_of_freedom);
        return 2.0 * distribution.cumulative_probability(-t);

    }

    /**
     * Computes p-value for 2-sided, 2-sample t-test, under the assumption
     * of equal subpopulation variances.
     * <p>
     * The sum of the sample sizes minus 2 is used as degrees of freedom.</p>
     *
     * @param m1 first sample mean
     * @param m2 second sample mean
     * @param v1 first sample variance
     * @param v2 second sample variance
     * @param n1 first sample n
     * @param n2 second sample n
     * @return p-value
     * @Math_Illegal_State_Exception if an error occurs computing the p-value
     * @ if the estimated degrees of freedom is not
     * strictly positive
     */
    protected double homoscedastic_t_test(double m1, double m2, double v1, double v2, double n1, double n2)
        , Math_Illegal_State_Exception 
        {

        const double t = std::abs(homoscedastic_t(m1, m2, v1, v2, n1, n2));
        const double degrees_of_freedom = n1 + n2 - 2;
        const T_Distribution distribution = T_Distribution(degrees_of_freedom);
        return 2.0 * distribution.cumulative_probability(-t);

    }

    /**
     * Check significance level.
     *
     * @param alpha significance level
     * @ if the significance level is out of bounds.
     */
    private void check_significance_level(const double& alpha)
         
        {

        if (alpha <= 0 || alpha > 0.5) 
        {
            throw (Localized_Stat_Formats.SIGNIFICANCE_LEVEL, alpha, 0.0, 0.5);
        }

    }

    /**
     * Check sample data.
     *
     * @param data Sample data.
     * @Null_Argument_Exception if {@code data} is {@code NULL}.
     * @ if there is not enough sample data.
     */
    private void check_sample_data(const std::vector<double> data)
        , Null_Argument_Exception 
        {

        //Math_Utils::check_not_null(data, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
        if (data.size() < 2) 
        {
            throw (
                    Localized_Stat_Formats.INSUFFICIENT_DATA_FOR_T_STATISTIC, data.size(), 2, true);
        }

    }

    /**
     * Check sample data.
     *
     * @param stat Statistical summary.
     * @Null_Argument_Exception if {@code data} is {@code NULL}.
     * @ if there is not enough sample data.
     */
    private void check_sample_data(const Statistical_Summary stat)
        , Null_Argument_Exception 
        {

        //Math_Utils::check_not_null(stat);
        if (stat.get_n() < 2) 
        {
            throw (
                    Localized_Stat_Formats.INSUFFICIENT_DATA_FOR_T_STATISTIC, stat.get_n(), 2, true);
        }

    }

}


