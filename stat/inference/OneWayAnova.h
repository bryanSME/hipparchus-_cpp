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
//import java.util.Collection;

//import org.hipparchus.distribution.continuous.F_Distribution;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.stat.Localized_Stat_Formats;
//import org.hipparchus.stat.descriptive.Streaming_Statistics;
//import org.hipparchus.util.Math_Utils;
#include <vector>
#include "../LocalizedStatFormats.h"

/**
 * Implements one-way ANOVA (analysis of variance) statistics.
 *
 * <p> Tests for differences between two or more categories of univariate data
 * (for example, the body mass index of accountants, lawyers, doctors and
 * computer programmers).  When two categories are given, this is equivalent to
 * the {@link org.hipparchus.stat.inference.T_Test}.
 * </p><p>
 * Uses the {@link org.hipparchus.distribution.continuous.F_Distribution
 * Hipparchus F Distribution implementation} to estimate exact p-values.</p>
 * <p>This implementation is based on a description at
 * http://faculty.vassar.edu/lowry/ch13pt1.html</p>
 * <pre>
 * Abbreviations: bg = between groups, *                wg = within groups, *                ss = sum squared deviations
 * </pre>
 *
 */
class OneWay_anova 
{
private:
    /**
     * This method calls the method that actually does the calculations (except
     * P-value).
     *
     * @param category_data
     *            <code>Collection</code> of <code>std::vector<double></code> arrays each
     *            containing data for one category
     * @return computed Anova_Stats
     * @Null_Argument_Exception
     *             if <code>category_data</code> is <code>null</code>
     * @
     *             if the length of the <code>category_data</code> array is less
     *             than 2 or a contained <code>std::vector<double></code> array does not
     *             contain at least two values
     */
    Anova_Stats anova_stats(const Collection<std::vector<double>> category_data)
    {
        //Math_Utils::check_not_null(category_data);

        const Collection<Streaming_Statistics> category_data_summary_statistics =
            Array_list<>(category_data.size());

        // convert arrays to Summary_Statistics
        for (const std::vector<double> data : category_data)
        {
            const Streaming_Statistics data_summary_statistics = Streaming_Statistics();
            category_data_summary_statistics.add(data_summary_statistics);
            for (const double val : data)
            {
                data_summary_statistics.add_value(val);
            }
        }

        return anova_stats(category_data_summary_statistics, false);

    }

    /**
 * This method actually does the calculations (except P-value).
 *
 * @param category_data <code>Collection</code> of <code>std::vector<double></code>
 * arrays each containing data for one category
 * @param allow_one_element_data if true, allow computation for one category
 * only or for one data element per category
 * @return computed Anova_Stats
 * @Null_Argument_Exception if <code>category_data</code> is <code>null</code>
 * @ if <code>allow_one_element_data</code> is false and the number of
 * categories is less than 2 or a contained Summary_Statistics does not contain
 * at least two values
 */
    Anova_Stats anova_stats(const Collection<Streaming_Statistics> category_data, const bool allow_one_element_data)
    {
        //Math_Utils::check_not_null(category_data);

        if (!allow_one_element_data)
        {
            // check if we have enough categories
            if (category_data.size() < 2)
            {
                throw (Localized_Stat_Formats.TWO_OR_MORE_CATEGORIES_REQUIRED, category_data.size(), 2);
            }

            // check if each category has enough data
            for (const Streaming_Statistics array : category_data)
            {
                if (array.get_n() <= 1)
                {
                    throw (Localized_Stat_Formats.TWO_OR_MORE_VALUES_IN_CATEGORY_REQUIRED, static_cast<int>(array.get_n(), 2);
                }
            }
        }

        int dfwg = 0;
        double sswg = 0;
        double totsum = 0;
        double totsumsq = 0;
        int totnum = 0;

        for (const Streaming_Statistics data : category_data)
        {

            const double sum = data.get_sum();
            const double sumsq = data.get_sum_of_squares();
            const int& num = static_cast<int>(data.get_n();
            totnum += num;
            totsum += sum;
            totsumsq += sumsq;

            dfwg += num - 1;
            const double ss = sumsq - ((sum * sum) / num);
            sswg += ss;
        }

        const double sst = totsumsq - ((totsum * totsum) / totnum);
        const double ssbg = sst - sswg;
        const int dfbg = category_data.size() - 1;
        const double msbg = ssbg / dfbg;
        const double mswg = sswg / dfwg;
        const double F = msbg / mswg;

        return Anova_Stats(dfbg, dfwg, F);

    }

    /**
     * Convenience class to pass dfbg,dfwg,F values around within OneWay_anova.
     * No get/set methods provided.
     */
    class Anova_Stats
    {
    private:
        /** Degrees of freedom in numerator (between groups). */
        const int my_dfbg;

        /** Degrees of freedom in denominator (within groups). */
        const int my_dfwg;

        /** Statistic. */
        const double my_F;

        /**
         * Constructor
         * @param dfbg degrees of freedom in numerator (between groups)
         * @param dfwg degrees of freedom in denominator (within groups)
         * @param F statistic
         */
        Anova_Stats(const int& dfbg, int dfwg, double F) 
            :
            my_dfbg{dfbg},
            my_dfwg{dfwg},
            my_F{F}
        {
        }

        int get_dfbg() const
        {
            return my_dfbg;
        }

        int get_dfwg() const
        {
            return my_dfwg;
        }

        double get_f() const
        {
            return my_F;
        }
    };

public:
    /**
     * Computes the ANOVA F-value for a collection of <code>std::vector<double></code>
     * arrays.
     *
     * <p><strong>Preconditions</strong>: <ul>
     * <li>The category_data <code>Collection</code> must contain
     * <code>std::vector<double></code> arrays.</li>
     * <li> There must be at least two <code>std::vector<double></code> arrays in the
     * <code>category_data</code> collection and each of these arrays must
     * contain at least two values.</li></ul></p><p>
     * This implementation computes the F statistic using the definitional
     * formula<pre>
     *   F = msbg/mswg</pre>
     * where<pre>
     *  msbg = between group mean square
     *  mswg = within group mean square</pre>
     * are as defined <a href="http://faculty.vassar.edu/lowry/ch13pt1.html">
     * here</a></p>
     *
     * @param category_data <code>Collection</code> of <code>std::vector<double></code>
     * arrays each containing data for one category
     * @return Fvalue
     * @Null_Argument_Exception if <code>category_data</code> is <code>null</code>
     * @ if the length of the <code>category_data</code>
     * array is less than 2 or a contained <code>std::vector<double></code> array does not have
     * at least two values
     */
    double anova_f_value(const Collection<std::vector<double>> category_data)
    {
        return anova_stats(category_data).F;

    }

    /**
     * Computes the ANOVA P-value for a collection of <code>std::vector<double></code>
     * arrays.
     *
     * <p><strong>Preconditions</strong>: <ul>
     * <li>The category_data <code>Collection</code> must contain
     * <code>std::vector<double></code> arrays.</li>
     * <li> There must be at least two <code>std::vector<double></code> arrays in the
     * <code>category_data</code> collection and each of these arrays must
     * contain at least two values.</li></ul></p><p>
     * This implementation uses the
     * {@link org.hipparchus.distribution.continuous.F_Distribution
     * Hipparchus F Distribution implementation} to estimate the exact
     * p-value, using the formula<pre>
     *   p = 1 - cumulative_probability(F)</pre>
     * where <code>F</code> is the F value and <code>cumulative_probability</code>
     * is the Hipparchus implementation of the F distribution.</p>
     *
     * @param category_data <code>Collection</code> of <code>std::vector<double></code>
     * arrays each containing data for one category
     * @return Pvalue
     * @Null_Argument_Exception if <code>category_data</code> is <code>null</code>
     * @ if the length of the <code>category_data</code>
     * array is less than 2 or a contained <code>std::vector<double></code> array does not have
     * at least two values
     * @Math_Illegal_State_Exception if the p-value can not be computed due to a convergence error
     * @Math_Illegal_State_Exception if the maximum number of iterations is exceeded
     */
    double anova_p_value(const Collection<std::vector<double>> category_data)
    {
        const Anova_Stats a = anova_stats(category_data);
        // No try-catch or advertised exception because args are valid
        const F_Distribution fdist = F_Distribution(a.dfbg, a.dfwg);
        return 1.0 - fdist.cumulative_probability(a.F);

    }

    /**
     * Computes the ANOVA P-value for a collection of {@link Streaming_Statistics}.
     *
     * <p><strong>Preconditions</strong>: <ul>
     * <li>The category_data <code>Collection</code> must contain
     * {@link Streaming_Statistics}.</li>
     * <li> There must be at least two {@link Streaming_Statistics} in the
     * <code>category_data</code> collection and each of these statistics must
     * contain at least two values.</li></ul></p><p>
     * This implementation uses the
     * {@link org.hipparchus.distribution.continuous.F_Distribution
     * Hipparchus F Distribution implementation} to estimate the exact
     * p-value, using the formula<pre>
     *   p = 1 - cumulative_probability(F)</pre>
     * where <code>F</code> is the F value and <code>cumulative_probability</code>
     * is the Hipparchus implementation of the F distribution.</p>
     *
     * @param category_data <code>Collection</code> of {@link Streaming_Statistics}
     * each containing data for one category
     * @param allow_one_element_data if true, allow computation for one catagory
     * only or for one data element per category
     * @return Pvalue
     * @Null_Argument_Exception if <code>category_data</code> is <code>null</code>
     * @ if the length of the <code>category_data</code>
     * array is less than 2 or a contained {@link Streaming_Statistics} does not have
     * at least two values
     * @Math_Illegal_State_Exception if the p-value can not be computed due to a convergence error
     * @Math_Illegal_State_Exception if the maximum number of iterations is exceeded
     */
     double anova_p_value(const Collection<Streaming_Statistics> category_data, const bool allow_one_element_data)
     {
        const Anova_Stats a = anova_stats(category_data, allow_one_element_data);
        const F_Distribution fdist = F_Distribution(a.get_dfbg(), a.get_dfwg());
        return 1.0 - fdist.cumulative_probability(a.get_f());
     }



    /**
     * Performs an ANOVA test, evaluating the NULL hypothesis that there
     * is no difference among the means of the data categories.
     *
     * <p><strong>Preconditions</strong>: <ul>
     * <li>The category_data <code>Collection</code> must contain
     * <code>std::vector<double></code> arrays.</li>
     * <li> There must be at least two <code>std::vector<double></code> arrays in the
     * <code>category_data</code> collection and each of these arrays must
     * contain at least two values.</li>
     * <li>alpha must be strictly greater than 0 and less than or equal to 0.5.
     * </li></ul></p><p>
     * This implementation uses the
     * {@link org.hipparchus.distribution.continuous.F_Distribution
     * Hipparchus F Distribution implementation} to estimate the exact
     * p-value, using the formula<pre>
     *   p = 1 - cumulative_probability(F)</pre>
     * where <code>F</code> is the F value and <code>cumulative_probability</code>
     * is the Hipparchus implementation of the F distribution.</p>
     * <p>True is returned iff the estimated p-value is less than alpha.</p>
     *
     * @param category_data <code>Collection</code> of <code>std::vector<double></code>
     * arrays each containing data for one category
     * @param alpha significance level of the test
     * @return true if the NULL hypothesis can be rejected with
     * confidence 1 - alpha
     * @Null_Argument_Exception if <code>category_data</code> is <code>null</code>
     * @ if the length of the <code>category_data</code>
     * array is less than 2 or a contained <code>std::vector<double></code> array does not have
     * at least two values
     * @ if <code>alpha</code> is not in the range (0, 0.5]
     * @Math_Illegal_State_Exception if the p-value can not be computed due to a convergence error
     * @Math_Illegal_State_Exception if the maximum number of iterations is exceeded
     */
    bool anova_test(const Collection<std::vector<double>> category_data, const double& alpha)
    {
        if ((alpha <= 0) || (alpha > 0.5)) 
        {
            throw std::exception("Not implemented");
            //throw (Localized_Stat_Formats::OUT_OF_BOUND_SIGNIFICANCE_LEVEL, alpha, 0, 0.5);
        }
        return anova_p_value(category_data) < alpha;
    }
};