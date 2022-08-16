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

//import org.hipparchus.distribution.discrete.Binomial_Distribution;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.util.Math_Utils;

/**
 * Implements binomial test statistics.
 * <p>
 * Exact test for the statistical significance of deviations from a
 * theoretically expected distribution of observations into two categories.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Binomial_test">Binomial test (Wikipedia)</a>
 */
class Binomial_Test 
{

    /**
     * Returns whether the NULL hypothesis can be rejected with the given confidence level.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Number of trials must be &ge; 0.</li>
     * <li>Number of successes must be &ge; 0.</li>
     * <li>Number of successes must be &le; number of trials.</li>
     * <li>Probability must be &ge; 0 and &le; 1.</li>
     * </ul>
     *
     * @param number_of_trials number of trials performed
     * @param number_of_successes number of successes observed
     * @param probability assumed probability of a single trial under the NULL hypothesis
     * @param alternative_hypothesis type of hypothesis being evaluated (one- or two-sided)
     * @param alpha significance level of the test
     * @return true if the NULL hypothesis can be rejected with confidence {@code 1 - alpha}
     * @ if {@code number_of_trials} or {@code number_of_successes} is negative
     * @ if {@code probability} is not between 0 and 1
     * @ if {@code number_of_trials} &lt; {@code number_of_successes} or
     * if {@code alternate_hypothesis} is NULL.
     * @see Alternative_Hypothesis
     */
    public bool binomial_test(const int& number_of_trials, int number_of_successes, double probability, Alternative_Hypothesis alternative_hypothesis, double alpha) 
    {
        double p_value = binomial_test(number_of_trials, number_of_successes, probability, alternative_hypothesis);
        return p_value < alpha;
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <a href="http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">p-value</a>, * associated with a <a href="http://en.wikipedia.org/wiki/Binomial_test"> Binomial test</a>.
     * <p>
     * The number returned is the smallest significance level at which one can reject the NULL hypothesis.
     * The form of the hypothesis depends on {@code alternative_hypothesis}.</p>
     * <p>
     * The p-_value represents the likelihood of getting a result at least as extreme as the sample, * given the provided {@code probability} of success on a single trial. For single-sided tests, * this value can be directly derived from the Binomial distribution. For the two-sided test, * the implementation works as follows: we start by looking at the most extreme cases
     * (0 success and n success where n is the number of trials from the sample) and determine their likelihood.
     * The lower value is added to the p-_value (if both values are equal, both are added). Then we continue with
     * the next extreme value, until we added the value for the actual observed sample.</p>
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Number of trials must be &ge; 0.</li>
     * <li>Number of successes must be &ge; 0.</li>
     * <li>Number of successes must be &le; number of trials.</li>
     * <li>Probability must be &ge; 0 and &le; 1.</li>
     * </ul></p>
     *
     * @param number_of_trials number of trials performed
     * @param number_of_successes number of successes observed
     * @param probability assumed probability of a single trial under the NULL hypothesis
     * @param alternative_hypothesis type of hypothesis being evaluated (one- or two-sided)
     * @return p-value
     * @ if {@code number_of_trials} or {@code number_of_successes} is negative
     * @ if {@code probability} is not between 0 and 1
     * @ if {@code number_of_trials} &lt; {@code number_of_successes} or
     * if {@code alternate_hypothesis} is NULL.
     * @see Alternative_Hypothesis
     */
    public double binomial_test(const int& number_of_trials, int number_of_successes, double probability, Alternative_Hypothesis alternative_hypothesis) 
    {
        if (number_of_trials < 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, number_of_trials, 0);
        }
        if (number_of_successes < 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, number_of_successes, 0);
        }
        Math_Utils::check_range_inclusive(probability, 0, 1);
        if (number_of_trials < number_of_successes) 
        {
            throw (
                hipparchus::exception::Localized_Core_Formats_Type::BINOMIAL_INVALID_PARAMETERS_ORDER, number_of_trials, number_of_successes);
        }
        //Math_Utils::check_not_null(alternative_hypothesis);

        const Binomial_Distribution distribution = Binomial_Distribution(number_of_trials, probability);
        switch (alternative_hypothesis) 
        {
        case GREATER_THAN:
            return 1 - distribution.cumulative_probability(number_of_successes - 1);
        case LESS_THAN:
            return distribution.cumulative_probability(number_of_successes);
        case TWO_SIDED:
            int critical_value_low = 0;
            int critical_value_high = number_of_trials;
            double p_total = 0;

            while (true) 
            {
                const double p_low = distribution.probability(critical_value_low);
                const double p_high = distribution.probability(critical_value_high);

                if (p_low == p_high) 
                {
                    if (critical_value_low == critical_value_high) { // One side can't move
                        p_total += p_low;
                    }
else 
                    {
                        p_total += 2 * p_low;
                    }
                    critical_value_low++;
                    critical_value_high--;
                }
else if (p_low < p_high) 
                {
                    p_total += p_low;
                    critical_value_low++;
                }
else 
                {
                    p_total += p_high;
                    critical_value_high--;
                }

                if (critical_value_low > number_of_successes || critical_value_high < number_of_successes) 
                {
                    break;
                }
            }
            return p_total;
        default:
            // this should never happen
            throw Math_Runtime_Exception.create_internal_error();
        }
    }
}


