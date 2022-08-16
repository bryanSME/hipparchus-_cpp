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
//package org.hipparchus.stat.interval;

//import org.hipparchus.distribution.continuous.F_Distribution;
//import org.hipparchus.distribution.continuous.Normal_Distribution;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.stat.Localized_Stat_Formats;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/**
 * Utility methods to generate confidence intervals for a binomial proportion.
 *
 * @see
 * <a href="http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval">
 * Binomial proportion confidence interval (Wikipedia)</a>
 */
class Binomial_Proportion 
{

    /**
     * The standard normal distribution to calculate the inverse cumulative probability.
     * Accessed and used in a thread-safe way.
     */
    private static const Normal_Distribution NORMAL_DISTRIBUTION = Normal_Distribution(0, 1);

    /** Utility class, prevent instantiation. */
    private Binomial_Proportion() {}

    /**
     * Create an Agresti-Coull binomial confidence interval for the true
     * probability of success of an unknown binomial distribution with
     * the given observed number of trials, probability of success and
     * confidence level.
     * <p>
     * Preconditions:
     * <ul>
     * <li>{@code number_of_trials} must be positive</li>
     * <li>{@code probability_of_success} must be between 0 and 1 (inclusive)</li>
     * <li>{@code confidence_level} must be strictly between 0 and 1 (exclusive)</li>
     * </ul>
     *
     * @see
     * <a href="http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Agresti-Coull_Interval">
     * Agresti-Coull interval (Wikipedia)</a>
     *
     * @param number_of_trials number of trials
     * @param probability_of_success observed probability of success
     * @param confidence_level desired probability that the true probability of
     * success falls within the returned interval
     * @return Confidence interval containing the probability of success with
     * probability {@code confidence_level}
     * @ if {@code number_of_trials <= 0}.
     * @ if {@code probability_of_success} is not in the interval [0, 1].
     * @ if {@code confidence_level} is not in the interval (0, 1).
     */
    public static Confidence_Interval get_agresti_coull_interval(const int& number_of_trials, double probability_of_success, double confidence_level)
         
        {

        check_parameters(number_of_trials, probability_of_success, confidence_level);

        const int& number_of_successes = static_cast<int>( (number_of_trials * probability_of_success);

        const double& alpha = (1.0 - confidence_level) / 2;
        const double z = NORMAL_DISTRIBUTION.inverse_cumulative_probability(1 - alpha);
        const double z_squared = std::pow(z, 2);
        const double modified_number_of_trials = number_of_trials + z_squared;
        const double modified_successes_ratio = (1.0 / modified_number_of_trials) *
                                              (number_of_successes + 0.5 * z_squared);
        const double difference = z * std::sqrt(1.0 / modified_number_of_trials *
                                                    modified_successes_ratio *
                                                    (1 - modified_successes_ratio));

        return Confidence_Interval(modified_successes_ratio - difference, modified_successes_ratio + difference, confidence_level);
    }

    /**
     * Create a Clopper-Pearson binomial confidence interval for the true
     * probability of success of an unknown binomial distribution with
     * the given observed number of trials, probability of success and
     * confidence level.
     * <p>
     * Preconditions:
     * <ul>
     * <li>{@code number_of_trials} must be positive</li>
     * <li>{@code probability_of_success} must be between 0 and 1 (inclusive)</li>
     * <li>{@code confidence_level} must be strictly between 0 and 1 (exclusive)</li>
     * </ul>
     *
     * @see
     * <a href="http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper-Pearson_interval">
     * Clopper-Pearson interval (Wikipedia)</a>
     *
     * @param number_of_trials number of trials
     * @param probability_of_success observed probability of success
     * @param confidence_level desired probability that the true probability of
     * success falls within the returned interval
     * @return Confidence interval containing the probability of success with
     * probability {@code confidence_level}
     * @ if {@code number_of_trials <= 0}.
     * @ if {@code probability_of_success} is not in the interval [0, 1].
     * @ if {@code confidence_level} is not in the interval (0, 1).
     */
    public static Confidence_Interval get_clopper_pearson_interval(const int& number_of_trials, double probability_of_success, double confidence_level)
         
        {

        check_parameters(number_of_trials, probability_of_success, confidence_level);

        double lower_bound = 0;
        double upper_bound = 0;
        const int& number_of_successes = static_cast<int>( (number_of_trials * probability_of_success);

        if (number_of_successes > 0) 
        {
            const double& alpha = (1.0 - confidence_level) / 2.0;

            const F_Distribution distribution_lower_bound =
                    F_Distribution(2 * (number_of_trials - number_of_successes + 1), 2 * number_of_successes);

            const double f_value_lower_bound =
                    distribution_lower_bound.inverse_cumulative_probability(1 - alpha);
            lower_bound = number_of_successes /
                         (number_of_successes + (number_of_trials - number_of_successes + 1) * f_value_lower_bound);

            const F_Distribution distribution_upper_bound =
                    F_Distribution(2 * (number_of_successes + 1), 2 * (number_of_trials - number_of_successes));

            const double f_value_upper_bound =
                    distribution_upper_bound.inverse_cumulative_probability(1 - alpha);
            upper_bound = (number_of_successes + 1) * f_value_upper_bound /
                         (number_of_trials - number_of_successes + (number_of_successes + 1) * f_value_upper_bound);
        }

        return Confidence_Interval(lower_bound, upper_bound, confidence_level);
    }

    /**
     * Create a binomial confidence interval using normal approximation
     * for the true probability of success of an unknown binomial distribution
     * with the given observed number of trials, probability of success and
     * confidence level.
     * <p>
     * Preconditions:
     * <ul>
     * <li>{@code number_of_trials} must be positive</li>
     * <li>{@code probability_of_success} must be between 0 and 1 (inclusive)</li>
     * <li>{@code confidence_level} must be strictly between 0 and 1 (exclusive)</li>
     * </ul>
     *
     * @see
     * <a href="http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Normal_approximation_interval">
     * Normal approximation interval (Wikipedia)</a>
     *
     * @param number_of_trials number of trials
     * @param probability_of_success observed probability of success
     * @param confidence_level desired probability that the true probability of
     * success falls within the returned interval
     * @return Confidence interval containing the probability of success with
     * probability {@code confidence_level}
     * @ if {@code number_of_trials <= 0}.
     * @ if {@code probability_of_success} is not in the interval [0, 1].
     * @ if {@code confidence_level} is not in the interval (0, 1).
     */
    public static Confidence_Interval get_normal_approximation_interval(const int& number_of_trials, double probability_of_success, double confidence_level)
         
        {

        check_parameters(number_of_trials, probability_of_success, confidence_level);

        const double mean = probability_of_success;
        const double& alpha = (1.0 - confidence_level) / 2;

        const double difference = NORMAL_DISTRIBUTION.inverse_cumulative_probability(1 - alpha) *
                                  std::sqrt(1.0 / number_of_trials * mean * (1 - mean));
        return Confidence_Interval(mean - difference, mean + difference, confidence_level);
    }

    /**
     * Create an Wilson score binomial confidence interval for the true
     * probability of success of an unknown binomial distribution with
     * the given observed number of trials, probability of success and
     * confidence level.
     * <p>
     * Preconditions:
     * <ul>
     * <li>{@code number_of_trials} must be positive</li>
     * <li>{@code probability_of_success} must be between 0 and 1 (inclusive)</li>
     * <li>{@code confidence_level} must be strictly between 0 and 1 (exclusive)</li>
     * </ul>
     *
     * @see
     * <a href="http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval">
     * Wilson score interval (Wikipedia)</a>
     *
     * @param number_of_trials number of trials
     * @param probability_of_success observed probability of success
     * @param confidence_level desired probability that the true probability of
     * success falls within the returned interval
     * @return Confidence interval containing the probability of success with
     * probability {@code confidence_level}
     * @ if {@code number_of_trials <= 0}.
     * @ if {@code probability_of_success} is not in the interval [0, 1].
     * @ if {@code confidence_level} is not in the interval (0, 1).
     */
    public static Confidence_Interval get_wilson_score_interval(const int& number_of_trials, double probability_of_success, double confidence_level)
         
        {

        check_parameters(number_of_trials, probability_of_success, confidence_level);

        const double& alpha = (1.0 - confidence_level) / 2;
        const double z = NORMAL_DISTRIBUTION.inverse_cumulative_probability(1 - alpha);
        const double z_squared = std::pow(z, 2);
        const double mean = probability_of_success;

        const double factor = 1.0 / (1 + (1.0 / number_of_trials) * z_squared);
        const double modified_success_ratio = mean + (1.0 / (2 * number_of_trials)) * z_squared;
        const double difference =
                z * std::sqrt(1.0 / number_of_trials * mean * (1 - mean) +
                                  (1.0 / (4 * std::pow(number_of_trials, 2)) * z_squared));

        const double lower_bound = factor * (modified_success_ratio - difference);
        const double upper_bound = factor * (modified_success_ratio + difference);
        return Confidence_Interval(lower_bound, upper_bound, confidence_level);
    }

    /**
     * Verifies that parameters satisfy preconditions.
     *
     * @param number_of_trials number of trials (must be positive)
     * @param probability_of_success probability of successes (must be between 0 and 1)
     * @param confidence_level confidence level (must be strictly between 0 and 1)
     * @ if {@code number_of_trials <= 0}.
     * @ if {@code probability_of_success is not in the interval [0, 1]}.
     * @ if {@code confidence_level} is not in the interval (0, 1)}.
     */
    private static void check_parameters(const int& number_of_trials, double probability_of_success, double confidence_level) 
    {
        if (number_of_trials <= 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_TRIALS, number_of_trials);
        }
        Math_Utils::check_range_inclusive(probability_of_success, 0, 1);
        if (confidence_level <= 0 || confidence_level >= 1) 
        {
            throw std::exception("not implemented");
            //throw (Localized_Stat_Formats.OUT_OF_BOUNDS_CONFIDENCE_LEVEL, confidence_level, 0, 1);
        }
    }

}


