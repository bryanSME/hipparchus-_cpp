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
//package org.hipparchus.distribution;

//import java.io.Serializable;
//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Pair;
//import org.hipparchus.util.Precision;
#include<vector>

/**
 * A generic implementation of a
 * <a href="http://en.wikipedia.org/wiki/Probability_distribution#Discrete_probability_distribution">
 * discrete probability distribution (Wikipedia)</a> over a finite sample space, * based on an enumerated list of &lt;value, probability&gt; pairs.
 * <p>
 * Input probabilities must all be non-negative, but zero values are allowed and
 * their sum does not have to equal one. Constructors will normalize input
 * probabilities to make them sum to one.
 * <p>
 * The list of &lt;value, probability&gt; pairs does not, strictly speaking, have
 * to be a function and it can contain NULL values.  The pmf created by the constructor
 * will combine probabilities of equal values and will treat NULL values as equal.
 * <p>
 * For example, if the list of pairs &lt;"dog", 0.2&gt;, &lt;null, 0.1&gt;, * &lt;"pig", 0.2&gt;, &lt;"dog", 0.1&gt;, &lt;null, 0.4&gt; is provided to the
 * constructor, the resulting pmf will assign mass of 0.5 to NULL, 0.3 to "dog"
 * and 0.2 to NULL.
 *
 * @param <T> type of the elements in the sample space.
 */
template<typename T>
class Enumerated_Distribution
{   
private:
    /**
     * List of random variable values.
     */
    const std::vector<T> my_singletons;

    /**
     * Probabilities of respective random variable values. For i = 0, ..., singletons.size() - 1, * probability[i] is the probability that a random variable following this distribution takes
     * the value singletons[i].
     */
    const std::vector<double> my_probabilities;

public:
    /**
     * Create an enumerated distribution using the given probability mass function
     * enumeration.
     *
     * @param pmf probability mass function enumerated as a list of &lt;T, probability&gt;
     * pairs.
     * @ of weights includes negative, NaN or infinite values or only 0's
     */
    Enumerated_Distribution(const std::vector<std::pair<T, double>>& pmf)     
    {

        my_singletons = std::vector<>(pmf.size());
        auto probs = std::vector<double>(pmf.size());

        for (int i{}; i < pmf.size(); i++) 
        {
            const auto& sample = pmf.at(i);
            my_singletons.add(sample.first);
            probs[i] = sample.second;
        }
        my_probabilities = check_and_normalize(probs);
    }

    /**
     * For a random variable {@code X} whose values are distributed according to
     * this distribution, this method returns {@code P(X = x)}. In other words, * this method represents the probability mass function (PMF) for the
     * distribution.
     * <p>
     * Note that if {@code x1} and {@code x2} satisfy {@code x1.equals(x2)}, * or both are NULL, then {@code probability(x1) = probability(x2)}.
     *
     * @param x the point at which the PMF is evaluated
     * @return the value of the probability mass function at {@code x}
     */
    double probability(const T& x) 
    {
        double probability{};

        for (int i{}; i < probabilities.size(); i++) 
        {
            if ((x == NULL && my_singletons.get(i) == NULL) ||
                (x != NULL && x.equals(my_singletons.get(i)))) 
                {
                probability += probabilities[i];
            }
        }

        return probability;
    }

    /**
     * Return the probability mass function as a list of (value, probability) pairs.
     * <p>
     * Note that if duplicate and / or NULL values were provided to the constructor
     * when creating this Enumerated_Distribution, the returned list will contain these
     * values.  If duplicates values exist, what is returned will not represent
     * a pmf (i.e., it is up to the caller to consolidate duplicate mass points).
     *
     * @return the probability mass function.
     */
    std::vector<std::pair<T, double>> get_pmf() 
    {
        auto samples = std::vector<>(probabilities.size());

        for (int i{}; i < probabilities.size(); i++) 
        {
            samples.push_back(std::pair<>(my_singletons.get(i), probabilities[i]));
        }

        return samples;
    }

    /**
     * Checks to make sure that weights is neither NULL nor empty and contains only non-negative, finite, * non-NaN values and if necessary normalizes it to sum to 1.
     *
     * @param weights input array to be used as the basis for the values of a PMF
     * @return a possibly rescaled copy of the array that sums to 1 and contains only valid probability values
     * @ of weights is NULL or empty or includes negative, NaN or
     *         infinite values or only 0's
     */
    static std::vector<double> check_and_normalize(const std::vector<double>& weights)
    {
        if (weights == NULL || weights.size() == 0)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::ARRAY_ZERO_LENGTH_OR_NULL_NOT_ALLOWED);
        }
        const int len = weights.size();
        double sum_wt{};
        bool pos_wt = false;
        for (int i{}; i < len; i++)
        {
            if (weights[i] < 0)
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, weights[i], 0);
            }
            if (weights[i] > 0)
            {
                pos_wt = true;
            }
            if (std::isnan(weights[i]))
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::NAN_ELEMENT_AT_INDEX, i);
            }
            if (std::isinf(weights[i]))
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::INFINITE_ARRAY_ELEMENT, weights[i], i);
            }
            sum_wt += weights[i];
        }
        if (!pos_wt)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::WEIGHT_AT_LEAST_ONE_NON_ZERO);
        }
        std::vector<double> norm_wt;
        if (Precision.equals(sum_wt, 1, 10)) // allow small error (10 ulps)
        {
            norm_wt = weights;
        }
        else
        {
            norm_wt = std::vector<double>(len);
            for (int i{}; i < len; i++)
            {
                norm_wt[i] = weights[i] / sum_wt;
            }
        }
        return norm_wt;
    };

};