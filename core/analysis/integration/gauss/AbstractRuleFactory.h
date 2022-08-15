#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
//package org.hipparchus.analysis.integration.gauss;

//import java.util.Arrays;
//import java.util.Sorted_Map;
//import java.util.Tree_Map;

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Incrementor;
//import org.hipparchus.util.Pair;

#include <algorithm>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include "RuleFactory.h"


/**
 * Base class for rules that determines the integration nodes and their
 * weights.
 * Subclasses must implement the {@link #compute_rulestatic_cast<int>( compute_rule} method.
 *
 * @since 2.0
 */
class AbstractRule_Factory : public Rule_Factory 
{
private:
    /** List of points and weights, indexed by the order of the rule. */
    const std::map<int, std:pair<std::vector<double>, std::vector<double>>> points_and_weights = Tree_Map<>();

public:
    /** {@inherit_doc} */
    //override
    std::pair<std::vector<double>, std::vector<double>> get_rule(const int& number_of_points)
         
        {

        if (number_of_points <= 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, number_of_points);
        }
        if (number_of_points > 1000) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, number_of_points, 1000);
        }

        std::pair<std::vector<double>, std::vector<double>> rule;
        synchronized (points_and_weights) 
        {
            // Try to obtain the rule from the cache.
            rule = points_and_weights.get(number_of_points);

            if (rule == NULL) 
            {
                // Rule not computed yet.

                // Compute the rule.
                rule = compute_rule(number_of_points);

                // Cache it.
                points_and_weights.put(number_of_points, rule);
            }
        }

        // Return a copy.
        return std::pair<>(rule.get_first().clone(), rule.get_second().clone());

    }

protected:
    /**
     * Computes the rule for the given order.
     *
     * @param number_of_points Order of the rule to be computed.
     * @return the computed rule.
     * @ if the elements of the pair do not
     * have the same length.
     */
    virtual std::pair<std::vector<double>, std::vector<double>> compute_rule(const int& number_of_points);

    /** Computes roots of the associated orthogonal polynomials.
     * <p>
     * The roots are found using the <a href="https://en.wikipedia.org/wiki/Aberth_method">Aberth method</a>.
     * The guess points for initializing search for degree n are fixed for degrees 1 and 2 and are
     * selected from n-1 roots of rule n-1 (the two extreme roots are used, plus the n-1 intermediate
     * points between all roots).
     * </p>
     * @param n number of roots to search for
     * @param ratio_evaluator function evaluating the ratio Pₙ(x)/Pₙ'(x)
     * @return sorted array of roots
     */
    std::vector<double> find_roots(const int& n, const Univariate_Function ratio_evaluator) 
    {

        auto roots  = std::vector<double>(n);

        // set up initial guess
        if (n == 1) 
        {
            // arbitrary guess
            roots[0] = 0;
        }
        else if (n == 2) 
        {
            // arbitrary guess
            roots[0] = -1;
            roots[1] = +1;
        }
        else 
        {

            // get roots from previous rule.
            // If it has not been computed yet it will trigger a recursive call
            auto previous_points = get_rule(n - 1).get_first();

            // first guess at previous first root
            roots[0] = previous_points[0];

            // intermediate guesses between previous roots
            for (int i{ 1 }; i < n - 1; ++i) 
            {
                roots[i] = (previous_points[i - 1] + previous_points[i]) * 0.5;
            }

            // last guess at previous last root
            roots[n - 1] = previous_points[n - 2];

        }

        // use Aberth method to find all roots simultaneously
        const std::vector<double>    ratio       = std::vector<double>(n];
        const Incrementor incrementor = Incrementor(1000);
        double            tol;
        double            max_offset;
        do 
        {

            // safety check that triggers an exception if too much iterations are made
            incrementor.increment();

            // find the ratio P(xᵢ)/P'(xᵢ) for all current roots approximations
            for (int i{}; i < n; ++i) 
            {
                ratio[i] = ratio_evaluator.value(roots[i]);
            }

            // move roots approximations all at once, using Aberth method
            max_offset = 0;
            for (int i{}; i < n; ++i) 
            {
                double sum{};
                for (int j{}; j < n; ++j) 
                {
                    if (j != i) 
                    {
                        sum += 1 / (roots[i] - roots[j]);
                    }
                }
                const double offset = ratio[i] / (1 - ratio[i] * sum);
                max_offset = std::max(max_offset, std::abs(offset));
                roots[i] -= offset;
            }

            // we set tolerance to 1 ulp of the largest root
            tol = 0;
            for (const double r : roots) 
            {
                tol = std::max(tol, FastMath.ulp(r));
            }

        } while (max_offset > tol);

        // sort the roots
        Arrays.sort(roots);

        return roots;
    }

    /** Enforce symmetry of roots.
     * @param roots roots to process in place
     */
    void enforce_symmetry(std::vector<double>& roots) 
    {

        const int n = roots.size();

        // enforce symmetry
        for (int i{}; i < n / 2; ++i) 
        {
            const int idx = n - i - 1;
            const double c = (roots[i] - roots[idx]) * 0.5;
            roots[i]   = +c;
            roots[idx] = -c;
        }

        // If n is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a Find_Bugs warning.
        if (n % 2 != 0) 
        {
            roots[n / 2] = 0;
        }
    }
};