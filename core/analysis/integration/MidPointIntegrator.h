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
//package org.hipparchus.analysis.integration;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.FastMath;
#include "BaseAbstractUnivariateIntegrator.h"
#include <cmath>
#include "../../exception/LocalizedCoreFormats.h"

/**
 * Implements the <a href="http://en.wikipedia.org/wiki/Midpoint_method">
 * Midpoint Rule</a> for integration of real univariate functions. For
 * reference, see <b>Numerical Mathematics</b>, ISBN 0387989595, * chapter 9.2.
 * <p>
 * The function should be integrable.</p>
 *
 */
class Mid_pointIntegrator : public Base_Abstract_Univariate_Integrator 
{
public:
    /** Maximum number of iterations for midpoint. */
    static constexpr int MIDPOINT_MAX_ITERATIONS_COUNT{ 64 };

    /**
     * Build a midpoint integrator with given accuracies and iterations counts.
     * @param relative_accuracy relative accuracy of the result
     * @param absolute_accuracy absolute accuracy of the result
     * @param minimal_iteration_count minimum number of iterations
     * @param maximal_iteration_count maximum number of iterations
     * (must be less than or equal to {@link #MIDPOINT_MAX_ITERATIONS_COUNT}
     * @exception  if minimal number of iterations
     * is not strictly positive
     * @exception  if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     * @exception  if maximal number of iterations
     * is greater than {@link #MIDPOINT_MAX_ITERATIONS_COUNT}
     */
    Mid_pointIntegrator(const double& relative_accuracy, const double& absolute_accuracy, const int minimal_iteration_count, const int maximal_iteration_count)
    {
        super(relative_accuracy, absolute_accuracy, minimal_iteration_count, maximal_iteration_count);
        if (maximal_iteration_count > MIDPOINT_MAX_ITERATIONS_COUNT) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, maximal_iteration_count, MIDPOINT_MAX_ITERATIONS_COUNT);
        }
    }

    /**
     * Build a midpoint integrator with given iteration counts.
     * @param minimal_iteration_count minimum number of iterations
     * @param maximal_iteration_count maximum number of iterations
     * (must be less than or equal to {@link #MIDPOINT_MAX_ITERATIONS_COUNT}
     * @exception  if minimal number of iterations
     * is not strictly positive
     * @exception  if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     * @exception  if maximal number of iterations
     * is greater than {@link #MIDPOINT_MAX_ITERATIONS_COUNT}
     */
    Mid_pointIntegrator(const int minimal_iteration_count, const int maximal_iteration_count)
    {
        super(minimal_iteration_count, maximal_iteration_count);
        if (maximal_iteration_count > MIDPOINT_MAX_ITERATIONS_COUNT) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, maximal_iteration_count, MIDPOINT_MAX_ITERATIONS_COUNT);
        }
    }

    /**
     * Construct a midpoint integrator with default settings.
     * (max iteration count set to {@link #MIDPOINT_MAX_ITERATIONS_COUNT})
     */
    Mid_pointIntegrator() 
    {
        super(DEFAULT_MIN_ITERATIONS_COUNT, MIDPOINT_MAX_ITERATIONS_COUNT);
    }

private:
    /**
     * Compute the n-th stage integral of midpoint rule.
     * This function should only be called by API <code>integrate()</code> in the //package.
     * To save time it does not verify arguments - caller does.
     * <p>
     * The interval is divided equally into 2^n sections rather than an
     * arbitrary m sections because this configuration can best utilize the
     * already computed values.</p>
     *
     * @param n the stage of 1/2 refinement. Must be larger than 0.
     * @param previous_stage_result Result from the previous call to the
     * {@code stage} method.
     * @param min Lower bound of the integration interval.
     * @param diff_max_min Difference between the lower bound and upper bound
     * of the integration interval.
     * @return the value of n-th stage integral
     * @Math_Illegal_State_Exception if the maximal number of evaluations
     * is exceeded.
     */
    double stage(const int& n, double previous_stage_result, const double& min,  double diff_max_min)
    {

        // number of points in this stage
        const long np = 1L << (n - 1);
        double sum{};

        // spacing between adjacent points
        const double spacing = diff_max_min / np;

        // the first point
        double x = min + 0.5 * spacing;
        for (long i = 0; i < np; i++) 
        {
            sum += compute_objective_value(x);
            x += spacing;
        }
        // add the sum to previously calculated result
        return 0.5 * (previous_stage_result + sum * spacing);
    }

protected:

    /** {@inherit_doc} */
    //override
    double do_integrate()
    {
        const double min = get_min();
        const double diff = get_max() - min;
        const double mid_point = min + 0.5 * diff;

        double oldt = diff * compute_objective_value(mid_point);

        while (true) 
        {
            iterations.increment();
            const int i = iterations.get_count();
            const double t = stage(i, oldt, min, diff);
            if (i >= get_minimal_iteration_count()) 
            {
                const double delta = std::abs(t - oldt);
                const double r_limit =
                        get_relative_accuracy() * (std::abs(oldt) + std::abs(t)) * 0.5;
                if ((delta <= r_limit) || (delta <= get_absolute_accuracy())) 
                {
                    return t;
                }
            }
            oldt = t;
        }

    }

};