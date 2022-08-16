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

#include <cmath>
#include "BaseAbstractUnivariateIntegrator.h"
#include "../../exception/LocalizedCoreFormats.h"
#include "TrapezoidIntegrator.h"
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.FastMath;

/**
 * Implements <a href="http://mathworld.wolfram.com/SimpsonsRule.html">
 * Simpson's Rule</a> for integration of real univariate functions. For
 * reference, see <b>Introduction to Numerical Analysis</b>, ISBN 038795452X, * chapter 3.
 * <p>
 * This implementation employs the basic trapezoid rule to calculate Simpson's
 * rule.</p>
 *
 */
class Simpson_Integrator : public Base_Abstract_Univariate_Integrator
{
public:
    /** Maximal number of iterations for Simpson. */
    static constexpr int SIMPSON_MAX_ITERATIONS_COUNT = 64;

    /**
     * Build a Simpson integrator with given accuracies and iterations counts.
     * @param relative_accuracy relative accuracy of the result
     * @param absolute_accuracy absolute accuracy of the result
     * @param minimal_iteration_count minimum number of iterations
     * @param maximal_iteration_count maximum number of iterations
     * (must be less than or equal to {@link #SIMPSON_MAX_ITERATIONS_COUNT})
     * @exception  if minimal number of iterations
     * is not strictly positive
     * @exception  if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     * @exception  if maximal number of iterations
     * is greater than {@link #SIMPSON_MAX_ITERATIONS_COUNT}
     */
    Simpson_Integrator(const double& relative_accuracy, const double& absolute_accuracy, const int minimal_iteration_count, const int maximal_iteration_count)
    {
        super(relative_accuracy, absolute_accuracy, minimal_iteration_count, maximal_iteration_count);
        if (maximal_iteration_count > SIMPSON_MAX_ITERATIONS_COUNT)
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, maximal_iteration_count, SIMPSON_MAX_ITERATIONS_COUNT);
        }
    }

    /**
     * Build a Simpson integrator with given iteration counts.
     * @param minimal_iteration_count minimum number of iterations
     * @param maximal_iteration_count maximum number of iterations
     * (must be less than or equal to {@link #SIMPSON_MAX_ITERATIONS_COUNT})
     * @exception  if minimal number of iterations
     * is not strictly positive
     * @exception  if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     * @exception  if maximal number of iterations
     * is greater than {@link #SIMPSON_MAX_ITERATIONS_COUNT}
     */
    Simpson_Integrator(const int minimal_iteration_count, const int maximal_iteration_count)
    {
        super(minimal_iteration_count, maximal_iteration_count);
        if (maximal_iteration_count > SIMPSON_MAX_ITERATIONS_COUNT)
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, maximal_iteration_count, SIMPSON_MAX_ITERATIONS_COUNT);
        }
    }

    /**
     * Construct an integrator with default settings.
     * (max iteration count set to {@link #SIMPSON_MAX_ITERATIONS_COUNT})
     */
    Simpson_Integrator()
    {
        super(DEFAULT_MIN_ITERATIONS_COUNT, SIMPSON_MAX_ITERATIONS_COUNT);
    }

protected:
    /** {@inherit_doc} */
    //override
    double do_integrate()
    {

        auto qtrap = Trapezoid_Integrator();
        if (get_minimal_iteration_count() == 1)
        {
            return (4 * qtrap.stage(this, 1) - qtrap.stage(this, 0)) / 3.0;
        }

        // Simpson's rule requires at least two trapezoid stages.
        double olds{};
        double oldt = qtrap.stage(this, 0);
        while (true)
        {
            const double t = qtrap.stage(this, iterations.get_count());
            iterations.increment();
            const double s = (4 * t - oldt) / 3.0;
            if (iterations.get_count() >= get_minimal_iteration_count())
            {
                const double delta = std::abs(s - olds);
                const double r_limit = get_relative_accuracy() * (std::abs(olds) + std::abs(s)) * 0.5;
                if ((delta <= r_limit) || (delta <= get_absolute_accuracy()))
                {
                    return s;
                }
            }
            olds = s;
            oldt = t;
        }

    }

};