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

#include <cmath>
#include "AbstractUnivariateSolver.h"

/**
 * Implements the <em>Secant</em> method for root-finding (approximating a
 * zero of a univariate real function). The solution that is maintained is
 * not bracketed, and as such convergence is not guaranteed.
 *
 * <p>Implementation based on the following article: M. Dowell and P. Jarratt, * <em>A modified regula falsi method for computing the root of an
 * equation</em>, BIT Numerical Mathematics, volume 11, number 2, * pages 168-174, Springer, 1971.</p>
 *
 * <p>Note that since release 3.0 this class : the actual
 * <em>Secant</em> algorithm, and not a modified one. As such, the 3.0 version
 * is not backwards compatible with previous versions. To use an algorithm
 * similar to the pre-3.0 releases, use the
 * {@link Illinois_Solver <em>Illinois</em>} algorithm or the
 * {@link Pegasus_Solver <em>Pegasus</em>} algorithm.</p>
 *
 */
class Secant_Solver : public Abstract_Univariate_Solver
{
protected:
    /** Default absolute accuracy. */
    static constexpr double DEFAULT_ABSOLUTE_ACCURACY{ 1e-6 };

public:
    /** Construct a solver with default accuracy (1e-6). */
    Secant_Solver()
    {
        super(DEFAULT_ABSOLUTE_ACCURACY);
    }

    /**
     * Construct a solver.
     *
     * @param absolute_accuracy absolute accuracy
     */
    Secant_Solver(const double& absolute_accuracy)
    {
        super(absolute_accuracy);
    }

    /**
     * Construct a solver.
     *
     * @param relative_accuracy relative accuracy
     * @param absolute_accuracy absolute accuracy
     */
    Secant_Solver(const double& relative_accuracy, const double& absolute_accuracy)
    {
        super(relative_accuracy, absolute_accuracy);
    }

protected:
    /** {@inherit_doc} */
    //override
    const double do_solve()
    {
        // Get initial solution
        double x0 = get_min();
        double x1 = get_max();
        double f0 = compute_objective_value(x0);
        double f1 = compute_objective_value(x1);

        // If one of the bounds is the exact root, return it. sin_ce these are
        // not under-approximations or over-approximations, we can return them
        // regardless of the allowed solutions.
        if (f0 == 0.0)
        {
            return x0;
        }
        if (f1 == 0.0)
        {
            return x1;
        }

        // Verify bracketing of initial solution.
        verify_bracketing(x0, x1);

        // Get accuracies.
        const double ftol = get_function_value_accuracy();
        const double atol = get_absolute_accuracy();
        const double rtol = get_relative_accuracy();

        // Keep finding better approximations.
        while (true)
        {
            // Calculate the next approximation.
            const double x = x1 - ((f1 * (x1 - x0)) / (f1 - f0));
            const double fx = compute_objective_value(x);

            // If the approximation is the exact root, return it. sin_ce
            // this is not an under-approximation or an over-approximation, // we can return it regardless of the allowed solutions.
            if (fx == 0.0)
            {
                return x;
            }

            // Update the bounds with the approximation.
            x0 = x1;
            f0 = f1;
            x1 = x;
            f1 = fx;

            // If the function value of the last approximation is too small, // given the function value accuracy, then we can't get closer to
            // the root than we already are.
            if (std::abs(f1) <= ftol)
            {
                return x1;
            }

            // If the current interval is within the given accuracies, we
            // are satisfied with the current approximation.
            if (std::abs(x1 - x0) < std::max(rtol * std::abs(x1), atol))
            {
                return x1;
            }
        }
    }

};