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

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.FastMath;
#include "AbstractUnivariateSolver.h"

/**
 * This class : the <a href="http://mathworld.wolfram.com/Mullers_method.html">
 * Muller's Method</a> for root finding of real univariate functions. For
 * reference, see <b>Elementary Numerical Analysis</b>, ISBN 0070124477, * chapter 3.
 * <p>
 * Muller's method applies to both real and complex functions, but here we
 * restrict ourselves to real functions.
 * This class differs from {@link Muller_Solver} in the way it avoids complex
 * operations.</p><p>
 * Muller's original method would have function evaluation at complex point.
 * sin_ce our f(x) is real, we have to find ways to avoid that. Bracketing
 * condition is one way to go: by requiring bracketing in every iteration, * the newly computed approximation is guaranteed to be real.</p>
 * <p>
 * Normally Muller's method converges quadratically in the vicinity of a
 * zero, however it may be very slow in regions far away from zeros. For
 * example, f(x) = exp(x) - 1, min = -50, max = 100. In such case we use
 * bisection as a safety backup if it performs very poorly.</p>
 * <p>
 * The formulas here use divided differences directly.</p>
 *
 * @see Muller_Solver2
 */
class Muller_Solver : public Abstract_Univariate_Solver 
{
private:
    /** Default absolute accuracy. */
    static constexpr double DEFAULT_ABSOLUTE_ACCURACY{ 1e-6 };

    template<typename T>
    static int sgn(const T& val)
    {
        return (T(0) < val) - (val < T(0));
    }

    /**
     * Find a real root in the given interval.
     *
     * @param min Lower bound for the interval.
     * @param max Upper bound for the interval.
     * @param f_min function value at the lower bound.
     * @param f_max function value at the upper bound.
     * @return the point at which the function value is zero.
     * @Math_Illegal_State_Exception if the allowed number of calls to
     * the function to be solved has been exhausted.
     */
    double solve(const double& min, const double& max, const double& f_min, const double& f_max)
    {
        const double relative_accuracy = get_relative_accuracy();
        const double absolute_accuracy = get_absolute_accuracy();
        const double function_value_accuracy = get_function_value_accuracy();

        // [x0, x2] is the bracketing interval in each iteration
        // x1 is the last approximation and an interpolation point in (x0, x2)
        // x is the root approximation and x1 for next round
        // d01, d12, d012 are divided differences

        double x0 = min;
        double y0 = f_min;
        double x2 = max;
        double y2 = f_max;
        double x1 = 0.5 * (x0 + x2);
        double y1 = compute_objective_value(x1);

        double oldx = INFINITY;
        while (true)
        {
            // Muller's method employs quadratic interpolation through
            // x0, x1, x2 and x is the zero of the interpolating parabola.
            // Due to bracketing condition, this parabola must have two
            // real roots and we choose one in [x0, x2] to be x.
            const double d01 = (y1 - y0) / (x1 - x0);
            const double d12 = (y2 - y1) / (x2 - x1);
            const double d012 = (d12 - d01) / (x2 - x0);
            const double c1 = d01 + (x1 - x0) * d012;
            const double delta = c1 * c1 - 4 * y1 * d012;
            const double xplus = x1 + (-2.0 * y1) / (c1 + std::sqrt(delta));
            const double xminus = x1 + (-2.0 * y1) / (c1 - std::sqrt(delta));
            // xplus and xminus are two roots of parabola and at least
            // one of them should lie in (x0, x2)
            const double x = is_sequence(x0, xplus, x2) ? xplus : xminus;
            const double y = compute_objective_value(x);

            // check for convergence
            const double& tolerance = std::max(relative_accuracy * std::abs(x), absolute_accuracy);
            if (std::abs(x - oldx) <= tolerance ||
                std::abs(y) <= function_value_accuracy)
            {
                return x;
            }

            // Bisect if convergence is too slow. Bisection would waste
            // our calculation of x, hopefully it won't happen often.
            // the real number equality test x == x1 is intentional and
            // completes the proximity tests above it
            bool bisect = (x < x1 && (x1 - x0) > 0.95 * (x2 - x0)) ||
                (x > x1 && (x2 - x1) > 0.95 * (x2 - x0)) ||
                (x == x1);
            // prepare the bracketing interval for next iteration
            if (!bisect)
            {
                x0 = x < x1 ? x0 : x1;
                y0 = x < x1 ? y0 : y1;
                x2 = x > x1 ? x2 : x1;
                y2 = x > x1 ? y2 : y1;
                x1 = x; y1 = y;
                oldx = x;
            }
            else
            {
                const double xm = 0.5 * (x0 + x2);
                const double ym = compute_objective_value(xm);
                if (sgn(y0) + sgn(ym) == 0.0)
                {
                    x2 = xm;
                    y2 = ym;
                }
                else
                {
                    x0 = xm;
                    y0 = ym;
                }
                x1 = 0.5 * (x0 + x2);
                y1 = compute_objective_value(x1);
                oldx = INFINITY;
            }
        }
    }

public:
    /**
     * Construct a solver with default accuracy (1e-6).
     */
    Muller_Solver() 
    {
        Muller_Solver(DEFAULT_ABSOLUTE_ACCURACY);
    }
    /**
     * Construct a solver.
     *
     * @param absolute_accuracy Absolute accuracy.
     */
    Muller_Solver(const double& absolute_accuracy) 
    {
        super(absolute_accuracy);
    }
    /**
     * Construct a solver.
     *
     * @param relative_accuracy Relative accuracy.
     * @param absolute_accuracy Absolute accuracy.
     */
    Muller_Solver(const double& relative_accuracy, const double& absolute_accuracy) 
    {
        super(relative_accuracy, absolute_accuracy);
    }

protected:
    /**
     * {@inherit_doc}
     */
    //override
    double do_solve()
    {
        const double min = get_min();
        const double max = get_max();
        const double initial = get_start_value();

        const double function_value_accuracy = get_function_value_accuracy();

        verify_sequence(min, initial, max);

        // check for zeros before verifying bracketing
        const double f_min = compute_objective_value(min);
        if (std::abs(f_min) < function_value_accuracy) 
        {
            return min;
        }
        const double f_max = compute_objective_value(max);
        if (std::abs(f_max) < function_value_accuracy) 
        {
            return max;
        }
        const double f_initial = compute_objective_value(initial);
        if (std::abs(f_initial) <  function_value_accuracy) 
        {
            return initial;
        }

        verify_bracketing(min, max);

        if (is_bracketing(min, initial)) 
        {
            return solve(min, initial, f_min, f_initial);
        }
        return solve(initial, max, f_initial, f_max);
    }  
};