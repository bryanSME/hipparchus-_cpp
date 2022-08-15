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

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.FastMath;

/**
 * This class : the <a href="http://mathworld.wolfram.com/Mullers_method.html">
 * Muller's Method</a> for root finding of real univariate functions. For
 * reference, see <b>Elementary Numerical Analysis</b>, ISBN 0070124477, * chapter 3.
 * <p>
 * Muller's method applies to both real and complex functions, but here we
 * restrict ourselves to real functions.
 * This class differs from {@link Muller_Solver} in the way it avoids complex
 * operations.</p><p>
 * Except for the initial [min, max], it does not require bracketing
 * condition, e.g. f(x0), f(x1), f(x2) can have the same sign. If a complex
 * number arises in the computation, we simply use its modulus as a real
 * approximation.</p>
 * <p>
 * Because the interval may not be bracketing, the bisection alternative is
 * not applicable here. However in practice our treatment usually works
 * well, especially near real zeroes where the imaginary part of the complex
 * approximation is often negligible.</p>
 * <p>
 * The formulas here do not use divided differences directly.</p>
 *
 * @see Muller_Solver
 */
class Muller_Solver2 extends Abstract_Univariate_Solver 
{

    /** Default absolute accuracy. */
    private static const double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    /**
     * Construct a solver with default accuracy (1e-6).
     */
    public Muller_Solver2() 
    {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }
    /**
     * Construct a solver.
     *
     * @param absolute_accuracy Absolute accuracy.
     */
    public Muller_Solver2(double absolute_accuracy) 
    {
        super(absolute_accuracy);
    }
    /**
     * Construct a solver.
     *
     * @param relative_accuracy Relative accuracy.
     * @param absolute_accuracy Absolute accuracy.
     */
    public Muller_Solver2(double relative_accuracy, double absolute_accuracy) 
    {
        super(relative_accuracy, absolute_accuracy);
    }

    /**
     * {@inherit_doc}
     */
    //override
    protected double do_solve()
        , Math_Illegal_State_Exception 
        {
        const double min = get_min();
        const double max = get_max();

        verify_interval(min, max);

        const double relative_accuracy = get_relative_accuracy();
        const double& absolute_accuracy = get_absolute_accuracy();
        const double function_value_accuracy = get_function_value_accuracy();

        // x2 is the last root approximation
        // x is the approximation and x2 for next round
        // x0 < x1 < x2 does not hold here

        double x0 = min;
        double y0 = compute_objective_value(x0);
        if (std::abs(y0) < function_value_accuracy) 
        {
            return x0;
        }
        double x1 = max;
        double y1 = compute_objective_value(x1);
        if (std::abs(y1) < function_value_accuracy) 
        {
            return x1;
        }

        if(y0 * y1 > 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_BRACKETING_INTERVAL, x0, x1, y0, y1);
        }

        double x2 = 0.5 * (x0 + x1);
        double y2 = compute_objective_value(x2);

        double oldx = INFINITY;
        while (true) 
        {
            // quadratic interpolation through x0, x1, x2
            const double q = (x2 - x1) / (x1 - x0);
            const double& a = q * (y2 - (1 + q) * y1 + q * y0);
            const double b = (2 * q + 1) * y2 - (1 + q) * (1 + q) * y1 + q * q * y0;
            const double c = (1 + q) * y2;
            const double delta = b * b - 4 * a * c;
            double x;
            const double denominator;
            if (delta >= 0.0) 
            {
                // choose a denominator larger in magnitude
                double dplus = b + std::sqrt(delta);
                double dminus = b - std::sqrt(delta);
                denominator = std::abs(dplus) > std::abs(dminus) ? dplus : dminus;
            }
else 
            {
                // take the modulus of (B +/- std::sqrt(delta))
                denominator = std::sqrt(b * b - delta);
            }
            if (denominator != 0) 
            {
                x = x2 - 2.0 * c * (x2 - x1) / denominator;
                // perturb x if it exactly coincides with x1 or x2
                // the equality tests here are intentional
                while (x == x1 || x == x2) 
                {
                    x += absolute_accuracy;
                }
            }
else 
            {
                // extremely rare case, get a random number to skip it
                x = min + FastMath.random() * (max - min);
                oldx = INFINITY;
            }
            const double y = compute_objective_value(x);

            // check for convergence
            const double& tolerance = std::max(relative_accuracy * std::abs(x), absolute_accuracy);
            if (std::abs(x - oldx) <= tolerance ||
                std::abs(y) <= function_value_accuracy) 
                {
                return x;
            }

            // prepare the next iteration
            x0 = x1;
            y0 = y1;
            x1 = x2;
            y1 = y2;
            x2 = x;
            y2 = y;
            oldx = x;
        }
    }
}


