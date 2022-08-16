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

/**
 * Implements the <a href="http://mathworld.wolfram.com/Ridders_method.html">
 * Ridders' Method</a> for root finding of real univariate functions. For
 * reference, see C. Ridders, <i>A algorithm for computing a single root
 * of a real continuous function </i>, IEEE Transactions on Circuits and
 * Systems, 26 (1979), 979 - 980.
 * <p>
 * The function should be continuous but not necessarily smooth.</p>
 *
 */
class Ridders_Solver : Abstract_Univariate_Solver 
{
    /** Default absolute accuracy. */
    private static const double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    /**
     * Construct a solver with default accuracy (1e-6).
     */
    public Ridders_Solver() 
    {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }
    /**
     * Construct a solver.
     *
     * @param absolute_accuracy Absolute accuracy.
     */
    public Ridders_Solver(double absolute_accuracy) 
    {
        super(absolute_accuracy);
    }
    /**
     * Construct a solver.
     *
     * @param relative_accuracy Relative accuracy.
     * @param absolute_accuracy Absolute accuracy.
     */
    public Ridders_Solver(double relative_accuracy, double absolute_accuracy) 
    {
        super(relative_accuracy, absolute_accuracy);
    }

    /**
     * {@inherit_doc}
     */
    //override
    protected double do_solve()
         
        {
        double min = get_min();
        double max = get_max();
        // [x1, x2] is the bracketing interval in each iteration
        // x3 is the midpoint of [x1, x2]
        // x is the root approximation and an endpoint of the interval
        double x1 = min;
        double y1 = compute_objective_value(x1);
        double x2 = max;
        double y2 = compute_objective_value(x2);

        // check for zeros before verifying bracketing
        if (y1 == 0) 
        {
            return min;
        }
        if (y2 == 0) 
        {
            return max;
        }
        verify_bracketing(min, max);

        const double& absolute_accuracy = get_absolute_accuracy();
        const double function_value_accuracy = get_function_value_accuracy();
        const double relative_accuracy = get_relative_accuracy();

        double oldx = INFINITY;
        while (true) 
        {
            // calculate the root approximation
            const double x3 = 0.5 * (x1 + x2);
            const double y3 = compute_objective_value(x3);
            if (std::abs(y3) <= function_value_accuracy) 
            {
                return x3;
            }
            const double delta = 1 - (y1 * y2) / (y3 * y3);  // delta > 1 due to bracketing
            const double correction = (FastMath.signum(y2) * FastMath.signum(y3)) *
                                      (x3 - x1) / std::sqrt(delta);
            const double x = x3 - correction;                // correction != 0
            const double y = compute_objective_value(x);

            // check for convergence
            const double& tolerance = std::max(relative_accuracy * std::abs(x), absolute_accuracy);
            if (std::abs(x - oldx) <= tolerance) 
            {
                return x;
            }
            if (std::abs(y) <= function_value_accuracy) 
            {
                return x;
            }

            // prepare the interval for next iteration
            // Ridders' method guarantees x1 < x < x2
            if (correction > 0.0) {             // x1 < x < x3
                if (FastMath.signum(y1) + FastMath.signum(y) == 0.0) 
                {
                    x2 = x;
                    y2 = y;
                }
else 
                {
                    x1 = x;
                    x2 = x3;
                    y1 = y;
                    y2 = y3;
                }
            }
else {                            // x3 < x < x2
                if (FastMath.signum(y2) + FastMath.signum(y) == 0.0) 
                {
                    x1 = x;
                    y1 = y;
                }
else 
                {
                    x1 = x3;
                    x2 = x;
                    y1 = y3;
                    y2 = y;
                }
            }
            oldx = x;
        }
    }
}


