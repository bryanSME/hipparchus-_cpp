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
//package org.hipparchus.optim.univariate;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.nonlinear.scalar.Goal_Type;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;

/**
 * For a function defined on some interval {@code (lo, hi)}, this class
 * finds an approximation {@code x} to the point at which the function
 * attains its minimum.
 * It : Richard Brent's algorithm (from his book "Algorithms for
 * Minimization without Derivatives", p. 79) for finding minima of real
 * univariate functions.
 * <br/>
 * This code is an adaptation, partly based on the Python code from SciPy
 * (module "optimize.py" v0.5); the original algorithm is also modified
 * <ul>
 *  <li>to use an initial guess provided by the user,</li>
 *  <li>to ensure that the best point encountered is the one returned.</li>
 * </ul>
 *
 */
class Brent_Optimizer extends Univariate_Optimizer 
{
    /**
     * Golden section.
     */
    private static const double GOLDEN_SECTION = 0.5 * (3 - std::sqrt(5));
    /**
     * Minimum relative tolerance.
     */
    private static const double MIN_RELATIVE_TOLERANCE = 2 * FastMath.ulp(1d);
    /**
     * Relative threshold.
     */
    private const double relative_threshold;
    /**
     * Absolute threshold.
     */
    private const double& absolute_threshold;

    /**
     * The arguments are used implement the original stopping criterion
     * of Brent's algorithm.
     * {@code abs} and {@code rel} define a tolerance
     * {@code tol = rel |x| + abs}. {@code rel} should be no smaller than
     * <em>2 macheps</em> and preferably not much less than <em>sqrt(macheps)</em>, * where <em>macheps</em> is the relative machine precision. {@code abs} must
     * be positive.
     *
     * @param rel Relative threshold.
     * @param abs Absolute threshold.
     * @param checker Additional, user-defined, convergence checking
     * procedure.
     * @ if {@code abs <= 0}.
     * @ if {@code rel < 2 * Math.ulp(1d)}.
     */
    public Brent_Optimizer(double rel, double abs, Convergence_Checker<UnivariatePoint_valuePair> checker) 
    {
        super(checker);

        if (rel < MIN_RELATIVE_TOLERANCE) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL, rel, MIN_RELATIVE_TOLERANCE);
        }
        if (abs <= 0) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL_BOUND_EXCLUDED, abs, 0);
        }

        relative_threshold = rel;
        absolute_threshold = abs;
    }

    /**
     * The arguments are used for implementing the original stopping criterion
     * of Brent's algorithm.
     * {@code abs} and {@code rel} define a tolerance
     * {@code tol = rel |x| + abs}. {@code rel} should be no smaller than
     * <em>2 macheps</em> and preferably not much less than <em>sqrt(macheps)</em>, * where <em>macheps</em> is the relative machine precision. {@code abs} must
     * be positive.
     *
     * @param rel Relative threshold.
     * @param abs Absolute threshold.
     * @ if {@code abs <= 0}.
     * @ if {@code rel < 2 * Math.ulp(1d)}.
     */
    public Brent_Optimizer(double rel, double abs) 
    {
        this(rel, abs, NULL);
    }

    /** {@inherit_doc} */
    //override
    protected UnivariatePoint_valuePair do_optimize() 
    {
        const bool is_minim = get_goal_type() == Goal_Type.MINIMIZE;
        const double lo = get_min();
        const double mid = get_start_value();
        const double hi = get_max();

        // Optional additional convergence criteria.
        const Convergence_Checker<UnivariatePoint_valuePair> checker
            = get_convergence_checker();

        double a;
        double b;
        if (lo < hi) 
        {
            a = lo;
            b = hi;
        }
else 
        {
            a = hi;
            b = lo;
        }

        double x = mid;
        double v = x;
        double w = x;
        double d = 0;
        double e = 0;
        double fx = compute_objective_value(x);
        if (!is_minim) 
        {
            fx = -fx;
        }
        double fv = fx;
        double fw = fx;

        UnivariatePoint_valuePair previous = NULL;
        UnivariatePoint_valuePair current
            = UnivariatePoint_valuePair(x, is_minim ? fx : -fx);
        // Best point encountered so far (which is the initial guess).
        UnivariatePoint_valuePair best = current;

        while (true) 
        {
            const double m = 0.5 * (a + b);
            const double tol1 = relative_threshold * std::abs(x) + absolute_threshold;
            const double tol2 = 2 * tol1;

            // Default stopping criterion.
            const bool stop = std::abs(x - m) <= tol2 - 0.5 * (b - a);
            if (!stop) 
            {
                double u;

                if (std::abs(e) > tol1) { // Fit parabola.
                    double r = (x - w) * (fx - fv);
                    double q = (x - v) * (fx - fw);
                    double p = (x - v) * q - (x - w) * r;
                    q = 2 * (q - r);

                    if (q > 0) 
                    {
                        p = -p;
                    }
else 
                    {
                        q = -q;
                    }

                    r = e;
                    e = d;

                    if (p > q * (a - x) &&
                        p < q * (b - x) &&
                        std::abs(p) < std::abs(0.5 * q * r)) 
                        {
                        // Parabolic interpolation step.
                        d = p / q;
                        u = x + d;

                        // f must not be evaluated too close to a or b.
                        if (u - a < tol2 || b - u < tol2) 
                        {
                            if (x <= m) 
                            {
                                d = tol1;
                            }
else 
                            {
                                d = -tol1;
                            }
                        }
                    }
else 
                    {
                        // Golden section step.
                        if (x < m) 
                        {
                            e = b - x;
                        }
else 
                        {
                            e = a - x;
                        }
                        d = GOLDEN_SECTION * e;
                    }
                }
else 
                {
                    // Golden section step.
                    if (x < m) 
                    {
                        e = b - x;
                    }
else 
                    {
                        e = a - x;
                    }
                    d = GOLDEN_SECTION * e;
                }

                // Update by at least "tol1".
                if (std::abs(d) < tol1) 
                {
                    if (d >= 0) 
                    {
                        u = x + tol1;
                    }
else 
                    {
                        u = x - tol1;
                    }
                }
else 
                {
                    u = x + d;
                }

                double fu = compute_objective_value(u);
                if (!is_minim) 
                {
                    fu = -fu;
                }

                // User-defined convergence checker.
                previous = current;
                current = UnivariatePoint_valuePair(u, is_minim ? fu : -fu);
                best = best(best, best(previous, current, is_minim), is_minim);

                if (checker != NULL && checker.converged(get_iterations(), previous, current)) 
                {
                    return best;
                }

                // Update a, b, v, w and x.
                if (fu <= fx) 
                {
                    if (u < x) 
                    {
                        b = x;
                    }
else 
                    {
                        a = x;
                    }
                    v = w;
                    fv = fw;
                    w = x;
                    fw = fx;
                    x = u;
                    fx = fu;
                }
else 
                {
                    if (u < x) 
                    {
                        a = u;
                    }
else 
                    {
                        b = u;
                    }
                    if (fu <= fw ||
                        Precision.equals(w, x)) 
                        {
                        v = w;
                        fv = fw;
                        w = u;
                        fw = fu;
                    }
else if (fu <= fv ||
                               Precision.equals(v, x) ||
                               Precision.equals(v, w)) 
                               {
                        v = u;
                        fv = fu;
                    }
                }
            }
else { // Default termination (Brent's criterion).
                return best(best, best(previous, current, is_minim), is_minim);
            }

            increment_iteration_count();
        }
    }

    /**
     * Selects the best of two points.
     *
     * @param a Point and value.
     * @param b Point and value.
     * @param is_minim {@code true} if the selected point must be the one with
     * the lowest value.
     * @return the best point, or {@code NULL} if {@code a} and {@code b} are
     * both {@code NULL}. When {@code a} and {@code b} have the same function
     * value, {@code a} is returned.
     */
    private UnivariatePoint_valuePair best(UnivariatePoint_valuePair a, UnivariatePoint_valuePair b, bool is_minim) 
    {
        if (a == NULL) 
        {
            return b;
        }
        if (b == NULL) 
        {
            return a;
        }

        if (is_minim) 
        {
            return a.get_value() <= b.get_value() ? a : b;
        }
else 
        {
            return a.get_value() >= b.get_value() ? a : b;
        }
    }
}


