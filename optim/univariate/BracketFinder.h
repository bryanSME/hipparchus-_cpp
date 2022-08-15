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

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.optim.nonlinear.scalar.Goal_Type;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Incrementor;

/**
 * Provide an interval that brackets a local optimum of a function.
 * This code is based on a Python implementation (from <em>SciPy</em>, * module {@code optimize.py} v0.5).
 *
 */
class Bracket_Finder 
{
    /** Tolerance to avoid division by zero. */
    private static const double EPS_MIN = 1e-21;
    /**
     * Golden section.
     */
    private static const double GOLD = 1.618034;
    /**
     * Factor for expanding the interval.
     */
    private const double grow_limit;
    /**
     * Number of allowed function evaluations.
     */
    private const int max_evaluations;
    /**
     * Number of function evaluations performed in the last search.
     */
    private int evaluations;
    /**
     * Lower bound of the bracket.
     */
    private double lo;
    /**
     * Higher bound of the bracket.
     */
    private double hi;
    /**
     * Point inside the bracket.
     */
    private double mid;
    /**
     * Function value at {@link #lo}.
     */
    private double f_lo;
    /**
     * Function value at {@link #hi}.
     */
    private double f_hi;
    /**
     * Function value at {@link #mid}.
     */
    private double f_mid;

    /**
     * Constructor with default values {@code 100, 500} (see the
     * {@link #Bracket_Finder(double,int) other constructor}).
     */
    public Bracket_Finder() 
    {
        this(100, 500);
    }

    /**
     * Create a bracketing interval finder.
     *
     * @param grow_limit Expanding factor.
     * @param max_evaluations Maximum number of evaluations allowed for finding
     * a bracketing interval.
     */
    public Bracket_Finder(double grow_limit, int max_evaluations) 
    {
        if (grow_limit <= 0) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL_BOUND_EXCLUDED, grow_limit, 0);
        }
        if (max_evaluations <= 0) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL_BOUND_EXCLUDED, max_evaluations, 0);
        }

        this.grow_limit = grow_limit;
        this.max_evaluations = max_evaluations;
    }

    /**
     * Search points that bracket a local optimum of the function.
     *
     * @param func Function whose optimum should be bracketed.
     * @param goal {@link Goal_Type Goal type}.
     * @param x_a Initial point.
     * @param x_b Initial point.
     * @org.hipparchus.exception.Math_Illegal_State_Exception if the maximum number of evaluations
     * is exceeded.
     */
    public void search(Univariate_Function func, Goal_Type goal, double x_a, double x_b) 
    {
        const Function_Evaluator eval = Function_Evaluator(func);
        const bool is_minim = goal == Goal_Type.MINIMIZE;

        double fA = eval.value(x_a);
        double fb = eval.value(x_b);
        if (is_minim ?
            fA < fb :
            fA > fb) 
            {

            double tmp = x_a;
            x_a = x_b;
            x_b = tmp;

            tmp = fA;
            fA = fb;
            fb = tmp;
        }

        double x_c = x_b + GOLD * (x_b - x_a);
        double fC = eval.value(x_c);

        while (is_minim ? fC < fb : fC > fb) 
        {
            double tmp1 = (x_b - x_a) * (fb - fC);
            double tmp2 = (x_b - x_c) * (fb - fA);

            double val = tmp2 - tmp1;
            double denom = std::abs(val) < EPS_MIN ? 2 * EPS_MIN : 2 * val;

            double w = x_b - ((x_b - x_c) * tmp2 - (x_b - x_a) * tmp1) / denom;
            double w_lim = x_b + grow_limit * (x_c - x_b);

            double fW;
            if ((w - x_c) * (x_b - w) > 0) 
            {
                fW = eval.value(w);
                if (is_minim ?
                    fW < fC :
                    fW > fC) 
                    {
                    x_a = x_b;
                    x_b = w;
                    fA = fb;
                    fb = fW;
                    break;
                }
else if (is_minim ?
                           fW > fb :
                           fW < fb) 
                           {
                    x_c = w;
                    fC = fW;
                    break;
                }
                w = x_c + GOLD * (x_c - x_b);
                fW = eval.value(w);
            }
else if ((w - w_lim) * (w_lim - x_c) >= 0) 
            {
                w = w_lim;
                fW = eval.value(w);
            }
else if ((w - w_lim) * (x_c - w) > 0) 
            {
                fW = eval.value(w);
                if (is_minim ?
                    fW < fC :
                    fW > fC) 
                    {
                    x_b = x_c;
                    x_c = w;
                    w = x_c + GOLD * (x_c - x_b);
                    fb = fC;
                    fC =fW;
                    fW = eval.value(w);
                }
            }
else 
            {
                w = x_c + GOLD * (x_c - x_b);
                fW = eval.value(w);
            }

            x_a = x_b;
            fA = fb;
            x_b = x_c;
            fb = fC;
            x_c = w;
            fC = fW;
        }

        lo = x_a;
        f_lo = fA;
        mid = x_b;
        f_mid = fb;
        hi = x_c;
        f_hi = fC;

        if (lo > hi) 
        {
            double tmp = lo;
            lo = hi;
            hi = tmp;

            tmp = f_lo;
            f_lo = f_hi;
            f_hi = tmp;
        }
    }

    /**
     * @return the number of evaluations.
     */
    public int get_max_evaluations() 
    {
        return max_evaluations;
    }

    /**
     * @return the number of evaluations.
     */
    public int get_evaluations() 
    {
        return evaluations;
    }

    /**
     * @return the lower bound of the bracket.
     * @see #get_f_lo()
     */
    public double get_lo() 
    {
        return lo;
    }

    /**
     * Get function value at {@link #get_lo()}.
     * @return function value at {@link #get_lo()}
     */
    public double get_f_lo() 
    {
        return f_lo;
    }

    /**
     * @return the higher bound of the bracket.
     * @see #get_f_hi()
     */
    public double get_hi() 
    {
        return hi;
    }

    /**
     * Get function value at {@link #get_hi()}.
     * @return function value at {@link #get_hi()}
     */
    public double get_f_hi() 
    {
        return f_hi;
    }

    /**
     * @return a point in the middle of the bracket.
     * @see #get_f_mid()
     */
    public double get_mid() 
    {
        return mid;
    }

    /**
     * Get function value at {@link #get_mid()}.
     * @return function value at {@link #get_mid()}
     */
    public double get_f_mid() 
    {
        return f_mid;
    }

    /**
     * Utility for incrementing a counter at each function evaluation.
     */
    private class Function_Evaluator 
    {
        /** Function. */
        private const Univariate_Function func;
        /** Counter. */
        private const Incrementor inc;

        /**
         * @param func Function.
         */
        Function_Evaluator(Univariate_Function func) 
        {
            this.func = func;
            inc = Incrementor(max_evaluations);
            evaluations = 0;
        }

        /**
         * @param x Argument.
         * @return {@code f(x)}
         * @org.hipparchus.exception.Math_Illegal_State_Exception if the maximal number of evaluations is
         * exceeded.
         */
        double value(double x) 
        {
            inc.increment();
            evaluations = inc.get_count();
            return func.value(x);
        }
    }
}


