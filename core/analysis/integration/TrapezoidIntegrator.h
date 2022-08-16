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

/**
 * Implements the <a href="http://mathworld.wolfram.com/TrapezoidalRule.html">
 * Trapezoid Rule</a> for integration of real univariate functions. For
 * reference, see <b>Introduction to Numerical Analysis</b>, ISBN 038795452X, * chapter 3.
 * <p>
 * The function should be integrable.</p>
 *
 */
class Trapezoid_Integrator : public Base_Abstract_Univariate_Integrator 
{
private:
    /** Intermediate result. */
    double my_s;

public:
    /** Maximum number of iterations for trapezoid. */
    static constexpr int TRAPEZOID_MAX_ITERATIONS_COUNT{ 64 };



    /**
     * Build a trapezoid integrator with given accuracies and iterations counts.
     * @param relative_accuracy relative accuracy of the result
     * @param absolute_accuracy absolute accuracy of the result
     * @param minimal_iteration_count minimum number of iterations
     * @param maximal_iteration_count maximum number of iterations
     * (must be less than or equal to {@link #TRAPEZOID_MAX_ITERATIONS_COUNT}
     * @exception  if minimal number of iterations
     * is not strictly positive
     * @exception  if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     * @exception  if maximal number of iterations
     * is greater than {@link #TRAPEZOID_MAX_ITERATIONS_COUNT}
     */
    Trapezoid_Integrator(const double& relative_accuracy, const double& absolute_accuracy, const int minimal_iteration_count, const int maximal_iteration_count)
    {
        super(relative_accuracy, absolute_accuracy, minimal_iteration_count, maximal_iteration_count);
        if (maximal_iteration_count > TRAPEZOID_MAX_ITERATIONS_COUNT) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, maximal_iteration_count, TRAPEZOID_MAX_ITERATIONS_COUNT);
        }
    }

    /**
     * Build a trapezoid integrator with given iteration counts.
     * @param minimal_iteration_count minimum number of iterations
     * @param maximal_iteration_count maximum number of iterations
     * (must be less than or equal to {@link #TRAPEZOID_MAX_ITERATIONS_COUNT}
     * @exception  if minimal number of iterations
     * is not strictly positive
     * @exception  if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     * @exception  if maximal number of iterations
     * is greater than {@link #TRAPEZOID_MAX_ITERATIONS_COUNT}
     */
    Trapezoid_Integrator(const int minimal_iteration_count, const int maximal_iteration_count)
    {
        super(minimal_iteration_count, maximal_iteration_count);
        if (maximal_iteration_count > TRAPEZOID_MAX_ITERATIONS_COUNT) 
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, maximal_iteration_count, TRAPEZOID_MAX_ITERATIONS_COUNT);
        }
    }

    /**
     * Construct a trapezoid integrator with default settings.
     * (max iteration count set to {@link #TRAPEZOID_MAX_ITERATIONS_COUNT})
     */
    Trapezoid_Integrator() 
    {
        super(DEFAULT_MIN_ITERATIONS_COUNT, TRAPEZOID_MAX_ITERATIONS_COUNT);
    }

    /**
     * Compute the n-th stage integral of trapezoid rule. This function
     * should only be called by API <code>integrate()</code> in the //package.
     * To save time it does not verify arguments - caller does.
     * <p>
     * The interval is divided equally into 2^n sections rather than an
     * arbitrary m sections because this configuration can best utilize the
     * already computed values.</p>
     *
     * @param base_integrator integrator holding integration parameters
     * @param n the stage of 1/2 refinement, n = 0 is no refinement
     * @return the value of n-th stage integral
     * @Math_Illegal_State_Exception if the maximal number of evaluations
     * is exceeded.
     */
    double stage(const Base_Abstract_Univariate_Integrator& base_integrator, const int& n)
    {
        if (n == 0) 
        {
            const double max = base_integrator.get_max();
            const double min = base_integrator.get_min();
            my_s = 0.5 * (max - min) * (base_integrator.compute_objective_value(min) + base_integrator.compute_objective_value(max));
            return my_s;
        }
        const long np = 1L << (n-1);           // number of points in this stage
        double sum{};
        const double max = base_integrator.get_max();
        const double min = base_integrator.get_min();
        // spacing between adjacent points
        const double spacing = (max - min) / np;
        double x = min + 0.5 * spacing;    // the first point
        for (long i = 0; i < np; i++) 
        {
            sum += base_integrator.compute_objective_value(x);
            x += spacing;
        }
        // add the sum to previously calculated result
        s = 0.5 * (s + sum * spacing);
        return s;
    }

protected:
    /** {@inherit_doc} */
    //override
    double do_integrate()
    {
        double oldt = stage(this, 0);
        iterations.increment();
        while (true) 
        {
            const int i = iterations.get_count();
            const double t = stage(this, i);
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
            iterations.increment();
        }

    }

};