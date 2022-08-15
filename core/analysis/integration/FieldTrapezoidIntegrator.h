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

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.FastMath;
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/**
 * Implements the <a href="http://mathworld.wolfram.com/TrapezoidalRule.html">
 * Trapezoid Rule</a> for integration of real univariate functions. For
 * reference, see <b>Introduction to Numerical Analysis</b>, ISBN 038795452X, * chapter 3.
 * <p>
 * The function should be integrable.</p>
 * @param <T> Type of the field elements.
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Trapezoid_Integrator : public BaseAbstractField_Univariate_Integrator<T> 
{

    /** Maximum number of iterations for trapezoid. */
    public static const int TRAPEZOID_MAX_ITERATIONS_COUNT = 64;

    /** Intermediate result. */
    private T s;

    /**
     * Build a trapezoid integrator with given accuracies and iterations counts.
     * @param field field to which function argument and value belong
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
    public Field_Trapezoid_Integrator(const Field<T> field, const double relative_accuracy, const double& absolute_accuracy, const int minimal_iteration_count, const int maximal_iteration_count)
         
        {
        super(field, relative_accuracy, absolute_accuracy, minimal_iteration_count, maximal_iteration_count);
        if (maximal_iteration_count > TRAPEZOID_MAX_ITERATIONS_COUNT) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, maximal_iteration_count, TRAPEZOID_MAX_ITERATIONS_COUNT);
        }
    }

    /**
     * Build a trapezoid integrator with given iteration counts.
     * @param field field to which function argument and value belong
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
    public Field_Trapezoid_Integrator(const Field<T> field, const int minimal_iteration_count, const int maximal_iteration_count)
         
        {
        super(field, minimal_iteration_count, maximal_iteration_count);
        if (maximal_iteration_count > TRAPEZOID_MAX_ITERATIONS_COUNT) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, maximal_iteration_count, TRAPEZOID_MAX_ITERATIONS_COUNT);
        }
    }

    /**
     * Construct a trapezoid integrator with default settings.
     * @param field field to which function argument and value belong
     * (max iteration count set to {@link #TRAPEZOID_MAX_ITERATIONS_COUNT})
     */
    public Field_Trapezoid_Integrator(const Field<T> field) 
    {
        super(field, DEFAULT_MIN_ITERATIONS_COUNT, TRAPEZOID_MAX_ITERATIONS_COUNT);
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
    T stage(const BaseAbstractField_Univariate_Integrator<T> base_integrator, const int& n)
        Math_Illegal_State_Exception 
        {

        if (n == 0) 
        {
            const T max = base_integrator.get_max();
            const T min = base_integrator.get_min();
            s = max.subtract(min).multiply(0.5).
                multiply(base_integrator.compute_objective_value(min).
                         add(base_integrator.compute_objective_value(max)));
            return s;
        }
else 
        {
            const long np = 1L << (n-1);           // number of points in this stage
            T sum = get_field().get_zero();
            const T max = base_integrator.get_max();
            const T min = base_integrator.get_min();
            // spacing between adjacent points
            const T spacing = max.subtract(min).divide(np);
            T x = min.add(spacing.multiply(0.5));    // the first point
            for (long i = 0; i < np; i++) 
            {
                sum = sum.add(base_integrator.compute_objective_value(x));
                x = x.add(spacing);
            }
            // add the sum to previously calculated result
            s = s.add(sum.multiply(spacing)).multiply(0.5);
            return s;
        }
    }

    /** {@inherit_doc} */
    //override
    protected T do_integrate()
        , Math_Illegal_State_Exception 
        {

        T oldt = stage(this, 0);
        iterations.increment();
        while (true) 
        {
            const int i = iterations.get_count();
            const T t = stage(this, i);
            if (i >= get_minimal_iteration_count()) 
            {
                const double delta  = std::abs(t.subtract(oldt)).get_real();
                const double rlimit = std::abs(oldt).add(std::abs(t)).multiply(0.5 * get_relative_accuracy()).get_real();
                if ((delta <= rlimit) || (delta <= get_absolute_accuracy())) 
                {
                    return t;
                }
            }
            oldt = t;
            iterations.increment();
        }

    }

}


