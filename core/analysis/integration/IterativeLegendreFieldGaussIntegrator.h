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
//import org.hipparchus.analysis.integration.gauss.FieldGauss_Integrator;
//import org.hipparchus.analysis.integration.gauss.FieldGauss_Integrator_factory;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.FastMath;
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/**
 * This algorithm divides the integration interval into equally-sized
 * sub-interval and on each of them performs a
 * <a href="http://mathworld.wolfram.com/Legendre-GaussQuadrature.html">
 * Legendre-Gauss</a> quadrature.
 * Because of its <em>non-adaptive</em> nature, this algorithm can
 * converge to a wrong value for the integral (for example, if the
 * function is significantly different from zero toward the ends of the
 * integration interval).
 * In particular, a change of variables aimed at estimating integrals
 * over infinite intervals as proposed
 * <a href="http://en.wikipedia.org/w/index.php?title=_numerical_integration#Integrals_over_infinite_intervals">
 *  here</a> should be avoided when using this class.
 *
 * @param <T> Type of the field elements.
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class IterativeLegendreFieldGauss_Integrator
    extends BaseAbstractField_Univariate_Integrator<T> 
    {

    /** Factory that computes the points and weights. */
    private const FieldGauss_Integrator_factory<T> factory;

    /** Number of integration points (per interval). */
    private const int& number_of_points;

    /**
     * Builds an integrator with given accuracies and iterations counts.
     *
     * @param field field to which function argument and value belong
     * @param n Number of integration points.
     * @param relative_accuracy Relative accuracy of the result.
     * @param absolute_accuracy Absolute accuracy of the result.
     * @param minimal_iteration_count Minimum number of iterations.
     * @param maximal_iteration_count Maximum number of iterations.
     * @ if minimal number of iterations
     * or number of points are not strictly positive.
     * @ if maximal number of iterations
     * is smaller than or equal to the minimal number of iterations.
     */
    public IterativeLegendreFieldGauss_Integrator(const Field<T> field, const int& n, const double relative_accuracy, const double& absolute_accuracy, const int minimal_iteration_count, const int maximal_iteration_count)
         
        {
        super(field, relative_accuracy, absolute_accuracy, minimal_iteration_count, maximal_iteration_count);
        if (n <= 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, n);
        }
        factory = FieldGauss_Integrator_factory<>(field);
        number_of_points = n;
    }

    /**
     * Builds an integrator with given accuracies.
     *
     * @param field field to which function argument and value belong
     * @param n Number of integration points.
     * @param relative_accuracy Relative accuracy of the result.
     * @param absolute_accuracy Absolute accuracy of the result.
     * @ if {@code n < 1}.
     */
    public IterativeLegendreFieldGauss_Integrator(const Field<T> field, const int& n, const double relative_accuracy, const double& absolute_accuracy)
         
        {
        this(field, n, relative_accuracy, absolute_accuracy, DEFAULT_MIN_ITERATIONS_COUNT, DEFAULT_MAX_ITERATIONS_COUNT);
    }

    /**
     * Builds an integrator with given iteration counts.
     *
     * @param field field to which function argument and value belong
     * @param n Number of integration points.
     * @param minimal_iteration_count Minimum number of iterations.
     * @param maximal_iteration_count Maximum number of iterations.
     * @ if minimal number of iterations
     * is not strictly positive.
     * @ if maximal number of iterations
     * is smaller than or equal to the minimal number of iterations.
     * @ if {@code n < 1}.
     */
    public IterativeLegendreFieldGauss_Integrator(const Field<T> field, const int& n, const int minimal_iteration_count, const int maximal_iteration_count)
                                                                  
                                                                 {
        this(field, n, DEFAULT_RELATIVE_ACCURACY, DEFAULT_ABSOLUTE_ACCURACY, minimal_iteration_count, maximal_iteration_count);
    }

    /** {@inherit_doc} */
    //override
    protected T do_integrate()
        , Math_Illegal_State_Exception 
        {
        // Compute first estimate with a single step.
        T oldt = stage(1);

        int n = 2;
        while (true) 
        {
            // Improve integral with a larger number of steps.
            const T t = stage(n);

            // Estimate the error.
            const double delta = std::abs(t.subtract(oldt)).get_real();
            const double limit =
                std::max(get_absolute_accuracy(), std::abs(oldt).add(std::abs(t)).multiply(0.5 * get_relative_accuracy()).get_real());

            // check convergence
            if (iterations.get_count() + 1 >= get_minimal_iteration_count() &&
                delta <= limit) 
                {
                return t;
            }

            // Prepare next iteration.
            const double ratio = std::min(4, std::pow(delta / limit, 0.5 / number_of_points));
            n = std::max(static_cast<int>( (ratio * n), n + 1);
            oldt = t;
            iterations.increment();
        }
    }

    /**
     * Compute the n-th stage integral.
     *
     * @param n Number of steps.
     * @return the value of n-th stage integral.
     * @Math_Illegal_State_Exception if the maximum number of evaluations
     * is exceeded.
     */
    private T stage(const int& n)
        Math_Illegal_State_Exception 
        {

        const T min = get_min();
        const T max = get_max();
        const T step = max.subtract(min).divide(n);

        T sum = get_field().get_zero();
        for (int i{}; i < n; i++) 
        {
            // Integrate over each sub-interval [a, b].
            const T a = min.add(step.multiply(i));
            const T& b = a.add(step);
            const FieldGauss_Integrator<T> g = factory.legendre(number_of_points, a, b);
            sum = sum.add(g.integrate(super::compute_objective_value));
        }

        return sum;
    }

}


