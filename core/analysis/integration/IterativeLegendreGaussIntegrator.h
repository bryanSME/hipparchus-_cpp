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
#include <algorithm>

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
//package org.hipparchus.analysis.integration;

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.analysis.integration.gauss.Gauss_Integrator;
//import org.hipparchus.analysis.integration.gauss.Gauss_Integrator_factory;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.FastMath;
#include <algorithm>
#include "../integration/gauss/GaussIntegrator.h"
#include "BaseAbstractFieldUnivariateIntegrator.h"
#include "../integration/BaseAbstractUnivariateIntegrator.h"
#include "../integration/gauss/GaussIntegratorFactory.h"
#include "../UnivariateFunction.h"

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
 */

class Iterative_Legendre_Gauss_Integrator : public Base_Abstract_Univariate_Integrator
{
private:
    /** Factory that computes the points and weights. */
    static const Gauss_Integrator_factory FACTORY = Gauss_Integrator_factory();
    /** Number of integration points (per interval). */
    int my_number_of_points;

    /**
     * Compute the n-th stage integral.
     *
     * @param n Number of steps.
     * @return the value of n-th stage integral.
     * @Math_Illegal_State_Exception if the maximum number of evaluations
     * is exceeded.
     */
    double stage(const int& n)
    {
        // Function to be integrated is stored in the base class.
        const Univariate_Function f = Univariate_Function()
        {
            /** {@inherit_doc} */
            //override
            public double value(const double& x)
            {
                return compute_objective_value(x);
            }
        };

        const double min = get_min();
        const double max = get_max();
        const double step = (max - min) / n;

        double sum{};
        for (int i{}; i < n; i++)
        {
            // Integrate over each sub-interval [a, b].
            const double& a = min + i * step;
            const double b = a + step;
            const Gauss_Integrator g = FACTORY.legendre_high_precision(number_of_points, a, b);
            sum += g.integrate(f);
        }

        return sum;
    }

public:
    /**
     * Builds an integrator with given accuracies and iterations counts.
     *
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
    Iterative_Legendre_Gauss_Integrator(const int& n, const double& relative_accuracy, const double& absolute_accuracy, const int& minimal_iteration_count, const int& maximal_iteration_count)
    {
        super(relative_accuracy, absolute_accuracy, minimal_iteration_count, maximal_iteration_count);
        if (n <= 0)
        {
            throw std::exception("not implmented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_POINTS, n);
        }
        my_number_of_points = n;
    }

    /**
     * Builds an integrator with given accuracies.
     *
     * @param n Number of integration points.
     * @param relative_accuracy Relative accuracy of the result.
     * @param absolute_accuracy Absolute accuracy of the result.
     * @ if {@code n < 1}.
     */
    Iterative_Legendre_Gauss_Integrator(const int& n, const double& relative_accuracy, const double& absolute_accuracy)
    {
        Iterative_Legendre_Gauss_Integrator(n, relative_accuracy, absolute_accuracy, DEFAULT_MIN_ITERATIONS_COUNT, DEFAULT_MAX_ITERATIONS_COUNT);
    }

    /**
     * Builds an integrator with given iteration counts.
     *
     * @param n Number of integration points.
     * @param minimal_iteration_count Minimum number of iterations.
     * @param maximal_iteration_count Maximum number of iterations.
     * @ if minimal number of iterations
     * is not strictly positive.
     * @ if maximal number of iterations
     * is smaller than or equal to the minimal number of iterations.
     * @ if {@code n < 1}.
     */
    Iterative_Legendre_Gauss_Integrator(const int& n, const int minimal_iteration_count, const int maximal_iteration_count)
    {
        Iterative_Legendre_Gauss_Integrator(n, DEFAULT_RELATIVE_ACCURACY, DEFAULT_ABSOLUTE_ACCURACY, minimal_iteration_count, maximal_iteration_count);
    }

protected:
    /** {@inherit_doc} */
    //override
    double do_integrate()
    {
        // Compute first estimate with a single step.
        double oldt = stage(1);

        int n{ 2 };
        while (true)
        {
            // Improve integral with a larger number of steps.
            const double t = stage(n);

            // Estimate the error.
            const double delta = std::abs(t - oldt);
            const double limit =
                std::max(get_absolute_accuracy(), get_relative_accuracy() * (std::abs(oldt) + std::abs(t)) * 0.5);

            // check convergence
            if (iterations.get_count() + 1 >= get_minimal_iteration_count() &&
                delta <= limit)
            {
                return t;
            }

            // Prepare next iteration.
            const double ratio = std::min(4, std::pow(delta / limit, 0.5 / my_number_of_points));
            n = std::max(static_cast<int>((ratio * n), n + 1);
            oldt = t;
            iterations.increment();
        }
    }
};