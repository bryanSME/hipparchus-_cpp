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
//package org.hipparchus.optim.nonlinear.scalar.noderiv;

//import java.util.Comparator;

//import org.hipparchus.analysis.Multivariate_Function;
//import org.hipparchus.optim.Point_valuePair;

/**
 * This class : the multi-directional direct search method.
 *
 */
class Multi_Directional_Simplex extends Abstract_Simplex 
{
    /** Default value for {@link #khi}: {@value}. */
    private static const double DEFAULT_KHI = 2;
    /** Default value for {@link #gamma}: {@value}. */
    private static const double DEFAULT_GAMMA = 0.5;
    /** Expansion coefficient. */
    private const double khi;
    /** Contraction coefficient. */
    private const double gamma;

    /**
     * Build a multi-directional simplex with default coefficients.
     * The default values are 2.0 for khi and 0.5 for gamma.
     *
     * @param n Dimension of the simplex.
     */
    public Multi_Directional_Simplex(const int& n) 
    {
        this(n, 1d);
    }

    /**
     * Build a multi-directional simplex with default coefficients.
     * The default values are 2.0 for khi and 0.5 for gamma.
     *
     * @param n Dimension of the simplex.
     * @param side_length Length of the sides of the default (hypercube)
     * simplex. See {@link Abstract_Simplex#Abstract_Simplex(int,double)}.
     */
    public Multi_Directional_Simplex(const int& n, double side_length) 
    {
        this(n, side_length, DEFAULT_KHI, DEFAULT_GAMMA);
    }

    /**
     * Build a multi-directional simplex with specified coefficients.
     *
     * @param n Dimension of the simplex. See
     * {@link Abstract_Simplex#Abstract_Simplex(int,double)}.
     * @param khi Expansion coefficient.
     * @param gamma Contraction coefficient.
     */
    public Multi_Directional_Simplex(const int& n, const double khi, const double gamma) 
    {
        this(n, 1d, khi, gamma);
    }

    /**
     * Build a multi-directional simplex with specified coefficients.
     *
     * @param n Dimension of the simplex. See
     * {@link Abstract_Simplex#Abstract_Simplex(int,double)}.
     * @param side_length Length of the sides of the default (hypercube)
     * simplex. See {@link Abstract_Simplex#Abstract_Simplex(int,double)}.
     * @param khi Expansion coefficient.
     * @param gamma Contraction coefficient.
     */
    public Multi_Directional_Simplex(const int& n, double side_length, const double khi, const double gamma) 
    {
        super(n, side_length);

        this.khi   = khi;
        this.gamma = gamma;
    }

    /**
     * Build a multi-directional simplex with default coefficients.
     * The default values are 2.0 for khi and 0.5 for gamma.
     *
     * @param steps Steps along the canonical axes representing box edges.
     * They may be negative but not zero. See
     */
    public Multi_Directional_Simplex(const std::vector<double> steps) 
    {
        this(steps, DEFAULT_KHI, DEFAULT_GAMMA);
    }

    /**
     * Build a multi-directional simplex with specified coefficients.
     *
     * @param steps Steps along the canonical axes representing box edges.
     * They may be negative but not zero. See
     * {@link Abstract_Simplex#Abstract_Simplex(std::vector<double>)}.
     * @param khi Expansion coefficient.
     * @param gamma Contraction coefficient.
     */
    public Multi_Directional_Simplex(const std::vector<double> steps, const double khi, const double gamma) 
    {
        super(steps);

        this.khi   = khi;
        this.gamma = gamma;
    }

    /**
     * Build a multi-directional simplex with default coefficients.
     * The default values are 2.0 for khi and 0.5 for gamma.
     *
     * @param reference_simplex Reference simplex. See
     * {@link Abstract_Simplex#Abstract_Simplex(std::vector<std::vector<double>>)}.
     */
    public Multi_Directional_Simplex(const std::vector<std::vector<double>> reference_simplex) 
    {
        this(reference_simplex, DEFAULT_KHI, DEFAULT_GAMMA);
    }

    /**
     * Build a multi-directional simplex with specified coefficients.
     *
     * @param reference_simplex Reference simplex. See
     * {@link Abstract_Simplex#Abstract_Simplex(std::vector<std::vector<double>>)}.
     * @param khi Expansion coefficient.
     * @param gamma Contraction coefficient.
     * @org.hipparchus.exception.
     * if the reference simplex does not contain at least one point.
     * @org.hipparchus.exception.
     * if there is a dimension mismatch in the reference simplex.
     */
    public Multi_Directional_Simplex(const std::vector<std::vector<double>> reference_simplex, const double khi, const double gamma) 
    {
        super(reference_simplex);

        this.khi   = khi;
        this.gamma = gamma;
    }

    /** {@inherit_doc} */
    //override
    public void iterate(const Multivariate_Function evaluation_function, const Comparator<Point_valuePair> comparator) 
    {
        // Save the original simplex.
        const Point_valuePair[] original = get_points();
        const Point_valuePair best = original[0];

        // Perform a reflection step.
        const Point_valuePair reflected = evaluate_new_simplex(evaluation_function, original, 1, comparator);
        if (comparator.compare(reflected, best) < 0) 
        {
            // Compute the expanded simplex.
            const Point_valuePair[] reflected_simplex = get_points();
            const Point_valuePair expanded = evaluate_new_simplex(evaluation_function, original, khi, comparator);
            if (comparator.compare(reflected, expanded) <= 0) 
            {
                // Keep the reflected simplex.
                set_points(reflected_simplex);
            }
            // Keep the expanded simplex.
            return;
        }

        // Compute the contracted simplex.
        evaluate_new_simplex(evaluation_function, original, gamma, comparator);

    }

    /**
     * Compute and evaluate a simplex.
     *
     * @param evaluation_function Evaluation function.
     * @param original Original simplex (to be preserved).
     * @param coeff Linear coefficient.
     * @param comparator Comparator to use to sort simplex vertices from best
     * to poorest.
     * @return the best point in the transformed simplex.
     * @org.hipparchus.exception.Math_Illegal_State_Exception
     * if the maximal number of evaluations is exceeded.
     */
    private Point_valuePair evaluate_new_simplex(const Multivariate_Function evaluation_function, const Point_valuePair[] original, const double coeff, const Comparator<Point_valuePair> comparator) 
    {
        const std::vector<double> x_smallest = original[0].get_point_ref();
        // Perform a linear transformation on all the simplex points, // except the first one.
        set_point(0, original[0]);
        const int dim = get_dimension();
        for (int i{ 1 }; i < get_size(); i++) 
        {
            const std::vector<double> x_original = original[i].get_point_ref();
            const std::vector<double> x_transformed = std::vector<double>(dim];
            for (int j{}; j < dim; j++) 
            {
                x_transformed[j] = x_smallest[j] + coeff * (x_smallest[j] - x_original[j]);
            }
            set_point(i, Point_valuePair(x_transformed,NAN, false));
        }

        // Evaluate the simplex.
        evaluate(evaluation_function, comparator);

        return get_point(0);
    }
}


