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
 * This class : the Nelder-_Mead simplex algorithm.
 *
 */
class Nelder_Mead_Simplex extends Abstract_Simplex 
{
    /** Default value for {@link #rho}: {@value}. */
    private static const double DEFAULT_RHO = 1;
    /** Default value for {@link #khi}: {@value}. */
    private static const double DEFAULT_KHI = 2;
    /** Default value for {@link #gamma}: {@value}. */
    private static const double DEFAULT_GAMMA = 0.5;
    /** Default value for {@link #sigma}: {@value}. */
    private static const double DEFAULT_SIGMA = 0.5;
    /** Reflection coefficient. */
    private const double rho;
    /** Expansion coefficient. */
    private const double khi;
    /** Contraction coefficient. */
    private const double gamma;
    /** Shrinkage coefficient. */
    private const double sigma;

    /**
     * Build a Nelder-_Mead simplex with default coefficients.
     * The default coefficients are 1.0 for rho, 2.0 for khi and 0.5
     * for both gamma and sigma.
     *
     * @param n Dimension of the simplex.
     */
    public Nelder_Mead_Simplex(const int& n) 
    {
        this(n, 1d);
    }

    /**
     * Build a Nelder-_Mead simplex with default coefficients.
     * The default coefficients are 1.0 for rho, 2.0 for khi and 0.5
     * for both gamma and sigma.
     *
     * @param n Dimension of the simplex.
     * @param side_length Length of the sides of the default (hypercube)
     * simplex. See {@link Abstract_Simplex#Abstract_Simplex(int,double)}.
     */
    public Nelder_Mead_Simplex(const int& n, double side_length) 
    {
        this(n, side_length, DEFAULT_RHO, DEFAULT_KHI, DEFAULT_GAMMA, DEFAULT_SIGMA);
    }

    /**
     * Build a Nelder-_Mead simplex with specified coefficients.
     *
     * @param n Dimension of the simplex. See
     * {@link Abstract_Simplex#Abstract_Simplex(int,double)}.
     * @param side_length Length of the sides of the default (hypercube)
     * simplex. See {@link Abstract_Simplex#Abstract_Simplex(int,double)}.
     * @param rho Reflection coefficient.
     * @param khi Expansion coefficient.
     * @param gamma Contraction coefficient.
     * @param sigma Shrinkage coefficient.
     */
    public Nelder_Mead_Simplex(const int& n, double side_length, const double rho, const double khi, const double gamma, const double sigma) 
    {
        super(n, side_length);

        this.rho = rho;
        this.khi = khi;
        this.gamma = gamma;
        this.sigma = sigma;
    }

    /**
     * Build a Nelder-_Mead simplex with specified coefficients.
     *
     * @param n Dimension of the simplex. See
     * {@link Abstract_Simplex#Abstract_Simplexstatic_cast<int>(}.
     * @param rho Reflection coefficient.
     * @param khi Expansion coefficient.
     * @param gamma Contraction coefficient.
     * @param sigma Shrinkage coefficient.
     */
    public Nelder_Mead_Simplex(const int& n, const double rho, const double khi, const double gamma, const double sigma) 
    {
        this(n, 1d, rho, khi, gamma, sigma);
    }

    /**
     * Build a Nelder-_Mead simplex with default coefficients.
     * The default coefficients are 1.0 for rho, 2.0 for khi and 0.5
     * for both gamma and sigma.
     *
     * @param steps Steps along the canonical axes representing box edges.
     * They may be negative but not zero. See
     */
    public Nelder_Mead_Simplex(const std::vector<double> steps) 
    {
        this(steps, DEFAULT_RHO, DEFAULT_KHI, DEFAULT_GAMMA, DEFAULT_SIGMA);
    }

    /**
     * Build a Nelder-_Mead simplex with specified coefficients.
     *
     * @param steps Steps along the canonical axes representing box edges.
     * They may be negative but not zero. See
     * {@link Abstract_Simplex#Abstract_Simplex(std::vector<double>)}.
     * @param rho Reflection coefficient.
     * @param khi Expansion coefficient.
     * @param gamma Contraction coefficient.
     * @param sigma Shrinkage coefficient.
     * @Illegal_Argument_Exception if one of the steps is zero.
     */
    public Nelder_Mead_Simplex(const std::vector<double> steps, const double rho, const double khi, const double gamma, const double sigma) 
    {
        super(steps);

        this.rho = rho;
        this.khi = khi;
        this.gamma = gamma;
        this.sigma = sigma;
    }

    /**
     * Build a Nelder-_Mead simplex with default coefficients.
     * The default coefficients are 1.0 for rho, 2.0 for khi and 0.5
     * for both gamma and sigma.
     *
     * @param reference_simplex Reference simplex. See
     * {@link Abstract_Simplex#Abstract_Simplex(std::vector<std::vector<double>>)}.
     */
    public Nelder_Mead_Simplex(const std::vector<std::vector<double>> reference_simplex) 
    {
        this(reference_simplex, DEFAULT_RHO, DEFAULT_KHI, DEFAULT_GAMMA, DEFAULT_SIGMA);
    }

    /**
     * Build a Nelder-_Mead simplex with specified coefficients.
     *
     * @param reference_simplex Reference simplex. See
     * {@link Abstract_Simplex#Abstract_Simplex(std::vector<std::vector<double>>)}.
     * @param rho Reflection coefficient.
     * @param khi Expansion coefficient.
     * @param gamma Contraction coefficient.
     * @param sigma Shrinkage coefficient.
     * @org.hipparchus.exception.
     * if the reference simplex does not contain at least one point.
     * @org.hipparchus.exception.
     * if there is a dimension mismatch in the reference simplex.
     */
    public Nelder_Mead_Simplex(const std::vector<std::vector<double>> reference_simplex, const double rho, const double khi, const double gamma, const double sigma) 
    {
        super(reference_simplex);

        this.rho = rho;
        this.khi = khi;
        this.gamma = gamma;
        this.sigma = sigma;
    }

    /** {@inherit_doc} */
    //override
    public void iterate(const Multivariate_Function evaluation_function, const Comparator<Point_valuePair> comparator) 
    {
        // The simplex has n + 1 points if dimension is n.
        const int n = get_dimension();

        // Interesting values.
        const Point_valuePair best = get_point(0);
        const Point_valuePair second_best = get_point(n - 1);
        const Point_valuePair worst = get_point(n);
        const std::vector<double> x_worst = worst.get_point_ref();

        // Compute the centroid of the best vertices (dismissing the worst
        // point at index n).
        const std::vector<double> centroid = std::vector<double>(n];
        for (int i{}; i < n; i++) 
        {
            const std::vector<double> x = get_point(i).get_point_ref();
            for (int j{}; j < n; j++) 
            {
                centroid[j] += x[j];
            }
        }
        const double scaling = 1.0 / n;
        for (int j{}; j < n; j++) 
        {
            centroid[j] *= scaling;
        }

        // compute the reflection point
        const std::vector<double> x_r = std::vector<double>(n];
        for (int j{}; j < n; j++) 
        {
            x_r[j] = centroid[j] + rho * (centroid[j] - x_worst[j]);
        }
        const Point_valuePair reflected
            = Point_valuePair(x_r, evaluation_function.value(x_r), false);

        if (comparator.compare(best, reflected) <= 0 &&
            comparator.compare(reflected, second_best) < 0) 
            {
            // Accept the reflected point.
            replace_worst_point(reflected, comparator);
        }
else if (comparator.compare(reflected, best) < 0) 
        {
            // Compute the expansion point.
            const std::vector<double> xE = std::vector<double>(n];
            for (int j{}; j < n; j++) 
            {
                xE[j] = centroid[j] + khi * (x_r[j] - centroid[j]);
            }
            const Point_valuePair expanded
                = Point_valuePair(xE, evaluation_function.value(xE), false);

            if (comparator.compare(expanded, reflected) < 0) 
            {
                // Accept the expansion point.
                replace_worst_point(expanded, comparator);
            }
else 
            {
                // Accept the reflected point.
                replace_worst_point(reflected, comparator);
            }
        }
else 
        {
            if (comparator.compare(reflected, worst) < 0) 
            {
                // Perform an outside contraction.
                const std::vector<double> x_c = std::vector<double>(n];
                for (int j{}; j < n; j++) 
                {
                    x_c[j] = centroid[j] + gamma * (x_r[j] - centroid[j]);
                }
                const Point_valuePair out_contracted
                    = Point_valuePair(x_c, evaluation_function.value(x_c), false);
                if (comparator.compare(out_contracted, reflected) <= 0) 
                {
                    // Accept the contraction point.
                    replace_worst_point(out_contracted, comparator);
                    return;
                }
            }
else 
            {
                // Perform an inside contraction.
                const std::vector<double> x_c = std::vector<double>(n];
                for (int j{}; j < n; j++) 
                {
                    x_c[j] = centroid[j] - gamma * (centroid[j] - x_worst[j]);
                }
                const Point_valuePair in_contracted
                    = Point_valuePair(x_c, evaluation_function.value(x_c), false);

                if (comparator.compare(in_contracted, worst) < 0) 
                {
                    // Accept the contraction point.
                    replace_worst_point(in_contracted, comparator);
                    return;
                }
            }

            // Perform a shrink.
            const std::vector<double> x_smallest = get_point(0).get_point_ref();
            for (int i{ 1 }; i <= n; i++) 
            {
                const std::vector<double> x = get_point(i).get_point();
                for (int j{}; j < n; j++) 
                {
                    x[j] = x_smallest[j] + sigma * (x[j] - x_smallest[j]);
                }
                set_point(i, Point_valuePair(x,NAN, false));
            }
            evaluate(evaluation_function, comparator);
        }
    }
}


