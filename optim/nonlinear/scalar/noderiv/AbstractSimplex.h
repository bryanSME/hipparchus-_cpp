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

//import java.util.Arrays;
//import java.util.Comparator;

//import org.hipparchus.analysis.Multivariate_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.optim.Localized_Optim_Formats;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.optim.Point_valuePair;
//import org.hipparchus.util.Math_Utils;

/**
 * This class : the simplex concept.
 * It is intended to be used in conjunction with {@link Simplex_Optimizer}.
 * <br/>
 * The initial configuration of the simplex is set by the constructors
 * {@link #Abstract_Simplex(std::vector<double>)} or {@link #Abstract_Simplex(std::vector<std::vector<double>>)}.
 * The other {@link #Abstract_Simplexstatic_cast<int>( constructor} will set all steps
 * to 1, thus building a default configuration from a unit hypercube.
 * <br/>
 * Users <em>must</em> call the {@link #build(std::vector<double>) build} method in order
 * to create the data structure that will be acted on by the other methods of
 * this class.
 *
 * @see Simplex_Optimizer
 */
class Abstract_Simplex : Optimization_data 
{
    /** Simplex. */
    private Point_valuePair[] simplex;
    /** Start simplex configuration. */
    private std::vector<std::vector<double>> start_configuration;
    /** Simplex dimension (must be equal to {@code simplex.size() - 1}). */
    private const int dimension;

    /**
     * Build a unit hypercube simplex.
     *
     * @param n Dimension of the simplex.
     */
    protected Abstract_Simplex(const int& n) 
    {
        this(n, 1d);
    }

    /**
     * Build a hypercube simplex with the given side length.
     *
     * @param n Dimension of the simplex.
     * @param side_length Length of the sides of the hypercube.
     */
    protected Abstract_Simplex(const int& n, double side_length) 
    {
        this(create_hypercube_steps(n, side_length));
    }

    /**
     * The start configuration for simplex is built from a box parallel to
     * the canonical axes of the space. The simplex is the subset of vertices
     * of a box parallel to the canonical axes. It is built as the path followed
     * while traveling from one vertex of the box to the diagonally opposite
     * vertex moving only along the box edges. The first vertex of the box will
     * be located at the start point of the optimization.
     * As an example, in dimension 3 a simplex has 4 vertices. Setting the
     * steps to (1, 10, 2) and the start point to (1, 1, 1) would imply the
     * start simplex would be: { (1, 1, 1), (2, 1, 1), (2, 11, 1), (2, 11, 3) }.
     * The first vertex would be set to the start point at (1, 1, 1) and the
     * last vertex would be set to the diagonally opposite vertex at (2, 11, 3).
     *
     * @param steps Steps along the canonical axes representing box edges. They
     * may be negative but not zero.
     * @Null_Argument_Exception if {@code steps} is {@code NULL}.
     * @ if one of the steps is zero.
     */
    protected Abstract_Simplex(const std::vector<double> steps) 
    {
        if (steps == NULL) 
        {
            throw Null_Argument_Exception();
        }
        if (steps.size() == 0) 
        {
            throw (Localized_Core_Formats.ZERO_NOT_ALLOWED);
        }
        dimension = steps.size();

        // Only the relative position of the n const vertices with respect
        // to the first one are stored.
        start_configuration = std::vector<double>(dimension][dimension];
        for (int i{}; i < dimension; i++) 
        {
            const std::vector<double>& vertex_i = start_configuration[i];
            for (int j{}; j < i + 1; j++) 
            {
                if (steps[j] == 0) 
                {
                    throw (Localized_Optim_Formats.EQUAL_VERTICES_IN_SIMPLEX);
                }
                System.arraycopy(steps, 0, vertex_i, 0, j + 1);
            }
        }
    }

    /**
     * The real initial simplex will be set up by moving the reference
     * simplex such that its first point is located at the start point of the
     * optimization.
     *
     * @param reference_simplex Reference simplex.
     * @ if the reference simplex does not
     * contain at least one point.
     * @ if there is a dimension mismatch
     * in the reference simplex.
     * @Illegal_Argument_Exception if one of its vertices is duplicated.
     */
    protected Abstract_Simplex(const std::vector<std::vector<double>> reference_simplex) 
    {
        if (reference_simplex.size() <= 0) 
        {
            throw (Localized_Optim_Formats.SIMPLEX_NEED_ONE_POINT, reference_simplex.size());
        }
        dimension = reference_simplex.size() - 1;

        // Only the relative position of the n const vertices with respect
        // to the first one are stored.
        start_configuration = std::vector<double>(dimension][dimension];
        const std::vector<double> ref0 = reference_simplex[0];

        // Loop over vertices.
        for (int i{}; i < reference_simplex.size(); i++) 
        {
            const std::vector<double> ref_i = reference_simplex[i];

            // Safety checks.
            if (ref_i.size() != dimension) 
            {
                throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, ref_i.size(), dimension);
            }
            for (int j{}; j < i; j++) 
            {
                const std::vector<double> ref_j = reference_simplex[j];
                bool all_equals = true;
                for (int k{}; k < dimension; k++) 
                {
                    if (ref_i[k] != ref_j[k]) 
                    {
                        all_equals = false;
                        break;
                    }
                }
                if (all_equals) 
                {
                    throw (Localized_Optim_Formats.EQUAL_VERTICES_IN_SIMPLEX, i, j);
                }
            }

            // Store vertex i position relative to vertex 0 position.
            if (i > 0) 
            {
                const std::vector<double> conf_i = start_configuration[i - 1];
                for (int k{}; k < dimension; k++) 
                {
                    conf_i[k] = ref_i[k] - ref0[k];
                }
            }
        }
    }

    /**
     * Get simplex dimension.
     *
     * @return the dimension of the simplex.
     */
    public int get_dimension() 
    {
        return dimension;
    }

    /**
     * Get simplex size.
     * After calling the {@link #build(std::vector<double>) build} method, this method will
     * will be equivalent to {@code get_dimension() + 1}.
     *
     * @return the size of the simplex.
     */
    public int get_size() 
    {
        return simplex.size();
    }

    /**
     * Compute the next simplex of the algorithm.
     *
     * @param evaluation_function Evaluation function.
     * @param comparator Comparator to use to sort simplex vertices from best
     * to worst.
     * @org.hipparchus.exception.Math_Illegal_State_Exception
     * if the algorithm fails to converge.
     */
    public virtual void iterate(Multivariate_Function evaluation_function, Comparator<Point_valuePair> comparator);

    /**
     * Build an initial simplex.
     *
     * @param start_point First point of the simplex.
     * @ if the start point does not match
     * simplex dimension.
     */
    public void build(const std::vector<double> start_point) 
    {
        if (dimension != start_point.size()) 
        {
            throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, dimension, start_point.size());
        }

        // Set first vertex.
        simplex = Point_valuePair[dimension + 1];
        simplex[0] = Point_valuePair(start_point,NAN);

        // Set remaining vertices.
        for (int i{}; i < dimension; i++) 
        {
            const std::vector<double> conf_i = start_configuration[i];
            const std::vector<double>& vertex_i = std::vector<double>(dimension];
            for (int k{}; k < dimension; k++) 
            {
                vertex_i[k] = start_point[k] + conf_i[k];
            }
            simplex[i + 1] = Point_valuePair(vertex_i,NAN);
        }
    }

    /**
     * Evaluate all the non-evaluated points of the simplex.
     *
     * @param evaluation_function Evaluation function.
     * @param comparator Comparator to use to sort simplex vertices from best to worst.
     * @org.hipparchus.exception.Math_Illegal_State_Exception
     * if the maximal number of evaluations is exceeded.
     */
    public void evaluate(const Multivariate_Function evaluation_function, const Comparator<Point_valuePair> comparator) 
    {
        // Evaluate the objective function at all non-evaluated simplex points.
        for (int i{}; i < simplex.size(); i++) 
        {
            const Point_valuePair vertex = simplex[i];
            const std::vector<double> point = vertex.get_point_ref();
            if (std::isnan(vertex.get_value())) 
            {
                simplex[i] = Point_valuePair(point, evaluation_function.value(point), false);
            }
        }

        // Sort the simplex from best to worst.
        Arrays.sort(simplex, comparator);
    }

    /**
     * Replace the worst point of the simplex by a point.
     *
     * @param point_value_pair Point to insert.
     * @param comparator Comparator to use for sorting the simplex vertices
     * from best to worst.
     */
    protected void replace_worst_point(Point_valuePair point_value_pair, const Comparator<Point_valuePair> comparator) 
    {
        for (int i{}; i < dimension; i++) 
        {
            if (comparator.compare(simplex[i], point_value_pair) > 0) 
            {
                Point_valuePair tmp = simplex[i];
                simplex[i] = point_value_pair;
                point_value_pair = tmp;
            }
        }
        simplex[dimension] = point_value_pair;
    }

    /**
     * Get the points of the simplex.
     *
     * @return all the simplex points.
     */
    public Point_valuePair[] get_points() 
    {
        const Point_valuePair[] copy = Point_valuePair[simplex.size()];
        System.arraycopy(simplex, 0, copy, 0, simplex.size());
        return copy;
    }

    /**
     * Get the simplex point stored at the requested {@code index}.
     *
     * @param index Location.
     * @return the point at location {@code index}.
     */
    public Point_valuePair get_point(const int& index) 
    {
        Math_Utils::check_range_inclusive(index, 0, simplex.size() - 1);
        return simplex[index];
    }

    /**
     * Store a point at location {@code index}.
     * Note that no deep-copy of {@code point} is performed.
     *
     * @param index Location.
     * @param point New value.
     */
    protected void set_point(const int& index, Point_valuePair point) 
    {
        Math_Utils::check_range_inclusive(index, 0, simplex.size() - 1);
        simplex[index] = point;
    }

    /**
     * Replace all points.
     * Note that no deep-copy of {@code points} is performed.
     *
     * @param points New Points.
     */
    protected void set_points(Point_valuePair[] points) 
    {
        if (points.size() != simplex.size()) 
        {
            throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, points.size(), simplex.size());
        }
        simplex = points.clone();
    }

    /**
     * Create steps for a unit hypercube.
     *
     * @param n Dimension of the hypercube.
     * @param side_length Length of the sides of the hypercube.
     * @return the steps.
     */
    private static std::vector<double> create_hypercube_steps(const int& n, double side_length) 
    {
        const std::vector<double> steps = std::vector<double>(n];
        for (int i{}; i < n; i++) 
        {
            steps[i] = side_length;
        }
        return steps;
    }
}


