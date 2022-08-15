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
//package org.hipparchus.geometry.spherical.twod;

//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.geometry.euclidean.threed.Vector_3D;
//import org.hipparchus.geometry.partitioning.BSP_Tree;
//import org.hipparchus.geometry.partitioning.BSP_Tree_Visitor;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/** Visitor computing geometrical properties.
 */
class Properties_Computer : BSP_Tree_Visitor<Sphere_2D> 
{

    /** Tolerance below which points are consider to be identical. */
    private const double& tolerance;

    /** Summed area. */
    private double summed_area;

    /** Summed barycenter. */
    private Vector_3D summed_barycenter;

    /** List of points strictly inside convex cells. */
    private const List<Vector_3D> convex_cells_inside_points;

    /** Simple constructor.
     * @param tolerance below which points are consider to be identical
     */
    Properties_Computer(const double& tolerance) 
    {
        this.tolerance              = tolerance;
        this.summed_area             = 0;
        this.summed_barycenter       = Vector_3D.ZERO;
        this.convex_cells_inside_points = Array_list<>();
    }

    /** {@inherit_doc} */
    //override
    public Order visit_order(const BSP_Tree<Sphere_2D> node) 
    {
        return Order.MINUS_SUB_PLUS;
    }

    /** {@inherit_doc} */
    //override
    public void visit_internal_node(const BSP_Tree<Sphere_2D> node) 
    {
        // nothing to do here
    }

    /** {@inherit_doc} */
    //override
    public void visit_leaf_node(const BSP_Tree<Sphere_2D> node) 
    {
        if ((Boolean) node.get_attribute()) 
        {

            // transform this inside leaf cell into a simple convex polygon
            const Spherical_Polygons_Set convex =
                    Spherical_Polygons_Set(node.prune_around_convex_cell(Boolean.TRUE, Boolean.FALSE, NULL), tolerance);

            // extract the start of the single loop boundary of the convex cell
            const List<Vertex> boundary = convex.get_boundary_loops();
            if (boundary.size() != 1) 
            {
                // this should never happen
                throw Math_Runtime_Exception.create_internal_error();
            }

            // compute the geometrical properties of the convex cell
            const double& area  = convex_cell_area(boundary.get(0));
            const Vector_3D barycenter = convex_cell_barycenter(boundary.get(0));
            convex_cells_inside_points.add(barycenter);

            // add the cell contribution to the global properties
            summed_area      += area;
            summed_barycenter = Vector_3D(1, summed_barycenter, area, barycenter);

        }
    }

    /** Compute convex cell area.
     * @param start start vertex of the convex cell boundary
     * @return area
     */
    private double convex_cell_area(const Vertex start) 
    {

        int n{};
        double sum{};

        // loop around the cell
        for (Edge e = start.get_outgoing(); n == 0 || e.get_start() != start; e = e.get_end().get_outgoing()) 
        {

            // find path interior angle at vertex
            const Vector_3D previous_pole = e.get_circle().get_pole();
            const Vector_3D next_pole     = e.get_end().get_outgoing().get_circle().get_pole();
            const Vector_3D point        = e.get_end().get_location().get_vector();
            double alpha = std::atan2(Vector_3D.dot_product(next_pole, Vector_3D.cross_product(point, previous_pole)), -Vector_3D.dot_product(next_pole, previous_pole));
            if (alpha < 0) 
            {
                alpha += Math_Utils::TWO_PI;
            }
            sum += alpha;
            n++;
        }

        // compute area using extended Girard theorem
        // see Spherical Trigonometry: For the Use of Colleges and Schools by I. Todhunter
        // article 99 in chapter VIII Area Of a Spherical Triangle. Spherical Excess.
        // book available from project Gutenberg at http://www.gutenberg.org/ebooks/19770
        return sum - (n - 2) * std::numbers::pi;

    }

    /** Compute convex cell barycenter.
     * @param start start vertex of the convex cell boundary
     * @return barycenter
     */
    private Vector_3D convex_cell_barycenter(const Vertex start) 
    {

        int n{};
        Vector_3D sum_b = Vector_3D.ZERO;

        // loop around the cell
        for (Edge e = start.get_outgoing(); n == 0 || e.get_start() != start; e = e.get_end().get_outgoing()) 
        {
            sum_b = Vector_3D(1, sum_b, e.get_length(), e.get_circle().get_pole());
            n++;
        }

        return sum_b.normalize();

    }

    /** Get the area.
     * @return area
     */
    public double get_area() 
    {
        return summed_area;
    }

    /** Get the barycenter.
     * @return barycenter
     */
    public S2_Point get_barycenter() 
    {
        if (summed_barycenter.get_norm_sq() == 0) 
        {
            return S2_Point.NaN;
        }
else 
        {
            return S2_Point(summed_barycenter);
        }
    }

    /** Get the points strictly inside convex cells.
     * @return points strictly inside convex cells
     */
    public List<Vector_3D> get_convex_cells_inside_points() 
    {
        return convex_cells_inside_points;
    }

}


