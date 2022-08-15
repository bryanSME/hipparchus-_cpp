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
//package org.hipparchus.geometry.euclidean.twod.hull;

//import java.util.Array_list;
//import java.util.Collection;
//import java.util.Collections;
//import java.util.Comparator;
//import java.util.List;

//import org.hipparchus.geometry.euclidean.twod.Line;
//import org.hipparchus.geometry.euclidean.twod.Vector_2D;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;

/**
 * Implements Andrew's monotone chain method to generate the convex hull of a finite set of
 * points in the two-dimensional euclidean space.
 * <p>
 * The runtime complexity is O(n log n), with n being the number of input points. If the
 * point set is already sorted (by x-coordinate), the runtime complexity is O(n).
 * <p>
 * The implementation is not sensitive to collinear points on the hull. The parameter
 * {@code include_collinear_points} allows to control the behavior with regard to collinear points.
 * If {@code true}, all points on the boundary of the hull will be added to the hull vertices, * otherwise only the extreme points will be present. By default, collinear points are not added
 * as hull vertices.
 * <p>
 * The {@code tolerance} parameter (default: 1e-10) is used as epsilon criteria to determine
 * identical and collinear points.
 *
 * @see <a href="http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain">
 * Andrew's monotone chain algorithm (Wikibooks)</a>
 */
class Monotone_Chain extends Abstract_Convex_Hull_Generator_2D 
{

    /**
     * Create a Monotone_Chain instance.
     */
    public Monotone_Chain() 
    {
        this(false);
    }

    /**
     * Create a Monotone_Chain instance.
     * @param include_collinear_points whether collinear points shall be added as hull vertices
     */
    public Monotone_Chain(const bool include_collinear_points) 
    {
        super(include_collinear_points);
    }

    /**
     * Create a Monotone_Chain instance.
     * @param include_collinear_points whether collinear points shall be added as hull vertices
     * @param tolerance tolerance below which points are considered identical
     */
    public Monotone_Chain(const bool include_collinear_points, const double& tolerance) 
    {
        super(include_collinear_points, tolerance);
    }

    /** {@inherit_doc} */
    //override
    public Collection<Vector_2D> find_hull_vertices(const Collection<Vector_2D> points) 
    {

        const List<Vector_2D> points_sorted_by_x_axis = Array_list<>(points);

        // sort the points in increasing order on the x-axis
        Collections.sort(points_sorted_by_x_axis, Comparator<Vector_2D>() 
        {
            /** {@inherit_doc} */
            //override
            public int compare(const Vector_2D o1, const Vector_2D o2) 
            {
                const double& tolerance = get_tolerance();
                // need to take the tolerance value into account, otherwise collinear points
                // will not be handled correctly when building the upper/lower hull
                const int diff = Precision.compare_to(o1.get_x(), o2.get_x(), tolerance);
                if (diff == 0) 
                {
                    return Precision.compare_to(o1.get_y(), o2.get_y(), tolerance);
                }
else 
                {
                    return diff;
                }
            }
        });

        // build lower hull
        const List<Vector_2D> lower_hull = Array_list<>();
        for (Vector_2D p : points_sorted_by_x_axis) 
        {
            update_hull(p, lower_hull);
        }

        // build upper hull
        const List<Vector_2D> upper_hull = Array_list<>();
        for (const int& idx = points_sorted_by_x_axis.size() - 1; idx >= 0; idx--) 
        {
            const Vector_2D p = points_sorted_by_x_axis.get(idx);
            update_hull(p, upper_hull);
        }

        // concatenate the lower and upper hulls
        // the last point of each list is omitted as it is repeated at the beginning of the other list
        const List<Vector_2D> hull_vertices = Array_list<>(lower_hull.size() + upper_hull.size() - 2);
        for (const int& idx = 0; idx < lower_hull.size() - 1; idx++) 
        {
            hull_vertices.add(lower_hull.get(idx));
        }
        for (const int& idx = 0; idx < upper_hull.size() - 1; idx++) 
        {
            hull_vertices.add(upper_hull.get(idx));
        }

        // special case: if the lower and upper hull may contain only 1 point if all are identical
        if (hull_vertices.is_empty() && ! lower_hull.is_empty()) 
        {
            hull_vertices.add(lower_hull.get(0));
        }

        return hull_vertices;
    }

    /**
     * Update the partial hull with the current point.
     *
     * @param point the current point
     * @param hull the partial hull
     */
    private void update_hull(const Vector_2D point, const List<Vector_2D> hull) 
    {
        const double& tolerance = get_tolerance();

        if (hull.size() == 1) 
        {
            // ensure that we do not add an identical point
            const Vector_2D p1 = hull.get(0);
            if (p1.distance(point) < tolerance) 
            {
                return;
            }
        }

        while (hull.size() >= 2) 
        {
            const int size = hull.size();
            const Vector_2D p1 = hull.get(size - 2);
            const Vector_2D p2 = hull.get(size - 1);

            const double offset = Line(p1, p2, tolerance).get_offset(point);
            if (std::abs(offset) < tolerance) 
            {
                // the point is collinear to the line (p1, p2)

                const double distance_to_current = p1.distance(point);
                if (distance_to_current < tolerance || p2.distance(point) < tolerance) 
                {
                    // the point is assumed to be identical to either p1 or p2
                    return;
                }

                const double distance_to_last = p1.distance(p2);
                if (is_include_collinear_points()) 
                {
                    const int index = distance_to_current < distance_to_last ? size - 1 : size;
                    hull.add(index, point);
                }
else 
                {
                    if (distance_to_current > distance_to_last) 
                    {
                        hull.remove(size - 1);
                        hull.add(point);
                    }
                }
                return;
            }
else if (offset > 0) 
            {
                hull.remove(size - 1);
            }
else 
            {
                break;
            }
        }
        hull.add(point);
    }

}


