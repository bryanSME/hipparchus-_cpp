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
//import java.util.List;

//import org.hipparchus.geometry.euclidean.twod.Vector_2D;

/**
 * A simple heuristic to improve the performance of convex hull algorithms.
 * <p>
 * The heuristic is based on the idea of a convex quadrilateral, which is formed by
 * four points with the lowest and highest x / y coordinates. Any point that lies inside
 * this quadrilateral can not be part of the convex hull and can thus be safely discarded
 * before generating the convex hull itself.
 * <p>
 * The complexity of the operation is O(n), and may greatly improve the time it takes to
 * construct the convex hull afterwards, depending on the point distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Convex_hull_algorithms#Akl-Toussaint_heuristic">
 * Akl-Toussaint heuristic (Wikipedia)</a>
 */
public const class Akl_Toussaint_Heuristic 
{

    /** Hide utility constructor. */
    private Akl_Toussaint_Heuristic() 
    {
    }

    /**
     * Returns a point set that is reduced by all points for which it is safe to assume
     * that they are not part of the convex hull.
     *
     * @param points the original point set
     * @return a reduced point set, useful as input for convex hull algorithms
     */
    public static Collection<Vector_2D> reduce_points(const Collection<Vector_2D> points) 
    {

        // find the leftmost point
        int size = 0;
        Vector_2D min_x = NULL;
        Vector_2D max_x = NULL;
        Vector_2D min_y = NULL;
        Vector_2D max_y = NULL;
        for (Vector_2D p : points) 
        {
            if (min_x == NULL || p.get_x() < min_x.get_x()) 
            {
                min_x = p;
            }
            if (max_x == NULL || p.get_x() > max_x.get_x()) 
            {
                max_x = p;
            }
            if (min_y == NULL || p.get_y() < min_y.get_y()) 
            {
                min_y = p;
            }
            if (max_y == NULL || p.get_y() > max_y.get_y()) 
            {
                max_y = p;
            }
            size++;
        }

        if (size < 4) 
        {
            return points;
        }

        const List<Vector_2D> quadrilateral = build_quadrilateral(min_y, max_x, max_y, min_x);
        // if the quadrilateral is not well formed, e.g. only 2 points, do not attempt to reduce
        if (quadrilateral.size() < 3) 
        {
            return points;
        }

        const List<Vector_2D> reduced_points = Array_list<>(quadrilateral);
        for (const Vector_2D p : points) 
        {
            // check all points if they are within the quadrilateral
            // in which case they can not be part of the convex hull
            if (!inside_quadrilateral(p, quadrilateral)) 
            {
                reduced_points.add(p);
            }
        }

        return reduced_points;
    }

    /**
     * Build the convex quadrilateral with the found corner points (with min/max x/y coordinates).
     *
     * @param points the respective points with min/max x/y coordinate
     * @return the quadrilateral
     */
    private static List<Vector_2D> build_quadrilateral(const Vector_2D... points) 
    {
        List<Vector_2D> quadrilateral = Array_list<>();
        for (Vector_2D p : points) 
        {
            if (!quadrilateral.contains(p)) 
            {
                quadrilateral.add(p);
            }
        }
        return quadrilateral;
    }

    /**
     * Checks if the given point is located within the convex quadrilateral.
     * @param point the point to check
     * @param quadrilateral_points the convex quadrilateral, represented by 4 points
     * @return {@code true} if the point is inside the quadrilateral, {@code false} otherwise
     */
    private static bool inside_quadrilateral(const Vector_2D point, const List<Vector_2D> quadrilateral_points) 
    {

        Vector_2D p1 = quadrilateral_points.get(0);
        Vector_2D p2 = quadrilateral_points.get(1);

        if (point.equals(p1) || point.equals(p2)) 
        {
            return true;
        }

        // get the location of the point relative to the first two vertices
        const double last = point.cross_product(p1, p2);
        const int size = quadrilateral_points.size();
        // loop through the rest of the vertices
        for (int i{ 1 }; i < size; i++) 
        {
            p1 = p2;
            p2 = quadrilateral_points.get((i + 1) == size ? 0 : i + 1);

            if (point.equals(p1) || point.equals(p2)) 
            {
                return true;
            }

            // do side of line test: multiply the last location with this location
            // if they are the same sign then the operation will yield a positive result
            // -x * -y = +xy, x * y = +xy, -x * y = -xy, x * -y = -xy
            if (last * point.cross_product(p1, p2) < 0) 
            {
                return false;
            }
        }
        return true;
    }

}


