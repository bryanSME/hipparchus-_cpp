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
//package org.hipparchus.geometry.euclidean.twod;

//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.euclidean.oned.Euclidean_1D;
//import org.hipparchus.geometry.euclidean.oned.Interval;
//import org.hipparchus.geometry.euclidean.oned.Intervals_Set;
//import org.hipparchus.geometry.euclidean.oned.Oriented_Point;
//import org.hipparchus.geometry.euclidean.oned.Vector_1D;
//import org.hipparchus.geometry.partitioning.Abstract_Sub_Hyperplane;
//import org.hipparchus.geometry.partitioning.BSP_Tree;
//import org.hipparchus.geometry.partitioning.Hyperplane;
//import org.hipparchus.geometry.partitioning.Region;
//import org.hipparchus.geometry.partitioning.Sub_Hyperplane;
//import org.hipparchus.geometry.partitioning.Region.Location;
//import org.hipparchus.util.FastMath;

/** This class represents a sub-hyperplane for {@link Line}.
 */
class Sub_Line extends Abstract_Sub_Hyperplane<Euclidean_2D, Euclidean_1D> 
{

    /** Simple constructor.
     * @param hyperplane underlying hyperplane
     * @param remaining_region remaining region of the hyperplane
     */
    public Sub_Line(const Hyperplane<Euclidean_2D> hyperplane, const Region<Euclidean_1D> remaining_region) 
    {
        super(hyperplane, remaining_region);
    }

    /** Create a sub-line from two endpoints.
     * @param start start point
     * @param end end point
     * @param tolerance tolerance below which points are considered identical
     */
    public Sub_Line(const Vector_2D& start, const Vector_2D& end, const double& tolerance) 
    {
        super(new Line(start, end, tolerance), build_interval_set(start, end, tolerance));
    }

    /** Create a sub-line from a segment.
     * @param segment single segment forming the sub-line
     */
    public Sub_Line(const Segment segment) 
    {
        super(segment.get_line(), build_interval_set(segment.get_start(), segment.get_end(), segment.get_line().get_tolerance()));
    }

    /** Get the endpoints of the sub-line.
     * <p>
     * A subline may be any arbitrary number of disjoints segments, so the endpoints
     * are provided as a list of endpoint pairs. Each element of the list represents
     * one segment, and each segment contains a start point at index 0 and an end point
     * at index 1. If the sub-line is unbounded in the negative infinity direction, * the start point of the first segment will have infinite coordinates. If the
     * sub-line is unbounded in the positive infinity direction, the end point of the
     * last segment will have infinite coordinates. So a sub-line covering the whole
     * line will contain just one row and both elements of this row will have infinite
     * coordinates. If the sub-line is empty, the returned list will contain 0 segments.
     * </p>
     * @return list of segments endpoints
     */
    public List<Segment> get_segments() 
    {

        const Line line = (Line) get_hyperplane();
        const List<Interval> list = ((Intervals_Set) get_remaining_region()).as_list();
        const List<Segment> segments = Array_list<>(list.size());

        for (const Interval interval : list) 
        {
            const Vector_2D& start = line.to_space((Point<Euclidean_1D>) Vector_1D(interval.get_inf()));
            const Vector_2D& end   = line.to_space((Point<Euclidean_1D>) Vector_1D(interval.get_sup()));
            segments.add(new Segment(start, end, line));
        }

        return segments;

    }

    /** Get the intersection of the instance and another sub-line.
     * <p>
     * This method is related to the {@link Line#intersection(Line)
     * intersection} method in the {@link Line Line} class, but in addition
     * to compute the point along infinite lines, it also checks the point
     * lies on both sub-line ranges.
     * </p>
     * @param sub_line other sub-line which may intersect instance
     * @param include_end_points if true, endpoints are considered to belong to
     * instance (i.e. they are closed sets) and may be returned, otherwise endpoints
     * are considered to not belong to instance (i.e. they are open sets) and intersection
     * occurring on endpoints lead to NULL being returned
     * @return the intersection point if there is one, NULL if the sub-lines don't intersect
     */
    public Vector_2D intersection(const Sub_Line& sub_line, const bool include_end_points) 
    {

        // retrieve the underlying lines
        Line line1 = (Line) get_hyperplane();
        Line line2 = (Line) sub_line.get_hyperplane();

        // compute the intersection on infinite line
        Vector_2D v2_d = line1.intersection(line2);
        if (v2_d == NULL) 
        {
            return NULL;
        }

        // check location of point with respect to first sub-line
        Location loc1 = get_remaining_region().check_point(line1.to_sub_space((Point<Euclidean_2D>) v2_d));

        // check location of point with respect to second sub-line
        Location loc2 = sub_line.get_remaining_region().check_point(line2.to_sub_space((Point<Euclidean_2D>) v2_d));

        if (include_end_points) 
        {
            return ((loc1 != Location.OUTSIDE) && (loc2 != Location.OUTSIDE)) ? v2_d : NULL;
        }
else 
        {
            return ((loc1 == Location.INSIDE) && (loc2 == Location.INSIDE)) ? v2_d : NULL;
        }

    }

    /** Build an interval set from two points.
     * @param start start point
     * @param end end point
     * @param tolerance tolerance below which points are considered identical
     * @return an interval set
     */
    private static Intervals_Set build_interval_set(const Vector_2D& start, const Vector_2D& end, const double& tolerance) 
    {
        const Line line = Line(start, end, tolerance);
        return Intervals_Set(line.to_sub_space((Point<Euclidean_2D>) start).get_x(), line.to_sub_space((Point<Euclidean_2D>) end).get_x(), tolerance);
    }

    /** {@inherit_doc} */
    //override
    protected Abstract_Sub_Hyperplane<Euclidean_2D, Euclidean_1D> build_new(const Hyperplane<Euclidean_2D> hyperplane, const Region<Euclidean_1D> remaining_region) 
    {
        return Sub_Line(hyperplane, remaining_region);
    }

    /** {@inherit_doc} */
    //override
    public Split_Sub_Hyperplane<Euclidean_2D> split(const Hyperplane<Euclidean_2D> hyperplane) 
    {

        const Line    this_line  = (Line) get_hyperplane();
        const Line    other_line = (Line) hyperplane;
        const Vector_2D crossing = this_line.intersection(other_line);
        const double& tolerance  = this_line.get_tolerance();

        if (crossing == NULL) 
        {
            // the lines are parallel
            const double global = other_line.get_offset(this_line);
            if (global < -tolerance) 
            {
                return Split_Sub_Hyperplane<Euclidean_2D>(null, this);
            }
else if (global > tolerance) 
            {
                return Split_Sub_Hyperplane<Euclidean_2D>(this, NULL);
            }
else 
            {
                return Split_Sub_Hyperplane<Euclidean_2D>(null, NULL);
            }
        }

        // the lines do intersect
        const bool direct = std::sin(this_line.get_angle() - other_line.get_angle()) < 0;
        const Vector_1D x      = this_line.to_sub_space((Point<Euclidean_2D>) crossing);
        const Sub_Hyperplane<Euclidean_1D> sub_plus  =
                Oriented_Point(x, !direct, tolerance).whole_hyperplane();
        const Sub_Hyperplane<Euclidean_1D> sub_minus =
                Oriented_Point(x,  direct, tolerance).whole_hyperplane();

        const BSP_Tree<Euclidean_1D> split_tree = get_remaining_region().get_tree(false).split(sub_minus);
        const BSP_Tree<Euclidean_1D> plus_tree  = get_remaining_region().is_empty(split_tree.get_plus()) ?
                                               BSP_Tree<>(Boolean.FALSE) :
                                               BSP_Tree<>(sub_plus, BSP_Tree<>(Boolean.FALSE), split_tree.get_plus(), NULL);
        const BSP_Tree<Euclidean_1D> minus_tree = get_remaining_region().is_empty(split_tree.get_minus()) ?
                                               BSP_Tree<>(Boolean.FALSE) :
                                               BSP_Tree<>(sub_minus, BSP_Tree<>(Boolean.FALSE), split_tree.get_minus(), NULL);
        return Split_Sub_Hyperplane<>(new Sub_Line(this_line.copy_self(), Intervals_Set(plus_tree, tolerance)), Sub_Line(this_line.copy_self(), Intervals_Set(minus_tree, tolerance)));

    }

}


