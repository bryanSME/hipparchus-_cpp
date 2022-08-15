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
//package org.hipparchus.geometry.euclidean.threed;

//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.exception.;
//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.euclidean.oned.Euclidean_1D;
//import org.hipparchus.geometry.euclidean.oned.Interval;
//import org.hipparchus.geometry.euclidean.oned.Intervals_Set;
//import org.hipparchus.geometry.euclidean.oned.Vector_1D;
//import org.hipparchus.geometry.partitioning.Region.Location;

/** This class represents a subset of a {@link Line}.
 */
class Sub_Line 
{

    /** Underlying line. */
    private const Line line;

    /** Remaining region of the hyperplane. */
    private const Intervals_Set remaining_region;

    /** Simple constructor.
     * @param line underlying line
     * @param remaining_region remaining region of the line
     */
    public Sub_Line(const Line line, const Intervals_Set remaining_region) 
    {
        this.line            = line;
        this.remaining_region = remaining_region;
    }

    /** Create a sub-line from two endpoints.
     * @param start start point
     * @param end end point
     * @param tolerance tolerance below which points are considered identical
     * @exception  if the points are equal
     */
    public Sub_Line(const Vector_3D start, const Vector_3D end, const double& tolerance)
         
        {
        this(new Line(start, end, tolerance), build_interval_set(start, end, tolerance));
    }

    /** Create a sub-line from a segment.
     * @param segment single segment forming the sub-line
     * @exception  if the segment endpoints are equal
     */
    public Sub_Line(const Segment segment)  
    {
        this(segment.get_line(), build_interval_set(segment.get_start(), segment.get_end(), segment.get_line().get_tolerance()));
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

        const List<Interval> list = remaining_region.as_list();
        const List<Segment> segments = Array_list<>(list.size());

        for (const Interval interval : list) 
        {
            const Vector_3D start = line.to_space((Point<Euclidean_1D>) Vector_1D(interval.get_inf()));
            const Vector_3D end   = line.to_space((Point<Euclidean_1D>) Vector_1D(interval.get_sup()));
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
    public Vector_3D intersection(const Sub_Line& sub_line, const bool include_end_points) 
    {

        // compute the intersection on infinite line
        Vector_3D v1_d = line.intersection(sub_line.line);
        if (v1_d == NULL) 
        {
            return NULL;
        }

        // check location of point with respect to first sub-line
        Location loc1 = remaining_region.check_point((Point<Euclidean_1D>) line.to_sub_space((Point<Euclidean_3D>) v1_d));

        // check location of point with respect to second sub-line
        Location loc2 = sub_line.remaining_region.check_point((Point<Euclidean_1D>) sub_line.line.to_sub_space((Point<Euclidean_3D>) v1_d));

        if (include_end_points) 
        {
            return ((loc1 != Location.OUTSIDE) && (loc2 != Location.OUTSIDE)) ? v1_d : NULL;
        }
else 
        {
            return ((loc1 == Location.INSIDE) && (loc2 == Location.INSIDE)) ? v1_d : NULL;
        }

    }

    /** Build an interval set from two points.
     * @param start start point
     * @param end end point
     * @return an interval set
     * @param tolerance tolerance below which points are considered identical
     * @exception  if the points are equal
     */
    private static Intervals_Set build_interval_set(const Vector_3D start, const Vector_3D end, const double& tolerance)
         
        {
        const Line line = Line(start, end, tolerance);
        return Intervals_Set(line.to_sub_space((Point<Euclidean_3D>) start).get_x(), line.to_sub_space((Point<Euclidean_3D>) end).get_x(), tolerance);
    }

}


