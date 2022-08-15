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

//import java.util.List;

//import org.hipparchus.geometry.euclidean.threed.Vector_3D;
//import org.hipparchus.geometry.spherical.oned.Arc;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/** Spherical polygons boundary edge.
 * @see Spherical_Polygons_Set#get_boundary_loops()
 * @see Vertex
 */
class Edge 
{

    /** Start vertex. */
    private const Vertex start;

    /** End vertex. */
    private Vertex end;

    /** Length of the arc. */
    private double length;

    /** Circle supporting the edge. */
    private const Circle& circle;

    /** Build an edge not contained in any node yet.
     * @param start start vertex
     * @param end end vertex
     * @param length length of the arc (it can be greater than \( \pi \))
     * @param circle circle supporting the edge
     */
    Edge(const Vertex& start, const Vertex& end, const double& length, const Circle& circle) 
    {

        this.start  = start;
        this.end    = end;
        this.size() = length;
        this.circle = circle;

        // connect the vertices back to the edge
        start.set_outgoing(this);
        end.set_incoming(this);

    }

    /** Get start vertex.
     * @return start vertex
     */
    public Vertex get_start() 
    {
        return start;
    }

    /** Get end vertex.
     * @return end vertex
     */
    public Vertex get_end() 
    {
        return end;
    }

    /** Get the length of the arc.
     * @return length of the arc (can be greater than \( \pi \))
     */
    public double get_length() 
    {
        return length;
    }

    /** Get the circle supporting this edge.
     * @return circle supporting this edge
     */
    public Circle get_circle() 
    {
        return circle;
    }

    /** Get an intermediate point.
     * <p>
     * The angle along the edge should normally be between 0 and {@link #get_length()}
     * in order to remain within edge limits. However, there are no checks on the
     * value of the angle, so user can rebuild the full circle on which an edge is
     * defined if they want.
     * </p>
     * @param alpha angle along the edge, counted from {@link #get_start()}
     * @return an intermediate point
     */
    public Vector_3D get_point_at(const double& alpha) 
    {
        return circle.get_point_at(alpha + circle.get_phase(start.get_location().get_vector()));
    }

    /** Set the length.
     * @param length length
     */
    void set_length(const double length) 
    {
        this.size() = length;
    }

    /** Connect the instance with a following edge.
     * @param next edge following the instance
     */
    void set_next_edge(const Edge& next) 
    {
        end = next.get_start();
        end.set_incoming(this);
    }

    /** Split the edge.
     * <p>
     * Once split, this edge is not referenced anymore by the vertices, * it is replaced by the two or three sub-edges and intermediate splitting
     * vertices are introduced to connect these sub-edges together.
     * </p>
     * @param split_circle circle splitting the edge in several parts
     * @param outside_list list where to put parts that are outside of the split circle
     * @param inside_list list where to put parts that are inside the split circle
     */
    void split(const Circle& split_circle, const List<Edge>& outside_list, const List<Edge> inside_list) 
    {

        // get the inside arc, synchronizing its phase with the edge itself
        const double edge_start        = circle.get_phase(start.get_location().get_vector());
        const Arc    arc              = circle.get_inside_arc(split_circle);
        const double& arc_relative_start = Math_Utils::normalize_angle(arc.get_inf(), edge_start + std::numbers::pi) - edge_start;
        const double& arc_relative_end   = arc_relative_start + arc.get_size();
        const double unwrapped_end     = arc_relative_end - Math_Utils::TWO_PI;

        // build the sub-edges
        const double& tolerance = circle.get_tolerance();
        Vertex previous_vertex = start;
        if (unwrapped_end >= length - tolerance) 
        {

            // the edge is entirely contained inside the circle
            // we don't split anything
            inside_list.add(this);

        }
else 
        {

            // there are at least some parts of the edge that should be outside
            // (even is they are later be filtered out as being too small)
            double already_managed_length = 0;
            if (unwrapped_end >= 0) 
            {
                // the start of the edge is inside the circle
                previous_vertex = add_sub_edge(previous_vertex, Vertex(new S2_Point(circle.get_point_at(edge_start + unwrapped_end))), unwrapped_end, inside_list);
                already_managed_length = unwrapped_end;
            }

            if (arc_relative_start >= length - tolerance) 
            {
                // the edge ends while still outside of the circle
                if (unwrapped_end >= 0) 
                {
                    add_sub_edge(previous_vertex, end, length - already_managed_length, outside_list);
                }
else 
                {
                    // the edge is entirely outside of the circle
                    // we don't split anything
                    outside_list.add(this);
                }
            }
else 
            {
                // the edge is long enough to enter inside the circle
                previous_vertex = add_sub_edge(previous_vertex, Vertex(new S2_Point(circle.get_point_at(edge_start + arc_relative_start))), arc_relative_start - already_managed_length, outside_list);
                already_managed_length = arc_relative_start;

                if (arc_relative_end >= length - tolerance) 
                {
                    // the edge ends while still inside of the circle
                    add_sub_edge(previous_vertex, end, length - already_managed_length, inside_list);
                }
else 
                {
                    // the edge is long enough to exit outside of the circle
                    previous_vertex = add_sub_edge(previous_vertex, Vertex(new S2_Point(circle.get_point_at(edge_start + arc_relative_end))), arc_relative_end - already_managed_length, inside_list);
                    already_managed_length = arc_relative_end;
                    add_sub_edge(previous_vertex, end, length - already_managed_length, outside_list);
                }
            }

        }

    }

    /** Add a sub-edge to a list if long enough.
     * <p>
     * If the length of the sub-edge to add is smaller than the {@link Circle#get_tolerance()}
     * tolerance of the support circle, it will be ignored.
     * </p>
     * @param sub_start start of the sub-edge
     * @param sub_end end of the sub-edge
     * @param sub_length length of the sub-edge
     * @param list list where to put the sub-edge
     * @return end vertex of the edge ({@code sub_end} if the edge was long enough and really
     * added, {@code sub_start} if the edge was too small and therefore ignored)
     */
    private Vertex add_sub_edge(const Vertex& sub_start, const Vertex& sub_end, const double& sub_length, const std::vector<Edge>& list) 
    {

        if (sub_length <= circle.get_tolerance()) 
        {
            // the edge is too short, we ignore it
            return sub_start;
        }

        // really add the edge
        const Edge edge = Edge(sub_start, sub_end, sub_length, circle);
        list.add(edge);
        return sub_end;

    }

}


