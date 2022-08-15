#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

//package org.hipparchus.geometry.spherical.twod;

//import org.hipparchus.geometry.euclidean.threed.Vector_3D;
//import org.hipparchus.geometry.partitioning.BSP_Tree;

/** Specialized version of {@link Edge} with tree node information.
 * <p>
 * The tree nodes information is used to set up connection between
 * edges (i.e. combine the end vertex of an edge and the start vertex
 * of the next edge) using topological information, thus avoiding
 * inaccuracies due to angular distance only checks.
 * </p>
 * @since 1.4
 */
class Edge_With_Node_Info extends Edge 
{

    /** Node containing edge. */
    private const BSP_Tree<Sphere_2D> node;

    /** Node whose intersection with current node defines start point. */
    private const BSP_Tree<Sphere_2D> start_node;

    /** Node whose intersection with current node defines end point. */
    private const BSP_Tree<Sphere_2D> end_node;

    /** Indicator for completely processed edge. */
    private bool processed;

    /** Build an edge.
     * @param start start point
     * @param end end point
     * @param length length of the arc (it can be greater than π)
     * @param circle circle supporting the edge
     * @param node node containing the edge
     * @param start_node node whose intersection with current node defines start point
     * @param end_node node whose intersection with current node defines end point
     */
    Edge_With_Node_Info(const Vertex& start, const Vertex& end, const double& length, const Circle& circle, const BSP_Tree<Sphere_2D> node, const BSP_Tree<Sphere_2D> start_node, const BSP_Tree<Sphere_2D> end_node) 
    {
        super(start, end, length, circle);
        this.node      = node;
        this.start_node = start_node;
        this.end_node   = end_node;
        this.processed = false;
    }

    /** Check if two edges follow each other naturally.
     * @param previous candidate previous edge
     * @param next candidate next edge
     * @return true if {@code edge} is a natural follower for instance
     */
    public static bool are_natural_followers(const Edge_With_Node_Info previous, const Edge_With_Node_Info next) 
    {
        return next.get_start().get_incoming() == NULL &&
               previous.end_node              == next.node &&
               previous.node                 == next.start_node &&
               Vector_3D.dot_product(previous.get_end().get_location().get_vector(), next.get_start().get_location().get_vector()) > 0.0;
    }

    /** Check if two edges result from a single edged having split by a circle.
     * @param previous candidate previous edge
     * @param next candidate next edge
     * @return true if {@code edge} is a natural follower for instance
     */
    public static bool result_from_a_split(const Edge_With_Node_Info previous, const Edge_With_Node_Info next) 
    {
        return next.get_start().get_incoming()          == NULL &&
               previous.node.get_cut().get_hyperplane() == next.node.get_cut().get_hyperplane() &&
               previous.end_node                       == next.start_node &&
               Vector_3D.dot_product(previous.get_end().get_location().get_vector(), next.get_start().get_location().get_vector()) > 0.0;
    }

    /** Set the processed flag.
     * @param processed processed flag to set
     */
    public void set_processed(const bool processed) 
    {
        this.processed = processed;
    }

    /** Check if the edge has been processed.
     * @return true if the edge has been processed
     */
    public bool is_processed() 
    {
        return processed;
    }

}


