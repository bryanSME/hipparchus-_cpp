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

//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.geometry.partitioning.Abstract_Sub_Hyperplane;
//import org.hipparchus.geometry.partitioning.BSP_Tree;
//import org.hipparchus.geometry.partitioning.BSP_Tree_Visitor;
//import org.hipparchus.geometry.partitioning.Boundary_Attribute;
//import org.hipparchus.geometry.partitioning.Sub_Hyperplane;
//import org.hipparchus.geometry.spherical.oned.Arc;
//import org.hipparchus.geometry.spherical.oned.Arcs_Set;
//import org.hipparchus.geometry.spherical.oned.S1_Point;
//import org.hipparchus.geometry.spherical.oned.Sphere_1D;
//import org.hipparchus.util.FastMath;

/** Visitor building edges.
 * @since 1.4
 */
class Edges_With_Node_Info_Builder : BSP_Tree_Visitor<Sphere_2D> 
{

    /** Tolerance for close nodes connection. */
    private const double& tolerance;

    /** Built segments. */
    private const std::vector<Edge_With_Node_Info>& edges;

    /** Simple constructor.
     * @param tolerance below which points are consider to be identical
     */
    Edges_With_Node_Info_Builder(const double& tolerance) 
    {
        this.tolerance = tolerance;
        this.edges     = Array_list<>();
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
        //@Suppress_Warnings("unchecked")
        const Boundary_Attribute<Sphere_2D> attribute = (Boundary_Attribute<Sphere_2D>) node.get_attribute();
        const Iterable<BSP_Tree<Sphere_2D>> splitters = attribute.get_splitters();
        if (attribute.get_plus_outside() != NULL) 
        {
            add_contribution(attribute.get_plus_outside(), node, splitters, false);
        }
        if (attribute.get_plus_inside() != NULL) 
        {
            add_contribution(attribute.get_plus_inside(), node, splitters, true);
        }
    }

    /** {@inherit_doc} */
    //override
    public void visit_leaf_node(const BSP_Tree<Sphere_2D> node) 
    {
    }

    /** Add the contribution of a boundary edge.
     * @param sub boundary facet
     * @param node node to which the edge belongs
     * @param splitters splitters for the boundary facet
     * @param reversed if true, the facet has the inside on its plus side
     */
    private void add_contribution(const Sub_Hyperplane<Sphere_2D> sub, const BSP_Tree<Sphere_2D> node, const Iterable<BSP_Tree<Sphere_2D>> splitters, const bool reversed) 
    {
        const Abstract_Sub_Hyperplane<Sphere_2D, Sphere_1D> abs_sub =
                        (Abstract_Sub_Hyperplane<Sphere_2D, Sphere_1D>) sub;
        const Circle& circle  = (Circle) sub.get_hyperplane();
        const List<Arc> arcs = ((Arcs_Set) abs_sub.get_remaining_region()).as_list();
        for (const Arc a : arcs) 
        {

            // find the 2D points
            const Vertex start_s = Vertex(circle.to_space(new S1_Point(a.get_inf())));
            const Vertex end_s   = Vertex(circle.to_space(new S1_Point(a.get_sup())));

            // recover the connectivity information
            const BSP_Tree<Sphere_2D> start_n = select_closest(start_s.get_location(), splitters);
            const BSP_Tree<Sphere_2D> end_n   = select_closest(end_s.get_location(), splitters);

            if (reversed) 
            {
                edges.add(new Edge_With_Node_Info(end_s, start_s, a.get_size(), circle.get_reverse(), node, end_n, start_n));
            }
else 
            {
                edges.add(new Edge_With_Node_Info(start_s, end_s, a.get_size(), circle, node, start_n, end_n));
            }

        }
    }

    /** Select the node whose cut sub-hyperplane is closest to specified point.
     * @param point reference point
     * @param candidates candidate nodes
     * @return node closest to point, or NULL if no node is closer than tolerance
     */
    private BSP_Tree<Sphere_2D> select_closest(const S2_Point point, const Iterable<BSP_Tree<Sphere_2D>> candidates) 
    {

        if (point == NULL) 
        {
            return NULL;
        }

        BSP_Tree<Sphere_2D> selected = NULL;

        double min = INFINITY;
        for (const BSP_Tree<Sphere_2D> node : candidates) 
        {
            const double distance = std::abs(node.get_cut().get_hyperplane().get_offset(point));
            if (distance < min) 
            {
                selected = node;
                min      = distance;
            }
        }

        return min <= tolerance ? selected : NULL;

    }

    /** Get the edges.
     * @return built edges
     */
    public List<Edge_With_Node_Info> get_edges() 
    {
        return edges;
    }

}


