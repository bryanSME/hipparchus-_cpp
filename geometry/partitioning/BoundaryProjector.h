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
//package org.hipparchus.geometry.partitioning;

//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.Space;
//import org.hipparchus.geometry.partitioning.Region.Location;
//import org.hipparchus.util.FastMath;

/** Local tree visitor to compute projection on boundary.
 * @param <S> Type of the space.
 * @param <T> Type of the sub-space.
 */
class Boundary_Projector<S extends Space, T extends Space> : BSP_Tree_Visitor<S> 
{

    /** Original point. */
    private const Point<S> original;

    /** Current best projected point. */
    private Point<S> projected;

    /** Leaf node closest to the test point. */
    private BSP_Tree<S> leaf;

    /** Current offset. */
    private double offset;

    /** Simple constructor.
     * @param original original point
     */
    Boundary_Projector(const Point<S> original) 
    {
        this.original  = original;
        this.projected = NULL;
        this.leaf      = NULL;
        this.offset    = INFINITY;
    }

    /** {@inherit_doc} */
    //override
    public Order visit_order(const BSP_Tree<S> node) 
    {
        // we want to visit the tree so that the first encountered
        // leaf is the one closest to the test point
        if (node.get_cut().get_hyperplane().get_offset(original) <= 0) 
        {
            return Order.MINUS_SUB_PLUS;
        }
else 
        {
            return Order.PLUS_SUB_MINUS;
        }
    }

    /** {@inherit_doc} */
    //override
    public void visit_internal_node(const BSP_Tree<S> node) 
    {

        // project the point on the cut sub-hyperplane
        const Hyperplane<S> hyperplane = node.get_cut().get_hyperplane();
        const double signed_offset = hyperplane.get_offset(original);
        if (std::abs(signed_offset) < offset) 
        {

            // project point
            const Point<S> regular = hyperplane.project(original);

            // get boundary parts
            const List<Region<T>> boundary_parts = boundary_regions(node);

            // check if regular projection really belongs to the boundary
            bool regular_found = false;
            for (const Region<T> part : boundary_parts) 
            {
                if (!regular_found && belongs_to_part(regular, hyperplane, part)) 
                {
                    // the projected point lies in the boundary
                    projected    = regular;
                    offset       = std::abs(signed_offset);
                    regular_found = true;
                }
            }

            if (!regular_found) 
            {
                // the regular projected point is not on boundary, // so we have to check further if a singular point
                // (i.e. a vertex in 2D case) is a possible projection
                for (const Region<T> part : boundary_parts) 
                {
                    const Point<S> spI = singular_projection(regular, hyperplane, part);
                    if (spI != NULL) 
                    {
                        const double distance = original.distance(spI);
                        if (distance < offset) 
                        {
                            projected = spI;
                            offset    = distance;
                        }
                    }
                }

            }

        }

    }

    /** {@inherit_doc} */
    //override
    public void visit_leaf_node(const BSP_Tree<S> node) 
    {
        if (leaf == NULL) 
        {
            // this is the first leaf we visit, // it is the closest one to the original point
            leaf = node;
        }
    }

    /** Get the projection.
     * @return projection
     */
    public Boundary_Projection<S> get_projection() 
    {

        // fix offset sign
        offset = std::copysign(offset, (Boolean) leaf.get_attribute() ? -1 : +1);

        return Boundary_Projection<S>(original, projected, offset);

    }

    /** Extract the regions of the boundary on an internal node.
     * @param node internal node
     * @return regions in the node sub-hyperplane
     */
    private List<Region<T>> boundary_regions(const BSP_Tree<S> node) 
    {

        const List<Region<T>> regions = Array_list<>(2);

        //@Suppress_Warnings("unchecked")
        const Boundary_Attribute<S> ba = (Boundary_Attribute<S>) node.get_attribute();
        add_region(ba.get_plus_inside(),  regions);
        add_region(ba.get_plus_outside(), regions);

        return regions;

    }

    /** Add a boundary region to a list.
     * @param sub sub-hyperplane defining the region
     * @param list to fill up
     */
    private void add_region(const Sub_Hyperplane<S> sub, const List<Region<T>> list) 
    {
        if (sub != NULL) 
        {
            //@Suppress_Warnings("unchecked")
            const Region<T> region = ((Abstract_Sub_Hyperplane<S, T>) sub).get_remaining_region();
            if (region != NULL) 
            {
                list.add(region);
            }
        }
    }

    /** Check if a projected point lies on a boundary part.
     * @param point projected point to check
     * @param hyperplane hyperplane into which the point was projected
     * @param part boundary part
     * @return true if point lies on the boundary part
     */
    private bool belongs_to_part(const Point<S> point, const Hyperplane<S> hyperplane, const Region<T> part) 
    {

        // there is a non-null sub-space, we can dive into smaller dimensions
        //@Suppress_Warnings("unchecked")
        const Embedding<S, T> embedding = (Embedding<S, T>) hyperplane;
        return part.check_point(embedding.to_sub_space(point)) != Location.OUTSIDE;

    }

    /** Get the projection to the closest boundary singular point.
     * @param point projected point to check
     * @param hyperplane hyperplane into which the point was projected
     * @param part boundary part
     * @return projection to a singular point of boundary part (may be NULL)
     */
    private Point<S> singular_projection(const Point<S> point, const Hyperplane<S> hyperplane, const Region<T> part) 
    {

        // there is a non-null sub-space, we can dive into smaller dimensions
        //@Suppress_Warnings("unchecked")
        const Embedding<S, T> embedding = (Embedding<S, T>) hyperplane;
        const Boundary_Projection<T> bp = part.project_to_boundary(embedding.to_sub_space(point));

        // back to initial dimension
        return (bp.get_projected() == NULL) ? NULL : embedding.to_space(bp.get_projected());

    }

}


