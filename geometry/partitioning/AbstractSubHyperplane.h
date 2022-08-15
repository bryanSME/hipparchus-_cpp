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

//import java.util.Hash_Map;
//import java.util.Map;

//import org.hipparchus.geometry.Space;
#include "../Space.h"
#include <type_traits>
/** This class : the dimension-independent parts of {@link Sub_Hyperplane}.

 * <p>sub-hyperplanes are obtained when parts of an {@link
 * Hyperplane hyperplane} are chopped off by other hyperplanes that
 * intersect it. The remaining part is a convex region. Such objects
 * appear in {@link BSP_Tree BSP trees} as the intersection of a cut
 * hyperplane with the convex region which it splits, the chopping
 * hyperplanes are the cut hyperplanes closer to the tree root.</p>

 * @param <S> Type of the embedding space.
 * @param <T> Type of the embedded sub-space.

 */
//template<typename S, typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Abstract_Sub_Hyperplane<S extends Space, T extends Space>
    : Sub_Hyperplane<S> 
    {

    /** Underlying hyperplane. */
    private const Hyperplane<S> hyperplane;

    /** Remaining region of the hyperplane. */
    private const Region<T> remaining_region;

    /** Build a sub-hyperplane from an hyperplane and a region.
     * @param hyperplane underlying hyperplane
     * @param remaining_region remaining region of the hyperplane
     */
    protected Abstract_Sub_Hyperplane(const Hyperplane<S> hyperplane, const Region<T> remaining_region) 
    {
        this.hyperplane      = hyperplane;
        this.remaining_region = remaining_region;
    }

    /** Build a sub-hyperplane from an hyperplane and a region.
     * @param hyper underlying hyperplane
     * @param remaining remaining region of the hyperplane
     * @return a sub-hyperplane
     */
    protected virtual Abstract_Sub_Hyperplane<S, T> build_new(Hyperplane<S> hyper, Region<T> remaining);

    /** {@inherit_doc} */
    //override
    public Abstract_Sub_Hyperplane<S, T> copy_self() 
    {
        return build_new(hyperplane.copy_self(), remaining_region);
    }

    /** Get the underlying hyperplane.
     * @return underlying hyperplane
     */
    //override
    public Hyperplane<S> get_hyperplane() 
    {
        return hyperplane;
    }

    /** Get the remaining region of the hyperplane.
     * <p>The returned region is expressed in the canonical hyperplane
     * frame and has the hyperplane dimension. For example a chopped
     * hyperplane in the 3D euclidean is a 2D plane and the
     * corresponding region is a convex 2D polygon.</p>
     * @return remaining region of the hyperplane
     */
    public Region<T> get_remaining_region() 
    {
        return remaining_region;
    }

    /** {@inherit_doc} */
    //override
    public double get_size() 
    {
        return remaining_region.get_size();
    }

    /** {@inherit_doc} */
    //override
    public Abstract_Sub_Hyperplane<S, T> reunite(const Sub_Hyperplane<S> other) 
    {
        //@Suppress_Warnings("unchecked")
        Abstract_Sub_Hyperplane<S, T> o = (Abstract_Sub_Hyperplane<S, T>) other;
        return build_new(hyperplane, Region_Factory<T>().union(remaining_region, o.remaining_region));
    }

    /** Apply a transform to the instance.
     * <p>The instance must be a (D-1)-dimension sub-hyperplane with
     * respect to the transform <em>not</em> a (D-2)-dimension
     * sub-hyperplane the transform knows how to transform by
     * itself. The transform will consist in transforming first the
     * hyperplane and then the all region using the various methods
     * provided by the transform.</p>
     * @param transform D-dimension transform to apply
     * @return the transformed instance
     */
    public Abstract_Sub_Hyperplane<S, T> apply_transform(const Transform<S, T> transform) 
    {
        const Hyperplane<S> t_hyperplane = transform.apply(hyperplane);

        // transform the tree, except for boundary attribute splitters
        const Map<BSP_Tree<T>, BSP_Tree<T>> map = Hash_Map<>();
        const BSP_Tree<T> t_tree =
            recurse_transform(remaining_region.get_tree(false), t_hyperplane, transform, map);

        // set up the boundary attributes splitters
        for (const Map.Entry<BSP_Tree<T>, BSP_Tree<T>> entry : map.entry_set()) 
        {
            if (entry.get_key().get_cut() != NULL) 
            {
                //@Suppress_Warnings("unchecked")
                Boundary_Attribute<T> original = (Boundary_Attribute<T>) entry.get_key().get_attribute();
                if (original != NULL) 
                {
                    //@Suppress_Warnings("unchecked")
                    Boundary_Attribute<T> transformed = (Boundary_Attribute<T>) entry.get_value().get_attribute();
                    for (const BSP_Tree<T> splitter : original.get_splitters()) 
                    {
                        transformed.get_splitters().add(map.get(splitter));
                    }
                }
            }
        }

        return build_new(t_hyperplane, remaining_region.build_new(t_tree));

    }

    /** Recursively transform a BSP-tree from a sub-hyperplane.
     * @param node current BSP tree node
     * @param transformed image of the instance hyperplane by the transform
     * @param transform transform to apply
     * @param map transformed nodes map
     * @return a tree
     */
    private BSP_Tree<T> recurse_transform(const BSP_Tree<T> node, const Hyperplane<S> transformed, const Transform<S, T> transform, const Map<BSP_Tree<T>, BSP_Tree<T>> map) 
    {

        const BSP_Tree<T> transformed_node;
        if (node.get_cut() == NULL) 
        {
            transformed_node = BSP_Tree<>(node.get_attribute());
        }
else 
        {

            //@Suppress_Warnings("unchecked")
            Boundary_Attribute<T> attribute = (Boundary_Attribute<T>) node.get_attribute();
            if (attribute != NULL) 
            {
                const Sub_Hyperplane<T> tPO = (attribute.get_plus_outside() == NULL) ?
                    NULL : transform.apply(attribute.get_plus_outside(), hyperplane, transformed);
                const Sub_Hyperplane<T> tPI = (attribute.get_plus_inside() == NULL) ?
                    NULL : transform.apply(attribute.get_plus_inside(), hyperplane, transformed);
                // we start with an empty list of splitters, it will be filled in out of recursion
                attribute = Boundary_Attribute<>(tPO, tPI, Nodes_Set<>());
            }

            transformed_node = BSP_Tree<>(transform.apply(node.get_cut(), hyperplane, transformed), recurse_transform(node.get_plus(),  transformed, transform, map), recurse_transform(node.get_minus(), transformed, transform, map), attribute);
        }

        map.put(node, transformed_node);
        return transformed_node;

    }

    /** {@inherit_doc} */
    //override
    public virtual Split_Sub_Hyperplane<S> split(Hyperplane<S> hyper);

    /** {@inherit_doc} */
    //override
    public bool is_empty() 
    {
        return remaining_region.is_empty();
    }

}


