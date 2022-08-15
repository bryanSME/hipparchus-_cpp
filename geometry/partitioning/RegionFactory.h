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
#include "../Space.h"
#include <type_traits>
//import org.hipparchus.exception.;
//import org.hipparchus.geometry.Localized_Geometry_Formats;
//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.Space;
//import org.hipparchus.geometry.partitioning.BSP_Tree.Vanishing_Cut_Handler;
//import org.hipparchus.geometry.partitioning.Region.Location;
//import org.hipparchus.geometry.partitioning.Sub_Hyperplane.Split_Sub_Hyperplane;

/** This class is a factory for {@link Region}.

 * @param <S> Type of the space.

 */
template<typename S, typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Region_Factory
{

    /** Visitor removing internal nodes attributes. */
    private const Nodes_Cleaner node_cleaner;

    /** Simple constructor.
     */
    public Region_Factory() 
    {
        node_cleaner = Nodes_Cleaner();
    }

    /** Build a convex region from a collection of bounding hyperplanes.
     * @param hyperplanes collection of bounding hyperplanes
     * @return a convex region, or NULL if the collection is empty
     */
    //@Safe_Varargs
    public const Region<S> build_convex(const Hyperplane<S> ... hyperplanes) 
    {
        if ((hyperplanes == NULL) || (hyperplanes.size() == 0)) 
        {
            return NULL;
        }

        // use the first hyperplane to build the right class
        const Region<S> region = hyperplanes[0].whole_space();

        // chop off parts of the space
        BSP_Tree<S> node = region.get_tree(false);
        node.set_attribute(Boolean.TRUE);
        for (const Hyperplane<S> hyperplane : hyperplanes) 
        {
            if (node.insert_cut(hyperplane)) 
            {
                node.set_attribute(null);
                node.get_plus().set_attribute(Boolean.FALSE);
                node = node.get_minus();
                node.set_attribute(Boolean.TRUE);
            }
else 
            {
                // the hyperplane could not be inserted in the current leaf node
                // either it is completely outside (which means the input hyperplanes
                // are wrong), or it is parallel to a previous hyperplane
                Sub_Hyperplane<S> s = hyperplane.whole_hyperplane();
                for (BSP_Tree<S> tree = node; tree.get_parent() != NULL && s != NULL; tree = tree.get_parent()) 
                {
                    const Hyperplane<S>         other = tree.get_parent().get_cut().get_hyperplane();
                    const Split_Sub_Hyperplane<S> split = s.split(other);
                    switch (split.get_side()) 
                    {
                        case HYPER :
                            // the hyperplane is parallel to a previous hyperplane
                            if (!hyperplane.same_orientation_as(other)) 
                            {
                                // this hyperplane is opposite to the other one, // the region is thinner than the tolerance, we consider it empty
                                return get_complement(hyperplanes[0].whole_space());
                            }
                            // the hyperplane is an extension of an already known hyperplane, we just ignore it
                            break;
                        case PLUS :
                        // the hyperplane is outside of the current convex zone, // the input hyperplanes are inconsistent
                        throw (Localized_Geometry_Formats.NOT_CONVEX_HYPERPLANES);
                        default :
                            s = split.get_minus();
                    }
                }
            }
        }

        return region;

    }

    /** Compute the union of two regions.
     * @param region1 first region (will be unusable after the operation as
     * parts of it will be reused in the region)
     * @param region2 second region (will be unusable after the operation as
     * parts of it will be reused in the region)
     * @return a region, result of {@code region1 union region2}
     */
    public Region<S> union(const Region<S> region1, const Region<S> region2) 
    {
        const BSP_Tree<S> tree =
            region1.get_tree(false).merge(region2.get_tree(false), Union_Merger());
        tree.visit(node_cleaner);
        return region1.build_new(tree);
    }

    /** Compute the intersection of two regions.
     * @param region1 first region (will be unusable after the operation as
     * parts of it will be reused in the region)
     * @param region2 second region (will be unusable after the operation as
     * parts of it will be reused in the region)
     * @return a region, result of {@code region1 intersection region2}
     */
    public Region<S> intersection(const Region<S> region1, const Region<S> region2) 
    {
        const BSP_Tree<S> tree =
            region1.get_tree(false).merge(region2.get_tree(false), Intersection_Merger());
        tree.visit(node_cleaner);
        return region1.build_new(tree);
    }

    /** Compute the symmetric difference (exclusive or) of two regions.
     * @param region1 first region (will be unusable after the operation as
     * parts of it will be reused in the region)
     * @param region2 second region (will be unusable after the operation as
     * parts of it will be reused in the region)
     * @return a region, result of {@code region1 xor region2}
     */
    public Region<S> xor(const Region<S> region1, const Region<S> region2) 
    {
        const BSP_Tree<S> tree =
            region1.get_tree(false).merge(region2.get_tree(false), Xor_Merger());
        tree.visit(node_cleaner);
        return region1.build_new(tree);
    }

    /** Compute the difference of two regions.
     * @param region1 first region (will be unusable after the operation as
     * parts of it will be reused in the region)
     * @param region2 second region (will be unusable after the operation as
     * parts of it will be reused in the region)
     * @return a region, result of {@code region1 minus region2}
     */
    public Region<S> difference(const Region<S> region1, const Region<S> region2) 
    {
        const BSP_Tree<S> tree =
            region1.get_tree(false).merge(region2.get_tree(false), Difference_Merger(region1, region2));
        tree.visit(node_cleaner);
        return region1.build_new(tree);
    }

    /** Get the complement of the region (exchanged interior/exterior).
     * @param region region to complement, it will not modified, a new
     * region independent region will be built
     * @return a region, complement of the specified one
     */
    /** Get the complement of the region (exchanged interior/exterior).
     * @param region region to complement, it will not modified, a new
     * region independent region will be built
     * @return a region, complement of the specified one
     */
    public Region<S> get_complement(const Region<S> region) 
    {
        return region.build_new(recurse_complement(region.get_tree(false)));
    }

    /** Recursively build the complement of a BSP tree.
     * @param node current node of the original tree
     * @return tree, complement of the node
     */
    private BSP_Tree<S> recurse_complement(const BSP_Tree<S> node) 
    {

        // transform the tree, except for boundary attribute splitters
        const Map<BSP_Tree<S>, BSP_Tree<S>> map = Hash_Map<>();
        const BSP_Tree<S> transformed_tree = recurse_complement(node, map);

        // set up the boundary attributes splitters
        for (const Map.Entry<BSP_Tree<S>, BSP_Tree<S>> entry : map.entry_set()) 
        {
            if (entry.get_key().get_cut() != NULL) 
            {
                //@Suppress_Warnings("unchecked")
                Boundary_Attribute<S> original = (Boundary_Attribute<S>) entry.get_key().get_attribute();
                if (original != NULL) 
                {
                    //@Suppress_Warnings("unchecked")
                    Boundary_Attribute<S> transformed = (Boundary_Attribute<S>) entry.get_value().get_attribute();
                    for (const BSP_Tree<S> splitter : original.get_splitters()) 
                    {
                        transformed.get_splitters().add(map.get(splitter));
                    }
                }
            }
        }

        return transformed_tree;

    }

    /** Recursively build the complement of a BSP tree.
     * @param node current node of the original tree
     * @param map transformed nodes map
     * @return tree, complement of the node
     */
    private BSP_Tree<S> recurse_complement(const BSP_Tree<S> node, const Map<BSP_Tree<S>, BSP_Tree<S>> map) 
    {

        const BSP_Tree<S> transformed_node;
        if (node.get_cut() == NULL) 
        {
            transformed_node = BSP_Tree<>(((Boolean) node.get_attribute()) ? Boolean.FALSE : Boolean.TRUE);
        }
else 
        {

            //@Suppress_Warnings("unchecked")
            Boundary_Attribute<S> attribute = (Boundary_Attribute<S>) node.get_attribute();
            if (attribute != NULL) 
            {
                const Sub_Hyperplane<S> plus_outside =
                        (attribute.get_plus_inside() == NULL) ? NULL : attribute.get_plus_inside().copy_self();
                const Sub_Hyperplane<S> plus_inside  =
                        (attribute.get_plus_outside() == NULL) ? NULL : attribute.get_plus_outside().copy_self();
                // we start with an empty list of splitters, it will be filled in out of recursion
                attribute = Boundary_Attribute<>(plus_outside, plus_inside, Nodes_Set<S>());
            }

            transformed_node = BSP_Tree<>(node.get_cut().copy_self(), recurse_complement(node.get_plus(),  map), recurse_complement(node.get_minus(), map), attribute);
        }

        map.put(node, transformed_node);
        return transformed_node;

    }

    /** BSP tree leaf merger computing union of two regions. */
    private class Union_Merger : BSP_Tree.Leaf_Merger<S> 
    {
        /** {@inherit_doc} */
        //override
        public BSP_Tree<S> merge(const BSP_Tree<S> leaf, const BSP_Tree<S> tree, const BSP_Tree<S> parent_tree, const bool is_plus_child, const bool leaf_from_instance) 
        {
            if ((Boolean) leaf.get_attribute()) 
            {
                // the leaf node represents an inside cell
                leaf.insert_in_tree(parent_tree, is_plus_child, Vanishing_To_Leaf(true));
                return leaf;
            }
            // the leaf node represents an outside cell
            tree.insert_in_tree(parent_tree, is_plus_child, Vanishing_To_Leaf(false));
            return tree;
        }
    }

    /** BSP tree leaf merger computing intersection of two regions. */
    private class Intersection_Merger : BSP_Tree.Leaf_Merger<S> 
    {
        /** {@inherit_doc} */
        //override
        public BSP_Tree<S> merge(const BSP_Tree<S> leaf, const BSP_Tree<S> tree, const BSP_Tree<S> parent_tree, const bool is_plus_child, const bool leaf_from_instance) 
        {
            if ((Boolean) leaf.get_attribute()) 
            {
                // the leaf node represents an inside cell
                tree.insert_in_tree(parent_tree, is_plus_child, Vanishing_To_Leaf(true));
                return tree;
            }
            // the leaf node represents an outside cell
            leaf.insert_in_tree(parent_tree, is_plus_child, Vanishing_To_Leaf(false));
            return leaf;
        }
    }

    /** BSP tree leaf merger computing symmetric difference (exclusive or) of two regions. */
    private class Xor_Merger : BSP_Tree.Leaf_Merger<S> 
    {
        /** {@inherit_doc} */
        //override
        public BSP_Tree<S> merge(const BSP_Tree<S> leaf, const BSP_Tree<S> tree, const BSP_Tree<S> parent_tree, const bool is_plus_child, const bool leaf_from_instance) 
        {
            BSP_Tree<S> t = tree;
            if ((Boolean) leaf.get_attribute()) 
            {
                // the leaf node represents an inside cell
                t = recurse_complement(t);
            }
            t.insert_in_tree(parent_tree, is_plus_child, Vanishing_To_Leaf(true));
            return t;
        }
    }

    /** BSP tree leaf merger computing difference of two regions. */
    private class Difference_Merger : BSP_Tree.Leaf_Merger<S>, Vanishing_Cut_Handler<S> 
    {

        /** Region to subtract from. */
        private const Region<S> region1;

        /** Region to subtract. */
        private const Region<S> region2;

        /** Simple constructor.
         * @param region1 region to subtract from
         * @param region2 region to subtract
         */
        Difference_Merger(const Region<S> region1, const Region<S> region2) 
        {
            this.region1 = region1.copy_self();
            this.region2 = region2.copy_self();
        }

        /** {@inherit_doc} */
        //override
        public BSP_Tree<S> merge(const BSP_Tree<S> leaf, const BSP_Tree<S> tree, const BSP_Tree<S> parent_tree, const bool is_plus_child, const bool leaf_from_instance) 
        {
            if ((Boolean) leaf.get_attribute()) 
            {
                // the leaf node represents an inside cell
                const BSP_Tree<S> arg_tree =
                    recurse_complement(leaf_from_instance ? tree : leaf);
                arg_tree.insert_in_tree(parent_tree, is_plus_child, this);
                return arg_tree;
            }
            // the leaf node represents an outside cell
            const BSP_Tree<S> instance_tree =
                leaf_from_instance ? leaf : tree;
            instance_tree.insert_in_tree(parent_tree, is_plus_child, this);
            return instance_tree;
        }

        /** {@inherit_doc} */
        //override
        public BSP_Tree<S> fix_node(const BSP_Tree<S> node) 
        {
            // get a representative point in the degenerate cell
            const BSP_Tree<S> cell = node.prune_around_convex_cell(Boolean.TRUE, Boolean.FALSE, NULL);
            const Region<S> r = region1.build_new(cell);
            const Point<S> p = r.get_barycenter();
            return BSP_Tree<S>(region1.check_point(p) == Location.INSIDE &&
                                  region2.check_point(p) == Location.OUTSIDE);
        }

    }

    /** Visitor removing internal nodes attributes. */
    private class Nodes_Cleaner :  BSP_Tree_Visitor<S> 
    {

        /** {@inherit_doc} */
        //override
        public Order visit_order(const BSP_Tree<S> node) 
        {
            return Order.PLUS_SUB_MINUS;
        }

        /** {@inherit_doc} */
        //override
        public void visit_internal_node(const BSP_Tree<S> node) 
        {
            node.set_attribute(null);
        }

        /** {@inherit_doc} */
        //override
        public void visit_leaf_node(const BSP_Tree<S> node) 
        {
        }

    }

    /** Handler replacing nodes with vanishing cuts with leaf nodes. */
    private class Vanishing_To_Leaf : Vanishing_Cut_Handler<S> 
    {

        /** Inside/outside indocator to use for ambiguous nodes. */
        private const bool inside;

        /** Simple constructor.
         * @param inside inside/outside indicator to use for ambiguous nodes
         */
        Vanishing_To_Leaf(const bool inside) 
        {
            this.inside = inside;
        }

        /** {@inherit_doc} */
        //override
        public BSP_Tree<S> fix_node(const BSP_Tree<S> node) 
        {
            if (node.get_plus().get_attribute().equals(node.get_minus().get_attribute())) 
            {
                // no ambiguity
                return BSP_Tree<S>(node.get_plus().get_attribute());
            }
else 
            {
                // ambiguous node
                return BSP_Tree<S>(inside);
            }
        }

    }

}


