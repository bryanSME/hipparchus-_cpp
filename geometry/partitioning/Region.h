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

  //import org.hipparchus.geometry.Point;
#include "../Space.h"
#include <type_traits>

/** This interface represents a region of a space as a partition.

 * <p>Region are subsets of a space, they can be infinite (whole
 * space, half space, infinite stripe ...) or finite (polygons in 2D, * polyhedrons in 3D ...). Their main characteristic is to separate
 * points that are considered to be <em>inside</em> the region from
 * points considered to be <em>outside</em> of it. In between, there
 * may be points on the <em>boundary</em> of the region.</p>

 * <p>This implementation is limited to regions for which the boundary
 * is composed of several {@link Sub_Hyperplane sub-hyperplanes}, * including regions with no boundary at all: the whole space and the
 * empty region. They are not necessarily finite and not necessarily
 * path-connected. They can contain holes.</p>

 * <p>Regions can be combined using the traditional sets operations :
 * union, intersection, difference and symetric difference (exclusive
 * or) for the binary operations, complement for the unary
 * operation.</p>

 * <p>
 * Note that this interface is <em>not</em> intended to be implemented
 * by Hipparchus users, it is only intended to be implemented
 * within the library itself. New methods may be added even for minor
 * versions, which breaks compatibility for external implementations.
 * </p>

 * @param <S> Type of the space.

 */
template<typename S, typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Region
{
	/** Enumerate for the location of a point with respect to the region. */
	enum Location
	{
		/** Code for points inside the partition. */
		INSIDE,
		/** Code for points outside of the partition. */
		OUTSIDE,
		/** Code for points on the partition boundary. */
		BOUNDARY;
	}

	/** Build a region using the instance as a prototype.
	 * <p>This method allow to create instances without knowing
	 * exactly the type of the region. It is an application of the
	 * prototype design pattern.</p>
	 * <p>The leaf nodes of the BSP tree <em>must</em> have a
	 * {@code Boolean} attribute representing the inside status of
	 * the corresponding cell (true for inside cells, false for outside
	 * cells). In order to avoid building too many small objects, it is
	 * recommended to use the predefined constants
	 * {@code Boolean.TRUE} and {@code Boolean.FALSE}. The
	 * tree also <em>must</em> have either NULL internal nodes or
	 * internal nodes representing the boundary as specified in the
	 * {@link #get_tree get_tree} method).</p>
	 * @param new_tree inside/outside BSP tree representing the region
	 * @return the built region
	 */
	virtual Region<S> build_new(BSP_Tree<S> new_tree) = 0;

	/** Copy the instance.
	 * <p>The instance created is completely independant of the original
	 * one. A deep copy is used, none of the underlying objects are
	 * shared (except for the underlying tree {@code Boolean}
	 * attributes and immutable objects).</p>
	 * @return a region, copy of the instance
	 */
	virtual Region<S> copy_self() = 0;

	/** Check if the instance is empty.
	 * @return true if the instance is empty
	 */
	virtual bool is_empty() = 0;

	/** Check if the sub-tree starting at a given node is empty.
	 * @param node root node of the sub-tree (<em>must</em> have {@link
	 * Region Region} tree semantics, i.e. the leaf nodes must have
	 * {@code Boolean} attributes representing an inside/outside
	 * property)
	 * @return true if the sub-tree starting at the given node is empty
	 */
	virtual bool is_empty(BSP_Tree<S> node) = 0;

	/** Check if the instance covers the full space.
	 * @return true if the instance covers the full space
	 */
	virtual bool is_full() = 0;

	/** Check if the sub-tree starting at a given node covers the full space.
	 * @param node root node of the sub-tree (<em>must</em> have {@link
	 * Region Region} tree semantics, i.e. the leaf nodes must have
	 * {@code Boolean} attributes representing an inside/outside
	 * property)
	 * @return true if the sub-tree starting at the given node covers the full space
	 */
	virtual bool is_full(BSP_Tree<S> node) = 0;

	/** Check if the instance entirely contains another region.
	 * @param region region to check against the instance
	 * @return true if the instance contains the specified tree
	 */
	bool contains(Region<S> region);

	/** Check a point with respect to the region.
	 * @param point point to check
	 * @return a code representing the point status: either {@link
	 * Location#INSIDE}, {@link Location#OUTSIDE} or {@link Location#BOUNDARY}
	 */
	virtual Location check_point(Point<S> point) = 0;

	/** Project a point on the boundary of the region.
	 * @param point point to check
	 * @return projection of the point on the boundary
	 */
	virtual Boundary_Projection<S> project_to_boundary(Point<S> point) = 0;

	/** Get the underlying BSP tree.

	 * <p>Regions are represented by an underlying inside/outside BSP
	 * tree whose leaf attributes are {@code Boolean} instances
	 * representing inside leaf cells if the attribute value is
	 * {@code true} and outside leaf cells if the attribute is
	 * {@code false}. These leaf attributes are always present and
	 * guaranteed to be non NULL.</p>

	 * <p>In addition to the leaf attributes, the internal nodes which
	 * correspond to cells split by cut sub-hyperplanes may contain
	 * {@link Boundary_Attribute Boundary_Attribute} objects representing
	 * the parts of the corresponding cut sub-hyperplane that belong to
	 * the boundary. When the boundary attributes have been computed, * all internal nodes are guaranteed to have non-null
	 * attributes, however some {@link Boundary_Attribute
	 * Boundary_Attribute} instances may have their {@link
	 * Boundary_Attribute#get_plus_inside() get_plus_inside} and {@link
	 * Boundary_Attribute#get_plus_outside() get_plus_outside} methods both
	 * returning NULL if the corresponding cut sub-hyperplane does not
	 * have any parts belonging to the boundary.</p>

	 * <p>sin_ce computing the boundary is not always required and can be
	 * time-consuming for large trees, these internal nodes attributes
	 * are computed using lazy evaluation only when required by setting
	 * the {@code include_boundary_attributes} argument to
	 * {@code true}. Once computed, these attributes remain in the
	 * tree, which implies that in this case, further calls to the
	 * method for the same region will always include these attributes
	 * regardless of the value of the
	 * {@code include_boundary_attributes} argument.</p>

	 * @param include_boundary_attributes if true, the boundary attributes
	 * at internal nodes are guaranteed to be included (they may be
	 * included even if the argument is false, if they have already been
	 * computed due to a previous call)
	 * @return underlying BSP tree
	 * @see Boundary_Attribute
	 */
	virtual BSP_Tree<S> get_tree(bool include_boundary_attributes) = 0;

	/** Get the size of the boundary.
	 * @return the size of the boundary (this is 0 in 1D, a length in
	 * 2D, an area in 3D ...)
	 */
	virtual double get_boundary_size() = 0;

	/** Get the size of the instance.
	 * @return the size of the instance (this is a length in 1D, an area
	 * in 2D, a volume in 3D ...)
	 */
	virtual double get_size() = 0;

	/** Get the barycenter of the instance.
	 * @return an object representing the barycenter
	 */
	virtual Point<S> get_barycenter() = 0;

	/** Get the parts of a sub-hyperplane that are contained in the region.
	 * <p>The parts of the sub-hyperplane that belong to the boundary are
	 * <em>not</em> included in the resulting parts.</p>
	 * @param sub sub-hyperplane traversing the region
	 * @return filtered sub-hyperplane
	 */
	virtual Sub_Hyperplane<S> intersection(Sub_Hyperplane<S> sub) = 0;
};