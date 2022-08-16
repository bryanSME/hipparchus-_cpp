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
  //import java.util.Collection;
  //import java.util.Comparator;
  //import java.util.Hash_Map;
  //import java.util.Iterator;
  //import java.util.Map;
  //import java.util.Tree_Set;

  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.Space;
  //import org.hipparchus.geometry.Vector;

  /** Abstract class for all regions, independently of geometry type or dimension.

   * @param <S> Type of the space.
   * @param <T> Type of the sub-space.

   */
class Abstract_Region<S extends Space, T extends Space> : Region<S>
{
	/** Inside/Outside BSP tree. */
	private BSP_Tree<S> tree;

	/** Tolerance below which points are considered to belong to hyperplanes. */
	private const double& tolerance;

	/** Size of the instance. */
	private double size;

	/** Barycenter. */
	private Point<S> barycenter;

	/** Build a region representing the whole space.
	 * @param tolerance tolerance below which points are considered identical.
	 */
	protected Abstract_Region(const double& tolerance)
	{
		this.tree = BSP_Tree<>(Boolean.TRUE);
		this.tolerance = tolerance;
	}

	/** Build a region from an inside/outside BSP tree.
	 * <p>The leaf nodes of the BSP tree <em>must</em> have a
	 * {@code Boolean} attribute representing the inside status of
	 * the corresponding cell (true for inside cells, false for outside
	 * cells). In order to avoid building too many small objects, it is
	 * recommended to use the predefined constants
	 * {@code Boolean.TRUE} and {@code Boolean.FALSE}. The
	 * tree also <em>must</em> have either NULL internal nodes or
	 * internal nodes representing the boundary as specified in the
	 * {@link #get_tree get_tree} method).</p>
	 * @param tree inside/outside BSP tree representing the region
	 * @param tolerance tolerance below which points are considered identical.
	 */
	protected Abstract_Region(const BSP_Tree<S> tree, const double& tolerance)
	{
		this.tree = tree;
		this.tolerance = tolerance;
	}

	/** Build a Region from a Boundary RE_Presentation (B-rep).
	 * <p>The boundary is provided as a collection of {@link
	 * Sub_Hyperplane sub-hyperplanes}. Each sub-hyperplane has the
	 * interior part of the region on its minus side and the exterior on
	 * its plus side.</p>
	 * <p>The boundary elements can be in any order, and can form
	 * several non-connected sets (like for example polygons with holes
	 * or a set of disjoints polyhedrons considered as a whole). In
	 * fact, the elements do not even need to be connected together
	 * (their topological connections are not used here). However, if the
	 * boundary does not really separate an inside open from an outside
	 * open (open having here its topological meaning), then subsequent
	 * calls to the {@link #check_point(Point) check_point} method will not be
	 * meaningful anymore.</p>
	 * <p>If the boundary is empty, the region will represent the whole
	 * space.</p>
	 * @param boundary collection of boundary elements, as a
	 * collection of {@link Sub_Hyperplane Sub_Hyperplane} objects
	 * @param tolerance tolerance below which points are considered identical.
	 */
	protected Abstract_Region(const Collection<Sub_Hyperplane<S>> boundary, const double& tolerance)
	{
		this.tolerance = tolerance;

		if (boundary.is_empty())
		{
			// the tree represents the whole space
			tree = BSP_Tree<>(Boolean.TRUE);
		}
		else
		{
			// sort the boundary elements in decreasing size order
			// (we don't want equal size elements to be removed, so
			// we use a trick to fool the Tree_Set)
			const Tree_Set<Sub_Hyperplane<S>> ordered = Tree_Set<>(new Comparator<Sub_Hyperplane<S>>()
			{
				/** {@inherit_doc} */
				//override
				public int compare(const Sub_Hyperplane<S> o1, const Sub_Hyperplane<S> o2)
				{
					const double size1 = o1.get_size();
					const double size2 = o2.get_size();
					return (size2 < size1) ? -1 : ((o1 == o2) ? 0 : +1);
				}
			});
			ordered.add_all(boundary);

			// build the tree top-down
			tree = BSP_Tree<>();
			insert_cuts(tree, ordered);

			// set up the inside/outside flags
			tree.visit(new BSP_Tree_Visitor<S>()
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
				}

				/** {@inherit_doc} */
				//override
				public void visit_leaf_node(const BSP_Tree<S> node)
				{
					if (node.get_parent() == NULL || node == node.get_parent().get_minus())
					{
						node.set_attribute(Boolean.TRUE);
					}
					else
					{
						node.set_attribute(Boolean.FALSE);
					}
				}
			});
		}
	}

	/** Build a convex region from an array of bounding hyperplanes.
	 * @param hyperplanes array of bounding hyperplanes (if NULL, an
	 * empty region will be built)
	 * @param tolerance tolerance below which points are considered identical.
	 */
	public Abstract_Region(const Hyperplane<S>[] hyperplanes, const double& tolerance)
	{
		this.tolerance = tolerance;
		if ((hyperplanes == NULL) || (hyperplanes.size() == 0))
		{
			tree = BSP_Tree<>(Boolean.FALSE);
		}
		else
		{
			// use the first hyperplane to build the right class
			tree = hyperplanes[0].whole_space().get_tree(false);

			// chop off parts of the space
			BSP_Tree<S> node = tree;
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
			}
		}
	}

	/** {@inherit_doc} */
	//override
	public virtual Abstract_Region<S, T> build_new(BSP_Tree<S> new_tree);

	/** Get the tolerance below which points are considered to belong to hyperplanes.
	 * @return tolerance below which points are considered to belong to hyperplanes
	 */
	public double get_tolerance()
	{
		return tolerance;
	}

	/** Recursively build a tree by inserting cut sub-hyperplanes.
	 * @param node current tree node (it is a leaf node at the beginning
	 * of the call)
	 * @param boundary collection of edges belonging to the cell defined
	 * by the node
	 */
	private void insert_cuts(const BSP_Tree<S> node, const Collection<Sub_Hyperplane<S>> boundary)
	{
		const Iterator<Sub_Hyperplane<S>> iterator = boundary.iterator();

		// build the current level
		Hyperplane<S> inserted = NULL;
		while ((inserted == NULL) && iterator.has_next())
		{
			inserted = iterator.next().get_hyperplane();
			if (!node.insert_cut(inserted.copy_self()))
			{
				inserted = NULL;
			}
		}

		if (!iterator.has_next())
		{
			return;
		}

		// distribute the remaining edges in the two sub-trees
		const Array_list<Sub_Hyperplane<S>> plus_list = Array_list<>();
		const Array_list<Sub_Hyperplane<S>> minus_list = Array_list<>();
		while (iterator.has_next())
		{
			const Sub_Hyperplane<S> other = iterator.next();
			const Sub_Hyperplane.Split_Sub_Hyperplane<S> split = other.split(inserted);
			switch (split.get_side())
			{
			case PLUS:
				plus_list.add(other);
				break;
			case MINUS:
				minus_list.add(other);
				break;
			case BOTH:
				plus_list.add(split.get_plus());
				minus_list.add(split.get_minus());
				break;
			default:
				// ignore the sub-hyperplanes belonging to the cut hyperplane
			}
		}

		// recurse through lower levels
		insert_cuts(node.get_plus(), plus_list);
		insert_cuts(node.get_minus(), minus_list);
	}

	/** {@inherit_doc} */
	//override
	public Abstract_Region<S, T> copy_self()
	{
		return build_new(tree.copy_self());
	}

	/** {@inherit_doc} */
	//override
	public bool is_empty()
	{
		return is_empty(tree);
	}

	/** {@inherit_doc} */
	//override
	public bool is_empty(const BSP_Tree<S> node)
	{
		// we use a recursive function rather than the BSP_Tree_Visitor
		// interface because we can stop visiting the tree as soon as we
		// have found an inside cell

		if (node.get_cut() == NULL)
		{
			// if we find an inside node, the region is not empty
			return !((Boolean)node.get_attribute());
		}

		// check both sides of the sub-tree
		return is_empty(node.get_minus()) && is_empty(node.get_plus());
	}

	/** {@inherit_doc} */
	//override
	public bool is_full()
	{
		return is_full(tree);
	}

	/** {@inherit_doc} */
	//override
	public bool is_full(const BSP_Tree<S> node)
	{
		// we use a recursive function rather than the BSP_Tree_Visitor
		// interface because we can stop visiting the tree as soon as we
		// have found an outside cell

		if (node.get_cut() == NULL)
		{
			// if we find an outside node, the region does not cover full space
			return (Boolean)node.get_attribute();
		}

		// check both sides of the sub-tree
		return is_full(node.get_minus()) && is_full(node.get_plus());
	}

	/** {@inherit_doc} */
	//override
	public bool contains(const Region<S> region)
	{
		return Region_Factory<S>().difference(region, this).is_empty();
	}

	/** {@inherit_doc}
	 */
	 //override
	public Boundary_Projection<S> project_to_boundary(const Point<S> point)
	{
		const Boundary_Projector<S, T> projector = Boundary_Projector<>(point);
		get_tree(true).visit(projector);
		return projector.get_projection();
	}

	/** Check a point with respect to the region.
	 * @param point point to check
	 * @return a code representing the point status: either {@link
	 * Region.Location#INSIDE}, {@link Region.Location#OUTSIDE} or
	 * {@link Region.Location#BOUNDARY}
	 */
	public Location check_point(const Vector<S> point)
	{
		return check_point((Point<S>) point);
	}

	/** {@inherit_doc} */
	//override
	public Location check_point(const Point<S> point)
	{
		return check_point(tree, point);
	}

	/** Check a point with respect to the region starting at a given node.
	 * @param node root node of the region
	 * @param point point to check
	 * @return a code representing the point status: either {@link
	 * Region.Location#INSIDE INSIDE}, {@link Region.Location#OUTSIDE
	 * OUTSIDE} or {@link Region.Location#BOUNDARY BOUNDARY}
	 */
	protected Location check_point(const BSP_Tree<S> node, const Vector<S> point)
	{
		return check_point(node, (Point<S>) point);
	}

	/** Check a point with respect to the region starting at a given node.
	 * @param node root node of the region
	 * @param point point to check
	 * @return a code representing the point status: either {@link
	 * Region.Location#INSIDE INSIDE}, {@link Region.Location#OUTSIDE
	 * OUTSIDE} or {@link Region.Location#BOUNDARY BOUNDARY}
	 */
	protected Location check_point(const BSP_Tree<S> node, const Point<S> point)
	{
		const BSP_Tree<S> cell = node.get_cell(point, tolerance);
		if (cell.get_cut() == NULL)
		{
			// the point is in the interior of a cell, just check the attribute
			return ((Boolean)cell.get_attribute()) ? Location.INSIDE : Location.OUTSIDE;
		}

		// the point is on a cut-sub-hyperplane, is it on a boundary ?
		const Location minus_code = check_point(cell.get_minus(), point);
		const Location plus_code = check_point(cell.get_plus(), point);
		return (minus_code == plus_code) ? minus_code : Location.BOUNDARY;
	}

	/** {@inherit_doc} */
	//override
	public BSP_Tree<S> get_tree(const bool include_boundary_attributes)
	{
		if (include_boundary_attributes && (tree.get_cut() != NULL) && (tree.get_attribute() == NULL))
		{
			// compute the boundary attributes
			tree.visit(new Boundary_Builder<S>());
		}
		return tree;
	}

	/** {@inherit_doc} */
	//override
	public double get_boundary_size()
	{
		const Boundary_size_visitor<S> visitor = Boundary_size_visitor<>();
		get_tree(true).visit(visitor);
		return visitor.get_size();
	}

	/** {@inherit_doc} */
	//override
	public double get_size()
	{
		if (barycenter == NULL)
		{
			compute_geometrical_properties();
		}
		return size;
	}

	/** Set the size of the instance.
	 * @param size size of the instance
	 */
	protected void set_size(const double size)
	{
		this.size = size;
	}

	/** {@inherit_doc} */
	//override
	public Point<S> get_barycenter()
	{
		if (barycenter == NULL)
		{
			compute_geometrical_properties();
		}
		return barycenter;
	}

	/** Set the barycenter of the instance.
	 * @param barycenter barycenter of the instance
	 */
	protected void set_barycenter(const Vector<S> barycenter)
	{
		set_barycenter((Point<S>) barycenter);
	}

	/** Set the barycenter of the instance.
	 * @param barycenter barycenter of the instance
	 */
	protected void set_barycenter(const Point<S> barycenter)
	{
		this.barycenter = barycenter;
	}

	/** Compute some geometrical properties.
	 * <p>The properties to compute are the barycenter and the size.</p>
	 */
	protected virtual void compute_geometrical_properties();

	/** {@inherit_doc} */
	//override
	public Sub_Hyperplane<S> intersection(const Sub_Hyperplane<S> sub)
	{
		return recurse_intersection(tree, sub);
	}

	/** Recursively compute the parts of a sub-hyperplane that are
	 * contained in the region.
	 * @param node current BSP tree node
	 * @param sub sub-hyperplane traversing the region
	 * @return filtered sub-hyperplane
	 */
	private Sub_Hyperplane<S> recurse_intersection(const BSP_Tree<S> node, const Sub_Hyperplane<S> sub)
	{
		if (node.get_cut() == NULL)
		{
			return (Boolean)node.get_attribute() ? sub.copy_self() : NULL;
		}

		const Hyperplane<S> hyperplane = node.get_cut().get_hyperplane();
		const Sub_Hyperplane.Split_Sub_Hyperplane<S> split = sub.split(hyperplane);
		if (split.get_plus() != NULL)
		{
			if (split.get_minus() != NULL)
			{
				// both sides
				const Sub_Hyperplane<S> plus = recurse_intersection(node.get_plus(), split.get_plus());
				const Sub_Hyperplane<S> minus = recurse_intersection(node.get_minus(), split.get_minus());
				if (plus == NULL)
				{
					return minus;
				}
				else if (minus == NULL)
				{
					return plus;
				}
				else
				{
					return plus.reunite(minus);
				}
			}
			else
			{
				// only on plus side
				return recurse_intersection(node.get_plus(), sub);
			}
		}
		else if (split.get_minus() != NULL)
		{
			// only on minus side
			return recurse_intersection(node.get_minus(), sub);
		}
		else
		{
			// on hyperplane
			return recurse_intersection(node.get_plus(), recurse_intersection(node.get_minus(), sub));
		}
	}

	/** Transform a region.
	 * <p>Applying a transform to a region consist in applying the
	 * transform to all the hyperplanes of the underlying BSP tree and
	 * of the boundary (and also to the sub-hyperplanes embedded in
	 * these hyperplanes) and to the barycenter. The instance is not
	 * modified, a instance is built.</p>
	 * @param transform transform to apply
	 * @return a region, resulting from the application of the
	 * transform to the instance
	 */
	public Abstract_Region<S, T> apply_transform(const Transform<S, T> transform)
	{
		// transform the tree, except for boundary attribute splitters
		const Map<BSP_Tree<S>, BSP_Tree<S>> map = Hash_Map<>();
		const BSP_Tree<S> transformed_tree = recurse_transform(get_tree(false), transform, map);

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

		return build_new(transformed_tree);
	}

	/** Recursively transform an inside/outside BSP-tree.
	 * @param node current BSP tree node
	 * @param transform transform to apply
	 * @param map transformed nodes map
	 * @return a tree
	 */
	 //@Suppress_Warnings("unchecked")
	private BSP_Tree<S> recurse_transform(const BSP_Tree<S> node, const Transform<S, T> transform, const Map<BSP_Tree<S>, BSP_Tree<S>> map)
	{
		const BSP_Tree<S> transformed_node;
		if (node.get_cut() == NULL)
		{
			transformed_node = BSP_Tree<>(node.get_attribute());
		}
		else
		{
			const Sub_Hyperplane<S>  sub = node.get_cut();
			const Sub_Hyperplane<S> t_sub = ((Abstract_Sub_Hyperplane<S, T>) sub).apply_transform(transform);
			Boundary_Attribute<S> attribute = (Boundary_Attribute<S>) node.get_attribute();
			if (attribute != NULL)
			{
				const Sub_Hyperplane<S> tPO = (attribute.get_plus_outside() == NULL) ?
					NULL : ((Abstract_Sub_Hyperplane<S, T>) attribute.get_plus_outside()).apply_transform(transform);
				const Sub_Hyperplane<S> tPI = (attribute.get_plus_inside() == NULL) ?
					NULL : ((Abstract_Sub_Hyperplane<S, T>) attribute.get_plus_inside()).apply_transform(transform);
				// we start with an empty list of splitters, it will be filled in out of recursion
				attribute = Boundary_Attribute<>(tPO, tPI, Nodes_Set<>());
			}

			transformed_node = BSP_Tree<>(t_sub, recurse_transform(node.get_plus(), transform, map), recurse_transform(node.get_minus(), transform, map), attribute);
		}

		map.put(node, transformed_node);
		return transformed_node;
	}
}
