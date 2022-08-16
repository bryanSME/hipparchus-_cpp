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

  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.Space;
  //import org.hipparchus.util.FastMath;

#include "../Space.h"
#include <type_traits>

/** This class represent a Binary Space Partition tree.

 * <p>BSP trees are an efficient way to represent space partitions and
 * to associate attributes with each cell. Each node in a BSP tree
 * represents a convex region which is partitioned in two convex
 * sub-regions at each side of a cut hyperplane. The root tree
 * contains the complete space.</p>

 * <p>The main use of such partitions is to use a bool attribute to
 * define an inside/outside property, hence representing arbitrary
 * polytopes (line segments in 1D, polygons in 2D and polyhedrons in
 * 3D) and to operate on them.</p>

 * <p>Another example would be to represent Voronoi tesselations, the
 * attribute of each cell holding the defining point of the cell.</p>

 * <p>The application-defined attributes are shared among copied
 * instances and propagated to split parts. These attributes are not
 * used by the BSP-tree algorithms themselves, so the application can
 * use them for any purpose. sin_ce the tree visiting method holds
 * internal and leaf nodes differently, it is possible to use
 * different classes for internal nodes attributes and leaf nodes
 * attributes. This should be used with care, though, because if the
 * tree is modified in any way after attributes have been set, some
 * internal nodes may become leaf nodes and some leaf nodes may become
 * internal nodes.</p>

 * <p>One of the main sources for the development of this //package was
 * Bruce Naylor, John Amanatides and William Thibault paper <a
 * href="http://www.cs.yorku.ca/~amana/research/bsptSetOp.pdf">Merging
 * BSP Trees Yields Polyhedral Set Operations</a> Proc. Siggraph '90, * Computer Graphics 24(4), August 1990, pp 115-124, published by the
 * Association for Computing Machinery (ACM).</p>

 * @param <S> Type of the space.

 */
template<typename S, typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class BSP_Tree
{
	/** Cut sub-hyperplane. */
	private Sub_Hyperplane<S> cut;

	/** Tree at the plus side of the cut hyperplane. */
	private BSP_Tree<S> plus;

	/** Tree at the minus side of the cut hyperplane. */
	private BSP_Tree<S> minus;

	/** Parent tree. */
	private BSP_Tree<S> parent;

	/** Application-defined attribute. */
	private Object attribute;

	/** Build a tree having only one root cell representing the whole space.
	 */
	public BSP_Tree()
	{
		cut = NULL;
		plus = NULL;
		minus = NULL;
		parent = NULL;
		attribute = NULL;
	}

	/** Build a tree having only one root cell representing the whole space.
	 * @param attribute attribute of the tree (may be NULL)
	 */
	public BSP_Tree(const Object attribute)
	{
		cut = NULL;
		plus = NULL;
		minus = NULL;
		parent = NULL;
		this.attribute = attribute;
	}

	/** Build a BSP_Tree from its underlying elements.
	 * <p>This method does <em>not</em> perform any verification on
	 * consistency of its arguments, it should therefore be used only
	 * when then caller knows what it is doing.</p>
	 * <p>This method is mainly useful to build trees
	 * bottom-up. Building trees top-down is realized with the help of
	 * method {@link #insert_cut insert_cut}.</p>
	 * @param cut cut sub-hyperplane for the tree
	 * @param plus plus side sub-tree
	 * @param minus minus side sub-tree
	 * @param attribute attribute associated with the node (may be NULL)
	 * @see #insert_cut
	 */
	public BSP_Tree(const Sub_Hyperplane<S> cut, const BSP_Tree<S> plus, const BSP_Tree<S> minus, const Object attribute)
	{
		this.cut = cut;
		this.plus = plus;
		this.minus = minus;
		this.parent = NULL;
		this.attribute = attribute;
		plus.parent = this;
		minus.parent = this;
	}

	/** Insert a cut sub-hyperplane in a node.
	 * <p>The sub-tree starting at this node will be completely
	 * overwritten. The cut sub-hyperplane will be built from the
	 * intersection of the provided hyperplane with the cell. If the
	 * hyperplane does intersect the cell, the cell will have two
	 * children cells with {@code NULL} attributes on each side of
	 * the inserted cut sub-hyperplane. If the hyperplane does not
	 * intersect the cell then <em>no</em> cut hyperplane will be
	 * inserted and the cell will be changed to a leaf cell. The
	 * attribute of the node is never changed.</p>
	 * <p>This method is mainly useful when called on leaf nodes
	 * (i.e. nodes for which {@link #get_cut get_cut} returns
	 * {@code NULL}), in this case it provides a way to build a
	 * tree top-down (whereas the {@link #BSP_Tree(Sub_Hyperplane, * BSP_Tree, BSP_Tree, Object) 4 arguments constructor} is devoted to
	 * build trees bottom-up).</p>
	 * @param hyperplane hyperplane to insert, it will be chopped in
	 * order to fit in the cell defined by the parent nodes of the
	 * instance
	 * @return true if a cut sub-hyperplane has been inserted (i.e. if
	 * the cell now has two leaf child nodes)
	 * @see #BSP_Tree(Sub_Hyperplane, BSP_Tree, BSP_Tree, Object)
	 */
	public bool insert_cut(const Hyperplane<S> hyperplane)
	{
		if (cut != NULL)
		{
			plus.parent = NULL;
			minus.parent = NULL;
		}

		const Sub_Hyperplane<S> chopped = fit_to_cell(hyperplane.whole_hyperplane());
		if (chopped == NULL || chopped.is_empty())
		{
			cut = NULL;
			plus = NULL;
			minus = NULL;
			return false;
		}

		cut = chopped;
		plus = BSP_Tree<>();
		plus.parent = this;
		minus = BSP_Tree<>();
		minus.parent = this;
		return true;
	}

	/** Copy the instance.
	 * <p>The instance created is completely independent of the original
	 * one. A deep copy is used, none of the underlying objects are
	 * shared (except for the nodes attributes and immutable
	 * objects).</p>
	 * @return a tree, copy of the instance
	 */
	public BSP_Tree<S> copy_self()
	{
		if (cut == NULL)
		{
			return BSP_Tree<S>(attribute);
		}

		return BSP_Tree<S>(cut.copy_self(), plus.copy_self(), minus.copy_self(), attribute);
	}

	/** Get the cut sub-hyperplane.
	 * @return cut sub-hyperplane, NULL if this is a leaf tree
	 */
	public Sub_Hyperplane<S> get_cut()
	{
		return cut;
	}

	/** Get the tree on the plus side of the cut hyperplane.
	 * @return tree on the plus side of the cut hyperplane, NULL if this
	 * is a leaf tree
	 */
	public BSP_Tree<S> get_plus()
	{
		return plus;
	}

	/** Get the tree on the minus side of the cut hyperplane.
	 * @return tree on the minus side of the cut hyperplane, NULL if this
	 * is a leaf tree
	 */
	public BSP_Tree<S> get_minus()
	{
		return minus;
	}

	/** Get the parent node.
	 * @return parent node, NULL if the node has no parents
	 */
	public BSP_Tree<S> get_parent()
	{
		return parent;
	}

	/** Associate an attribute with the instance.
	 * @param attribute attribute to associate with the node
	 * @see #get_attribute
	 */
	public void set_attribute(const Object attribute)
	{
		this.attribute = attribute;
	}

	/** Get the attribute associated with the instance.
	 * @return attribute associated with the node or NULL if no
	 * attribute has been explicitly set using the {@link #set_attribute
	 * set_attribute} method
	 * @see #set_attribute
	 */
	public Object get_attribute()
	{
		return attribute;
	}

	/** Visit the BSP tree nodes.
	 * @param visitor object visiting the tree nodes
	 */
	public void visit(const BSP_Tree_Visitor<S> visitor)
	{
		if (cut == NULL)
		{
			visitor.visit_leaf_node(this);
		}
		else
		{
			switch (visitor.visit_order(this))
			{
			case PLUS_MINUS_SUB:
				plus.visit(visitor);
				minus.visit(visitor);
				visitor.visit_internal_node(this);
				break;
			case PLUS_SUB_MINUS:
				plus.visit(visitor);
				visitor.visit_internal_node(this);
				minus.visit(visitor);
				break;
			case MINUS_PLUS_SUB:
				minus.visit(visitor);
				plus.visit(visitor);
				visitor.visit_internal_node(this);
				break;
			case MINUS_SUB_PLUS:
				minus.visit(visitor);
				visitor.visit_internal_node(this);
				plus.visit(visitor);
				break;
			case SUB_PLUS_MINUS:
				visitor.visit_internal_node(this);
				plus.visit(visitor);
				minus.visit(visitor);
				break;
			case SUB_MINUS_PLUS:
				visitor.visit_internal_node(this);
				minus.visit(visitor);
				plus.visit(visitor);
				break;
			default:
				throw Math_Runtime_Exception.create_internal_error();
			}
		}
	}

	/** Fit a sub-hyperplane inside the cell defined by the instance.
	 * <p>Fitting is done by chopping off the parts of the
	 * sub-hyperplane that lie outside of the cell using the
	 * cut-hyperplanes of the parent nodes of the instance.</p>
	 * @param sub sub-hyperplane to fit
	 * @return a sub-hyperplane, guaranteed to have no part outside
	 * of the instance cell
	 */
	private Sub_Hyperplane<S> fit_to_cell(const Sub_Hyperplane<S> sub)
	{
		Sub_Hyperplane<S> s = sub;
		for (BSP_Tree<S> tree = this; tree.parent != NULL && s != NULL; tree = tree.parent)
		{
			if (tree == tree.parent.plus)
			{
				s = s.split(tree.parent.cut.get_hyperplane()).get_plus();
			}
			else
			{
				s = s.split(tree.parent.cut.get_hyperplane()).get_minus();
			}
		}
		return s;
	}

	/** Get the cell to which a point belongs.
	 * <p>If the returned cell is a leaf node the points belongs to the
	 * interior of the node, if the cell is an internal node the points
	 * belongs to the node cut sub-hyperplane.</p>
	 * @param point point to check
	 * @param tolerance tolerance below which points close to a cut hyperplane
	 * are considered to belong to the hyperplane itself
	 * @return the tree cell to which the point belongs
	 */
	public BSP_Tree<S> get_cell(const Point<S> point, const double& tolerance)
	{
		if (cut == NULL)
		{
			return this;
		}

		// position of the point with respect to the cut hyperplane
		const double offset = cut.get_hyperplane().get_offset(point);

		if (std::abs(offset) < tolerance)
		{
			return this;
		}
		else if (offset <= 0)
		{
			// point is on the minus side of the cut hyperplane
			return minus.get_cell(point, tolerance);
		}
		else
		{
			// point is on the plus side of the cut hyperplane
			return plus.get_cell(point, tolerance);
		}
	}

	/** Get the cells whose cut sub-hyperplanes are close to the point.
	 * @param point point to check
	 * @param max_offset offset below which a cut sub-hyperplane is considered
	 * close to the point (in absolute value)
	 * @return close cells (may be empty if all cut sub-hyperplanes are farther
	 * than max_offset from the point)
	 */
	public List<BSP_Tree<S>> get_close_cuts(const Point<S> point, const double max_offset)
	{
		const List<BSP_Tree<S>> close = Array_list<>();
		recurse_close_cuts(point, max_offset, close);
		return close;
	}

	/** Get the cells whose cut sub-hyperplanes are close to the point.
	 * @param point point to check
	 * @param max_offset offset below which a cut sub-hyperplane is considered
	 * close to the point (in absolute value)
	 * @param close list to fill
	 */
	private void recurse_close_cuts(const Point<S> point, const double max_offset, const List<BSP_Tree<S>> close)
	{
		if (cut != NULL)
		{
			// position of the point with respect to the cut hyperplane
			const double offset = cut.get_hyperplane().get_offset(point);

			if (offset < -max_offset)
			{
				// point is on the minus side of the cut hyperplane
				minus.recurse_close_cuts(point, max_offset, close);
			}
			else if (offset > max_offset)
			{
				// point is on the plus side of the cut hyperplane
				plus.recurse_close_cuts(point, max_offset, close);
			}
			else
			{
				// point is close to the cut hyperplane
				close.add(this);
				minus.recurse_close_cuts(point, max_offset, close);
				plus.recurse_close_cuts(point, max_offset, close);
			}
		}
	}

	/** Perform condensation on a tree.
	 * <p>The condensation operation is not recursive, it must be called
	 * explicitly from leaves to root.</p>
	 */
	private void condense()
	{
		if ((cut != NULL) && (plus.cut == NULL) && (minus.cut == NULL) &&
			(((plus.attribute == NULL) && (minus.attribute == NULL)) ||
				((plus.attribute != NULL) && plus.attribute.equals(minus.attribute))))
		{
			attribute = (plus.attribute == NULL) ? minus.attribute : plus.attribute;
			cut = NULL;
			plus = NULL;
			minus = NULL;
		}
	}

	/** Merge a BSP tree with the instance.
	 * <p>All trees are modified (parts of them are reused in the new
	 * tree), it is the responsibility of the caller to ensure a copy
	 * has been done before if any of the former tree should be
	 * preserved, <em>no</em> such copy is done here!</p>
	 * <p>The algorithm used here is directly derived from the one
	 * described in the Naylor, Amanatides and Thibault paper (section
	 * III, Binary Partitioning of a BSP Tree).</p>
	 * @param tree other tree to merge with the instance (will be
	 * <em>unusable</em> after the operation, as well as the
	 * instance itself)
	 * @param leaf_merger object implementing the const merging phase
	 * (this is where the semantic of the operation occurs, generally
	 * depending on the attribute of the leaf node)
	 * @return a tree, result of <code>instance &lt;op&gt;
	 * tree</code>, this value can be ignored if parent_tree is not NULL
	 * since all connections have already been established
	 */
	public BSP_Tree<S> merge(const BSP_Tree<S> tree, const Leaf_Merger<S> leaf_merger)
	{
		return merge(tree, leaf_merger, NULL, false);
	}

	/** Merge a BSP tree with the instance.
	 * @param tree other tree to merge with the instance (will be
	 * <em>unusable</em> after the operation, as well as the
	 * instance itself)
	 * @param leaf_merger object implementing the const merging phase
	 * (this is where the semantic of the operation occurs, generally
	 * depending on the attribute of the leaf node)
	 * @param parent_tree parent tree to connect to (may be NULL)
	 * @param is_plus_child if true and if parent_tree is not NULL, the
	 * resulting tree should be the plus child of its parent, ignored if
	 * parent_tree is NULL
	 * @return a tree, result of <code>instance &lt;op&gt;
	 * tree</code>, this value can be ignored if parent_tree is not NULL
	 * since all connections have already been established
	 */
	private BSP_Tree<S> merge(const BSP_Tree<S> tree, const Leaf_Merger<S> leaf_merger, const BSP_Tree<S> parent_tree, const bool is_plus_child)
	{
		if (cut == NULL)
		{
			// cell/tree operation
			return leaf_merger.merge(this, tree, parent_tree, is_plus_child, true);
		}
		else if (tree.cut == NULL)
		{
			// tree/cell operation
			return leaf_merger.merge(tree, this, parent_tree, is_plus_child, false);
		}
		else
		{
			// tree/tree operation
			const BSP_Tree<S> merged = tree.split(cut);
			if (parent_tree != NULL)
			{
				merged.parent = parent_tree;
				if (is_plus_child)
				{
					parent_tree.plus = merged;
				}
				else
				{
					parent_tree.minus = merged;
				}
			}

			// merging phase
			plus.merge(merged.plus, leaf_merger, merged, true);
			minus.merge(merged.minus, leaf_merger, merged, false);
			merged.condense();
			if (merged.cut != NULL)
			{
				merged.cut = merged.fit_to_cell(merged.cut.get_hyperplane().whole_hyperplane());
			}

			return merged;
		}
	}

	/** This interface gather the merging operations between a BSP tree
	 * leaf and another BSP tree.
	 * <p>As explained in Bruce Naylor, John Amanatides and William
	 * Thibault paper <a
	 * href="http://www.cs.yorku.ca/~amana/research/bsptSetOp.pdf">Merging
	 * BSP Trees Yields Polyhedral Set Operations</a>, * the operations on {@link BSP_Tree BSP trees} can be expressed as a
	 * generic recursive merging operation where only the const part, * when one of the operand is a leaf, is specific to the real
	 * operation semantics. For example, a tree representing a region
	 * using a bool attribute to identify inside cells and outside
	 * cells would use four different objects to implement the const
	 * merging phase of the four set operations union, intersection, * difference and symmetric difference (exclusive or).</p>
	 * @param <S> Type of the space.
	 */
	class Leaf_Merger<S extends Space>
	{
		/** Merge a leaf node and a tree node.
		 * <p>This method is called at the end of a recursive merging
		 * resulting from a {@code tree1.merge(tree2, leaf_merger)}
		 * call, when one of the sub-trees involved is a leaf (i.e. when
		 * its cut-hyperplane is NULL). This is the only place where the
		 * precise semantics of the operation are required. For all upper
		 * level nodes in the tree, the merging operation is only a
		 * generic partitioning algorithm.</p>
		 * <p>sin_ce the const operation may be non-commutative, it is
		 * important to know if the leaf node comes from the instance tree
		 * ({@code tree1}) or the argument tree
		 * ({@code tree2}). The third argument of the method is
		 * devoted to this. It can be ignored for commutative
		 * operations.</p>
		 * <p>The {@link BSP_Tree#insert_in_tree BSP_Tree.insert_in_tree} method
		 * may be useful to implement this method.</p>
		 * @param leaf leaf node (its cut hyperplane is guaranteed to be
		 * NULL)
		 * @param tree tree node (its cut hyperplane may be NULL or not)
		 * @param parent_tree parent tree to connect to (may be NULL)
		 * @param is_plus_child if true and if parent_tree is not NULL, the
		 * resulting tree should be the plus child of its parent, ignored if
		 * parent_tree is NULL
		 * @param leaf_from_instance if true, the leaf node comes from the
		 * instance tree ({@code tree1}) and the tree node comes from
		 * the argument tree ({@code tree2})
		 * @return the BSP tree resulting from the merging (may be one of
		 * the arguments)
		 */
		BSP_Tree<S> merge(BSP_Tree<S> leaf, BSP_Tree<S> tree, BSP_Tree<S> parent_tree, bool is_plus_child, bool leaf_from_instance);
	}

	/** This interface handles the corner cases when an internal node cut sub-hyperplane vanishes.
	 * <p>
	 * Such cases happens for example when a cut sub-hyperplane is inserted into
	 * another tree (during a merge operation), and is split in several parts, * some of which becomes smaller than the tolerance. The corresponding node
	 * as then no cut sub-hyperplane anymore, but does have children. This interface
	 * specifies how to handle this situation.
	 * setting
	 * </p>
	 * @param <S> Type of the space.
	 */
	class Vanishing_Cut_Handler<S extends Space>
	{
		/** Fix a node with both vanished cut and children.
		 * @param node node to fix
		 * @return fixed node
		 */
		BSP_Tree<S> fix_node(BSP_Tree<S> node);
	}

	/** Split a BSP tree by an external sub-hyperplane.
	 * <p>Split a tree in two halves, on each side of the
	 * sub-hyperplane. The instance is not modified.</p>
	 * <p>The tree returned is not upward-consistent: despite all of its
	 * sub-trees cut sub-hyperplanes (including its own cut
	 * sub-hyperplane) are bounded to the current cell, it is <em>not</em>
	 * attached to any parent tree yet. This tree is intended to be
	 * later inserted into an higher level tree.</p>
	 * <p>The algorithm used here is the one given in Naylor, Amanatides
	 * and Thibault paper (section III, Binary Partitioning of a BSP
	 * Tree).</p>
	 * @param sub partitioning sub-hyperplane, must be already clipped
	 * to the convex region represented by the instance, will be used as
	 * the cut sub-hyperplane of the returned tree
	 * @return a tree having the specified sub-hyperplane as its cut
	 * sub-hyperplane, the two parts of the split instance as its two
	 * sub-trees and a NULL parent
	 */
	public BSP_Tree<S> split(const Sub_Hyperplane<S> sub)
	{
		if (cut == NULL)
		{
			return BSP_Tree<S>(sub, copy_self(), BSP_Tree<S>(attribute), NULL);
		}

		const Hyperplane<S> c_hyperplane = cut.get_hyperplane();
		const Hyperplane<S> s_hyperplane = sub.get_hyperplane();
		const Sub_Hyperplane.Split_Sub_Hyperplane<S> sub_parts = sub.split(c_hyperplane);
		switch (sub_parts.get_side())
		{
		case PLUS:
		{ // the partitioning sub-hyperplane is entirely in the plus sub-tree
			const BSP_Tree<S> split = plus.split(sub);
			if (cut.split(s_hyperplane).get_side() == Side.PLUS)
			{
				split.plus =
					BSP_Tree<>(cut.copy_self(), split.plus, minus.copy_self(), attribute);
				split.plus.condense();
				split.plus.parent = split;
			}
			else
			{
				split.minus =
					BSP_Tree<>(cut.copy_self(), split.minus, minus.copy_self(), attribute);
				split.minus.condense();
				split.minus.parent = split;
			}
			return split;
		}
		case MINUS:
		{ // the partitioning sub-hyperplane is entirely in the minus sub-tree
			const BSP_Tree<S> split = minus.split(sub);
			if (cut.split(s_hyperplane).get_side() == Side.PLUS)
			{
				split.plus =
					BSP_Tree<>(cut.copy_self(), plus.copy_self(), split.plus, attribute);
				split.plus.condense();
				split.plus.parent = split;
			}
			else
			{
				split.minus =
					BSP_Tree<>(cut.copy_self(), plus.copy_self(), split.minus, attribute);
				split.minus.condense();
				split.minus.parent = split;
			}
			return split;
		}
		case BOTH:

		{
			const Sub_Hyperplane.Split_Sub_Hyperplane<S> cut_parts = cut.split(s_hyperplane);
			const BSP_Tree<S> split =
				BSP_Tree<>(sub, plus.split(sub_parts.get_plus()), minus.split(sub_parts.get_minus()), NULL);
			const BSP_Tree<S> tmp = split.plus.minus;
			split.plus.minus = split.minus.plus;
			split.plus.minus.parent = split.plus;
			split.minus.plus = tmp;
			split.minus.plus.parent = split.minus;
			if (cut_parts.get_plus() == NULL)
			{
				split.plus.cut = cut.get_hyperplane().empty_hyperplane();
			}
			else
			{
				split.plus.cut = cut_parts.get_plus();
			}
			if (cut_parts.get_minus() == NULL)
			{
				split.minus.cut = cut.get_hyperplane().empty_hyperplane();
			}
			else
			{
				split.minus.cut = cut_parts.get_minus();
			}
			split.plus.condense();
			split.minus.condense();
			return split;
		}
		default:
			return c_hyperplane.same_orientation_as(s_hyperplane) ?
				BSP_Tree<>(sub, plus.copy_self(), minus.copy_self(), attribute) :
				BSP_Tree<>(sub, minus.copy_self(), plus.copy_self(), attribute);
		}
	}

	/** Insert the instance into another tree.
	 * <p>The instance itself is modified so its former parent should
	 * not be used anymore.</p>
	 * @param parent_tree parent tree to connect to (may be NULL)
	 * @param is_plus_child if true and if parent_tree is not NULL, the
	 * resulting tree should be the plus child of its parent, ignored if
	 * parent_tree is NULL
	 * @param vanishing_handler handler to use for handling very rare corner
	 * cases of vanishing cut sub-hyperplanes in internal nodes during merging
	 * @see Leaf_Merger
	 */
	public void insert_in_tree(const BSP_Tree<S> parent_tree, const bool is_plus_child, const Vanishing_Cut_Handler<S> vanishing_handler)
	{
		// set up parent/child links
		parent = parent_tree;
		if (parent_tree != NULL)
		{
			if (is_plus_child)
			{
				parent_tree.plus = this;
			}
			else
			{
				parent_tree.minus = this;
			}
		}

		// make sure the inserted tree lies in the cell defined by its parent nodes
		if (cut != NULL)
		{
			// explore the parent nodes from here towards tree root
			for (BSP_Tree<S> tree = this; tree.parent != NULL; tree = tree.parent)
			{
				// this is an hyperplane of some parent node
				const Hyperplane<S> hyperplane = tree.parent.cut.get_hyperplane();

				// chop off the parts of the inserted tree that extend
				// on the wrong side of this parent hyperplane
				if (tree == tree.parent.plus)
				{
					cut = cut.split(hyperplane).get_plus();
					plus.chop_off_minus(hyperplane, vanishing_handler);
					minus.chop_off_minus(hyperplane, vanishing_handler);
				}
				else
				{
					cut = cut.split(hyperplane).get_minus();
					plus.chop_off_plus(hyperplane, vanishing_handler);
					minus.chop_off_plus(hyperplane, vanishing_handler);
				}

				if (cut == NULL)
				{
					// the cut sub-hyperplane has vanished
					const BSP_Tree<S> fixed = vanishing_handler.fix_node(this);
					cut = fixed.cut;
					plus = fixed.plus;
					minus = fixed.minus;
					attribute = fixed.attribute;
					if (cut == NULL)
					{
						break;
					}
				}
			}

			// since we may have drop some parts of the inserted tree, // perform a condensation pass to keep the tree structure simple
			condense();
		}
	}

	/** Prune a tree around a cell.
	 * <p>
	 * This method can be used to extract a convex cell from a tree.
	 * The original cell may either be a leaf node or an internal node.
	 * If it is an internal node, it's subtree will be ignored (i.e. the
	 * extracted cell will be a leaf node in all cases). The original
	 * tree to which the original cell belongs is not touched at all, * a independent tree will be built.
	 * </p>
	 * @param cell_attribute attribute to set for the leaf node
	 * corresponding to the initial instance cell
	 * @param other_leafs_attributes attribute to set for the other leaf
	 * nodes
	 * @param internal_attributes attribute to set for the internal nodes
	 * @return a tree (the original tree is left untouched) containing
	 * a single branch with the cell as a leaf node, and other leaf nodes
	 * as the remnants of the pruned branches
	 */
	public BSP_Tree<S> prune_around_convex_cell(const Object cell_attribute, const Object& other_leafs_attributes, const Object internal_attributes)
	{
		// build the current cell leaf
		BSP_Tree<S> tree = BSP_Tree<>(cell_attribute);

		// build the pruned tree bottom-up
		for (BSP_Tree<S> current = this; current.parent != NULL; current = current.parent)
		{
			const Sub_Hyperplane<S> parent_cut = current.parent.cut.copy_self();
			const BSP_Tree<S>       sibling = BSP_Tree<>(other_leafs_attributes);
			if (current == current.parent.plus)
			{
				tree = BSP_Tree<>(parent_cut, tree, sibling, internal_attributes);
			}
			else
			{
				tree = BSP_Tree<>(parent_cut, sibling, tree, internal_attributes);
			}
		}

		return tree;
	}

	/** Chop off parts of the tree.
	 * <p>The instance is modified in place, all the parts that are on
	 * the minus side of the chopping hyperplane are discarded, only the
	 * parts on the plus side remain.</p>
	 * @param hyperplane chopping hyperplane
	 * @param vanishing_handler handler to use for handling very rare corner
	 * cases of vanishing cut sub-hyperplanes in internal nodes during merging
	 */
	private void chop_off_minus(const Hyperplane<S> hyperplane, const Vanishing_Cut_Handler<S> vanishing_handler)
	{
		if (cut != NULL)
		{
			cut = cut.split(hyperplane).get_plus();
			plus.chop_off_minus(hyperplane, vanishing_handler);
			minus.chop_off_minus(hyperplane, vanishing_handler);

			if (cut == NULL)
			{
				// the cut sub-hyperplane has vanished
				const BSP_Tree<S> fixed = vanishing_handler.fix_node(this);
				cut = fixed.cut;
				plus = fixed.plus;
				minus = fixed.minus;
				attribute = fixed.attribute;
			}
		}
	}

	/** Chop off parts of the tree.
	 * <p>The instance is modified in place, all the parts that are on
	 * the plus side of the chopping hyperplane are discarded, only the
	 * parts on the minus side remain.</p>
	 * @param hyperplane chopping hyperplane
	 * @param vanishing_handler handler to use for handling very rare corner
	 * cases of vanishing cut sub-hyperplanes in internal nodes during merging
	 */
	private void chop_off_plus(const Hyperplane<S> hyperplane, const Vanishing_Cut_Handler<S> vanishing_handler)
	{
		if (cut != NULL)
		{
			cut = cut.split(hyperplane).get_minus();
			plus.chop_off_plus(hyperplane, vanishing_handler);
			minus.chop_off_plus(hyperplane, vanishing_handler);

			if (cut == NULL)
			{
				// the cut sub-hyperplane has vanished
				const BSP_Tree<S> fixed = vanishing_handler.fix_node(this);
				cut = fixed.cut;
				plus = fixed.plus;
				minus = fixed.minus;
				attribute = fixed.attribute;
			}
		}
	}
}
