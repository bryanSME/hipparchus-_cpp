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
  //package org.hipparchus.geometry.euclidean.oned;

  //import java.util.Array_list;
  //import java.util.Collection;
  //import java.util.Iterator;
  //import java.util.List;
  //import java.util.No_Such_Element_Exception;

  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.partitioning.Abstract_Region;
  //import org.hipparchus.geometry.partitioning.BSP_Tree;
  //import org.hipparchus.geometry.partitioning.Boundary_Projection;
  //import org.hipparchus.geometry.partitioning.Sub_Hyperplane;
  //import org.hipparchus.util.Precision;

  /** This class represents a 1D region: a set of intervals.
   */
class Intervals_Set extends Abstract_Region<Euclidean_1D, Euclidean_1D> : Iterable<std::vector<double>>
{
public:
	/** Build an intervals set representing the whole real line.
	 * @param tolerance tolerance below which points are considered identical.
	 */
	Intervals_Set(const double& tolerance)
	{
		Iterable<std::vector<double>>(tolerance);
	}

	/** Build an intervals set corresponding to a single interval.
	 * @param lower lower bound of the interval, must be lesser or equal
	 * to {@code upper} (may be {@code -std::numeric_limits<double>::infinity()})
	 * @param upper upper bound of the interval, must be greater or equal
	 * to {@code lower} (may be {@code std::numeric_limits<double>::infinity()})
	 * @param tolerance tolerance below which points are considered identical.
	 */
	Intervals_Set(const double& lower, const double& upper, const double& tolerance)
	{
		Iterable<std::vector<double>>(build_tree(lower, upper, tolerance), tolerance);
	}

	/** Build an intervals set from an inside/outside BSP tree.
	 * <p>The leaf nodes of the BSP tree <em>must</em> have a
	 * {@code Boolean} attribute representing the inside status of
	 * the corresponding cell (true for inside cells, false for outside
	 * cells). In order to avoid building too many small objects, it is
	 * recommended to use the predefined constants
	 * {@code Boolean.TRUE} and {@code Boolean.FALSE}</p>
	 * @param tree inside/outside BSP tree representing the intervals set
	 * @param tolerance tolerance below which points are considered identical.
	 */
	Intervals_Set(const BSP_Tree<Euclidean_1D> tree, const double& tolerance)
	{
		Iterable<std::vector<double>>(tree, tolerance);
	}

	/** Build an intervals set from a Boundary RE_Presentation (B-rep).
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
	 * calls to the {@link
	 * org.hipparchus.geometry.partitioning.Region#check_point(org.hipparchus.geometry.Point)
	 * check_point} method will not be meaningful anymore.</p>
	 * <p>If the boundary is empty, the region will represent the whole
	 * space.</p>
	 * @param boundary collection of boundary elements
	 * @param tolerance tolerance below which points are considered identical.
	 */
	Intervals_Set(const Collection<Sub_Hyperplane<Euclidean_1D>> boundary, const double& tolerance)
	{
		Iterable<std::vector<double>>(boundary, tolerance);
	}

	/** {@inherit_doc} */
	//override
	Intervals_Set build_new(const BSP_Tree<Euclidean_1D>& tree)
	{
		return Intervals_Set(tree, get_tolerance());
	}

	/** Get the lowest value belonging to the instance.
	 * @return lowest value belonging to the instance
	 * ({@code -std::numeric_limits<double>::infinity()} if the instance doesn't
	 * have any low bound, {@code std::numeric_limits<double>::infinity()} if the
	 * instance is empty)
	 */
	double get_inf()
	{
		BSP_Tree<Euclidean_1D> node = get_tree(false);
		double inf = std::numeric_limits<double>::infinity();
		while (node.get_cut() != NULL)
		{
			const Oriented_Point op = (Oriented_Point)node.get_cut().get_hyperplane();
			inf = op.get_location().get_x();
			node = op.is_direct() ? node.get_minus() : node.get_plus();
		}
		return ((Boolean)node.get_attribute()) ? -std::numeric_limits<double>::infinity() : inf;
	}

	/** Get the highest value belonging to the instance.
	 * @return highest value belonging to the instance
	 * ({@code std::numeric_limits<double>::infinity()} if the instance doesn't
	 * have any high bound, {@code -std::numeric_limits<double>::infinity()} if the
	 * instance is empty)
	 */
	double get_sup()
	{
		BSP_Tree<Euclidean_1D> node = get_tree(false);
		double  sup = -std::numeric_limits<double>::infinity();
		while (node.get_cut() != NULL)
		{
			const Oriented_Point op = (Oriented_Point)node.get_cut().get_hyperplane();
			sup = op.get_location().get_x();
			node = op.is_direct() ? node.get_plus() : node.get_minus();
		}
		return ((Boolean)node.get_attribute()) ? std::numeric_limits<double>::infinity() : sup;
	}

	/** {@inherit_doc}
	 */
	 //override
	Boundary_Projection<Euclidean_1D> project_to_boundary(const Point<Euclidean_1D>& point)
	{
		// get position of test point
		const double x = ((Vector_1D)point).get_x();

		double previous = -std::numeric_limits<double>::infinity();
		for (const std::vector<double>& a : *this)
		{
			if (x < a[0])
			{
				// the test point lies between the previous and the current intervals
				// offset will be positive
				const double previous_offset = x - previous;
				const double current_offset = a[0] - x;
				return previous_offset < current_offset
					? Boundary_Projection<Euclidean_1D>(point, finite_or_null_point(previous), previous_offset)
					: Boundary_Projection<Euclidean_1D>(point, finite_or_null_point(a[0]), current_offset);
			}
			if (x <= a[1])
			{
				// the test point lies within the current interval
				// offset will be negative
				const double offset0 = a[0] - x;
				const double offset1 = x - a[1];
				return offset0 < offset1
					? Boundary_Projection<Euclidean_1D>(point, finite_or_null_point(a[1]), offset1)
					: Boundary_Projection<Euclidean_1D>(point, finite_or_null_point(a[0]), offset0);
			}
			previous = a[1];
		}
		// the test point if past the last sub-interval
		return Boundary_Projection<Euclidean_1D>(point, finite_or_null_point(previous), x - previous);
	}

	/** Build an ordered list of intervals representing the instance.
	 * <p>This method builds this intervals set as an ordered list of
	 * {@link Interval Interval} elements. If the intervals set has no
	 * lower limit, the first interval will have its low bound equal to
	 * {@code -std::numeric_limits<double>::infinity()}. If the intervals set has
	 * no upper limit, the last interval will have its upper bound equal
	 * to {@code std::numeric_limits<double>::infinity()}. An empty tree will
	 * build an empty list while a tree representing the whole real line
	 * will build a one element list with both bounds being
	 * infinite.</p>
	 * @return a ordered list containing {@link Interval Interval}
	 * elements
	 */
	std::vector<Interval> as_list()
	{
		std::vector<Interval> list;
		for (const std::vector<double>& a : *this)
		{
			list.add(Interval(a[0], a[1]));
		}
		return list;
	}

	/** {@inherit_doc}
	 * <p>
	 * The iterator returns the limit values of sub-intervals in ascending order.
	 * </p>
	 * <p>
	 * The iterator does <em>not</em> support the optional {@code remove} operation.
	 * </p>
	 */
	 //override
	std::vector<double>::Iterator iterator()
	{
		return Sub_Intervals_Iterator();
	}

private:
	/** Build an inside/outside tree representing a single interval.
	 * @param lower lower bound of the interval, must be lesser or equal
	 * to {@code upper} (may be {@code -std::numeric_limits<double>::infinity()})
	 * @param upper upper bound of the interval, must be greater or equal
	 * to {@code lower} (may be {@code std::numeric_limits<double>::infinity()})
	 * @param tolerance tolerance below which points are considered identical.
	 * @return the built tree
	 */
	static BSP_Tree<Euclidean_1D> build_tree(const double& lower, const double& upper, const double& tolerance)
	{
		if (std::isinf(lower) && (lower < 0))
		{
			if (std::isinf(upper) && (upper > 0))
			{
				// the tree must cover the whole real line
				return BSP_Tree<Euclidean_1D>(Boolean.TRUE);
			}
			// the tree must be open on the negative infinity side
			const Sub_Hyperplane<Euclidean_1D> upper_cut = Oriented_Point(Vector_1D(upper), true, tolerance).whole_hyperplane();
			return BSP_Tree<Euclidean_1D>(upper_cut, BSP_Tree<Euclidean_1D>(Boolean.FALSE), BSP_Tree<Euclidean_1D>(Boolean.TRUE), NULL);
		}
		const Sub_Hyperplane<Euclidean_1D> lower_cut = Oriented_Point(Vector_1D(lower), false, tolerance).whole_hyperplane();
		if (std::isinf(upper) && (upper > 0))
		{
			// the tree must be open on the positive infinity side
			return BSP_Tree<Euclidean_1D>(lower_cut, BSP_Tree<Euclidean_1D>(Boolean.FALSE), BSP_Tree<Euclidean_1D>(Boolean.TRUE), NULL);
		}

		// the tree must be bounded on the two sides
		const Sub_Hyperplane<Euclidean_1D> upper_cut = Oriented_Point(Vector_1D(upper), true, tolerance).whole_hyperplane();
		return BSP_Tree<Euclidean_1D>(lower_cut, BSP_Tree<Euclidean_1D>(Boolean.FALSE), BSP_Tree<Euclidean_1D>(upper_cut, BSP_Tree<Euclidean_1D>(Boolean.FALSE), BSP_Tree<Euclidean_1D>(Boolean.TRUE), NULL), NULL);
	}

	/** Build a finite point.
	 * @param x abscissa of the point
	 * @return a point for finite abscissa, NULL otherwise
	 */
	Vector_1D finite_or_null_point(const double& x)
	{
		return std::isinf(x) ? NULL : Vector_1D(x);
	}

	/** Get the first leaf node of a tree.
	 * @param root tree root
	 * @return first leaf node
	 */
	BSP_Tree<Euclidean_1D> get_first_leaf(const BSP_Tree<Euclidean_1D>& root)
	{
		if (root.get_cut() == NULL)
		{
			return root;
		}

		// find the smallest internal node
		BSP_Tree<Euclidean_1D> smallest = NULL;
		for (BSP_Tree<Euclidean_1D> n = root; n != NULL; n = previous_internal_node(n))
		{
			smallest = n;
		}

		return leaf_before(smallest);
	}

	/** Get the node corresponding to the first interval boundary.
	 * @return smallest internal node, * or NULL if there are no internal nodes (i.e. the set is either empty or covers the real line)
	 */
	BSP_Tree<Euclidean_1D> get_first_interval_boundary()
	{
		// start search at the tree root
		BSP_Tree<Euclidean_1D> node = get_tree(false);
		if (node.get_cut() == NULL)
		{
			return NULL;
		}

		// walk tree until we find the smallest internal node
		node = get_first_leaf(node).get_parent();

		// walk tree until we find an interval boundary
		while (node != NULL && !(is_interval_start(node) || is_interval_end(node)))
		{
			node = next_internal_node(node);
		}

		return node;
	}

	/** Check if an internal node corresponds to the start abscissa of an interval.
	 * @param node internal node to check
	 * @return true if the node corresponds to the start abscissa of an interval
	 */
	bool is_interval_start(const BSP_Tree<Euclidean_1D>& node)
	{
		if ((Boolean)leaf_before(node).get_attribute())
		{
			// it has an inside cell before it, it may end an interval but not start it
			return false;
		}

		if (!(Boolean)leaf_after(node).get_attribute())
		{
			// it has an outside cell after it, it is a dummy cut away from real intervals
			return false;
		}

		// the cell has an outside before and an inside after it
		// it is the start of an interval
		return true;
	}

	/** Check if an internal node corresponds to the end abscissa of an interval.
	 * @param node internal node to check
	 * @return true if the node corresponds to the end abscissa of an interval
	 */
	bool is_interval_end(const BSP_Tree<Euclidean_1D>& node)
	{
		if (!(Boolean)leaf_before(node).get_attribute())
		{
			// it has an outside cell before it, it may start an interval but not end it
			return false;
		}

		if ((Boolean)leaf_after(node).get_attribute())
		{
			// it has an inside cell after it, it is a dummy cut in the middle of an interval
			return false;
		}

		// the cell has an inside before and an outside after it
		// it is the end of an interval
		return true;
	}

	/** Get the next internal node.
	 * @param node current internal node
	 * @return next internal node in ascending order, or NULL
	 * if this is the last internal node
	 */
	BSP_Tree<Euclidean_1D> next_internal_node(BSP_Tree<Euclidean_1D>& node)
	{
		if (child_after(node).get_cut() != NULL)
		{
			// the next node is in the sub-tree
			return leaf_after(node).get_parent();
		}

		// there is nothing left deeper in the tree, we backtrack
		while (is_after_parent(node))
		{
			node = node.get_parent();
		}
		return node.get_parent();
	}

	/** Get the previous internal node.
	 * @param node current internal node
	 * @return previous internal node in ascending order, or NULL
	 * if this is the first internal node
	 */
	BSP_Tree<Euclidean_1D> previous_internal_node(BSP_Tree<Euclidean_1D>& node)
	{
		if (child_before(node).get_cut() != NULL)
		{
			// the next node is in the sub-tree
			return leaf_before(node).get_parent();
		}

		// there is nothing left deeper in the tree, we backtrack
		while (is_before_parent(node))
		{
			node = node.get_parent();
		}
		return node.get_parent();
	}

	/** Find the leaf node just before an internal node.
	 * @param node internal node at which the sub-tree starts
	 * @return leaf node just before the internal node
	 */
	BSP_Tree<Euclidean_1D> leaf_before(BSP_Tree<Euclidean_1D>& node)
	{
		node = child_before(node);
		while (node.get_cut() != NULL)
		{
			node = child_after(node);
		}

		return node;
	}

	/** Find the leaf node just after an internal node.
	 * @param node internal node at which the sub-tree starts
	 * @return leaf node just after the internal node
	 */
	BSP_Tree<Euclidean_1D> leaf_after(BSP_Tree<Euclidean_1D>& node)
	{
		node = child_after(node);
		while (node.get_cut() != NULL)
		{
			node = child_before(node);
		}

		return node;
	}

	/** Check if a node is the child before its parent in ascending order.
	 * @param node child node considered
	 * @return true is the node has a parent end is before it in ascending order
	 */
	bool is_before_parent(const BSP_Tree<Euclidean_1D>& node)
	{
		const BSP_Tree<Euclidean_1D> parent = node.get_parent();
		return parent == NULL
			? false
			? node == child_before(parent);
	}

	/** Check if a node is the child after its parent in ascending order.
	 * @param node child node considered
	 * @return true is the node has a parent end is after it in ascending order
	 */
	bool is_after_parent(const BSP_Tree<Euclidean_1D>& node)
	{
		const BSP_Tree<Euclidean_1D> parent = node.get_parent();
		return parent == NULL
			? false
			: node == child_after(parent);
	}

	/** Find the child node just before an internal node.
	 * @param node internal node at which the sub-tree starts
	 * @return child node just before the internal node
	 */
	BSP_Tree<Euclidean_1D> child_before(BSP_Tree<Euclidean_1D>& node)
	{
		return (is_direct(node))
			? node.get_minus() // smaller abscissas are on minus side, larger abscissas are on plus side
			: node.get_plus(); // smaller abscissas are on plus side, larger abscissas are on minus side
	}

	/** Find the child node just after an internal node.
	 * @param node internal node at which the sub-tree starts
	 * @return child node just after the internal node
	 */
	BSP_Tree<Euclidean_1D> child_after(BSP_Tree<Euclidean_1D>& node)
	{
		return is_direct(node)
			? node.get_plus() // smaller abscissas are on minus side, larger abscissas are on plus side
			: node.get_minus(); // smaller abscissas are on plus side, larger abscissas are on minus side
	}

	/** Check if an internal node has a direct oriented point.
	 * @param node internal node to check
	 * @return true if the oriented point is direct
	 */
	bool is_direct(const BSP_Tree<Euclidean_1D>& node)
	{
		return ((Oriented_Point)node.get_cut().get_hyperplane()).is_direct();
	}

	/** Get the abscissa of an internal node.
	 * @param node internal node to check
	 * @return abscissa
	 */
	double get_angle(const BSP_Tree<Euclidean_1D>& node)
	{
		return ((Oriented_Point)node.get_cut().get_hyperplane()).get_location().get_x();
	}

protected:
	/** {@inherit_doc} */
	//override
	void compute_geometrical_properties()
	{
		if (get_tree(false).get_cut() == NULL)
		{
			set_barycenter((Point<Euclidean_1D>) Vector_1D.NaN);
			set_size(((Boolean)get_tree(false).get_attribute()) ? std::numeric_limits<double>::infinity() : 0);
		}
		else
		{
			double size{};
			double sum{};
			for (const Interval interval : as_list())
			{
				size += interval.get_size();
				sum += interval.get_size() * interval.get_barycenter();
			}
			set_size(size);
			if (std::isinf(size))
			{
				set_barycenter((Point<Euclidean_1D>) Vector_1D.NaN);
			}
			else if (size >= Precision.SAFE_MIN)
			{
				set_barycenter((Point<Euclidean_1D>) Vector_1D(sum / size));
			}
			else
			{
				set_barycenter((Point<Euclidean_1D>) ((Oriented_Point)get_tree(false).get_cut().get_hyperplane()).get_location());
			}
		}
	}
};
