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
  //package org.hipparchus.geometry.spherical.oned;

  //import java.util.Array_list;
  //import java.util.Collection;
  //import java.util.Iterator;
  //import java.util.List;
  //import java.util.No_Such_Element_Exception;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.geometry.Localized_Geometry_Formats;
  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.partitioning.Abstract_Region;
  //import org.hipparchus.geometry.partitioning.BSP_Tree;
  //import org.hipparchus.geometry.partitioning.Boundary_Projection;
  //import org.hipparchus.geometry.partitioning.Side;
  //import org.hipparchus.geometry.partitioning.Sub_Hyperplane;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Precision;

  /** This class represents a region of a circle: a set of arcs.
   * <p>
   * Note that due to the wrapping around \(2 \pi\), barycenter is
   * ill-defined here. It was defined only in order to fulfill
   * the requirements of the {@link
   * org.hipparchus.geometry.partitioning.Region Region}
   * interface, but its use is discouraged.
   * </p>
   */
class Arcs_Set : public Abstract_Region<Sphere_1D, Sphere_1D>, Iterable<std::vector<double>>
{
public:
	/** Build an arcs set representing the whole circle.
	 * @param tolerance tolerance below which close sub-arcs are merged together
	 * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	Arcs_Set(const double& tolerance)
	{
		super(tolerance);
		Sphere_1D.check_tolerance(tolerance);
	}

	/** Build an arcs set corresponding to a single arc.
	 * <p>
	 * If either {@code lower} is equals to {@code upper} or
	 * the interval exceeds \( 2 \pi \), the arc is considered
	 * to be the full circle and its initial defining boundaries
	 * will be forgotten. {@code lower} is not allowed to be greater
	 * than {@code upper} (an exception is thrown in this case).
	 * </p>
	 * @param lower lower bound of the arc
	 * @param upper upper bound of the arc
	 * @param tolerance tolerance below which close sub-arcs are merged together
	 * @exception  if lower is greater than upper
	 * or tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	Arcs_Set(const double lower, const double upper, const double& tolerance)
	{
		super(build_tree(lower, upper, tolerance), tolerance);
		Sphere_1D.check_tolerance(tolerance);
	}

	/** Build an arcs set from an inside/outside BSP tree.
	 * <p>The leaf nodes of the BSP tree <em>must</em> have a
	 * {@code Boolean} attribute representing the inside status of
	 * the corresponding cell (true for inside cells, false for outside
	 * cells). In order to avoid building too many small objects, it is
	 * recommended to use the predefined constants
	 * {@code Boolean.TRUE} and {@code Boolean.FALSE}</p>
	 * @param tree inside/outside BSP tree representing the arcs set
	 * @param tolerance tolerance below which close sub-arcs are merged together
	 * @exception Inconsistent_State_At_2_Pi_Wrapping if the tree leaf nodes are not
	 * consistent across the \( 0, 2 \pi \) crossing
	 * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	Arcs_Set(const BSP_Tree<Sphere_1D>& tree, const double& tolerance)
	{
		super(tree, tolerance);
		Sphere_1D.check_tolerance(tolerance);
		check_2_pi_consistency();
	}

	/** Build an arcs set from a Boundary RE_Presentation (B-rep).
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
	 * @param tolerance tolerance below which close sub-arcs are merged together
	 * @exception Inconsistent_State_At_2_Pi_Wrapping if the tree leaf nodes are not
	 * consistent across the \( 0, 2 \pi \) crossing
	 * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	Arcs_Set(const Collection<Sub_Hyperplane<Sphere_1D>>& boundary, const double& tolerance)
	{
		super(boundary, tolerance);
		Sphere_1D.check_tolerance(tolerance);
		check_2_pi_consistency();
	}

	/** {@inherit_doc} */
	//override
	Arcs_Set build_new(const BSP_Tree<Sphere_1D>& tree)
	{
		return Arcs_Set(tree, get_tolerance());
	}

	/** {@inherit_doc}*/
	//override
	Boundary_Projection<Sphere_1D> project_to_boundary(const Point<Sphere_1D>& point)
	{
		// get position of test point
		const double& alpha = ((S1_Point)point).get_alpha();

		bool wrap_first{};
		double first = std::numeric_limits<double>::quiet_NaN();
		double previous = std::numeric_limits<double>::quiet_NaN();
		for (const std::vector<double> a : *this)
		{
			if (std::isnan(first))
			{
				// remember the first angle in case we need it later
				first = a[0];
			}

			if (!wrap_first)
			{
				if (alpha < a[0])
				{
					// the test point lies between the previous and the current arcs
					// offset will be positive
					if (std::isnan(previous))
					{
						// we need to wrap around the circle
						wrap_first = true;
					}
					else
					{
						const double previous_offset = alpha - previous;
						const double current_offset = a[0] - alpha;
						return (previous_offset < current_offset)
							? Boundary_Projection<Sphere_1D>(point, S1_Point(previous), previous_offset)
							: Boundary_Projection<Sphere_1D>(point, S1_Point(a[0]), current_offset);
					}
				}
				else if (alpha <= a[1])
				{
					// the test point lies within the current arc
					// offset will be negative
					const double offset0 = a[0] - alpha;
					const double offset1 = alpha - a[1];
					return offset0 < offset1
						? Boundary_Projection<Sphere_1D>(point, S1_Point(a[1]), offset1)
						: Boundary_Projection<Sphere_1D>(point, S1_Point(a[0]), offset0);
				}
			}
			previous = a[1];
		}

		if (std::isnan(previous))
		{
			// there are no points at all in the arcs set
			return Boundary_Projection<Sphere_1D>(point, NULL, Math_Utils::TWO_PI);
		}

		// the test point if before first arc and after last arc, // somewhere around the 0/2 \pi crossing
		if (wrap_first)
		{
			// the test point is between 0 and first
			const double previous_offset = alpha - (previous - Math_Utils::TWO_PI);
			const double current_offset = first - alpha;
			return previous_offset < current_offset
				? Boundary_Projection<Sphere_1D>(point, S1_Point(previous), previous_offset)
				: Boundary_Projection<Sphere_1D>(point, S1_Point(first), current_offset);
		}
		// the test point is between last and 2\pi
		const double previous_offset = alpha - previous;
		const double current_offset = first + Math_Utils::TWO_PI - alpha;
		return (previous_offset < current_offset)
			? Boundary_Projection<Sphere_1D>(point, S1_Point(previous), previous_offset)
			: Boundary_Projection<Sphere_1D>(point, S1_Point(first), current_offset);
		
	}

	/** Build an ordered list of arcs representing the instance.
	 * <p>This method builds this arcs set as an ordered list of
	 * {@link Arc Arc} elements. An empty tree will build an empty list
	 * while a tree representing the whole circle will build a one
	 * element list with bounds set to \( 0 and 2 \pi \).</p>
	 * @return a ordered list containing {@link Arc Arc} elements
	 */
	std::vector<Arc> as_list()
	{
		std::vector<Arc> list;
		for (const std::vector<double> a : *this)
		{
			list.add(Arc(a[0], a[1], get_tolerance()));
		}
		return list;
	}

	/** {@inherit_doc}
	 * <p>
	 * The iterator returns the limit angles pairs of sub-arcs in trigonometric order.
	 * </p>
	 * <p>
	 * The iterator does <em>not</em> support the optional {@code remove} operation.
	 * </p>
	 */
	 //override
	std::vector<double>::Iterator iterator()
	{
		return Sub_Arcs_Iterator();
	}

	/** Split the instance in two parts by an arc.
	 * @param arc splitting arc
	 * @return an object containing both the part of the instance
	 * on the plus side of the arc and the part of the
	 * instance on the minus side of the arc
	 */
	Split split(const Arc& arc)
	{
		const List<Double> minus = Array_list<>();
		const List<Double>  plus = Array_list<>();

		const double reference = std::numbers::pi + arc.get_inf();
		const double& arc_length = arc.get_sup() - arc.get_inf();

		for (const std::vector<double> a : this)
		{
			const double synced_start = Math_Utils::normalize_angle(a[0], reference) - arc.get_inf();
			const double& arc_offset = a[0] - synced_start;
			const double synced_end = a[1] - arc_offset;
			if (synced_start < arc_length)
			{
				// the start point a[0] is in the minus part of the arc
				minus.add(a[0]);
				if (synced_end > arc_length)
				{
					// the end point a[1] is past the end of the arc
					// so we leave the minus part and enter the plus part
					const double minus_to_plus = arc_length + arc_offset;
					minus.add(minus_to_plus);
					plus.add(minus_to_plus);
					if (synced_end > Math_Utils::TWO_PI)
					{
						// in fact the end point a[1] goes far enough that we
						// leave the plus part of the arc and enter the minus part again
						const double plus_to_minus = Math_Utils::TWO_PI + arc_offset;
						plus.add(plus_to_minus);
						minus.add(plus_to_minus);
						minus.add(a[1]);
					}
					else
					{
						// the end point a[1] is in the plus part of the arc
						plus.add(a[1]);
					}
				}
				else
				{
					// the end point a[1] is in the minus part of the arc
					minus.add(a[1]);
				}
			}
			else
			{
				// the start point a[0] is in the plus part of the arc
				plus.add(a[0]);
				if (synced_end > Math_Utils::TWO_PI)
				{
					// the end point a[1] wraps around to the start of the arc
					// so we leave the plus part and enter the minus part
					const double plus_to_minus = Math_Utils::TWO_PI + arc_offset;
					plus.add(plus_to_minus);
					minus.add(plus_to_minus);
					if (synced_end > Math_Utils::TWO_PI + arc_length)
					{
						// in fact the end point a[1] goes far enough that we
						// leave the minus part of the arc and enter the plus part again
						const double minus_to_plus = Math_Utils::TWO_PI + arc_length + arc_offset;
						minus.add(minus_to_plus);
						plus.add(minus_to_plus);
						plus.add(a[1]);
					}
					else
					{
						// the end point a[1] is in the minus part of the arc
						minus.add(a[1]);
					}
				}
				else
				{
					// the end point a[1] is in the plus part of the arc
					plus.add(a[1]);
				}
			}
		}

		return Split(create_split_part(plus), create_split_part(minus));
	}

protected:
	/** {@inherit_doc} */
	//override
	void compute_geometrical_properties()
	{
		if (get_tree(false).get_cut() == NULL)
		{
			set_barycenter(S1_Point.NaN);
			set_size(
				((bool)get_tree(false).get_attribute())
					? Math_Utils::TWO_PI
					: 0
			);
		}
		else
		{
			double size{};
			double sum{};
			for (const std::vector<double> a : *this)
			{
				const double length = a[1] - a[0];
				size += length;
				sum += length * (a[0] + a[1]);
			}
			set_size(size);
			if (Precision::equals(size, Math_Utils::TWO_PI, 0))
			{
				set_barycenter(S1_Point.NaN);
			}
			else if (size >= Precision.SAFE_MIN)
			{
				set_barycenter(new S1_Point(sum / (2 * size)));
			}
			else
			{
				const Limit_Angle limit = (Limit_Angle)get_tree(false).get_cut().get_hyperplane();
				set_barycenter(limit.get_location());
			}
		}
	}


private:
	/** Build an inside/outside tree representing a single arc.
	 * @param lower lower angular bound of the arc
	 * @param upper upper angular bound of the arc
	 * @param tolerance tolerance below which close sub-arcs are merged together
	 * @return the built tree
	 * @exception  if lower is greater than upper
	 * or tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	static BSP_Tree<Sphere_1D> build_tree(const double& lower, const double& upper, const double& tolerance)
	{
		Sphere_1D.check_tolerance(tolerance);
		if (Precision::equals(lower, upper, 0) || (upper - lower) >= Math_Utils::TWO_PI)
		{
			// the tree must cover the whole circle
			return BSP_Tree<Sphere_1D>(Boolean.TRUE);
		}
		if (lower > upper)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::ENDPOINTS_NOT_AN_INTERVAL, lower, upper, true);
		}

		// this is a regular arc, covering only part of the circle
		const double normalized_lower = Math_Utils::normalize_angle(lower, std::numbers::pi);
		const double normalized_upper = normalized_lower + (upper - lower);
		const Sub_Hyperplane<Sphere_1D> lower_cut = Limit_Angle(new S1_Point(normalized_lower), false, tolerance).whole_hyperplane();

		if (normalized_upper <= Math_Utils::TWO_PI)
		{
			// simple arc starting after 0 and ending before 2 \pi
			const Sub_Hyperplane<Sphere_1D> upper_cut = Limit_Angle(new S1_Point(normalized_upper), true, tolerance).whole_hyperplane();
			return BSP_Tree<Sphere_1D>(lower_cut, BSP_Tree<Sphere_1D>(Boolean.FALSE), BSP_Tree<Sphere_1D>(upper_cut, BSP_Tree<Sphere_1D>(Boolean.FALSE), BSP_Tree<Sphere_1D>(Boolean.TRUE), NULL), NULL);
		}
		
		// arc wrapping around 2 \pi
		const Sub_Hyperplane<Sphere_1D> upper_cut = Limit_Angle(new S1_Point(normalized_upper - Math_Utils::TWO_PI), true, tolerance).whole_hyperplane();
		return BSP_Tree<Sphere_1D>(lower_cut, BSP_Tree<Sphere_1D>(upper_cut, BSP_Tree<Sphere_1D>(Boolean.FALSE), BSP_Tree<Sphere_1D>(Boolean.TRUE), NULL), BSP_Tree<Sphere_1D>(Boolean.TRUE), NULL);
	}

	/** Check consistency.
	* @exception Inconsistent_State_At_2_Pi_Wrapping if the tree leaf nodes are not
	* consistent across the \( 0, 2 \pi \) crossing
	*/
	void check_2_pi_consistency()
	{
		// start search at the tree root
		BSP_Tree<Sphere_1D> root = get_tree(false);
		if (root.get_cut() == NULL)
		{
			return;
		}

		// find the inside/outside state before the smallest internal node
		const Boolean state_before = (bool)get_first_leaf(root).get_attribute();

		// find the inside/outside state after the largest internal node
		const Boolean state_after = (bool)get_last_leaf(root).get_attribute();

		if (state_before ^ state_after)
		{
			throw Inconsistent_State_At_2_Pi_Wrapping();
		}
	}

	/** Get the first leaf node of a tree.
	 * @param root tree root
	 * @return first leaf node (i.e. node corresponding to the region just after 0.0 radians)
	 */
	BSP_Tree<Sphere_1D> get_first_leaf(const BSP_Tree<Sphere_1D>& root)
	{
		if (root.get_cut() == NULL)
		{
			return root;
		}

		// find the smallest internal node
		BSP_Tree<Sphere_1D> smallest;
		for (BSP_Tree<Sphere_1D> n = root; n != NULL; n = previous_internal_node(n))
		{
			smallest = n;
		}

		return leaf_before(smallest);
	}

	/** Get the last leaf node of a tree.
	 * @param root tree root
	 * @return last leaf node (i.e. node corresponding to the region just before \( 2 \pi \) radians)
	 */
	BSP_Tree<Sphere_1D> get_last_leaf(const BSP_Tree<Sphere_1D>& root)
	{
		if (root.get_cut() == NULL)
		{
			return root;
		}

		// find the largest internal node
		BSP_Tree<Sphere_1D> largest = NULL;
		for (BSP_Tree<Sphere_1D> n = root; n != NULL; n = next_internal_node(n))
		{
			largest = n;
		}

		return leaf_after(largest);
	}

	/** Get the node corresponding to the first arc start.
	 * @return smallest internal node (i.e. first after 0.0 radians, in trigonometric direction), * or NULL if there are no internal nodes (i.e. the set is either empty or covers the full circle)
	 */
	BSP_Tree<Sphere_1D> get_first_arc_start()
	{
		// start search at the tree root
		BSP_Tree<Sphere_1D> node = get_tree(false);
		if (node.get_cut() == NULL)
		{
			return NULL;
		}

		// walk tree until we find the smallest internal node
		node = get_first_leaf(node).get_parent();

		// walk tree until we find an arc start
		while (node != NULL && !is_arc_start(node))
		{
			node = next_internal_node(node);
		}

		return node;
	}

	/** Check if an internal node corresponds to the start angle of an arc.
	 * @param node internal node to check
	 * @return true if the node corresponds to the start angle of an arc
	 */
	bool is_arc_start(const BSP_Tree<Sphere_1D>& node)
	{
		if ((Boolean)leaf_before(node).get_attribute())
		{
			// it has an inside cell before it, it may end an arc but not start it
			return false;
		}

		if (!(Boolean)leaf_after(node).get_attribute())
		{
			// it has an outside cell after it, it is a dummy cut away from real arcs
			return false;
		}

		// the cell has an outside before and an inside after it
		// it is the start of an arc
		return true;
	}

	/** Check if an internal node corresponds to the end angle of an arc.
	 * @param node internal node to check
	 * @return true if the node corresponds to the end angle of an arc
	 */
	bool is_arc_end(const BSP_Tree<Sphere_1D>& node)
	{
		if (!(Boolean)leaf_before(node).get_attribute())
		{
			// it has an outside cell before it, it may start an arc but not end it
			return false;
		}

		if ((Boolean)leaf_after(node).get_attribute())
		{
			// it has an inside cell after it, it is a dummy cut in the middle of an arc
			return false;
		}

		// the cell has an inside before and an outside after it
		// it is the end of an arc
		return true;
	}

	/** Get the next internal node.
	 * @param node current internal node
	 * @return next internal node in trigonometric order, or NULL
	 * if this is the last internal node
	 */
	BSP_Tree<Sphere_1D> next_internal_node(BSP_Tree<Sphere_1D>& node)
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
	 * @return previous internal node in trigonometric order, or NULL
	 * if this is the first internal node
	 */
	BSP_Tree<Sphere_1D> previous_internal_node(BSP_Tree<Sphere_1D>& node)
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
	BSP_Tree<Sphere_1D> leaf_before(BSP_Tree<Sphere_1D> node)
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
	BSP_Tree<Sphere_1D> leaf_after(BSP_Tree<Sphere_1D> node)
	{
		node = child_after(node);
		while (node.get_cut() != NULL)
		{
			node = child_before(node);
		}

		return node;
	}

	/** Check if a node is the child before its parent in trigonometric order.
	 * @param node child node considered
	 * @return true is the node has a parent end is before it in trigonometric order
	 */
	bool is_before_parent(const BSP_Tree<Sphere_1D> node)
	{
		const BSP_Tree<Sphere_1D> parent = node.get_parent();
		if (parent == NULL)
		{
			return false;
		}
		else
		{
			return node == child_before(parent);
		}
	}

	/** Check if a node is the child after its parent in trigonometric order.
	 * @param node child node considered
	 * @return true is the node has a parent end is after it in trigonometric order
	 */
	bool is_after_parent(const BSP_Tree<Sphere_1D> node)
	{
		const BSP_Tree<Sphere_1D> parent = node.get_parent();
		if (parent == NULL)
		{
			return false;
		}
		else
		{
			return node == child_after(parent);
		}
	}

	/** Find the child node just before an internal node.
	 * @param node internal node at which the sub-tree starts
	 * @return child node just before the internal node
	 */
	BSP_Tree<Sphere_1D> child_before(BSP_Tree<Sphere_1D> node)
	{
		if (is_direct(node))
		{
			// smaller angles are on minus side, larger angles are on plus side
			return node.get_minus();
		}
		else
		{
			// smaller angles are on plus side, larger angles are on minus side
			return node.get_plus();
		}
	}

	/** Find the child node just after an internal node.
	 * @param node internal node at which the sub-tree starts
	 * @return child node just after the internal node
	 */
	BSP_Tree<Sphere_1D> child_after(BSP_Tree<Sphere_1D> node)
	{
		if (is_direct(node))
		{
			// smaller angles are on minus side, larger angles are on plus side
			return node.get_plus();
		}
		else
		{
			// smaller angles are on plus side, larger angles are on minus side
			return node.get_minus();
		}
	}

	/** Check if an internal node has a direct limit angle.
	 * @param node internal node to check
	 * @return true if the limit angle is direct
	 */
	bool is_direct(const BSP_Tree<Sphere_1D> node)
	{
		return ((Limit_Angle)node.get_cut().get_hyperplane()).is_direct();
	}

	/** Get the limit angle of an internal node.
	 * @param node internal node to check
	 * @return limit angle
	 */
	double get_angle(const BSP_Tree<Sphere_1D> node)
	{
		return ((Limit_Angle)node.get_cut().get_hyperplane()).get_location().get_alpha();
	}

	/** Add an arc limit to a BSP tree under construction.
	 * @param tree BSP tree under construction
	 * @param alpha arc limit
	 * @param is_start if true, the limit is the start of an arc
	 */
	void add_arc_limit(const BSP_Tree<Sphere_1D>& tree, const double& alpha, const bool is_start)
	{
		const Limit_Angle limit = Limit_Angle(S1_Point(alpha), !is_start, get_tolerance());
		const BSP_Tree<Sphere_1D> node = tree.get_cell(limit.get_location(), get_tolerance());
		if (node.get_cut() != NULL)
		{
			// this should never happen
			throw Math_Runtime_Exception.create_internal_error();
		}

		node.insert_cut(limit);
		node.set_attribute(null);
		node.get_plus().set_attribute(Boolean.FALSE);
		node.get_minus().set_attribute(Boolean.TRUE);
	}

	/** Create a split part.
	 * <p>
	 * As per construction, the list of limit angles is known to have
	 * an even number of entries, with start angles at even indices and
	 * end angles at odd indices.
	 * </p>
	 * @param limits limit angles of the split part
	 * @return split part (may be NULL)
	 */
	Arcs_Set create_split_part(const std::vector<double>& limits)
	{
		if (limits.is_empty())
		{
			return NULL;
		}
		else
		{
			// collapse close limit angles
			for (int i{}; i < limits.size(); ++i)
			{
				const int j = (i + 1) % limits.size();
				const double lA = limits.get(i);
				const double lB = Math_Utils::normalize_angle(limits.get(j), lA);
				if (std::abs(lB - lA) <= get_tolerance())
				{
					// the two limits are too close to each other, we remove both of them
					if (j > 0)
					{
						// regular case, the two entries are consecutive ones
						limits.remove(j);
						limits.remove(i);
						i = i - 1;
					}
					else
					{
						// special case, i the the last entry and j is the first entry
						// we have wrapped around list end
						const double l_end = limits.remove(limits.size() - 1);
						const double l_start = limits.remove(0);
						if (limits.is_empty())
						{
							// the ends were the only limits, is it a full circle or an empty circle?
							if (l_end - l_start > std::numbers::pi)
							{
								// it was full circle
								return Arcs_Set(new BSP_Tree<Sphere_1D>(Boolean.TRUE), get_tolerance());
							}
							else
							{
								// it was an empty circle
								return NULL;
							}
						}
						else
						{
							// we have removed the first interval start, so our list
							// currently starts with an interval end, which is wrong
							// we need to move this interval end to the end of the list
							limits.add(limits.remove(0) + Math_Utils::TWO_PI);
						}
					}
				}
			}

			// build the tree by adding all angular sectors
			BSP_Tree<Sphere_1D> tree = BSP_Tree<>(false);
			for (int i{}; i < limits.size() - 1; i += 2)
			{
				add_arc_limit(tree, limits.get(i), true);
				add_arc_limit(tree, limits.get(i + 1), false);
			}

			if (tree.get_cut() == NULL)
			{
				// we did not insert anything
				return NULL;
			}

			return Arcs_Set(tree, get_tolerance());
		}
	}
};