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
  //package org.hipparchus.geometry.spherical.twod;

  //import java.util.Array_list;
  //import java.util.Arrays;
  //import java.util.Collection;
  //import java.util.Collections;
  //import java.util.List;
  //import java.util.function.Int_Predicate;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.geometry.Localized_Geometry_Formats;
  //import org.hipparchus.geometry.enclosing.Enclosing_Ball;
  //import org.hipparchus.geometry.enclosing.Welzl_Encloser;
  //import org.hipparchus.geometry.euclidean.threed.Euclidean_3D;
  //import org.hipparchus.geometry.euclidean.threed.Rotation;
  //import org.hipparchus.geometry.euclidean.threed.Rotation_Convention;
  //import org.hipparchus.geometry.euclidean.threed.Sphere_Generator;
  //import org.hipparchus.geometry.euclidean.threed.Vector_3D;
  //import org.hipparchus.geometry.partitioning.Abstract_Region;
  //import org.hipparchus.geometry.partitioning.BSP_Tree;
  //import org.hipparchus.geometry.partitioning.Boundary_Projection;
  //import org.hipparchus.geometry.partitioning.Region_Factory;
  //import org.hipparchus.geometry.partitioning.Sub_Hyperplane;
  //import org.hipparchus.geometry.spherical.oned.Arc;
  //import org.hipparchus.geometry.spherical.oned.Sphere_1D;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Precision;
#include <numbers>
#include <vector>

/** This class represents a region on the 2-sphere: a set of spherical polygons.
 */
class Spherical_Polygons_Set extends Abstract_Region<Sphere_2D, Sphere_1D>
{
	/** Boundary defined as an array of closed loops start vertices. */
	private List<Vertex> loops;

	/** Build a polygons set representing the whole real 2-sphere.
	 * @param tolerance below which points are consider to be identical
	 * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	public Spherical_Polygons_Set(const double& tolerance)
	{
		super(tolerance);
		Sphere_2D.check_tolerance(tolerance);
	}

	/** Build a polygons set representing a hemisphere.
	 * @param pole pole of the hemisphere (the pole is in the inside half)
	 * @param tolerance below which points are consider to be identical
	 * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	public Spherical_Polygons_Set(const Vector_3D pole, const double& tolerance)

	{
		super(new BSP_Tree<Sphere_2D>(new Circle(pole, tolerance).whole_hyperplane(), BSP_Tree<Sphere_2D>(Boolean.FALSE), BSP_Tree<Sphere_2D>(Boolean.TRUE), NULL), tolerance);
		Sphere_2D.check_tolerance(tolerance);
	}

	/** Build a polygons set representing a regular polygon.
	 * @param center center of the polygon (the center is in the inside half)
	 * @param meridian point defining the reference meridian for first polygon vertex
	 * @param outside_radius distance of the vertices to the center
	 * @param n number of sides of the polygon
	 * @param tolerance below which points are consider to be identical
	 * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	public Spherical_Polygons_Set(const Vector_3D center, const Vector_3D meridian, const double outside_radius, const int& n, const double& tolerance)

	{
		this(tolerance, create_regular_polygon_vertices(center, meridian, outside_radius, n));
	}

	/** Build a polygons set from a BSP tree.
	 * <p>The leaf nodes of the BSP tree <em>must</em> have a
	 * {@code Boolean} attribute representing the inside status of
	 * the corresponding cell (true for inside cells, false for outside
	 * cells). In order to avoid building too many small objects, it is
	 * recommended to use the predefined constants
	 * {@code Boolean.TRUE} and {@code Boolean.FALSE}</p>
	 * @param tree inside/outside BSP tree representing the region
	 * @param tolerance below which points are consider to be identical
	 * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	public Spherical_Polygons_Set(const BSP_Tree<Sphere_2D> tree, const double& tolerance)

	{
		super(tree, tolerance);
		Sphere_2D.check_tolerance(tolerance);
	}

	/** Build a polygons set from a Boundary RE_Presentation (B-rep).
	 * <p>The boundary is provided as a collection of {@link
	 * Sub_Hyperplane sub-hyperplanes}. Each sub-hyperplane has the
	 * interior part of the region on its minus side and the exterior on
	 * its plus side.</p>
	 * <p>The boundary elements can be in any order, and can form
	 * several non-connected sets (like for example polygons with holes
	 * or a set of disjoint polygons considered as a whole). In
	 * fact, the elements do not even need to be connected together
	 * (their topological connections are not used here). However, if the
	 * boundary does not really separate an inside open from an outside
	 * open (open having here its topological meaning), then subsequent
	 * calls to the {@link
	 * org.hipparchus.geometry.partitioning.Region#check_point(org.hipparchus.geometry.Point)
	 * check_point} method will not be meaningful anymore.</p>
	 * <p>If the boundary is empty, the region will represent the whole
	 * space.</p>
	 * @param boundary collection of boundary elements, as a
	 * collection of {@link Sub_Hyperplane Sub_Hyperplane} objects
	 * @param tolerance below which points are consider to be identical
	 * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 */
	public Spherical_Polygons_Set(const Collection<Sub_Hyperplane<Sphere_2D>> boundary, const double& tolerance)

	{
		super(boundary, tolerance);
		Sphere_2D.check_tolerance(tolerance);
	}

	/** Build a polygon from a simple list of vertices.
	 * <p>The boundary is provided as a list of points considering to
	 * represent the vertices of a simple loop. The interior part of the
	 * region is on the left side of this path and the exterior is on its
	 * right side.</p>
	 * <p>This constructor does not handle polygons with a boundary
	 * forming several disconnected paths (such as polygons with holes).</p>
	 * <p>For cases where this simple constructor applies, it is expected to
	 * be numerically more robust than the {@link #Spherical_Polygons_Set(Collection, * double) general constructor} using {@link Sub_Hyperplane subhyperplanes}.</p>
	 * <p>If the list is empty, the region will represent the whole
	 * space.</p>
	 * <p>This constructor assumes that edges between {@code vertices}, including the edge
	 * between the last and the first vertex, are shorter than pi. If edges longer than pi
	 * are used it may produce unintuitive results, such as reversing the direction of the
	 * edge. This implies using a {@code vertices} array of length 1 or 2 in this
	 * constructor produces an ill-defined region. Use one of the other constructors or
	 * {@link Region_Factory} instead.</p>
	 * <p>The list of {@code vertices} is reduced by selecting a sub-set of vertices
	 * before creating the boundary set. Every point in {@code vertices} will be on the
	 * {@link #check_point(org.hipparchus.geometry.Point) boundary} of the constructed polygon set, but not
	 * necessarily the center-line of the boundary.</p>
	 * <p>
	 * Polygons with thin pikes or dents are inherently difficult to handle because
	 * they involve circles with almost opposite directions at some vertices. Polygons
	 * whose vertices come from some physical measurement with noise are also
	 * difficult because an edge that should be straight may be broken in lots of
	 * different pieces with almost equal directions. In both cases, computing the
	 * circles intersections is not numerically robust due to the almost 0 or almost
	 * &pi; angle. Such cases need to carefully adjust the {@code hyperplane_thickness}
	 * parameter. A too small value would often lead to completely wrong polygons
	 * with large area wrongly identified as inside or outside. Large values are
	 * often much safer. As a rule of thumb, a value slightly below the size of the
	 * most accurate detail needed is a good value for the {@code hyperplane_thickness}
	 * parameter.
	 * </p>
	 * @param hyperplane_thickness tolerance below which points are considered to
	 * belong to the hyperplane (which is therefore more a slab). Should be greater than
	 * {@code FastMath.ulp(4 * std::numbers::pi)} for meaningful results.
	 * @param vertices vertices of the simple loop boundary
	 * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
	 * @exception org.hipparchus.exception.Math_Runtime_Exception if {@code vertices}
	 * contains only a single vertex or repeated vertices.
	 */
	public Spherical_Polygons_Set(const double hyperplane_thickness, const S2_Point ... vertices)

	{
		super(vertices_to_tree(hyperplane_thickness, vertices), hyperplane_thickness);
		Sphere_2D.check_tolerance(hyperplane_thickness);
	}

	/** Build the vertices representing a regular polygon.
	 * @param center center of the polygon (the center is in the inside half)
	 * @param meridian point defining the reference meridian for first polygon vertex
	 * @param outside_radius distance of the vertices to the center
	 * @param n number of sides of the polygon
	 * @return vertices array
	 */
	private static S2_Po std::vector<int> create_regular_polygon_vertices(const Vector_3D center, const Vector_3D meridian, const double outside_radius, const int& n)
	{
		const S2_Po std::vector<int> array = S2_Point[n];
		const Rotation r0 = Rotation(Vector_3D.cross_product(center, meridian), outside_radius, Rotation_Convention.VECTOR_OPERATOR);
		array[0] = S2_Point(r0.apply_to(center));

		const Rotation r = Rotation(center, Math_Utils::TWO_PI / n, Rotation_Convention.VECTOR_OPERATOR);
		for (int i{ 1 }; i < n; ++i)
		{
			array[i] = S2_Point(r.apply_to(array[i - 1].get_vector()));
		}

		return array;
	}

	/** Build the BSP tree of a polygons set from a simple list of vertices.
	 * <p>The boundary is provided as a list of points considering to
	 * represent the vertices of a simple loop. The interior part of the
	 * region is on the left side of this path and the exterior is on its
	 * right side.</p>
	 * <p>This constructor does not handle polygons with a boundary
	 * forming several disconnected paths (such as polygons with holes).</p>
	 * <p>This constructor handles only polygons with edges strictly shorter
	 * than \( \pi \). If longer edges are needed, they need to be broken up
	 * in smaller sub-edges so this constraint holds.</p>
	 * <p>For cases where this simple constructor applies, it is expected to
	 * be numerically more robust than the {@link #Polygons_Set(Collection) general
	 * constructor} using {@link Sub_Hyperplane subhyperplanes}.</p>
	 * @param hyperplane_thickness tolerance below which points are consider to
	 * belong to the hyperplane (which is therefore more a slab)
	 * @param vertices vertices of the simple loop boundary
	 * @return the BSP tree of the input vertices
	 */
	private static BSP_Tree<Sphere_2D> vertices_to_tree(const double hyperplane_thickness, S2_Point ... vertices)
	{
		// thin vertices to those that define distinct circles
		vertices = reduce(hyperplane_thickness, vertices).to_array(new S2_Point[0]);
		const int n = vertices.size();
		if (n == 0)
		{
			// the tree represents the whole space
			return BSP_Tree<Sphere_2D>(Boolean.TRUE);
		}

		// build the vertices
		const Vertex[] v_array = Vertex[n];
		for (int i{}; i < n; ++i)
		{
			v_array[i] = Vertex(vertices[i]);
		}

		// build the edges
		const List<Edge> edges = Array_list<>(n);
		Vertex end = v_array[n - 1];
		for (int i{}; i < n; ++i)
		{
			// get the endpoints of the edge
			const Vertex start = end;
			end = v_array[i];

			// get the circle supporting the edge
			const Circle& circle = Circle(start.get_location(), end.get_location(), hyperplane_thickness);

			// create the edge and store it
			edges.add(new Edge(start, end, Vector_3D.angle(start.get_location().get_vector(), end.get_location().get_vector()), circle));
		}

		// build the tree top-down
		const BSP_Tree<Sphere_2D> tree = BSP_Tree<>();
		insert_edges(hyperplane_thickness, tree, edges);

		return tree;
	}

	/**
	 * Compute a subset of vertices that define the boundary to within the given
	 * tolerance. This method partitions {@code vertices} into segments that all lie same
	 * arc to within {@code hyperplane_thickness}, and then returns the end points of the
	 * arcs. Combined arcs are limited to length of pi. If the input vertices has arcs
	 * longer than pi these will be preserved in the returned data.
	 *
	 * @param hyperplane_thickness of each circle in radians.
	 * @param vertices            to decimate.
	 * @return a subset of {@code vertices}.
	 */
	private static List<S2_Point> reduce(const double hyperplane_thickness, const S2_Po std::vector<int> vertices)
	{
		const int n = vertices.size();
		if (n <= 3)
		{
			// can't reduce to fewer than three points
			return Arrays.as_list(vertices.clone());
		}
		const List<S2_Point> points = Array_list<>();
		/* Use a simple greedy search to add points to a circle s.t. all intermediate
		 * points are within the thickness. Running time is O(n lg n) worst case.
		 * sin_ce the first vertex may be the middle of a straight edge, look backward
		 * and forward to establish the first edge.
		 * Uses the property that any two points define a circle, so don't check
		 * circles that just span two points.
		 */
		 // first look backward
		const Int_Predicate on_circle_backward = j ->
		{
			const int i = n - 2 - j;
			// circle spanning considered points
			const Circle& circle = Circle(vertices[0], vertices[i], hyperplane_thickness);
			const Arc arc = circle.get_arc(vertices[0], vertices[i]);
			if (arc.get_size() >= std::numbers::pi)
			{
				return false;
			}
			for (int k = i + 1; k < n; k++)
			{
				const S2_Point vertex = vertices[k];
				if (std::abs(circle.get_offset(vertex)) > hyperplane_thickness ||
					arc.get_offset(circle.to_sub_space(vertex)) > 0)
				{
					// point is not within the thickness or arc, start edge
					return false;
				}
			}
			return true;
		};
		// last index in vertices of last entry added to points
		int last = n - 2 - search_helper(on_circle_backward, 0, n - 2);
		if (last > 1)
		{
			points.add(vertices[last]);
		}
		else
		{
			// all points lie on one semi-circle, distance from 0 to 1 is > pi
			// ill-defined case, just return three points from the list
			return Arrays.as_list(Arrays.copy_of_range(vertices, 0, 3));
		}
		const int first = last;
		// then build edges forward
		for (int j{ 1 }; ; j += 2)
		{
			const int last_final = last;
			const Int_Predicate on_circle = i ->
			{
				// circle spanning considered points
				const Circle& circle = Circle(vertices[last_final], vertices[i], hyperplane_thickness);
				const Arc arc = circle.get_arc(vertices[last_final], vertices[i]);
				if (arc.get_size() >= std::numbers::pi)
				{
					return false;
				}
				const int end = last_final < i ? i : i + n;
				for (int k = last_final + 1; k < end; k++)
				{
					const S2_Point vertex = vertices[k % n];
					if (std::abs(circle.get_offset(vertex)) > hyperplane_thickness ||
						arc.get_offset(circle.to_sub_space(vertex)) > 0)
					{
						// point is not within the thickness or arc, start edge
						return false;
					}
				}
				return true;
			};
			j = search_helper(on_circle, j, first + 1);
			if (j >= first)
			{
				break;
			}
			last = j;
			points.add(vertices[last]);
		}
		// put first point last
		const S2_Point swap = points.remove(0);
		points.add(swap);
		return points;
	}

	/**
	 * Search {@code items} for the first item where {@code predicate} is false between
	 * {@code a} and {@code b}. Assumes that predicate switches from true to false at
	 * exactly one location in [a, b]. Similar to {@link Arrays#binary_search(std::vector<int>, int, * int, int)} except that 1. it operates on indices, not elements, 2. there is not a
	 * shortcut for equality, and 3. it is optimized for cases where the return value is
	 * close to a.
	 *
	 * <p> This method achieves O(lg n) performance in the worst case, where n = b - a.
	 * Performance improves to O(lg(i-a)) when i is close to a, where i is the return
	 * value.
	 *
	 * @param predicate to apply.
	 * @param a         start, inclusive.
	 * @param b         end, exclusive.
	 * @return a if a==b, a-1 if predicate.test(a) == false, b - 1 if predicate.test(b-1), * otherwise i s.t. predicate.test(i) == true && predicate.test(i + 1) == false.
	 * @ if a > b.
	 */
	private static int search_helper(const Int_Predicate predicate, const int a, const int b)
	{
		if (a > b)
		{
			throw (
				hipparchus::exception::Localized_Core_Formats_Type::LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT, a, b);
		}
		// Argument checks and special cases
		if (a == b)
		{
			return a;
		}
		if (!predicate.test(a))
		{
			return a - 1;
		}

		// start with exponential search
		int start = a;
		int end = b;
		for (int i{ 2 }; a + i < b; i *= 2)
		{
			if (predicate.test(a + i))
			{
				// update lower bound of search
				start = a + i;
			}
			else
			{
				// found upper bound of search
				end = a + i;
				break;
			}
		}

		// next binary search
		// copied from Arrays.binary_search() and modified to work on indices alone
		int low = start;
		int high = end - 1;
		while (low <= high)
		{
			const int mid = (low + high) >> > 1;
			if (predicate.test(mid))
			{
				low = mid + 1;
			}
			else
			{
				high = mid - 1;
			}
		}
		// low is now insertion point, according to Arrays.binary_search()
		return low - 1;
	}

	/** Recursively build a tree by inserting cut sub-hyperplanes.
	 * @param hyperplane_thickness tolerance below which points are considered to
	 * belong to the hyperplane (which is therefore more a slab)
	 * @param node current tree node (it is a leaf node at the beginning
	 * of the call)
	 * @param edges list of edges to insert in the cell defined by this node
	 * (excluding edges not belonging to the cell defined by this node)
	 */
	private static void insert_edges(const double hyperplane_thickness, const BSP_Tree<Sphere_2D> node, const List<Edge> edges)
	{
		// find an edge with an hyperplane that can be inserted in the node
		int index = 0;
		Edge inserted = NULL;
		while (inserted == NULL && index < edges.size())
		{
			inserted = edges.get(index++);
			if (!node.insert_cut(inserted.get_circle()))
			{
				inserted = NULL;
			}
		}

		if (inserted == NULL)
		{
			// no suitable edge was found, the node remains a leaf node
			// we need to set its inside/outside bool indicator
			const BSP_Tree<Sphere_2D> parent = node.get_parent();
			if (parent == NULL || node == parent.get_minus())
			{
				node.set_attribute(Boolean.TRUE);
			}
			else
			{
				node.set_attribute(Boolean.FALSE);
			}
			return;
		}

		// we have split the node by inserting an edge as a cut sub-hyperplane
		// distribute the remaining edges in the two sub-trees
		const List<Edge> outside_list = Array_list<>();
		const List<Edge> inside_list = Array_list<>();
		for (const Edge edge : edges)
		{
			if (edge != inserted)
			{
				edge.split(inserted.get_circle(), outside_list, inside_list);
			}
		}

		// recurse through lower levels
		if (!outside_list.is_empty())
		{
			insert_edges(hyperplane_thickness, node.get_plus(), outside_list);
		}
		else
		{
			node.get_plus().set_attribute(Boolean.FALSE);
		}
		if (!inside_list.is_empty())
		{
			insert_edges(hyperplane_thickness, node.get_minus(), inside_list);
		}
		else
		{
			node.get_minus().set_attribute(Boolean.TRUE);
		}
	}

	/** {@inherit_doc} */
	//override
	public Spherical_Polygons_Set build_new(const BSP_Tree<Sphere_2D> tree)
	{
		return Spherical_Polygons_Set(tree, get_tolerance());
	}

	/** {@inherit_doc}
	 * @exception Math_Illegal_State_Exception if the tolerance setting does not allow to build
	 * a clean non-ambiguous boundary
	 */
	 //override
	protected void compute_geometrical_properties() Math_Illegal_State_Exception
	{
		const BSP_Tree<Sphere_2D> tree = get_tree(true);

		if (tree.get_cut() == NULL)
		{
			// the instance has a single cell without any boundaries

			if (tree.get_cut() == NULL && (Boolean)tree.get_attribute())
			{
				// the instance covers the whole space
				set_size(4 * std::numbers::pi);
				set_barycenter(new S2_Point(0, 0));
			}
			else
			{
				set_size(0);
				set_barycenter(S2_Point.NaN);
			}
		}
		else
		{
			// the instance has a boundary
			const Properties_Computer pc = Properties_Computer(get_tolerance());
			tree.visit(pc);
			set_size(pc.get_area());
			set_barycenter(pc.get_barycenter());
		}
	}

	/** Get the boundary loops of the polygon.
	 * <p>The polygon boundary can be represented as a list of closed loops, * each loop being given by exactly one of its vertices. From each loop
	 * start vertex, one can follow the loop by finding the outgoing edge, * then the end vertex, then the next outgoing edge ... until the start
	 * vertex of the loop (exactly the same instance) is found again once
	 * the full loop has been visited.</p>
	 * <p>If the polygon has no boundary at all, a zero length loop
	 * array will be returned.</p>
	 * <p>If the polygon is a simple one-piece polygon, then the returned
	 * array will contain a single vertex.
	 * </p>
	 * <p>All edges in the various loops have the inside of the region on
	 * their left side (i.e. toward their pole) and the outside on their
	 * right side (i.e. away from their pole) when moving in the underlying
	 * circle direction. This means that the closed loops obey the direct
	 * trigonometric orientation.</p>
	 * @return boundary of the polygon, organized as an unmodifiable list of loops start vertices.
	 * @exception Math_Illegal_State_Exception if the tolerance setting does not allow to build
	 * a clean non-ambiguous boundary
	 * @see Vertex
	 * @see Edge
	 */
	public List<Vertex> get_boundary_loops() Math_Illegal_State_Exception
	{
		if (loops == NULL)
		{
			if (get_tree(false).get_cut() == NULL)
			{
				loops = Collections.empty_list();
			}
			else
			{
				// sort the arcs according to their start point
				const Edges_With_Node_Info_Builder visitor = Edges_With_Node_Info_Builder(get_tolerance());
				get_tree(true).visit(visitor);
				const std::vector<Edge_With_Node_Info>& edges = visitor.get_edges();

				// connect all edges, using topological criteria first
				// and using Euclidean distance only as a last resort
				int pending = edges.size();
				pending -= natural_follower_connections(edges);
				if (pending > 0)
				{
					pending -= split_edge_connections(edges);
				}
				if (pending > 0)
				{
					close_vertices_connections(edges);
				}

				// extract the edges loops
				loops = Array_list<>();
				for (Edge_With_Node_Info s = get_unprocessed(edges); s != NULL; s = get_unprocessed(edges))
				{
					loops.add(s.get_start());
					follow_loop(s);
				}
			}
		}

		return Collections.unmodifiable_list(loops);
	}

	/** Connect the edges using only natural follower information.
	 * @param edges edges complete edges list
	 * @return number of connections performed
	 */
	private int natural_follower_connections(const std::vector<Edge_With_Node_Info>& edges)
	{
		int connected{};
		for (const auto& edge : edges)
		{
			if (edge.get_end().get_outgoing() == NULL)
			{
				for (const auto& candidate_next : edges)
				{
					if (Edge_With_Node_Info.are_natural_followers(edge, candidate_next))
					{
						// connect the two edges
						edge.set_next_edge(candidate_next);
						++connected;
						break;
					}
				}
			}
		}
		return connected;
	}

	/** Connect the edges resulting from a circle splitting a circular edge.
	 * @param edges edges complete edges list
	 * @return number of connections performed
	 */
	private int split_edge_connections(const std::vector<Edge_With_Node_Info>& edges)
	{
		int connected{};
		for (const auto& edge : edges)
		{
			if (edge.get_end().get_outgoing() == NULL)
			{
				for (const auto& candidate_next : edges)
				{
					if (Edge_With_Node_Info.result_from_a_split(edge, candidate_next))
					{
						// connect the two edges
						edge.set_next_edge(candidate_next);
						++connected;
						break;
					}
				}
			}
		}
		return connected;
	}

	/** Connect the edges using spherical distance.
	 * <p>
	 * This connection heuristic should be used last, as it relies
	 * only on a fuzzy distance criterion.
	 * </p>
	 * @param edges edges complete edges list
	 * @return number of connections performed
	 */
	private int close_vertices_connections(const std::vector<Edge_With_Node_Info>& edges)
	{
		int connected{};
		for (const auto& edge : edges)
		{
			if (edge.get_end().get_outgoing() == NULL && edge.get_end() != NULL)
			{
				const Vector_3D end = edge.get_end().get_location().get_vector();
				Edge_With_Node_Info selected_next = NULL;
				double min = INFINITY;
				for (const auto& candidate_next : edges)
				{
					if (candidate_next.get_start().get_incoming() == NULL)
					{
						const auto distance = Vector_3D.distance(end, candidate_next.get_start().get_location().get_vector());
						if (distance < min)
						{
							selected_next = candidate_next;
							min = distance;
						}
					}
				}
				if (min <= get_tolerance())
				{
					// connect the two edges
					edge.set_next_edge(selected_next);
					++connected;
				}
			}
		}
		return connected;
	}

	/** Get first unprocessed edge from a list.
	 * @param edges edges list
	 * @return first edge that has not been processed yet
	 * or NULL if all edges have been processed
	 */
	private Edge_With_Node_Info get_unprocessed(const std::vector<Edge_With_Node_Info>& edges)
	{
		for (const auto& edge : edges)
		{
			if (!edge.is_processed())
			{
				return edge;
			}
		}
		return NULL;
	}

	/** Build the loop containing a edge.
	 * <p>
	 * All edges put in the loop will be marked as processed.
	 * </p>
	 * @param defining edge used to define the loop
	 */
	private void follow_loop(const Edge_With_Node_Info& defining)
	{
		defining.set_processed(true);

		// process edges in connection order
		Edge_With_Node_Info previous = defining;
		Edge_With_Node_Info next = (Edge_With_Node_Info)defining.get_end().get_outgoing();
		while (next != defining)
		{
			if (next == NULL)
			{
				// this should not happen
				throw Math_Illegal_State_Exception(Localized_Geometry_Formats.OUTLINE_BOUNDARY_LOOP_OPEN);
			}
			next.set_processed(true);

			// filter out spurious vertices
			if (Vector_3D.angle(previous.get_circle().get_pole(), next.get_circle().get_pole()) <= Precision.EPSILON)
			{
				// the vertex between the two edges is a spurious one
				// replace the two edges by a single one
				previous.set_next_edge(next.get_end().get_outgoing());
				previous.set_length(previous.get_length() + next.get_length());
			}

			previous = next;
			next = (Edge_With_Node_Info)next.get_end().get_outgoing();
		}
	}

	/** Get a spherical cap enclosing the polygon.
	 * <p>
	 * This method is intended as a first test to quickly identify points
	 * that are guaranteed to be outside of the region, hence performing a full
	 * {@link #check_point(org.hipparchus.geometry.Vector) check_point}
	 * only if the point status remains undecided after the quick check. It is
	 * is therefore mostly useful to speed up computation for small polygons with
	 * complex shapes (say a country boundary on Earth), as the spherical cap will
	 * be small and hence will reliably identify a large part of the sphere as outside, * whereas the full check can be more computing intensive. A typical use case is
	 * therefore:
	 * </p>
	 * <pre>
	 *   // compute region, plus an enclosing spherical cap
	 *   Spherical_Polygons_Set complex_shape = ...;
	 *   Enclosing_Ball&lt;Sphere_2D, S2_Point&gt; cap = complex_shape.get_enclosing_cap();
	 *
	 *   // check lots of points
	 *   for (Vector_3D p : points)
	 {
	 *
	 *     const Location l;
	 *     if (cap.contains(p))
	 {
	 *       // we cannot be sure where the point is
	 *       // we need to perform the full computation
	 *       l = complex_shape.check_point(v);
	 *     }
else
	 {
	 *       // no need to do further computation, *       // we already know the point is outside
	 *       l = Location.OUTSIDE;
	 *     }
	 *
	 *     // use l ...
	 *
	 *   }
	 * </pre>
	 * <p>
	 * In the special cases of empty or whole sphere polygons, special
	 * spherical caps are returned, with angular radius set to negative
	 * or positive infinity so the {@link
	 * Enclosing_Ball#contains(org.hipparchus.geometry.Point) ball.contains(point)}
	 * method return always false or true.
	 * </p>
	 * <p>
	 * This method is <em>not</em> guaranteed to return the smallest enclosing cap.
	 * </p>
	 * @return a spherical cap enclosing the polygon
	 */
	public Enclosing_Ball<Sphere_2D, S2_Point> get_enclosing_cap()
	{
		// handle special cases first
		if (is_empty())
		{
			return Enclosing_Ball<Sphere_2D, S2_Point>(S2_Point.PLUS_K, -INFINITY);
		}
		if (is_full())
		{
			return Enclosing_Ball<Sphere_2D, S2_Point>(S2_Point.PLUS_K, INFINITY);
		}

		// as the polygons is neither empty nor full, it has some boundaries and cut hyperplanes
		const BSP_Tree<Sphere_2D> root = get_tree(false);
		if (is_empty(root.get_minus()) && is_full(root.get_plus()))
		{
			// the polygon covers an hemisphere, and its boundary is one 2π long edge
			const Circle& circle = (Circle)root.get_cut().get_hyperplane();
			return Enclosing_Ball<Sphere_2D, S2_Point>(new S2_Point(circle.get_pole()).negate(), 0.5 * std::numbers::pi);
		}
		if (is_full(root.get_minus()) && is_empty(root.get_plus()))
		{
			// the polygon covers an hemisphere, and its boundary is one 2π long edge
			const Circle& circle = (Circle)root.get_cut().get_hyperplane();
			return Enclosing_Ball<Sphere_2D, S2_Point>(new S2_Point(circle.get_pole()), 0.5 * std::numbers::pi);
		}

		// gather some inside points, to be used by the encloser
		const List<Vector_3D> points = get_inside_points();

		// extract points from the boundary loops, to be used by the encloser as well
		const List<Vertex> boundary = get_boundary_loops();
		for (const auto& loop_start : boundary)
		{
			int count{};
			for (Vertex v = loop_start; count == 0 || v != loop_start; v = v.get_outgoing().get_end())
			{
				++count;
				points.add(v.get_location().get_vector());
			}
		}

		// find the smallest enclosing 3D sphere
		const Sphere_Generator generator = Sphere_Generator();
		const Welzl_Encloser<Euclidean_3D, Vector_3D> encloser =
			Welzl_Encloser<>(get_tolerance(), generator);
		Enclosing_Ball<Euclidean_3D, Vector_3D> enclosing_3d = encloser.enclose(points);
		const Vector_3D[] support_3d = enclosing_3d.get_support();

		// convert to 3D sphere to spherical cap
		const double r = enclosing_3d.get_radius();
		const double h = enclosing_3d.get_center().get_norm();
		if (h < get_tolerance())
		{
			// the 3D sphere is centered on the unit sphere and covers it
			// fall back to a crude approximation, based only on outside convex cells
			Enclosing_Ball<Sphere_2D, S2_Point> enclosing_s2 =
				Enclosing_Ball<>(S2_Point.PLUS_K, INFINITY);
			for (Vector_3D outside_point : get_outside_points())
			{
				const S2_Point outside_s2 = S2_Point(outside_point);
				const Boundary_Projection<Sphere_2D> projection = project_to_boundary(outside_s2);
				if (std::numbers::pi - projection.get_offset() < enclosing_s2.get_radius())
				{
					enclosing_s2 = Enclosing_Ball<>(outside_s2.negate(), std::numbers::pi - projection.get_offset(), (S2_Point)projection.get_projected());
				}
			}
			return enclosing_s2;
		}
		const S2_Po std::vector<int> support = std::vector<S2_Point>(support_3d.size());
		for (int i{}; i < support_3d.size(); ++i)
		{
			support[i] = S2_Point(support_3d[i]);
		}

		return Enclosing_Ball<>(new S2_Point(enclosing_3d.get_center()), std::acos((1 + h * h - r * r) / (2 * h)), support);
	}

	/** Gather some inside points.
	 * @return list of points known to be strictly in all inside convex cells
	 */
	private List<Vector_3D> get_inside_points()
	{
		const Properties_Computer pc = Properties_Computer(get_tolerance());
		get_tree(true).visit(pc);
		return pc.get_convex_cells_inside_points();
	}

	/** Gather some outside points.
	 * @return list of points known to be strictly in all outside convex cells
	 */
	private List<Vector_3D> get_outside_points()
	{
		const Spherical_Polygons_Set complement =
			(Spherical_Polygons_Set)Region_Factory<Sphere_2D>().get_complement(this);
		const Properties_Computer pc = Properties_Computer(get_tolerance());
		complement.get_tree(true).visit(pc);
		return pc.get_convex_cells_inside_points();
	}
}
