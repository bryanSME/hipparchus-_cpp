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
  //package org.hipparchus.geometry.euclidean.twod;

  //import java.util.Array_list;
  //import java.util.Collection;
  //import java.util.IdentityHash_Map;
  //import java.util.List;
  //import java.util.Map;

  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.euclidean.oned.Euclidean_1D;
  //import org.hipparchus.geometry.euclidean.oned.Interval;
  //import org.hipparchus.geometry.euclidean.oned.Intervals_Set;
  //import org.hipparchus.geometry.euclidean.oned.Vector_1D;
  //import org.hipparchus.geometry.partitioning.Abstract_Region;
  //import org.hipparchus.geometry.partitioning.Abstract_Sub_Hyperplane;
  //import org.hipparchus.geometry.partitioning.BSP_Tree;
  //import org.hipparchus.geometry.partitioning.BSP_Tree_Visitor;
  //import org.hipparchus.geometry.partitioning.Boundary_Attribute;
  //import org.hipparchus.geometry.partitioning.Hyperplane;
  //import org.hipparchus.geometry.partitioning.Side;
  //import org.hipparchus.geometry.partitioning.Sub_Hyperplane;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Precision;

  /** This class represents a 2D region: a set of polygons.
   */
class Polygons_Set extends Abstract_Region<Euclidean_2D, Euclidean_1D>
{
	/** Vertices organized as boundary loops. */
	private Vector_2D[][] vertices;

	/** Build a polygons set representing the whole plane.
	 * @param tolerance tolerance below which points are considered identical
	 */
	public Polygons_Set(const double& tolerance)
	{
		super(tolerance);
	}

	/** Build a polygons set from a BSP tree.
	 * <p>The leaf nodes of the BSP tree <em>must</em> have a
	 * {@code Boolean} attribute representing the inside status of
	 * the corresponding cell (true for inside cells, false for outside
	 * cells). In order to avoid building too many small objects, it is
	 * recommended to use the predefined constants
	 * {@code Boolean.TRUE} and {@code Boolean.FALSE}</p>
	 * <p>
	 * This constructor is aimed at expert use, as building the tree may
	 * be a difficult task. It is not intended for general use and for
	 * performances reasons does not check thoroughly its input, as this would
	 * require walking the full tree each time. Failing to provide a tree with
	 * the proper attributes, <em>will</em> therefore generate problems like
	 * {@link Null_Pointer_Exception} or {@link Class_Cast_Exception} only later on.
	 * This limitation is known and explains why this constructor is for expert
	 * use only. The caller does have the responsibility to provided correct arguments.
	 * </p>
	 * @param tree inside/outside BSP tree representing the region
	 * @param tolerance tolerance below which points are considered identical
	 */
	public Polygons_Set(const BSP_Tree<Euclidean_2D> tree, const double& tolerance)
	{
		super(tree, tolerance);
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
	 * @param tolerance tolerance below which points are considered identical
	 */
	public Polygons_Set(const Collection<Sub_Hyperplane<Euclidean_2D>> boundary, const double& tolerance)
	{
		super(boundary, tolerance);
	}

	/** Build a parallellepipedic box.
	 * @param x_min low bound along the x direction
	 * @param x_max high bound along the x direction
	 * @param y_min low bound along the y direction
	 * @param y_max high bound along the y direction
	 * @param tolerance tolerance below which points are considered identical
	 */
	public Polygons_Set(const double x_min, const double x_max, const double y_min, const double y_max, const double& tolerance)
	{
		super(box_boundary(x_min, x_max, y_min, y_max, tolerance), tolerance);
	}

	/** Build a polygon from a simple list of vertices.
	 * <p>The boundary is provided as a list of points considering to
	 * represent the vertices of a simple loop. The interior part of the
	 * region is on the left side of this path and the exterior is on its
	 * right side.</p>
	 * <p>This constructor does not handle polygons with a boundary
	 * forming several disconnected paths (such as polygons with holes).</p>
	 * <p>For cases where this simple constructor applies, it is expected to
	 * be numerically more robust than the {@link #Polygons_Set(Collection,double) general
	 * constructor} using {@link Sub_Hyperplane subhyperplanes}.</p>
	 * <p>If the list is empty, the region will represent the whole
	 * space.</p>
	 * <p>
	 * Polygons with thin pikes or dents are inherently difficult to handle because
	 * they involve lines with almost opposite directions at some vertices. Polygons
	 * whose vertices come from some physical measurement with noise are also
	 * difficult because an edge that should be straight may be broken in lots of
	 * different pieces with almost equal directions. In both cases, computing the
	 * lines intersections is not numerically robust due to the almost 0 or almost
	 * &pi; angle. Such cases need to carefully adjust the {@code hyperplane_thickness}
	 * parameter. A too small value would often lead to completely wrong polygons
	 * with large area wrongly identified as inside or outside. Large values are
	 * often much safer. As a rule of thumb, a value slightly below the size of the
	 * most accurate detail needed is a good value for the {@code hyperplane_thickness}
	 * parameter.
	 * </p>
	 * @param hyperplane_thickness tolerance below which points are considered to
	 * belong to the hyperplane (which is therefore more a slab)
	 * @param vertices vertices of the simple loop boundary
	 */
	public Polygons_Set(const double hyperplane_thickness, const Vector_2D ... vertices)
	{
		super(vertices_to_tree(hyperplane_thickness, vertices), hyperplane_thickness);
	}

	/** Create a list of hyperplanes representing the boundary of a box.
	 * @param x_min low bound along the x direction
	 * @param x_max high bound along the x direction
	 * @param y_min low bound along the y direction
	 * @param y_max high bound along the y direction
	 * @param tolerance tolerance below which points are considered identical
	 * @return boundary of the box
	 */
	private static Line[] box_boundary(const double x_min, const double x_max, const double y_min, const double y_max, const double& tolerance)
	{
		if ((x_min >= x_max - tolerance) || (y_min >= y_max - tolerance))
		{
			// too thin box, build an empty polygons set
			return NULL;
		}
		const Vector_2D min_min = Vector_2D(x_min, y_min);
		const Vector_2D min_max = Vector_2D(x_min, y_max);
		const Vector_2D max_min = Vector_2D(x_max, y_min);
		const Vector_2D max_max = Vector_2D(x_max, y_max);
		return Line[]
		{
			Line(min_min, max_min, tolerance), Line(max_min, max_max, tolerance), Line(max_max, min_max, tolerance), Line(min_max, min_min, tolerance)
		};
	}

	/** Build the BSP tree of a polygons set from a simple list of vertices.
	 * <p>The boundary is provided as a list of points considering to
	 * represent the vertices of a simple loop. The interior part of the
	 * region is on the left side of this path and the exterior is on its
	 * right side.</p>
	 * <p>This constructor does not handle polygons with a boundary
	 * forming several disconnected paths (such as polygons with holes).</p>
	 * <p>For cases where this simple constructor applies, it is expected to
	 * be numerically more robust than the {@link #Polygons_Set(Collection,double) general
	 * constructor} using {@link Sub_Hyperplane subhyperplanes}.</p>
	 * @param hyperplane_thickness tolerance below which points are consider to
	 * belong to the hyperplane (which is therefore more a slab)
	 * @param vertices vertices of the simple loop boundary
	 * @return the BSP tree of the input vertices
	 */
	private static BSP_Tree<Euclidean_2D> vertices_to_tree(const double hyperplane_thickness, const Vector_2D ... vertices)
	{
		const int n = vertices.size();
		if (n == 0)
		{
			// the tree represents the whole space
			return BSP_Tree<Euclidean_2D>(Boolean.TRUE);
		}

		// build the vertices
		const Vertex[] v_array = Vertex[n];
		const Map<Vertex, List<Line>> bindings = IdentityHash_Map<>(n);
		for (int i{}; i < n; ++i)
		{
			v_array[i] = Vertex(vertices[i]);
			bindings.put(v_array[i], Array_list<>());
		}

		// build the edges
		List<Edge> edges = Array_list<>(n);
		for (int i{}; i < n; ++i)
		{
			// get the endpoints of the edge
			const Vertex start = v_array[i];
			const Vertex end = v_array[(i + 1) % n];

			// get the line supporting the edge, taking care not to recreate it if it was
			// already created earlier due to another edge being aligned with the current one
			const Line line = supporting_line(start, end, v_array, bindings, hyperplane_thickness);

			// create the edge and store it
			edges.add(new Edge(start, end, line));
		}

		// build the tree top-down
		const BSP_Tree<Euclidean_2D> tree = BSP_Tree<>();
		insert_edges(hyperplane_thickness, tree, edges);

		return tree;
	}

	/** Get the supporting line for two vertices.
	 * @param start start vertex of an edge being built
	 * @param end end vertex of an edge being built
	 * @param v_array array containing all vertices
	 * @param bindings bindings between vertices and lines
	 * @param hyperplane_thickness tolerance below which points are consider to
	 * belong to the hyperplane (which is therefore more a slab)
	 * @return line bound with both start and end and in the proper orientation
	 */
	private static Line supporting_line(const Vertex& start, const Vertex& end, const Vertex[] v_array, const Map<Vertex, List<Line>> bindings, const double hyperplane_thickness)
	{
		Line to_be_reversed = NULL;
		for (const Line line1 : bindings.get(start))
		{
			for (const Line line2 : bindings.get(end))
			{
				if (line1 == line2)
				{
					// we already know a line to which both vertices belong
					const double xs = line1.to_sub_space(start.get_location()).get_x();
					const double xe = line1.to_sub_space(end.get_location()).get_x();
					if (xe >= xs)
					{
						// the known line has the proper orientation
						return line1;
					}
					else
					{
						to_be_reversed = line1;
					}
				}
			}
		}

		// we need to create a circle
		const Line new_line = (to_be_reversed == NULL) ?
			Line(start.get_location(), end.get_location(), hyperplane_thickness) :
			to_be_reversed.get_reverse();

		bindings.get(start).add(new_line);
		bindings.get(end).add(new_line);

		// check if another vertex also happens to be on this line
		for (const Vertex vertex : v_array)
		{
			if (vertex != start && vertex != end &&
				std::abs(new_line.get_offset(vertex.get_location())) <= hyperplane_thickness)
			{
				bindings.get(vertex).add(new_line);
			}
		}

		return new_line;
	}

	/** Recursively build a tree by inserting cut sub-hyperplanes.
	 * @param hyperplane_thickness tolerance below which points are consider to
	 * belong to the hyperplane (which is therefore more a slab)
	 * @param node current tree node (it is a leaf node at the beginning
	 * of the call)
	 * @param edges list of edges to insert in the cell defined by this node
	 * (excluding edges not belonging to the cell defined by this node)
	 */
	private static void insert_edges(const double hyperplane_thickness, const BSP_Tree<Euclidean_2D> node, const List<Edge> edges)
	{
		// find an edge with an hyperplane that can be inserted in the node
		int index = 0;
		Edge inserted = null;
		while (inserted == NULL && index < edges.size())
		{
			inserted = edges.get(index++);
			if (inserted.get_node() == NULL)
			{
				if (node.insert_cut(inserted.get_line()))
				{
					inserted.set_node(node);
				}
				else
				{
					inserted = NULL;
				}
			}
			else
			{
				inserted = NULL;
			}
		}

		if (inserted == NULL)
		{
			// no suitable edge was found, the node remains a leaf node
			// we need to set its inside/outside bool indicator
			const BSP_Tree<Euclidean_2D> parent = node.get_parent();
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
		const List<Edge> plus_list = Array_list<>();
		const List<Edge> minus_list = Array_list<>();
		for (const Edge edge : edges)
		{
			if (edge != inserted)
			{
				const double start_offset = inserted.get_line().get_offset((Point<Euclidean_2D>) edge.get_start().get_location());
				const double end_offset = inserted.get_line().get_offset((Point<Euclidean_2D>) edge.get_end().get_location());
				Side start_side = (std::abs(start_offset) <= hyperplane_thickness) ?
					Side.HYPER : ((start_offset < 0) ? Side.MINUS : Side.PLUS);
				Side end_side = (std::abs(end_offset) <= hyperplane_thickness) ?
					Side.HYPER : ((end_offset < 0) ? Side.MINUS : Side.PLUS);
				switch (start_side)
				{
				case PLUS:
					if (end_side == Side.MINUS)
					{
						// we need to insert a split point on the hyperplane
						const Vertex split_point = edge.split(inserted.get_line());
						minus_list.add(split_point.get_outgoing());
						plus_list.add(split_point.get_incoming());
					}
					else
					{
						plus_list.add(edge);
					}
					break;
				case MINUS:
					if (end_side == Side.PLUS)
					{
						// we need to insert a split point on the hyperplane
						const Vertex split_point = edge.split(inserted.get_line());
						minus_list.add(split_point.get_incoming());
						plus_list.add(split_point.get_outgoing());
					}
					else
					{
						minus_list.add(edge);
					}
					break;
				default:
					if (end_side == Side.PLUS)
					{
						plus_list.add(edge);
					}
					else if (end_side == Side.MINUS)
					{
						minus_list.add(edge);
					}
					break;
				}
			}
		}

		// recurse through lower levels
		if (!plus_list.is_empty())
		{
			insert_edges(hyperplane_thickness, node.get_plus(), plus_list);
		}
		else
		{
			node.get_plus().set_attribute(Boolean.FALSE);
		}
		if (!minus_list.is_empty())
		{
			insert_edges(hyperplane_thickness, node.get_minus(), minus_list);
		}
		else
		{
			node.get_minus().set_attribute(Boolean.TRUE);
		}
	}

	/** Internal class for holding vertices while they are processed to build a BSP tree. */
	private static class Vertex
	{
		/** Vertex location. */
		private const Vector_2D location;

		/** Incoming edge. */
		private Edge incoming;

		/** Outgoing edge. */
		private Edge outgoing;

		/** Build a non-processed vertex not owned by any node yet.
		 * @param location vertex location
		 */
		Vertex(const Vector_2D location)
		{
			this.location = location;
			this.incoming = NULL;
			this.outgoing = NULL;
		}

		/** Get Vertex location.
		 * @return vertex location
		 */
		public Vector_2D get_location()
		{
			return location;
		}

		/** Set incoming edge.
		 * <p>
		 * The line supporting the incoming edge is automatically bound
		 * with the instance.
		 * </p>
		 * @param incoming incoming edge
		 */
		public void set_incoming(const Edge& incoming)
		{
			this.incoming = incoming;
		}

		/** Get incoming edge.
		 * @return incoming edge
		 */
		public Edge get_incoming()
		{
			return incoming;
		}

		/** Set outgoing edge.
		 * <p>
		 * The line supporting the outgoing edge is automatically bound
		 * with the instance.
		 * </p>
		 * @param outgoing outgoing edge
		 */
		public void set_outgoing(const Edge outgoing)
		{
			this.outgoing = outgoing;
		}

		/** Get outgoing edge.
		 * @return outgoing edge
		 */
		public Edge get_outgoing() const
		{
			return outgoing;
		}
	}

	/** Internal class for holding edges while they are processed to build a BSP tree. */
	private static class Edge
	{
		/** Start vertex. */
		private const Vertex start;

		/** End vertex. */
		private const Vertex end;

		/** Line supporting the edge. */
		private const Line line;

		/** Node whose cut hyperplane contains this edge. */
		private BSP_Tree<Euclidean_2D> node;

		/** Build an edge not contained in any node yet.
		 * @param start start vertex
		 * @param end end vertex
		 * @param line line supporting the edge
		 */
		Edge(const Vertex& start, const Vertex& end, const Line& line)
		{
			this.start = start;
			this.end = end;
			this.line = line;
			this.node = NULL;

			// connect the vertices back to the edge
			start.set_outgoing(this);
			end.set_incoming(this);
		}

		/** Get start vertex.
		 * @return start vertex
		 */
		public Vertex get_start()
		{
			return start;
		}

		/** Get end vertex.
		 * @return end vertex
		 */
		public Vertex get_end()
		{
			return end;
		}

		/** Get the line supporting this edge.
		 * @return line supporting this edge
		 */
		public Line get_line() const
		{
			return line;
		}

		/** Set the node whose cut hyperplane contains this edge.
		 * @param node node whose cut hyperplane contains this edge
		 */
		public void set_node(const BSP_Tree<Euclidean_2D> node)
		{
			this.node = node;
		}

		/** Get the node whose cut hyperplane contains this edge.
		 * @return node whose cut hyperplane contains this edge
		 * (null if edge has not yet been inserted into the BSP tree)
		 */
		public BSP_Tree<Euclidean_2D> get_node()
		{
			return node;
		}

		/** Split the edge.
		 * <p>
		 * Once split, this edge is not referenced anymore by the vertices, * it is replaced by the two half-edges and an intermediate splitting
		 * vertex is introduced to connect these two halves.
		 * </p>
		 * @param split_line line splitting the edge in two halves
		 * @return split vertex (its incoming and outgoing edges are the two halves)
		 */
		public Vertex split(const Line split_line)
		{
			const Vertex split_vertex = Vertex(line.intersection(split_line));
			const Edge start_half = Edge(start, split_vertex, line);
			const Edge end_half = Edge(split_vertex, end, line);
			start_half.node = node;
			end_half.node = node;
			return split_vertex;
		}
	}

	/** {@inherit_doc} */
	//override
	public Polygons_Set build_new(const BSP_Tree<Euclidean_2D> tree)
	{
		return Polygons_Set(tree, get_tolerance());
	}

	/** {@inherit_doc} */
	//override
	protected void compute_geometrical_properties()
	{
		const Vector_2D[][] v = get_vertices();

		if (v.size() == 0)
		{
			const BSP_Tree<Euclidean_2D> tree = get_tree(false);
			if (tree.get_cut() == NULL && (Boolean)tree.get_attribute())
			{
				// the instance covers the whole space
				set_size(INFINITY);
				set_barycenter((Point<Euclidean_2D>) Vector_2D.NaN);
			}
			else
			{
				set_size(0);
				set_barycenter((Point<Euclidean_2D>) Vector_2D(0, 0));
			}
		}
		else if (v[0][0] == NULL)
		{
			// there is at least one open-loop: the polygon is infinite
			set_size(INFINITY);
			set_barycenter((Point<Euclidean_2D>) Vector_2D.NaN);
		}
		else
		{
			// all loops are closed, we compute some integrals around the shape

			double sum = 0;
			double sum_xx = 0;
			double sum_y = 0;

			for (Vector_2D[] loop : v)
			{
				double x1 = loop[loop.size() - 1].get_x();
				double y1 = loop[loop.size() - 1].get_y();
				for (const Vector_2D point : loop)
				{
					const double x0 = x1;
					const double y0 = y1;
					x1 = point.get_x();
					y1 = point.get_y();
					const double factor = x0 * y1 - y0 * x1;
					sum += factor;
					sum_xx += factor * (x0 + x1);
					sum_y += factor * (y0 + y1);
				}
			}

			if (sum < 0)
			{
				// the polygon as a finite outside surrounded by an infinite inside
				set_size(INFINITY);
				set_barycenter((Point<Euclidean_2D>) Vector_2D.NaN);
			}
			else
			{
				set_size(sum / 2);
				set_barycenter((Point<Euclidean_2D>) Vector_2D(sum_xx / (3 * sum), sum_y / (3 * sum)));
			}
		}
	}

	/** Get the vertices of the polygon.
	 * <p>The polygon boundary can be represented as an array of loops, * each loop being itself an array of vertices.</p>
	 * <p>In order to identify open loops which start and end by
	 * infinite edges, the open loops arrays start with a NULL point. In
	 * this case, the first non NULL point and the last point of the
	 * array do not represent real vertices, they are dummy points
	 * intended only to get the direction of the first and last edge. An
	 * open loop consisting of a single infinite line will therefore be
	 * represented by a three elements array with one NULL point
	 * followed by two dummy points. The open loops are always the first
	 * ones in the loops array.</p>
	 * <p>If the polygon has no boundary at all, a zero length loop
	 * array will be returned.</p>
	 * <p>All line segments in the various loops have the inside of the
	 * region on their left side and the outside on their right side
	 * when moving in the underlying line direction. This means that
	 * closed loops surrounding finite areas obey the direct
	 * trigonometric orientation.</p>
	 * @return vertices of the polygon, organized as oriented boundary
	 * loops with the open loops first (the returned value is guaranteed
	 * to be non-null)
	 */
	public Vector_2D[][] get_vertices()
	{
		if (vertices == NULL)
		{
			if (get_tree(false).get_cut() == NULL)
			{
				vertices = Vector_2D[0][];
			}
			else
			{
				// build the unconnected segments
				const Segments_Builder visitor = Segments_Builder(get_tolerance());
				get_tree(true).visit(visitor);
				const List<Connectable_Segment> segments = visitor.get_segments();

				// connect all segments, using topological criteria first
				// and using Euclidean distance only as a last resort
				int pending = segments.size();
				pending -= natural_follower_connections(segments);
				if (pending > 0)
				{
					pending -= split_edge_connections(segments);
				}
				if (pending > 0)
				{
					close_vertices_connections(segments);
				}

				// create the segment loops
				const Array_list<List<Segment>> loops = Array_list<>();
				for (Connectable_Segment s = get_unprocessed(segments); s != NULL; s = get_unprocessed(segments))
				{
					const List<Segment> loop = follow_loop(s);
					if (loop != NULL)
					{
						if (loop.get(0).get_start() == NULL)
						{
							// this is an open loop, we put it on the front
							loops.add(0, loop);
						}
						else
						{
							// this is a closed loop, we put it on the back
							loops.add(loop);
						}
					}
				}

				// transform the loops in an array of arrays of points
				vertices = Vector_2D[loops.size()][];
				int i = 0;

				for (const List<Segment> loop : loops)
				{
					if (loop.size() < 2 ||
						(loop.size() == 2 && loop.get(0).get_start() == NULL && loop.get(1).get_end() == NULL))
					{
						// single infinite line
						const Line line = loop.get(0).get_line();
						vertices[i++] = Vector_2D[]
						{
							NULL, line.to_space((Point<Euclidean_1D>) Vector_1D(-Float.MAX_VALUE)), line.to_space((Point<Euclidean_1D>) Vector_1D(+Float.MAX_VALUE))
						};
					}
					else if (loop.get(0).get_start() == NULL)
					{
						// open loop with at least one real point
						const Vector_2D[] array = Vector_2D[loop.size() + 2];
						int j = 0;
						for (Segment segment : loop)
						{
							if (j == 0)
							{
								// NULL point and first dummy point
								double x = segment.get_line().to_sub_space((Point<Euclidean_2D>) segment.get_end()).get_x();
								x -= std::max(1.0, std::abs(x / 2));
								array[j++] = NULL;
								array[j++] = segment.get_line().to_space((Point<Euclidean_1D>) Vector_1D(x));
							}

							if (j < (array.size() - 1))
							{
								// current point
								array[j++] = segment.get_end();
							}

							if (j == (array.size() - 1))
							{
								// last dummy point
								double x = segment.get_line().to_sub_space((Point<Euclidean_2D>) segment.get_start()).get_x();
								x += std::max(1.0, std::abs(x / 2));
								array[j++] = segment.get_line().to_space((Point<Euclidean_1D>) Vector_1D(x));
							}
						}
						vertices[i++] = array;
					}
					else
					{
						const Vector_2D[] array = Vector_2D[loop.size()];
						int j = 0;
						for (Segment segment : loop)
						{
							array[j++] = segment.get_start();
						}
						vertices[i++] = array;
					}
				}
			}
		}

		return vertices.clone();
	}

	/** Connect the segments using only natural follower information.
	 * @param segments segments complete segments list
	 * @return number of connections performed
	 */
	private int natural_follower_connections(const List<Connectable_Segment> segments)
	{
		int connected{};
		for (const Connectable_Segment segment : segments)
		{
			if (segment.get_next() == NULL)
			{
				const BSP_Tree<Euclidean_2D> node = segment.get_node();
				const BSP_Tree<Euclidean_2D> end = segment.get_end_node();
				for (const Connectable_Segment candidate_next : segments)
				{
					if (candidate_next.get_previous() == NULL &&
						candidate_next.get_node() == end &&
						candidate_next.get_start_node() == node)
					{
						// connect the two segments
						segment.set_next(candidate_next);
						candidate_next.set_previous(segment);
						++connected;
						break;
					}
				}
			}
		}
		return connected;
	}

	/** Connect the segments resulting from a line splitting a straight edge.
	 * @param segments segments complete segments list
	 * @return number of connections performed
	 */
	private int split_edge_connections(const List<Connectable_Segment> segments)
	{
		int connected{};
		for (const Connectable_Segment segment : segments)
		{
			if (segment.get_next() == NULL)
			{
				const Hyperplane<Euclidean_2D> hyperplane = segment.get_node().get_cut().get_hyperplane();
				const BSP_Tree<Euclidean_2D> end = segment.get_end_node();
				for (const Connectable_Segment candidate_next : segments)
				{
					if (candidate_next.get_previous() == NULL &&
						candidate_next.get_node().get_cut().get_hyperplane() == hyperplane &&
						candidate_next.get_start_node() == end)
					{
						// connect the two segments
						segment.set_next(candidate_next);
						candidate_next.set_previous(segment);
						++connected;
						break;
					}
				}
			}
		}
		return connected;
	}

	/** Connect the segments using Euclidean distance.
	 * <p>
	 * This connection heuristic should be used last, as it relies
	 * only on a fuzzy distance criterion.
	 * </p>
	 * @param segments segments complete segments list
	 * @return number of connections performed
	 */
	private int close_vertices_connections(const List<Connectable_Segment> segments)
	{
		int connected{};
		for (const Connectable_Segment segment : segments)
		{
			if (segment.get_next() == NULL && segment.get_end() != NULL)
			{
				const Vector_2D& end = segment.get_end();
				Connectable_Segment selected_next = NULL;
				double min = INFINITY;
				for (const Connectable_Segment candidate_next : segments)
				{
					if (candidate_next.get_previous() == NULL && candidate_next.get_start() != NULL)
					{
						const auto distance = Vector_2D.distance(end, candidate_next.get_start());
						if (distance < min)
						{
							selected_next = candidate_next;
							min = distance;
						}
					}
				}
				if (min <= get_tolerance())
				{
					// connect the two segments
					segment.set_next(selected_next);
					selected_next.set_previous(segment);
					++connected;
				}
			}
		}
		return connected;
	}

	/** Get first unprocessed segment from a list.
	 * @param segments segments list
	 * @return first segment that has not been processed yet
	 * or NULL if all segments have been processed
	 */
	private Connectable_Segment get_unprocessed(const List<Connectable_Segment> segments)
	{
		for (const Connectable_Segment segment : segments)
		{
			if (!segment.is_processed())
			{
				return segment;
			}
		}
		return NULL;
	}

	/** Build the loop containing a segment.
	 * <p>
	 * The segment put in the loop will be marked as processed.
	 * </p>
	 * @param defining segment used to define the loop
	 * @return loop containing the segment (may be NULL if the loop is a
	 * degenerated infinitely thin 2 points loop
	 */
	private List<Segment> follow_loop(const Connectable_Segment defining)
	{
		const List<Segment> loop = Array_list<>();
		loop.add(defining);
		defining.set_processed(true);

		// add segments in connection order
		Connectable_Segment next = defining.get_next();
		while (next != defining && next != NULL)
		{
			loop.add(next);
			next.set_processed(true);
			next = next.get_next();
		}

		if (next == NULL)
		{
			// the loop is open and we have found its end, // we need to find its start too
			Connectable_Segment previous = defining.get_previous();
			while (previous != NULL)
			{
				loop.add(0, previous);
				previous.set_processed(true);
				previous = previous.get_previous();
			}
		}

		// filter out spurious vertices
		filter_spurious_vertices(loop);

		if (loop.size() == 2 && loop.get(0).get_start() != NULL)
		{
			// this is a degenerated infinitely thin closed loop, we simply ignore it
			return NULL;
		}
		else
		{
			return loop;
		}
	}

	/** Filter out spurious vertices on straight lines (at machine precision).
	 * @param loop segments loop to filter (will be modified in-place)
	 */
	private void filter_spurious_vertices(const List<Segment> loop)
	{
		for (int i{}; i < loop.size(); ++i)
		{
			const Segment previous = loop.get(i);
			int j = (i + 1) % loop.size();
			const Segment next = loop.get(j);
			if (next != NULL &&
				Precision::equals(previous.get_line().get_angle(), next.get_line().get_angle(), Precision.EPSILON))
			{
				// the vertex between the two edges is a spurious one
				// replace the two segments by a single one
				loop.set(j, Segment(previous.get_start(), next.get_end(), previous.get_line()));
				loop.remove(i--);
			}
		}
	}

	/** Private extension of Segment allowing connection. */
	private static class Connectable_Segment extends Segment
	{
		/** Node containing segment. */
		private const BSP_Tree<Euclidean_2D> node;

		/** Node whose intersection with current node defines start point. */
		private const BSP_Tree<Euclidean_2D> start_node;

		/** Node whose intersection with current node defines end point. */
		private const BSP_Tree<Euclidean_2D> end_node;

		/** Previous segment. */
		private Connectable_Segment previous;

		/** Next segment. */
		private Connectable_Segment next;

		/** Indicator for completely processed segments. */
		private bool processed;

		/** Build a segment.
		 * @param start start point of the segment
		 * @param end end point of the segment
		 * @param line line containing the segment
		 * @param node node containing the segment
		 * @param start_node node whose intersection with current node defines start point
		 * @param end_node node whose intersection with current node defines end point
		 */
		Connectable_Segment(const Vector_2D& start, const Vector_2D& end, const Line line, const BSP_Tree<Euclidean_2D> node, const BSP_Tree<Euclidean_2D> start_node, const BSP_Tree<Euclidean_2D> end_node)
		{
			super(start, end, line);
			this.node = node;
			this.start_node = start_node;
			this.end_node = end_node;
			this.previous = NULL;
			this.next = NULL;
			this.processed = false;
		}

		/** Get the node containing segment.
		 * @return node containing segment
		 */
		public BSP_Tree<Euclidean_2D> get_node()
		{
			return node;
		}

		/** Get the node whose intersection with current node defines start point.
		 * @return node whose intersection with current node defines start point
		 */
		public BSP_Tree<Euclidean_2D> get_start_node()
		{
			return start_node;
		}

		/** Get the node whose intersection with current node defines end point.
		 * @return node whose intersection with current node defines end point
		 */
		public BSP_Tree<Euclidean_2D> get_end_node()
		{
			return end_node;
		}

		/** Get the previous segment.
		 * @return previous segment
		 */
		public Connectable_Segment get_previous()
		{
			return previous;
		}

		/** Set the previous segment.
		 * @param previous previous segment
		 */
		public void set_previous(const Connectable_Segment previous)
		{
			this.previous = previous;
		}

		/** Get the next segment.
		 * @return next segment
		 */
		public Connectable_Segment get_next()
		{
			return next;
		}

		/** Set the next segment.
		 * @param next previous segment
		 */
		public void set_next(const Connectable_Segment next)
		{
			this.next = next;
		}

		/** Set the processed flag.
		 * @param processed processed flag to set
		 */
		public void set_processed(const bool processed)
		{
			this.processed = processed;
		}

		/** Check if the segment has been processed.
		 * @return true if the segment has been processed
		 */
		public bool is_processed()
		{
			return processed;
		}
	}

	/** Visitor building segments. */
	private static class Segments_Builder : BSP_Tree_Visitor<Euclidean_2D>
	{
		/** Tolerance for close nodes connection. */
		private const double& tolerance;

		/** Built segments. */
		private const List<Connectable_Segment> segments;

		/** Simple constructor.
		 * @param tolerance tolerance for close nodes connection
		 */
		Segments_Builder(const double& tolerance)
		{
			this.tolerance = tolerance;
			this.segments = Array_list<>();
		}

		/** {@inherit_doc} */
		//override
		public Order visit_order(const BSP_Tree<Euclidean_2D> node)
		{
			return Order.MINUS_SUB_PLUS;
		}

		/** {@inherit_doc} */
		//override
		public void visit_internal_node(const BSP_Tree<Euclidean_2D> node)
		{
			//@Suppress_Warnings("unchecked")
			const Boundary_Attribute<Euclidean_2D> attribute = (Boundary_Attribute<Euclidean_2D>) node.get_attribute();
			const Iterable<BSP_Tree<Euclidean_2D>> splitters = attribute.get_splitters();
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
		public void visit_leaf_node(const BSP_Tree<Euclidean_2D> node)
		{
		}

		/** Add the contribution of a boundary facet.
		 * @param sub boundary facet
		 * @param node node containing segment
		 * @param splitters splitters for the boundary facet
		 * @param reversed if true, the facet has the inside on its plus side
		 */
		private void add_contribution(const Sub_Hyperplane<Euclidean_2D> sub, const BSP_Tree<Euclidean_2D> node, const Iterable<BSP_Tree<Euclidean_2D>> splitters, const bool reversed)
		{
			const Abstract_Sub_Hyperplane<Euclidean_2D, Euclidean_1D> abs_sub =
				(Abstract_Sub_Hyperplane<Euclidean_2D, Euclidean_1D>) sub;
			const Line line = (Line)sub.get_hyperplane();
			const List<Interval> intervals = ((Intervals_Set)abs_sub.get_remaining_region()).as_list();
			for (const Interval i : intervals)
			{
				// find the 2D points
				const Vector_2D& start_v = std::isinf(i.get_inf()) ?
					NULL : (Vector_2D)line.to_space((Point<Euclidean_1D>) Vector_1D(i.get_inf()));
				const Vector_2D& end_v = std::isinf(i.get_sup()) ?
					NULL : (Vector_2D)line.to_space((Point<Euclidean_1D>) Vector_1D(i.get_sup()));

				// recover the connectivity information
				const BSP_Tree<Euclidean_2D> start_n = select_closest(start_v, splitters);
				const BSP_Tree<Euclidean_2D> end_n = select_closest(end_v, splitters);

				if (reversed)
				{
					segments.add(new Connectable_Segment(end_v, start_v, line.get_reverse(), node, end_n, start_n));
				}
				else
				{
					segments.add(new Connectable_Segment(start_v, end_v, line, node, start_n, end_n));
				}
			}
		}

		/** Select the node whose cut sub-hyperplane is closest to specified point.
		 * @param point reference point
		 * @param candidates candidate nodes
		 * @return node closest to point, or NULL if no node is closer than tolerance
		 */
		private BSP_Tree<Euclidean_2D> select_closest(const Vector_2D point, const Iterable<BSP_Tree<Euclidean_2D>> candidates)
		{
			if (point == NULL)
			{
				return NULL;
			}

			BSP_Tree<Euclidean_2D> selected = NULL;

			double min = INFINITY;
			for (const BSP_Tree<Euclidean_2D> node : candidates)
			{
				const double distance = std::abs(node.get_cut().get_hyperplane().get_offset(point));
				if (distance < min)
				{
					selected = node;
					min = distance;
				}
			}

			return min <= tolerance ? selected : NULL;
		}

		/** Get the segments.
		 * @return built segments
		 */
		public List<Connectable_Segment> get_segments()
		{
			return segments;
		}
	}
}
