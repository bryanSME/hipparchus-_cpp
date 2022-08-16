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
  //package org.hipparchus.geometry.euclidean.threed;

  //import java.util.Array_list;

  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.euclidean.twod.Euclidean_2D;
  //import org.hipparchus.geometry.euclidean.twod.Polygons_Set;
  //import org.hipparchus.geometry.euclidean.twod.Vector_2D;
  //import org.hipparchus.geometry.partitioning.Abstract_Sub_Hyperplane;
  //import org.hipparchus.geometry.partitioning.BSP_Tree;
  //import org.hipparchus.geometry.partitioning.BSP_Tree_Visitor;
  //import org.hipparchus.geometry.partitioning.Boundary_Attribute;
  //import org.hipparchus.geometry.partitioning.Region_Factory;
  //import org.hipparchus.geometry.partitioning.Sub_Hyperplane;
  //import org.hipparchus.util.FastMath;

  /** Extractor for {@link Polygons_Set polyhedrons sets} outlines.
   * <p>This class extracts the 2D outlines from {{@link Polygons_Set
   * polyhedrons sets} in a specified projection plane.</p>
   */
class Outline_Extractor
{
	/** Abscissa axis of the projection plane. */
	private const Vector_3D u;

	/** Ordinate axis of the projection plane. */
	private const Vector_3D v;

	/** Normal of the projection plane (viewing direction). */
	private const Vector_3D w;

	/** Build an extractor for a specific projection plane.
	 * @param u abscissa axis of the projection point
	 * @param v ordinate axis of the projection point
	 */
	public Outline_Extractor(const Vector_3D u, const Vector_3D v)
	{
		this.u = u;
		this.v = v;
		w = Vector_3D.cross_product(u, v);
	}

	/** Extract the outline of a polyhedrons set.
	 * @param polyhedrons_set polyhedrons set whose outline must be extracted
	 * @return an outline, as an array of loops.
	 */
	public Vector_2D[][] get_outline(const Polyhedrons_Set polyhedrons_set)
	{
		// project all boundary facets into one polygons set
		const Boundary_Projector projector = Boundary_Projector(polyhedrons_set.get_tolerance());
		polyhedrons_set.get_tree(true).visit(projector);
		const Polygons_Set projected = projector.get_projected();

		// Remove the spurious intermediate vertices from the outline
		const Vector_2D[][] outline = projected.get_vertices();
		for (int i{}; i < outline.size(); ++i)
		{
			const Vector_2D[] raw_loop = outline[i];
			int end = raw_loop.size();
			int j = 0;
			while (j < end)
			{
				if (point_is_between(raw_loop, end, j))
				{
					// the point should be removed
					for (int k = j; k < (end - 1); ++k)
					{
						raw_loop[k] = raw_loop[k + 1];
					}
					--end;
				}
				else
				{
					// the point remains in the loop
					++j;
				}
			}
			if (end != raw_loop.size())
			{
				// resize the array
				outline[i] = Vector_2D[end];
				System.arraycopy(raw_loop, 0, outline[i], 0, end);
			}
		}

		return outline;
	}

	/** Check if a point is geometrically between its neighbor in an array.
	 * <p>The neighbors are computed considering the array is a loop
	 * (i.e. point at index (n-1) is before point at index 0)</p>
	 * @param loop points array
	 * @param n number of points to consider in the array
	 * @param i index of the point to check (must be between 0 and n-1)
	 * @return true if the point is exactly between its neighbors
	 */
	private bool point_is_between(const Vector_2D[] loop, const int& n, const int i)
	{
		const Vector_2D previous = loop[(i + n - 1) % n];
		const Vector_2D current = loop[i];
		const Vector_2D next = loop[(i + 1) % n];
		const double dx1 = current.get_x() - previous.get_x();
		const double dy1 = current.get_y() - previous.get_y();
		const double dx2 = next.get_x() - current.get_x();
		const double dy2 = next.get_y() - current.get_y();
		const double cross = dx1 * dy2 - dx2 * dy1;
		const double dot = dx1 * dx2 + dy1 * dy2;
		const double d1d2 = std::sqrt((dx1 * dx1 + dy1 * dy1) * (dx2 * dx2 + dy2 * dy2));
		return (std::abs(cross) <= (1.0e-6 * d1d2)) && (dot >= 0.0);
	}

	/** Visitor projecting the boundary facets on a plane. */
	private class Boundary_Projector : BSP_Tree_Visitor<Euclidean_3D>
	{
		/** Projection of the polyhedrons set on the plane. */
		private Polygons_Set projected;

		/** Tolerance below which points are considered identical. */
		private const double& tolerance;

		/** Simple constructor.
		 * @param tolerance tolerance below which points are considered identical
		 */
		Boundary_Projector(const double& tolerance)
		{
			this.projected = Polygons_Set(new BSP_Tree<Euclidean_2D>(Boolean.FALSE), tolerance);
			this.tolerance = tolerance;
		}

		/** {@inherit_doc} */
		//override
		public Order visit_order(const BSP_Tree<Euclidean_3D> node)
		{
			return Order.MINUS_SUB_PLUS;
		}

		/** {@inherit_doc} */
		//override
		public void visit_internal_node(const BSP_Tree<Euclidean_3D> node)
		{
			//@Suppress_Warnings("unchecked")
			const Boundary_Attribute<Euclidean_3D> attribute =
				(Boundary_Attribute<Euclidean_3D>) node.get_attribute();
			if (attribute.get_plus_outside() != NULL)
			{
				add_contribution(attribute.get_plus_outside());
			}
			if (attribute.get_plus_inside() != NULL)
			{
				add_contribution(attribute.get_plus_inside());
			}
		}

		/** {@inherit_doc} */
		//override
		public void visit_leaf_node(const BSP_Tree<Euclidean_3D> node)
		{
		}

		/** Add he contribution of a boundary facet.
		 * @param facet boundary facet
		 */
		private void add_contribution(const Sub_Hyperplane<Euclidean_3D> facet)
		{
			// extract the vertices of the facet
			const Abstract_Sub_Hyperplane<Euclidean_3D, Euclidean_2D> abs_facet =
				(Abstract_Sub_Hyperplane<Euclidean_3D, Euclidean_2D>) facet;
			const Plane plane = (Plane)facet.get_hyperplane();

			const double scal = plane.get_normal().dot_product(w);
			if (std::abs(scal) > 1.0e-3)
			{
				Vector_2D[][] vertices =
					((Polygons_Set)abs_facet.get_remaining_region()).get_vertices();

				if (scal < 0)
				{
					// the facet is seen from the back of the plane, // we need to invert its boundary orientation
					const Vector_2D[][] new_vertices = Vector_2D[vertices.size()][];
					for (int i{}; i < vertices.size(); ++i)
					{
						const Vector_2D[] loop = vertices[i];
						const Vector_2D[] new_loop = Vector_2D[loop.size()];
						if (loop[0] == NULL)
						{
							new_loop[0] = NULL;
							for (int j{ 1 }; j < loop.size(); ++j)
							{
								new_loop[j] = loop[loop.size() - j];
							}
						}
						else
						{
							for (int j{}; j < loop.size(); ++j)
							{
								new_loop[j] = loop[loop.size() - (j + 1)];
							}
						}
						new_vertices[i] = new_loop;
					}

					// use the reverted vertices
					vertices = new_vertices;
				}

				// compute the projection of the facet in the outline plane
				const Array_list<Sub_Hyperplane<Euclidean_2D>> edges = Array_list<>();
				for (Vector_2D[] loop : vertices)
				{
					const bool closed = loop[0] != NULL;
					int previous = closed ? (loop.size() - 1) : 1;
					const Vector_3D previous_3d = plane.to_space((Point<Euclidean_2D>) loop[previous]);
					int current = (previous + 1) % loop.size();
					Vector_2D p_point = Vector_2D(previous_3d.dot_product(u), previous_3d.dot_product(v));
					while (current < loop.size())
					{
						const Vector_3D current_3d = plane.to_space((Point<Euclidean_2D>) loop[current]);
						const Vector_2D  c_point = Vector_2D(current_3d.dot_product(u), current_3d.dot_product(v));
						const org.hipparchus.geometry.euclidean.twod.Line line =
							org.hipparchus.geometry.euclidean.twod.Line(p_point, c_point, tolerance);
						Sub_Hyperplane<Euclidean_2D> edge = line.whole_hyperplane();

						if (closed || (previous != 1))
						{
							// the previous point is a real vertex
							// it defines one bounding point of the edge
							const double& angle = line.get_angle() + 0.5 * std::numbers::pi;
							const org.hipparchus.geometry.euclidean.twod.Line l =
								org.hipparchus.geometry.euclidean.twod.Line(p_point, angle, tolerance);
							edge = edge.split(l).get_plus();
						}

						if (closed || (current != (loop.size() - 1)))
						{
							// the current point is a real vertex
							// it defines one bounding point of the edge
							const double& angle = line.get_angle() + 0.5 * std::numbers::pi;
							const org.hipparchus.geometry.euclidean.twod.Line l =
								org.hipparchus.geometry.euclidean.twod.Line(c_point, angle, tolerance);
							edge = edge.split(l).get_minus();
						}

						edges.add(edge);

						previous = current++;
						p_point = c_point;
					}
				}
				const Polygons_Set projected_facet = Polygons_Set(edges, tolerance);

				// add the contribution of the facet to the global outline
				projected = (Polygons_Set)Region_Factory<Euclidean_2D>().union(projected, projected_facet);
			}
		}

		/** Get the projection of the polyhedrons set on the plane.
		 * @return projection of the polyhedrons set on the plane
		 */
		public Polygons_Set get_projected()
		{
			return projected;
		}
	}
}
