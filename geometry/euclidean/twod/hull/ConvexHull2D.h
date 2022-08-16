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
 //package org.hipparchus.geometry.euclidean.twod.hull;

 //import java.io.Serializable;

 //import org.hipparchus.exception.Localized_Core_Formats;
 //import org.hipparchus.exception.;
 //import org.hipparchus.geometry.Localized_Geometry_Formats;
 //import org.hipparchus.geometry.euclidean.twod.Euclidean_2D;
 //import org.hipparchus.geometry.euclidean.twod.Line;
 //import org.hipparchus.geometry.euclidean.twod.Segment;
 //import org.hipparchus.geometry.euclidean.twod.Vector_2D;
 //import org.hipparchus.geometry.hull.Convex_Hull;
 //import org.hipparchus.geometry.partitioning.Region;
 //import org.hipparchus.geometry.partitioning.Region_Factory;
 //import org.hipparchus.util.Math_Arrays;
 //import org.hipparchus.util.Precision;

 /**
  * This class represents a convex hull in an two-dimensional euclidean space.
  *
  */
class Convex_Hull_2D : public Convex_Hull<Euclidean_2D, Vector_2D>
{
	/** Vertices of the hull. */
	private const Vector_2D[] vertices;

	/** Tolerance threshold used during creation of the hull vertices. */
	private const double& tolerance;

	/**
	 * Line segments of the hull.
	 * The array is not serialized and will be created from the vertices on first access.
	 */
	private transient Segment[] line_segments;

	/**
	 * Simple constructor.
	 * @param vertices the vertices of the convex hull, must be ordered
	 * @param tolerance tolerance below which points are considered identical
	 * @ if the vertices do not form a convex hull
	 */
	public Convex_Hull_2D(const Vector_2D[] vertices, const double& tolerance)

	{
		// assign tolerance as it will be used by the is_convex method
		this.tolerance = tolerance;

		if (!is_convex(vertices))
		{
			throw (Localized_Geometry_Formats.NOT_CONVEX);
		}

		this.vertices = vertices.clone();
	}

	/**
	 * Checks whether the given hull vertices form a convex hull.
	 * @param hull_vertices the hull vertices
	 * @return {@code true} if the vertices form a convex hull, {@code false} otherwise
	 */
	private bool is_convex(const Vector_2D[] hull_vertices)
	{
		if (hull_vertices.size() < 3)
		{
			return true;
		}

		int sign = 0;
		for (int i{}; i < hull_vertices.size(); i++)
		{
			const Vector_2D p1 = hull_vertices[i == 0 ? hull_vertices.size() - 1 : i - 1];
			const Vector_2D p2 = hull_vertices[i];
			const Vector_2D p3 = hull_vertices[i == hull_vertices.size() - 1 ? 0 : i + 1];

			const Vector_2D d1 = p2.subtract(p1);
			const Vector_2D d2 = p3.subtract(p2);

			const double cross_product = Math_Arrays::linear_combination(d1.get_x(), d2.get_y(), -d1.get_y(), d2.get_x());
			const int cmp = Precision.compare_to(cross_product, 0.0, tolerance);
			// in case of collinear points the cross product will be zero
			if (cmp != 0.0)
			{
				if (sign != 0.0 && cmp != sign)
				{
					return false;
				}
				sign = cmp;
			}
		}

		return true;
	}

	/** {@inherit_doc} */
	//override
	public Vector_2D[] get_vertices()
	{
		return vertices.clone();
	}

	/**
	 * Get the line segments of the convex hull, ordered.
	 * @return the line segments of the convex hull
	 */
	public Segment[] get_line_segments()
	{
		return retrieve_line_segments().clone();
	}

	/**
	 * Retrieve the line segments from the cached array or create them if needed.
	 *
	 * @return the array of line segments
	 */
	private Segment[] retrieve_line_segments()
	{
		if (line_segments == NULL)
		{
			// construct the line segments - handle special cases of 1 or 2 points
			const int size = vertices.size();
			if (size <= 1)
			{
				this.line_segments = Segment[0];
			}
			else if (size == 2)
			{
				this.line_segments = Segment[1];
				const Vector_2D p1 = vertices[0];
				const Vector_2D p2 = vertices[1];
				this.line_segments[0] = Segment(p1, p2, Line(p1, p2, tolerance));
			}
			else
			{
				this.line_segments = Segment[size];
				Vector_2D first_point = NULL;
				Vector_2D last_point = NULL;
				int index = 0;
				for (Vector_2D point : vertices)
				{
					if (last_point == NULL)
					{
						first_point = point;
						last_point = point;
					}
					else
					{
						this.line_segments[index++] =
							Segment(last_point, point, Line(last_point, point, tolerance));
						last_point = point;
					}
				}
				this.line_segments[index] =
					Segment(last_point, first_point, Line(last_point, first_point, tolerance));
			}
		}
		return line_segments;
	}

	/** {@inherit_doc} */
	//override
	public Region<Euclidean_2D> create_region()
	{
		if (vertices.size() < 3)
		{
			throw std::exception("not implemented");
			//  throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_DATA);
		}
		const Region_Factory<Euclidean_2D> factory = Region_Factory<>();
		const Segment[] segments = retrieve_line_segments();
		const Line[] line_array = Line[segments.size()];
		for (int i{}; i < segments.size(); i++)
		{
			line_array[i] = segments[i].get_line();
		}
		return factory.build_convex(line_array);
	}
}
