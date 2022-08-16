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

 //import java.util.Collection;

 //import org.hipparchus.exception.Localized_Core_Formats;
 //import org.hipparchus.exception.;
 //import org.hipparchus.exception.Math_Illegal_State_Exception;
 //import org.hipparchus.geometry.euclidean.twod.Vector_2D;
 //import org.hipparchus.util.Math_Utils;
#include <vector>
#include "ConvexHullGenerator2D.h"

/**
 * Abstract base class for convex hull generators in the two-dimensional euclidean space.
 *
 */
class Abstract_Convex_Hull_Generator_2D : public Convex_Hull_Generator_2D
{
private:
	/** Default value for tolerance. */
	static constexpr double DEFAULT_TOLERANCE{ 1e-10 };

	/** Tolerance below which points are considered identical. */
	const double my_tolerance;

	/**
	 * Indicates if collinear points on the hull shall be present in the output.
	 * If {@code false}, only the extreme points are added to the hull.
	 */
	const bool my_include_collinear_points;

protected:
	/**
	 * Simple constructor.
	 * <p>
	 * The default tolerance (1e-10) will be used to determine identical points.
	 *
	 * @param include_collinear_points indicates if collinear points on the hull shall be
	 * added as hull vertices
	 */
	Abstract_Convex_Hull_Generator_2D(const bool include_collinear_points)
	{
		Abstract_Convex_Hull_Generator_2D(include_collinear_points, DEFAULT_TOLERANCE);
	}

	/**
	 * Simple constructor.
	 *
	 * @param include_collinear_points indicates if collinear points on the hull shall be
	 * added as hull vertices
	 * @param tolerance tolerance below which points are considered identical
	 */
	Abstract_Convex_Hull_Generator_2D(const bool include_collinear_points, const double& tolerance)
		my_include_collinear_points {
		include_collinear_points
	},
		my_tolerance{ tolerance }
		{}

		/**
		 * Find the convex hull vertices from the set of input points.
		 * @param points the set of input points
		 * @return the convex hull vertices in CCW winding
		 */
		virtual Collection<Vector_2D> find_hull_vertices(const std::vector<Vector_2D>& points);

public:
	/**
	 * Get the tolerance below which points are considered identical.
	 * @return the tolerance below which points are considered identical
	 */
	double get_tolerance() const
	{
		return my_tolerance;
	}

	/**
	 * Returns if collinear points on the hull will be added as hull vertices.
	 * @return {@code true} if collinear points are added as hull vertices, or {@code false}
	 * if only extreme points are present.
	 */
	bool is_include_collinear_points() const
	{
		return my_include_collinear_points;
	}

	/** {@inherit_doc} */
	//override
	Convex_Hull_2D generate(const std::vector<Vector_2D> points)
	{
		// check for NULL points
		//Math_Utils::check_not_null(points);

		std::vector<Vector_2D> hull_vertices;
		if (points.size() < 2)
		{
			hull_vertices = points;
		}
		else
		{
			hull_vertices = find_hull_vertices(points);
		}

		try
		{
			return Convex_Hull_2D(hull_vertices.to_array(std::vector<Vector_2D>{}), my_tolerance);
		}
		catch (e)
		{
			// the hull vertices may not form a convex hull if the tolerance value is to large
			throw std::exception("not implemented");
			//throw Math_Illegal_State_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::CONVERGENCE_FAILED);
		}
	}
};