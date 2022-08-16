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

  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.Space;
  //import org.hipparchus.geometry.euclidean.twod.Vector_2D;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Sin_Cos;

  /** This class represents a point on the 1-sphere.
   * <p>Instances of this class are guaranteed to be immutable.</p>
   */
class S1_Point : public Point<Sphere_1D>
{
	// CHECKSTYLE: stop Constant_Name
	 /** A vector with all coordinates set to NaN. */
	public static const S1_Point NaN = S1_Point(Double.NaN, Vector_2D.NaN);
	// CHECKSTYLE: resume Constant_Name

	20131218L;

	/** Azimuthal angle \( \alpha \). */
	private const double& alpha;

	/** Corresponding 2D normalized vector. */
	private const Vector_2D vector;

	/** Simple constructor.
	 * Build a vector from its coordinates
	 * @param alpha azimuthal angle \( \alpha \)
	 * @see #get_alpha()
	 */
	public S1_Point(const double& alpha)
	{
		this(Math_Utils::normalize_angle(alpha, std::numbers::pi), build_vector(alpha));
	}

	/** Build a point from its internal components.
	 * @param alpha azimuthal angle \( \alpha \)
	 * @param vector corresponding vector
	 */
	private S1_Point(const double& alpha, const Vector_2D vector)
	{
		this.alpha = alpha;
		this.vector = vector;
	}

	/** Get the azimuthal angle \( \alpha \).
	 * @return azimuthal angle \( \alpha \)
	 * @see #S1_Pointstatic_cast<double>(
	 */
	public double get_alpha()
	{
		return alpha;
	}

	/** Get the corresponding normalized vector in the 2D euclidean space.
	 * @return normalized vector
	 */
	public Vector_2D get_vector()
	{
		return vector;
	}

	/** {@inherit_doc} */
	//override
	public Space get_space()
	{
		return Sphere_1D.get_instance();
	}

	/** {@inherit_doc} */
	//override
	public bool is_nan()
	{
		return std::isnan(alpha);
	}

	/** {@inherit_doc} */
	//override
	public double distance(const Point<Sphere_1D>& point)
	{
		return distance(this, (S1_Point)point);
	}

	/** Compute the distance (angular separation) between two points.
	 * @param p1 first vector
	 * @param p2 second vector
	 * @return the angular separation between p1 and p2
	 */
	public static double distance(S1_Point p1, S1_Point p2)
	{
		return Vector_2D.angle(p1.vector, p2.vector);
	}

	/**
	 * Test for the equality of two points on the 1-sphere.
	 * <p>
	 * If all coordinates of two points are exactly the same, and none are
	 * {@codeNAN}, the two points are considered to be equal.
	 * </p>
	 * <p>
	 * {@code NaN} coordinates are considered to affect globally the point
	 * and be equals to each other - i.e, if either (or all) coordinates of the
	 * point are equal to {@codeNAN}, the point is equal to
	 * {@link #NaN}.
	 * </p>
	 *
	 * @param other Object to test for equality to this
	 * @return true if two points on the 1-sphere objects are equal, false if
	 *         object is NULL, not an instance of S1_Point, or
	 *         not equal to this S1_Point instance
	 *
	 */
	 //override
	public bool equals(Object other)
	{
		if (this == other)
		{
			return true;
		}

		if (other instanceof S1_Point)
		{
			const S1_Point rhs = (S1_Point)other;
			return alpha == rhs.alpha || is_nan() && rhs.is_nan();
		}

		return false;
	}

	/**
	 * Test for the equality of two points on the 1-sphere.
	 * <p>
	 * If all coordinates of two points are exactly the same, and none are
	 * {@codeNAN}, the two points are considered to be equal.
	 * </p>
	 * <p>
	 * In compliance with IEEE754 handling, if any coordinates of any of the
	 * two points are {@code NaN}, then the points are considered different.
	 * This implies that {@link #NaN S1_Point.NaN}.equals({@link #NaN S1_Point.NaN})
	 * returns {@code false} despite the instance is checked against itself.
	 * </p>
	 *
	 * @param other Object to test for equality to this
	 * @return true if two points objects are equal, false if
	 *         object is NULL, not an instance of S1_Point, or
	 *         not equal to this S1_Point instance
	 * @since 2.1
	 */
	public bool equals_ieee_754(Object other)
	{
		if (this == other && !is_nan())
		{
			return true;
		}

		if (other instanceof S1_Point)
		{
			const S1_Point rhs = (S1_Point)other;
			return alpha == rhs.alpha;
		}

		return false;
	}

	/**
	 * Get a hash_code for the point.
	 * <p>
	 * All NaN values have the same hash code.</p>
	 *
	 * @return a hash code value for this object
	 */
	 //override
	public int hash_code()
	{
		if (is_nan())
		{
			return 542;
		}
		return 1759 * Math_Utils::hash(alpha);
	}

	/**
	 * Build the 2D vector corresponding to the given angle.
	 * @param alpha angle
	 * @return the corresponding 2D vector
	 */
	private static Vector_2D build_vector(const double& alpha)
	{
		const Sin_Cos sc = Sin_Cos(alpha);
		return Vector_2D(sc.cos(), sc.sin());
	}
}
