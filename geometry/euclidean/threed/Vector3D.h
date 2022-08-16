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

  //import java.io.Serializable;
  //import java.text.Number_Format;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.geometry.Localized_Geometry_Formats;
  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.Space;
  //import org.hipparchus.geometry.Vector;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Sin_Cos;

  /**
   * This class : vectors in a three-dimensional space.
   * <p>Instance of this class are guaranteed to be immutable.</p>
   */
class Vector_3D, Vector<Euclidean_3D>
{
	/** Null vector (coordinates: 0, 0, 0). */
	public static const Vector_3D ZERO = Vector_3D(0, 0, 0);

	/** First canonical vector (coordinates: 1, 0, 0). */
	public static const Vector_3D PLUS_I = Vector_3D(1, 0, 0);

	/** Opposite of the first canonical vector (coordinates: -1, 0, 0). */
	public static const Vector_3D MINUS_I = Vector_3D(-1, 0, 0);

	/** Second canonical vector (coordinates: 0, 1, 0). */
	public static const Vector_3D PLUS_J = Vector_3D(0, 1, 0);

	/** Opposite of the second canonical vector (coordinates: 0, -1, 0). */
	public static const Vector_3D MINUS_J = Vector_3D(0, -1, 0);

	/** Third canonical vector (coordinates: 0, 0, 1). */
	public static const Vector_3D PLUS_K = Vector_3D(0, 0, 1);

	/** Opposite of the third canonical vector (coordinates: 0, 0, -1).  */
	public static const Vector_3D MINUS_K = Vector_3D(0, 0, -1);

	// CHECKSTYLE: stop Constant_Name
	/** A vector with all coordinates set to NaN. */
	public static const Vector_3D NaN = Vector_3D(Double.NaN, NAN, NAN);
	// CHECKSTYLE: resume Constant_Name

	/** A vector with all coordinates set to positive infinity. */
	public static const Vector_3D POSITIVE_INFINITY =
		Vector_3D(INFINITY, INFINITY, INFINITY);

	/** A vector with all coordinates set to negative infinity. */
	public static const Vector_3D NEGATIVE_INFINITY =
		Vector_3D(-INFINITY, -INFINITY, -INFINITY);

	1313493323784566947L;

	/** Abscissa. */
	private const double x;

	/** Ordinate. */
	private const double y;

	/** Height. */
	private const double z;

	/** Simple constructor.
	 * Build a vector from its coordinates
	 * @param x abscissa
	 * @param y ordinate
	 * @param z height
	 * @see #get_x()
	 * @see #get_y()
	 * @see #get_z()
	 */
	public Vector_3D(const double& x, const double& y, const double& z)
	{
		this.x = x;
		this.y = y;
		this.z = z;
	}

	/** Simple constructor.
	 * Build a vector from its coordinates
	 * @param v coordinates array
	 * @exception  if array does not have 3 elements
	 * @see #to_array()
	 */
	public Vector_3D(std::vector<double> v)
	{
		if (v.size() != 3)
		{
			throw std::exception("not implemented");
			// throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), 3);
		}
		this.x = v[0];
		this.y = v[1];
		this.z = v[2];
	}

	/** Simple constructor.
	 * Build a vector from its azimuthal coordinates
	 * @param alpha azimuth (&alpha;) around Z
	 *              (0 is +X, &pi;/2 is +Y, &pi; is -X and 3&pi;/2 is -Y)
	 * @param delta elevation (&delta;) above (XY) plane, from -&pi;/2 to +&pi;/2
	 * @see #get_alpha()
	 * @see #get_delta()
	 */
	public Vector_3D(double alpha, double delta)
	{
		Sin_Cos sin_cos_alpha = Sin_Cos(alpha);
		Sin_Cos sin_cos_delta = Sin_Cos(delta);
		this.x = sin_cos_alpha.cos() * sin_cos_delta.cos();
		this.y = sin_cos_alpha.sin() * sin_cos_delta.cos();
		this.z = sin_cos_delta.sin();
	}

	/** Multiplicative constructor
	 * Build a vector from another one and a scale factor.
	 * The vector built will be a * u
	 * @param a scale factor
	 * @param u base (unscaled) vector
	 */
	public Vector_3D(const double& a, Vector_3D u)
	{
		this.x = a * u.x;
		this.y = a * u.y;
		this.z = a * u.z;
	}

	/** Linear constructor
	 * Build a vector from two other ones and corresponding scale factors.
	 * The vector built will be a1 * u1 + a2 * u2
	 * @param a1 first scale factor
	 * @param u1 first base (unscaled) vector
	 * @param a2 second scale factor
	 * @param u2 second base (unscaled) vector
	 */
	public Vector_3D(double a1, Vector_3D u1, double a2, Vector_3D u2)
	{
		this.x = Math_Arrays::linear_combination(a1, u1.x, a2, u2.x);
		this.y = Math_Arrays::linear_combination(a1, u1.y, a2, u2.y);
		this.z = Math_Arrays::linear_combination(a1, u1.z, a2, u2.z);
	}

	/** Linear constructor
	 * Build a vector from three other ones and corresponding scale factors.
	 * The vector built will be a1 * u1 + a2 * u2 + a3 * u3
	 * @param a1 first scale factor
	 * @param u1 first base (unscaled) vector
	 * @param a2 second scale factor
	 * @param u2 second base (unscaled) vector
	 * @param a3 third scale factor
	 * @param u3 third base (unscaled) vector
	 */
	public Vector_3D(double a1, Vector_3D u1, double a2, Vector_3D u2, double a3, Vector_3D u3)
	{
		this.x = Math_Arrays::linear_combination(a1, u1.x, a2, u2.x, a3, u3.x);
		this.y = Math_Arrays::linear_combination(a1, u1.y, a2, u2.y, a3, u3.y);
		this.z = Math_Arrays::linear_combination(a1, u1.z, a2, u2.z, a3, u3.z);
	}

	/** Linear constructor
	 * Build a vector from four other ones and corresponding scale factors.
	 * The vector built will be a1 * u1 + a2 * u2 + a3 * u3 + a4 * u4
	 * @param a1 first scale factor
	 * @param u1 first base (unscaled) vector
	 * @param a2 second scale factor
	 * @param u2 second base (unscaled) vector
	 * @param a3 third scale factor
	 * @param u3 third base (unscaled) vector
	 * @param a4 fourth scale factor
	 * @param u4 fourth base (unscaled) vector
	 */
	public Vector_3D(double a1, Vector_3D u1, double a2, Vector_3D u2, double a3, Vector_3D u3, double a4, Vector_3D u4)
	{
		this.x = Math_Arrays::linear_combination(a1, u1.x, a2, u2.x, a3, u3.x, a4, u4.x);
		this.y = Math_Arrays::linear_combination(a1, u1.y, a2, u2.y, a3, u3.y, a4, u4.y);
		this.z = Math_Arrays::linear_combination(a1, u1.z, a2, u2.z, a3, u3.z, a4, u4.z);
	}

	/** Get the abscissa of the vector.
	 * @return abscissa of the vector
	 * @see #Vector_3D(double, double, double)
	 */
	public double get_x()
	{
		return x;
	}

	/** Get the ordinate of the vector.
	 * @return ordinate of the vector
	 * @see #Vector_3D(double, double, double)
	 */
	public double get_y()
	{
		return y;
	}

	/** Get the height of the vector.
	 * @return height of the vector
	 * @see #Vector_3D(double, double, double)
	 */
	public double get_z()
	{
		return z;
	}

	/** Get the vector coordinates as a dimension 3 array.
	 * @return vector coordinates
	 * @see #Vector_3D(std::vector<double>)
	 */
	public std::vector<double> to_array()
	{
		return std::vector<double> { x, y, z };
	}

	/** {@inherit_doc} */
	//override
	public Space get_space()
	{
		return Euclidean_3D.get_instance();
	}

	/** {@inherit_doc} */
	//override
	public Vector_3D get_zero()
	{
		return ZERO;
	}

	/** {@inherit_doc} */
	//override
	public double get_norm1()
	{
		return std::abs(x) + std::abs(y) + std::abs(z);
	}

	/** {@inherit_doc} */
	//override
	public double get_norm()
	{
		// there are no cancellation problems here, so we use the straightforward formula
		return std::sqrt(x * x + y * y + z * z);
	}

	/** {@inherit_doc} */
	//override
	public double get_norm_sq()
	{
		// there are no cancellation problems here, so we use the straightforward formula
		return x * x + y * y + z * z;
	}

	/** {@inherit_doc} */
	//override
	public double get_norm_inf()
	{
		return std::max(std::max(std::abs(x), std::abs(y)), std::abs(z));
	}

	/** Get the azimuth of the vector.
	 * @return azimuth (&alpha;) of the vector, between -&pi; and +&pi;
	 * @see #Vector_3D(double, double)
	 */
	public double get_alpha()
	{
		return std::atan2(y, x);
	}

	/** Get the elevation of the vector.
	 * @return elevation (&delta;) of the vector, between -&pi;/2 and +&pi;/2
	 * @see #Vector_3D(double, double)
	 */
	public double get_delta()
	{
		return std::asin(z / get_norm());
	}

	/** {@inherit_doc} */
	//override
	public Vector_3D add(const Vector<Euclidean_3D> v)
	{
		const Vector_3D v3 = (Vector_3D)v;
		return Vector_3D(x + v3.x, y + v3.y, z + v3.z);
	}

	/** {@inherit_doc} */
	//override
	public Vector_3D add(double factor, const Vector<Euclidean_3D> v)
	{
		return Vector_3D(1, this, factor, (Vector_3D)v);
	}

	/** {@inherit_doc} */
	//override
	public Vector_3D subtract(const Vector<Euclidean_3D> v)
	{
		const Vector_3D v3 = (Vector_3D)v;
		return Vector_3D(x - v3.x, y - v3.y, z - v3.z);
	}

	/** {@inherit_doc} */
	//override
	public Vector_3D subtract(const double factor, const Vector<Euclidean_3D> v)
	{
		return Vector_3D(1, this, -factor, (Vector_3D)v);
	}

	/** {@inherit_doc} */
	//override
	public Vector_3D normalize() Math_Runtime_Exception
	{
		double s = get_norm();
		if (s == 0)
		{
			throw Math_Runtime_Exception(Localized_Geometry_Formats.CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR);
		}
		return scalar_multiply(1 / s);
	}

	/** Get a vector orthogonal to the instance.
	 * <p>There are an infinite number of normalized vectors orthogonal
	 * to the instance. This method picks up one of them almost
	 * arbitrarily. It is useful when one needs to compute a reference
	 * frame with one of the axes in a predefined direction. The
	 * following example shows how to build a frame having the k axis
	 * aligned with the known vector u :
	 * <pre><code>
	 *   Vector_3D k = u.normalize();
	 *   Vector_3D i = k.orthogonal();
	 *   Vector_3D j = Vector_3D.cross_product(k, i);
	 * </code></pre></p>
	 * @return a normalized vector orthogonal to the instance
	 * @exception Math_Runtime_Exception if the norm of the instance is NULL
	 */
	public Vector_3D orthogonal() Math_Runtime_Exception
	{
		double threshold = 0.6 * get_norm();
		if (threshold == 0)
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
		}

		if (std::abs(x) <= threshold)
		{
			double inverse = 1 / std::sqrt(y * y + z * z);
			return Vector_3D(0, inverse * z, -inverse * y);
		}
		else if (std::abs(y) <= threshold)
		{
			double inverse = 1 / std::sqrt(x * x + z * z);
			return Vector_3D(-inverse * z, 0, inverse * x);
		}
		double inverse = 1 / std::sqrt(x * x + y * y);
		return Vector_3D(inverse * y, -inverse * x, 0);
	}

	/** Compute the angular separation between two vectors.
	 * <p>This method computes the angular separation between two
	 * vectors using the dot product for well separated vectors and the
	 * cross product for almost aligned vectors. This allows to have a
	 * good accuracy in all cases, even for vectors very close to each
	 * other.</p>
	 * @param v1 first vector
	 * @param v2 second vector
	 * @return angular separation between v1 and v2
	 * @exception Math_Runtime_Exception if either vector has a NULL norm
	 */
	public static double angle(Vector_3D v1, Vector_3D v2) Math_Runtime_Exception
	{
		double norm_product = v1.get_norm() * v2.get_norm();
		if (norm_product == 0)
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
		}

		double dot = v1.dot_product(v2);
		double threshold = norm_product * 0.9999;
		if ((dot < -threshold) || (dot > threshold))
		{
			// the vectors are almost aligned, compute using the sine
			Vector_3D v3 = cross_product(v1, v2);
			if (dot >= 0)
			{
				return std::asin(v3.get_norm() / norm_product);
			}
			return std::numbers::pi - std::asin(v3.get_norm() / norm_product);
		}

		// the vectors are sufficiently separated to use the cosine
		return std::acos(dot / norm_product);
	}

	/** {@inherit_doc} */
	//override
	public Vector_3D negate()
	{
		return Vector_3D(-x, -y, -z);
	}

	/** {@inherit_doc} */
	//override
	public Vector_3D scalar_multiply(double a)
	{
		return Vector_3D(a * x, a * y, a * z);
	}

	/** {@inherit_doc} */
	//override
	public bool is_nan()
	{
		return std::isnan(x) || std::isnan(y) || std::isnan(z);
	}

	/** {@inherit_doc} */
	//override
	public bool is_infinite()
	{
		return !is_nan() && (std::isinf(x) || std::isinf(y) || std::isinf(z));
	}

	/**
	 * Test for the equality of two 3D vectors.
	 * <p>
	 * If all coordinates of two 3D vectors are exactly the same, and none are
	 * {@codeNAN}, the two 3D vectors are considered to be equal.
	 * </p>
	 * <p>
	 * {@code NaN} coordinates are considered to affect globally the vector
	 * and be equals to each other - i.e, if either (or all) coordinates of the
	 * 3D vector are equal to {@codeNAN}, the 3D vector is equal to
	 * {@link #NaN}.
	 * </p>
	 *
	 * @param other Object to test for equality to this
	 * @return true if two 3D vector objects are equal, false if
	 *         object is NULL, not an instance of Vector_3D, or
	 *         not equal to this Vector_3D instance
	 *
	 */
	 //override
	public bool equals(Object other)
	{
		if (this == other)
		{
			return true;
		}
		if (dynamic_cast<const Vector_3D*>(*other) != nullptr)
		{
			const Vector_3D rhs = (Vector_3D)other;
			return x == rhs.x && y == rhs.y && z == rhs.z || is_nan() && rhs.is_nan();
		}

		return false;
	}

	/**
	 * Test for the equality of two 3D vectors.
	 * <p>
	 * If all coordinates of two 3D vectors are exactly the same, and none are
	 * {@code NaN}, the two 3D vectors are considered to be equal.
	 * </p>
	 * <p>
	 * In compliance with IEEE754 handling, if any coordinates of any of the
	 * two vectors are {@code NaN}, then the vectors are considered different.
	 * This implies that {@link #NaN Vector_3D.NaN}.equals({@link #NaN Vector_3D.NaN})
	 * returns {@code false} despite the instance is checked against itself.
	 * </p>
	 *
	 * @param other Object to test for equality to this
	 * @return true if two 3D vector objects are equal, false if
	 *         object is NULL, not an instance of Vector_3D, or
	 *         not equal to this Vector_3D instance
	 * @since 2.1
	 */
	public bool equals_ieee_754(Object other)
	{
		if (this == other && !is_nan())
		{
			return true;
		}

		if (dynamic_cast<const Vector_3D*>(*other) != nullptr)
		{
			const Vector_3D rhs = (Vector_3D)other;
			return x == rhs.x && y == rhs.y && z == rhs.z;
		}

		return false;
	}

	/**
	 * Get a hash_code for the 3D vector.
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
			return 642;
		}
		return 643 * (164 * Math_Utils::hash(x) + 3 * Math_Utils::hash(y) + Math_Utils::hash(z));
	}

	/** {@inherit_doc}
	 * <p>
	 * The implementation uses specific multiplication and addition
	 * algorithms to preserve accuracy and reduce cancellation effects.
	 * It should be very accurate even for nearly orthogonal vectors.
	 * </p>
	 * @see Math_Arrays#linear_combination(double, double, double, double, double, double)
	 */
	 //override
	public double dot_product(const Vector<Euclidean_3D> v)
	{
		const Vector_3D v3 = (Vector_3D)v;
		return Math_Arrays::linear_combination(x, v3.x, y, v3.y, z, v3.z);
	}

	/** Compute the cross-product of the instance with another vector.
	 * @param v other vector
	 * @return the cross product this ^ v as a Vector_3D
	 */
	public Vector_3D cross_product(const Vector<Euclidean_3D> v)
	{
		const Vector_3D v3 = (Vector_3D)v;
		return Vector_3D(Math_Arrays::linear_combination(y, v3.z, -z, v3.y), Math_Arrays::linear_combination(z, v3.x, -x, v3.z), Math_Arrays::linear_combination(x, v3.y, -y, v3.x));
	}

	/** {@inherit_doc} */
	//override
	public double distance1(Vector<Euclidean_3D> v)
	{
		const Vector_3D v3 = (Vector_3D)v;
		const double dx = std::abs(v3.x - x);
		const double dy = std::abs(v3.y - y);
		const double dz = std::abs(v3.z - z);
		return dx + dy + dz;
	}

	/** {@inherit_doc} */
	//override
	public double distance(Point<Euclidean_3D> v)
	{
		const Vector_3D v3 = (Vector_3D)v;
		const double dx = v3.x - x;
		const double dy = v3.y - y;
		const double dz = v3.z - z;
		return std::sqrt(dx * dx + dy * dy + dz * dz);
	}

	/** {@inherit_doc} */
	//override
	public double distance_inf(Vector<Euclidean_3D> v)
	{
		const Vector_3D v3 = (Vector_3D)v;
		const double dx = std::abs(v3.x - x);
		const double dy = std::abs(v3.y - y);
		const double dz = std::abs(v3.z - z);
		return std::max(std::max(dx, dy), dz);
	}

	/** {@inherit_doc} */
	//override
	public double distance_sq(Vector<Euclidean_3D> v)
	{
		const Vector_3D v3 = (Vector_3D)v;
		const double dx = v3.x - x;
		const double dy = v3.y - y;
		const double dz = v3.z - z;
		return dx * dx + dy * dy + dz * dz;
	}

	/** Compute the dot-product of two vectors.
	 * @param v1 first vector
	 * @param v2 second vector
	 * @return the dot product v1.v2
	 */
	public static double dot_product(Vector_3D v1, Vector_3D v2)
	{
		return v1.dot_product(v2);
	}

	/** Compute the cross-product of two vectors.
	 * @param v1 first vector
	 * @param v2 second vector
	 * @return the cross product v1 ^ v2 as a Vector
	 */
	public static Vector_3D cross_product(const Vector_3D v1, const Vector_3D v2)
	{
		return v1.cross_product(v2);
	}

	/** Compute the distance between two vectors according to the L<sub>1</sub> norm.
	 * <p>Calling this method is equivalent to calling:
	 * <code>v1.subtract(v2).get_norm1()</code> except that no intermediate
	 * vector is built</p>
	 * @param v1 first vector
	 * @param v2 second vector
	 * @return the distance between v1 and v2 according to the L<sub>1</sub> norm
	 */
	public static double distance1(Vector_3D v1, Vector_3D v2)
	{
		return v1.distance1(v2);
	}

	/** Compute the distance between two vectors according to the L<sub>2</sub> norm.
	 * <p>Calling this method is equivalent to calling:
	 * <code>v1.subtract(v2).get_norm()</code> except that no intermediate
	 * vector is built</p>
	 * @param v1 first vector
	 * @param v2 second vector
	 * @return the distance between v1 and v2 according to the L<sub>2</sub> norm
	 */
	public static double distance(Vector_3D v1, Vector_3D v2)
	{
		return v1.distance(v2);
	}

	/** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
	 * <p>Calling this method is equivalent to calling:
	 * <code>v1.subtract(v2).get_norm_inf()</code> except that no intermediate
	 * vector is built</p>
	 * @param v1 first vector
	 * @param v2 second vector
	 * @return the distance between v1 and v2 according to the L<sub>&infin;</sub> norm
	 */
	public static double distance_inf(Vector_3D v1, Vector_3D v2)
	{
		return v1.distance_inf(v2);
	}

	/** Compute the square of the distance between two vectors.
	 * <p>Calling this method is equivalent to calling:
	 * <code>v1.subtract(v2).get_norm_sq()</code> except that no intermediate
	 * vector is built</p>
	 * @param v1 first vector
	 * @param v2 second vector
	 * @return the square of the distance between v1 and v2
	 */
	public static double distance_sq(Vector_3D v1, Vector_3D v2)
	{
		return v1.distance_sq(v2);
	}

	/** Get a string representation of this vector.
	 * @return a string representation of this vector
	 */
	 //override
	public std::string to_string() const
	{
		return Vector_3DFormat.get_vector_3d_format().format(this);
	}

	/** {@inherit_doc} */
	//override
	public std::string to_string(const Number_Format format)
	{
		return Vector_3DFormat(format).format(this);
	}
}
