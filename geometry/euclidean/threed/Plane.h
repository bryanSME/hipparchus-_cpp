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

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.Vector;
  //import org.hipparchus.geometry.euclidean.oned.Euclidean_1D;
  //import org.hipparchus.geometry.euclidean.oned.Vector_1D;
  //import org.hipparchus.geometry.euclidean.twod.Euclidean_2D;
  //import org.hipparchus.geometry.euclidean.twod.Polygons_Set;
  //import org.hipparchus.geometry.euclidean.twod.Vector_2D;
  //import org.hipparchus.geometry.partitioning.Embedding;
  //import org.hipparchus.geometry.partitioning.Hyperplane;
  //import org.hipparchus.geometry.partitioning.Region_Factory;
  //import org.hipparchus.util.FastMath;
#include "../../partitioning/Hyperplane.h"
#include "../../partitioning/Embedding.h"
#include "Euclidean3D.h"
#include "../twod/Euclidean2D.h"
#include "Vector3D.h"
#include "../twod/Vector2D.h"
#include <exception>
#include "../../Point.hpp"

  /** The class represent planes in a three dimensional space.
   */
class Plane : Hyperplane<Euclidean_3D>, Embedding<Euclidean_3D, Euclidean_2D>
{
private:
	/** Offset of the origin with respect to the plane. */
	double my_origin_offset;

	/** Origin of the plane frame. */
	Vector_3D my_origin;

	/** First vector of the plane frame (in plane). */
	Vector_3D my_u;

	/** Second vector of the plane frame (in plane). */
	Vector_3D my_v;

	/** Third vector of the plane frame (plane normal). */
	Vector_3D my_w;

	/** Tolerance below which points are considered identical. */
	double my_tolerance;

	/** Set the normal vactor.
	 * @param normal normal direction to the plane (will be copied)
	 * @exception Math_Runtime_Exception if the normal norm is too small
	 */
	void set_normal(const Vector_3D& normal)
	{
		const double norm = normal.get_norm();
		if (norm < 1.0e-10)
		{
			throw std::exception("not implemented");
			//throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
		}
		my_w = Vector_3D(1.0 / norm, normal);
	}

	/** Reset the plane frame.
	 */
	void set_frame()
	{
		my_origin = Vector_3D(-my_origin_offset, my_w);
		my_u = my_w.orthogonal();
		my_v = Vector_3D.cross_product(my_w, my_u);
	}

public:
	/** Build a plane normal to a given direction and containing the origin.
	 * @param normal normal direction to the plane
	 * @param tolerance tolerance below which points are considered identical
	 * @exception Math_Runtime_Exception if the normal norm is too small
	 */
	Plane(const Vector_3D& normal, const double& tolerance)
	{
		set_normal(normal);
		my_tolerance = tolerance;
		my_origin_offset = 0;
		set_frame();
	}

	/** Build a plane from a point and a normal.
	 * @param p point belonging to the plane
	 * @param normal normal direction to the plane
	 * @param tolerance tolerance below which points are considered identical
	 * @exception Math_Runtime_Exception if the normal norm is too small
	 */
	Plane(const Vector_3D& p, const Vector_3D& normal, const double& tolerance)
	{
		set_normal(normal);
		my_tolerance = tolerance;
		my_origin_offset = -p.dot_product(w);
		set_frame();
	}

	/** Build a plane from three points.
	 * <p>The plane is oriented in the direction of
	 * {@code (p2-p1) ^ (p3-p1)}</p>
	 * @param p1 first point belonging to the plane
	 * @param p2 second point belonging to the plane
	 * @param p3 third point belonging to the plane
	 * @param tolerance tolerance below which points are considered identical
	 * @exception Math_Runtime_Exception if the points do not constitute a plane
	 */
	Plane(const Vector_3D& p1, const Vector_3D& p2, const Vector_3D& p3, const double& tolerance)
	{
		Plane(p1, p2.subtract(p1).cross_product(p3.subtract(p1)), tolerance);
	}

	/** Copy constructor.
	 * <p>The instance created is completely independent of the original
	 * one. A deep copy is used, none of the underlying object are
	 * shared.</p>
	 * @param plane plane to copy
	 */
	Plane(const Plane& plane)
	{
		my_origin_offset = plane.origin_offset;
		my_origin = plane.origin;
		my_u = plane.u;
		my_v = plane.v;
		my_w = plane.w;
		my_tolerance = plane.tolerance;
	}

	/** Copy the instance.
	 * <p>The instance created is completely independant of the original
	 * one. A deep copy is used, none of the underlying objects are
	 * shared (except for immutable objects).</p>
	 * @return a hyperplane, copy of the instance
	 */
	 //override
	Plane copy_self()
	{
		return Plane(*this);
	}

	/** Reset the instance as if built from a point and a normal.
	 * @param p point belonging to the plane
	 * @param normal normal direction to the plane
	 * @exception Math_Runtime_Exception if the normal norm is too small
	 */
	void reset(const Vector_3D& p, const Vector_3D& normal)
	{
		set_normal(normal);
		my_origin_offset = -p.dot_product(my_w);
		set_frame();
	}

	/** Reset the instance from another one.
	 * <p>The updated instance is completely independant of the original
	 * one. A deep reset is used none of the underlying object is
	 * shared.</p>
	 * @param original plane to reset from
	 */
	void reset(const Plane& original)
	{
		my_origin_offset = original.get_origin_offset();
		my_origin = original.get_origin();
		my_u = original.get_u();
		my_v = original.get_v();
		my_w = original.get_w();
	}

	/** Get the origin point of the plane frame.
	 * <p>The point returned is the orthogonal projection of the
	 * 3D-space origin in the plane.</p>
	 * @return the origin point of the plane frame (point closest to the
	 * 3D-space origin)
	 */
	Vector_3D get_origin() const
	{
		return my_origin;
	}

	/** Get the normalized normal vector.
	 * <p>The frame defined by ({@link #get_u get_u}, {@link #get_v get_v}, * {@link #get_normal get_normal}) is a rigth-handed orthonormalized
	 * frame).</p>
	 * @return normalized normal vector
	 * @see #get_u
	 * @see #get_v
	 */
	Vector_3D get_normal() const
	{
		return my_w;
	}

	/** Get the plane first canonical vector.
	 * <p>The frame defined by ({@link #get_u get_u}, {@link #get_v get_v}, * {@link #get_normal get_normal}) is a rigth-handed orthonormalized
	 * frame).</p>
	 * @return normalized first canonical vector
	 * @see #get_v
	 * @see #get_normal
	 */
	Vector_3D get_u() const
	{
		return my_u;
	}

	/** Get the plane second canonical vector.
	 * <p>The frame defined by ({@link #get_u get_u}, {@link #get_v get_v}, * {@link #get_normal get_normal}) is a rigth-handed orthonormalized
	 * frame).</p>
	 * @return normalized second canonical vector
	 * @see #get_u
	 * @see #get_normal
	 */
	Vector_3D get_v() const
	{
		return my_v;
	}

	/** {@inherit_doc}
	 */
	 //override
	Point<Euclidean_3D> project(Point<Euclidean_3D>& point)
	{
		return to_space(to_sub_space(point));
	}

	/** {@inherit_doc}
	 */
	 //override
	double get_tolerance() const
	{
		return my_tolerance;
	}

	/** Revert the plane.
	 * <p>Replace the instance by a similar plane with opposite orientation.</p>
	 * <p>The plane frame is chosen in such a way that a 3D point that had
	 * {@code (x, y)} in-plane coordinates and {@code z} offset with
	 * respect to the plane and is unaffected by the change will have
	 * {@code (y, x)} in-plane coordinates and {@code -z} offset with
	 * respect to the plane. This means that the {@code u} and {@code v}
	 * vectors returned by the {@link #get_u} and {@link #get_v} methods are exchanged, * and the {@code w} vector returned by the {@link #get_normal} method is
	 * reversed.</p>
	 */
	void revert_self()
	{
		const Vector_3D tmp = my_u;
		my_u = my_v;
		my_v = tmp;
		my_w = my_w.negate();
		my_origin_offset = -my_origin_offset;
	}

	/** Transform a space point into a sub-space point.
	 * @param vector n-dimension point of the space
	 * @return (n-1)-dimension point of the sub-space corresponding to
	 * the specified space point
	 */
	Vector_2D to_sub_space(Vector<Euclidean_3D>& vector)
	{
		return to_sub_space((Point<Euclidean_3D>) vector);
	}

	/** Transform a sub-space point into a space point.
	 * @param vector (n-1)-dimension point of the sub-space
	 * @return n-dimension point of the space corresponding to the
	 * specified sub-space point
	 */
	Vector_3D to_space(const Vector<Euclidean_2D>& vector)
	{
		return to_space((Point<Euclidean_2D>) vector);
	}

	/** Transform a 3D space point into an in-plane point.
	 * @param point point of the space (must be a {@link Vector_3D
	 * Vector_3D} instance)
	 * @return in-plane point (really a {@link
	 * org.hipparchus.geometry.euclidean.twod.Vector_2D Vector_2D} instance)
	 * @see #to_space
	 */
	 //override
	Vector_2D to_sub_space(const Point<Euclidean_3D>& point)
	{
		const Vector_3D p3D = (Vector_3D)point;
		return Vector_2D(p3D.dot_product(u), p3D.dot_product(my_v));
	}

	/** Transform an in-plane point into a 3D space point.
	 * @param point in-plane point (must be a {@link
	 * org.hipparchus.geometry.euclidean.twod.Vector_2D Vector_2D} instance)
	 * @return 3D space point (really a {@link Vector_3D Vector_3D} instance)
	 * @see #to_sub_space
	 */
	 //override
	Vector_3D to_space(const Point<Euclidean_2D>& point)
	{
		const Vector_2D p2_d = (Vector_2D)point;
		return Vector_3D(p2_d.get_x(), my_u, p2_d.get_y(), my_v, -origin_offset, my_w);
	}

	/** Get one point from the 3D-space.
	 * @param in_plane desired in-plane coordinates for the point in the
	 * plane
	 * @param offset desired offset for the point
	 * @return one point in the 3D-space, with given coordinates and offset
	 * relative to the plane
	 */
	Vector_3D get_point_at(const Vector_2D in_plane, const double& offset)
	{
		return Vector_3D(in_plane.get_x(), my_u, in_plane.get_y(), my_v, offset - my_origin_offset, my_w);
	}

	/** Check if the instance is similar to another plane.
	 * <p>Planes are considered similar if they contain the same
	 * points. This does not mean they are equal since they can have
	 * opposite normals.</p>
	 * @param plane plane to which the instance is compared
	 * @return true if the planes are similar
	 */
	bool is_similar_to(const Plane plane)
	{
		const double& angle = Vector_3D.angle(my_w, plane.get_w());
		return ((angle < 1.0e-10) && (std::abs(origin_offset - plane.get_origin_offset()) < my_tolerance)) ||
			((angle > (std::numbers::pi - 1.0e-10)) && (std::abs(origin_offset + plane.get_origin_offset()) < tolerance));
	}

	/** Rotate the plane around the specified point.
	 * <p>The instance is not modified, a instance is created.</p>
	 * @param center rotation center
	 * @param rotation vectorial rotation operator
	 * @return a plane
	 */
	Plane rotate(const Vector_3D& center, const Rotation& rotation)
	{
		const Vector_3D delta = origin.subtract(center);
		const Plane plane = Plane(center.add(rotation.apply_to(delta)), rotation.apply_to(my_w), my_tolerance);

		// make sure the frame is transformed as desired
		plane.get_u() = rotation.apply_to(my_u);
		plane.get_v() = rotation.apply_to(my_v);

		return plane;
	}

	/** Translate the plane by the specified amount.
	 * <p>The instance is not modified, a instance is created.</p>
	 * @param translation translation to apply
	 * @return a plane
	 */
	Plane translate(const Vector_3D& translation)
	{
		const Plane plane = Plane(my_origin.add(translation), my_w, my_tolerance);

		// make sure the frame is transformed as desired
		plane.get_u() = my_u;
		plane.get_v() = my_v;

		return plane;
	}

	/** Get the intersection of a line with the instance.
	 * @param line line intersecting the instance
	 * @return intersection point between between the line and the
	 * instance (null if the line is parallel to the instance)
	 */
	Vector_3D intersection(const Line& line)
	{
		const Vector_3D& direction = line.get_direction();
		const double   dot = my_w.dot_product(direction);
		if (std::abs(dot) < 1.0e-10)
		{
			return NULL;
		}
		const Vector_3D point = line.to_space((Point<Euclidean_1D>) Vector_1D.ZERO);
		const double   k = -(origin_offset + my_w.dot_product(point)) / dot;
		return Vector_3D(1.0, point, k, direction);
	}

	/** Build the line shared by the instance and another plane.
	 * @param other other plane
	 * @return line at the intersection of the instance and the
	 * other plane (really a {@link Line Line} instance)
	 */
	Line intersection(const Plane other)
	{
		const Vector_3D direction = Vector_3D.cross_product(my_w, other.get_w());
		if (direction.get_norm() < my_tolerance)
		{
			return NULL;
		}
		const Vector_3D point = intersection(*this, other, Plane(direction, my_tolerance));
		return Line(point, point.add(direction), my_tolerance);
	}

	/** Get the intersection point of three planes.
	 * @param plane1 first plane1
	 * @param plane2 second plane2
	 * @param plane3 third plane2
	 * @return intersection point of three planes, NULL if some planes are parallel
	 */
	static Vector_3D intersection(const Plane& plane1, const Plane& plane2, const Plane& plane3)
	{
		// coefficients of the three planes linear equations
		const double a1 = plane1.get_w().get_x();
		const double b1 = plane1.get_w().get_y();
		const double c1 = plane1.get_w().get_z();
		const double d1 = plane1.get_origin_offset();

		const double a2 = plane2.get_w().get_x();
		const double b2 = plane2.get_w().get_y();
		const double c2 = plane2.get_w().get_z();
		const double d2 = plane2.get_origin_offset();

		const double a3 = plane3.get_w().get_x();
		const double b3 = plane3.get_w().get_y();
		const double c3 = plane3.get_w().get_z();
		const double d3 = plane3.get_origin_offset();

		// direct Cramer resolution of the linear system
		// (this is still feasible for a 3x3 system)
		const double a23 = b2 * c3 - b3 * c2;
		const double b23 = c2 * a3 - c3 * a2;
		const double c23 = a2 * b3 - a3 * b2;
		const double determinant = a1 * a23 + b1 * b23 + c1 * c23;
		if (std::abs(determinant) < 1.0e-10)
		{
			return NULL;
		}

		const double r = 1.0 / determinant;
		return Vector_3D(
			(-a23 * d1 - (c1 * b3 - c3 * b1) * d2 - (c2 * b1 - c1 * b2) * d3) * r, (-b23 * d1 - (c3 * a1 - c1 * a3) * d2 - (c1 * a2 - c2 * a1) * d3) * r, (-c23 * d1 - (b1 * a3 - b3 * a1) * d2 - (b2 * a1 - b1 * a2) * d3) * r);
	}

	/** Build a region covering the whole hyperplane.
	 * @return a region covering the whole hyperplane
	 */
	 //override
	Sub_Plane whole_hyperplane()
	{
		return Sub_Plane(*this, Polygons_Set(my_tolerance));
	}

	/** {@inherit_doc} */
	//override
	Sub_Plane empty_hyperplane()
	{
		return Sub_Plane(*this, Region_Factory<Euclidean_2D>().get_complement(Polygons_Set(my_tolerance)));
	}

	/** Build a region covering the whole space.
	 * @return a region containing the instance (really a {@link
	 * Polyhedrons_Set Polyhedrons_Set} instance)
	 */
	 //override
	Polyhedrons_Set whole_space()
	{
		return Polyhedrons_Set(my_tolerance);
	}

	/** Check if the instance contains a point.
	 * @param p point to check
	 * @return true if p belongs to the plane
	 */
	bool contains(const Vector_3D& p)
	{
		return std::abs(get_offset(p)) < my_tolerance;
	}

	/** Get the offset (oriented distance) of a parallel plane.
	 * <p>This method should be called only for parallel planes otherwise
	 * the result is not meaningful.</p>
	 * <p>The offset is 0 if both planes are the same, it is
	 * positive if the plane is on the plus side of the instance and
	 * negative if it is on the minus side, according to its natural
	 * orientation.</p>
	 * @param plane plane to check
	 * @return offset of the plane
	 */
	double get_offset(const Plane& plane)
	{
		return my_origin_offset + (same_orientation_as(plane) ? -plane.get_origin_offset() : plane.get_origin_offset());
	}

	/** Get the offset (oriented distance) of a vector.
	 * @param vector vector to check
	 * @return offset of the vector
	 */
	double get_offset(Vector<Euclidean_3D>& vector)
	{
		return get_offset((Point<Euclidean_3D>) vector);
	}

	/** Get the offset (oriented distance) of a point.
	 * <p>The offset is 0 if the point is on the underlying hyperplane, * it is positive if the point is on one particular side of the
	 * hyperplane, and it is negative if the point is on the other side, * according to the hyperplane natural orientation.</p>
	 * @param point point to check
	 * @return offset of the point
	 */
	 //override
	double get_offset(const Point<Euclidean_3D>& point)
	{
		return ((Vector_3D)point).dot_product(my_w) + my_origin_offset;
	}

	/** Check if the instance has the same orientation as another hyperplane.
	 * @param other other hyperplane to check against the instance
	 * @return true if the instance and the other hyperplane have
	 * the same orientation
	 */
	 //override
	bool same_orientation_as(const Hyperplane<Euclidean_3D>& other)
	{
		return (((Plane)other).get_w()).dot_product(my_w) > 0.0;
	}
};