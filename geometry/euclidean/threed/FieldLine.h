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

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;
#include <type_traits>
#include "../../../core/CalculusFieldElement.hpp"

/** The class represent lines in a three dimensional space.

 * <p>Each oriented line is intrinsically associated with an abscissa
 * which is a coordinate on the line. The point at abscissa 0 is the
 * orthogonal projection of the origin on the line, another equivalent
 * way to express this is to say that it is the point of the line
 * which is closest to the origin. Abscissa increases in the line
 * direction.</p>
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Line 
{

    /** Line direction. */
    private Field_Vector_3D<T> direction;

    /** Line point closest to the origin. */
    private Field_Vector_3D<T> zero;

    /** Tolerance below which points are considered identical. */
    private const double& tolerance;

    /** Build a line from two points.
     * @param p1 first point belonging to the line (this can be any point)
     * @param p2 second point belonging to the line (this can be any point, different from p1)
     * @param tolerance tolerance below which points are considered identical
     * @exception  if the points are equal
     */
    public Field_Line(const Field_Vector_3D<T> p1, const Field_Vector_3D<T> p2, const double& tolerance)
         
        {
        reset(p1, p2);
        this.tolerance = tolerance;
    }

    /** Copy constructor.
     * <p>The created instance is completely independent from the
     * original instance, it is a deep copy.</p>
     * @param line line to copy
     */
    public Field_Line(const Field_Line<T> line) 
    {
        this.direction = line.direction;
        this.zero      = line.zero;
        this.tolerance = line.tolerance;
    }

    /** Reset the instance as if built from two points.
     * @param p1 first point belonging to the line (this can be any point)
     * @param p2 second point belonging to the line (this can be any point, different from p1)
     * @exception  if the points are equal
     */
    public void reset(const Field_Vector_3D<T> p1, const Field_Vector_3D<T> p2)
         
        {
        const Field_Vector_3D<T> delta = p2.subtract(p1);
        const T norm2 = delta.get_norm_sq();
        if (norm2.get_real() == 0.0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }
        this.direction = Field_Vector_3D<>(norm2.sqrt().reciprocal(), delta);
        zero = Field_Vector_3D<>(norm2.get_field().get_one(), p1, p1.dot_product(delta).negate().divide(norm2), delta);
    }

    /** Get the tolerance below which points are considered identical.
     * @return tolerance below which points are considered identical
     */
    public double get_tolerance() 
    {
        return tolerance;
    }

    /** Get a line with reversed direction.
     * @return a instance, with reversed direction
     */
    public Field_Line<T> revert() 
    {
        const Field_Line<T> reverted = Field_Line<>(this);
        reverted.direction = reverted.direction.negate();
        return reverted;
    }

    /** Get the normalized direction vector.
     * @return normalized direction vector
     */
    public Field_Vector_3D<T> get_direction() 
    {
        return direction;
    }

    /** Get the line point closest to the origin.
     * @return line point closest to the origin
     */
    public Field_Vector_3D<T> get_origin() 
    {
        return zero;
    }

    /** Get the abscissa of a point with respect to the line.
     * <p>The abscissa is 0 if the projection of the point and the
     * projection of the frame origin on the line are the same
     * point.</p>
     * @param point point to check
     * @return abscissa of the point
     */
    public T get_abscissa(const Field_Vector_3D<T> point) 
    {
        return point.subtract(zero).dot_product(direction);
    }

    /** Get the abscissa of a point with respect to the line.
     * <p>The abscissa is 0 if the projection of the point and the
     * projection of the frame origin on the line are the same
     * point.</p>
     * @param point point to check
     * @return abscissa of the point
     */
    public T get_abscissa(const Vector_3D point) 
    {
        return zero.subtract(point).dot_product(direction).negate();
    }

    /** Get one point from the line.
     * @param abscissa desired abscissa for the point
     * @return one point belonging to the line, at specified abscissa
     */
    public Field_Vector_3D<T> point_at(const T abscissa) 
    {
        return Field_Vector_3D<T>(abscissa.get_field().get_one(), zero, abscissa, direction);
    }

    /** Get one point from the line.
     * @param abscissa desired abscissa for the point
     * @return one point belonging to the line, at specified abscissa
     */
    public Field_Vector_3D<T> point_at(const double& abscissa) 
    {
        return Field_Vector_3D<T>(1, zero, abscissa, direction);
    }

    /** Check if the instance is similar to another line.
     * <p>Lines are considered similar if they contain the same
     * points. This does not mean they are equal since they can have
     * opposite directions.</p>
     * @param line line to which instance should be compared
     * @return true if the lines are similar
     */
    public bool is_similar_to(const Field_Line<T> line) 
    {
        const double& angle = Field_Vector_3D.angle(direction, line.direction).get_real();
        return ((angle < tolerance) || (angle > (std::numbers::pi - tolerance))) && contains(line.zero);
    }

    /** Check if the instance contains a point.
     * @param p point to check
     * @return true if p belongs to the line
     */
    public bool contains(const Field_Vector_3D<T> p) 
    {
        return distance(p).get_real() < tolerance;
    }

    /** Check if the instance contains a point.
     * @param p point to check
     * @return true if p belongs to the line
     */
    public bool contains(const Vector_3D p) 
    {
        return distance(p).get_real() < tolerance;
    }

    /** Compute the distance between the instance and a point.
     * @param p to check
     * @return distance between the instance and the point
     */
    public T distance(const Field_Vector_3D<T> p) 
    {
        const Field_Vector_3D<T> d = p.subtract(zero);
        const Field_Vector_3D<T> n = Field_Vector_3D<>(zero.get_x().get_field().get_one(), d, d.dot_product(direction).negate(), direction);
        return n.get_norm();
    }

    /** Compute the distance between the instance and a point.
     * @param p to check
     * @return distance between the instance and the point
     */
    public T distance(const Vector_3D p) 
    {
        const Field_Vector_3D<T> d = zero.subtract(p).negate();
        const Field_Vector_3D<T> n = Field_Vector_3D<>(zero.get_x().get_field().get_one(), d, d.dot_product(direction).negate(), direction);
        return n.get_norm();
    }

    /** Compute the shortest distance between the instance and another line.
     * @param line line to check against the instance
     * @return shortest distance between the instance and the line
     */
    public T distance(const Field_Line<T> line) 
    {

        const Field_Vector_3D<T> normal = Field_Vector_3D.cross_product(direction, line.direction);
        const T n = normal.get_norm();
        if (n.get_real() < Precision.SAFE_MIN) 
        {
            // lines are parallel
            return distance(line.zero);
        }

        // signed separation of the two parallel planes that contains the lines
        const T offset = line.zero.subtract(zero).dot_product(normal).divide(n);

        return offset.abs();

    }

    /** Compute the point of the instance closest to another line.
     * @param line line to check against the instance
     * @return point of the instance closest to another line
     */
    public Field_Vector_3D<T> closest_point(const Field_Line<T> line) 
    {

        const T cos = direction.dot_product(line.direction);
        const T n = cos.multiply(cos).subtract(1).negate();
        if (n.get_real() < Precision.EPSILON) 
        {
            // the lines are parallel
            return zero;
        }

        const Field_Vector_3D<T> delta0 = line.zero.subtract(zero);
        const T a                     = delta0.dot_product(direction);
        const T& b                     = delta0.dot_product(line.direction);

        return Field_Vector_3D<T>(a.get_field().get_one(), zero, a.subtract(b.multiply(cos)).divide(n), direction);

    }

    /** Get the intersection point of the instance and another line.
     * @param line other line
     * @return intersection point of the instance and the other line
     * or NULL if there are no intersection points
     */
    public Field_Vector_3D<T> intersection(const Field_Line<T> line) 
    {
        const Field_Vector_3D<T> closest = closest_point(line);
        return line.contains(closest) ? closest : NULL;
    }

};