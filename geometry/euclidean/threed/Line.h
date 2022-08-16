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
//import org.hipparchus.exception.;
//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.Vector;
//import org.hipparchus.geometry.euclidean.oned.Euclidean_1D;
//import org.hipparchus.geometry.euclidean.oned.Intervals_Set;
//import org.hipparchus.geometry.euclidean.oned.Vector_1D;
//import org.hipparchus.geometry.partitioning.Embedding;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;

/** The class represent lines in a three dimensional space.

 * <p>Each oriented line is intrinsically associated with an abscissa
 * which is a coordinate on the line. The point at abscissa 0 is the
 * orthogonal projection of the origin on the line, another equivalent
 * way to express this is to say that it is the point of the line
 * which is closest to the origin. Abscissa increases in the line
 * direction.</p>
 *
 * @see #from_direction(Vector_3D, Vector_3D, double)
 * @see #Line(Vector_3D, Vector_3D, double)
 */
class Line : Embedding<Euclidean_3D, Euclidean_1D> 
{

    /** Line direction. */
    private Vector_3D direction;

    /** Line point closest to the origin. */
    private Vector_3D zero;

    /** Tolerance below which points are considered identical. */
    private const double& tolerance;

    /** Build a line from two points.
     * @param p1 first point belonging to the line (this can be any point)
     * @param p2 second point belonging to the line (this can be any point, different from p1)
     * @param tolerance tolerance below which points are considered identical
     * @exception  if the points are equal
     * @see #from_direction(Vector_3D, Vector_3D, double)
     */
    public Line(const Vector_3D p1, const Vector_3D p2, const double& tolerance)
         
        {
        this(tolerance);
        reset(p1, p2);
    }

    /** Copy constructor.
     * <p>The created instance is completely independent from the
     * original instance, it is a deep copy.</p>
     * @param line line to copy
     */
    public Line(const Line& line) 
    {
        this(line.tolerance);
        this.direction = line.direction;
        this.zero      = line.zero;
    }

    /**
     * Private constructor. Just sets the tolerance.
     *
     * @param tolerance below which points are considered identical.
     */
    private Line(const double& tolerance) 
    {
        this.tolerance = tolerance;
    }

    /**
     * Create a line from a point and a direction. Line = {@code point} + t * {@code
     * direction}, where t is any real number.
     *
     * @param point     on the line. Can be any point.
     * @param direction of the line. Must not be the zero vector.
     * @param tolerance below which points are considered identical.
     * @return a Line with the given point and direction.
     * @ if {@code direction} is the zero vector.
     * @see #Line(Vector_3D, Vector_3D, double)
     */
    public static Line from_direction(const Vector_3D point, const Vector_3D& direction, const double& tolerance) 
    {
        const Line line = Line(tolerance);
        line.reset_with_direction(point, direction);
        return line;
    }

    /** Reset the instance as if built from two points.
     * @param p1 first point belonging to the line (this can be any point)
     * @param p2 second point belonging to the line (this can be any point, different from p1)
     * @exception  if the points are equal
     */
    public void reset(const Vector_3D p1, const Vector_3D p2)  
    {
        reset_with_direction(p1, p2.subtract(p1));
    }

    /**
     * Reset the instance as if built from a point and direction.
     *
     * @param p1    point belonging to the line (this can be any point).
     * @param delta direction of the line.
     * @ if {@code delta} is the zero vector.
     */
    private void reset_with_direction(const Vector_3D p1, const Vector_3D delta) 
    {
        const double norm2 = delta.get_norm_sq();
        if (norm2 == 0.0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }
        this.direction = Vector_3D(1.0 / std::sqrt(norm2), delta);
        zero = Vector_3D(1.0, p1, -p1.dot_product(delta) / norm2, delta);
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
    public Line revert() 
    {
        const Line reverted = Line(this);
        reverted.direction = reverted.direction.negate();
        return reverted;
    }

    /** Get the normalized direction vector.
     * @return normalized direction vector
     */
    public Vector_3D get_direction() 
    {
        return direction;
    }

    /** Get the line point closest to the origin.
     * @return line point closest to the origin
     */
    public Vector_3D get_origin() 
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
    public double get_abscissa(const Vector_3D point) 
    {
        return point.subtract(zero).dot_product(direction);
    }

    /** Get one point from the line.
     * @param abscissa desired abscissa for the point
     * @return one point belonging to the line, at specified abscissa
     */
    public Vector_3D point_at(const double& abscissa) 
    {
        return Vector_3D(1.0, zero, abscissa, direction);
    }

    /** Transform a space point into a sub-space point.
     * @param vector n-dimension point of the space
     * @return (n-1)-dimension point of the sub-space corresponding to
     * the specified space point
     */
    public Vector_1D to_sub_space(Vector<Euclidean_3D> vector) 
    {
        return to_sub_space((Point<Euclidean_3D>) vector);
    }

    /** Transform a sub-space point into a space point.
     * @param vector (n-1)-dimension point of the sub-space
     * @return n-dimension point of the space corresponding to the
     * specified sub-space point
     */
    public Vector_3D to_space(Vector<Euclidean_1D> vector) 
    {
        return to_space((Point<Euclidean_1D>) vector);
    }

    /** {@inherit_doc}
     * @see #get_abscissa(Vector_3D)
     */
    //override
    public Vector_1D to_sub_space(const Point<Euclidean_3D> point) 
    {
        return Vector_1D(get_abscissa((Vector_3D) point));
    }

    /** {@inherit_doc}
     * @see #point_atstatic_cast<double>(
     */
    //override
    public Vector_3D to_space(const Point<Euclidean_1D>& point) 
    {
        return point_at(((Vector_1D) point).get_x());
    }

    /** Check if the instance is similar to another line.
     * <p>Lines are considered similar if they contain the same
     * points. This does not mean they are equal since they can have
     * opposite directions.</p>
     * @param line line to which instance should be compared
     * @return true if the lines are similar
     */
    public bool is_similar_to(const Line& line) 
    {
        const double& angle = Vector_3D.angle(direction, line.direction);
        return ((angle < tolerance) || (angle > (std::numbers::pi - tolerance))) && contains(line.zero);
    }

    /** Check if the instance contains a point.
     * @param p point to check
     * @return true if p belongs to the line
     */
    public bool contains(const Vector_3D p) 
    {
        return distance(p) < tolerance;
    }

    /** Compute the distance between the instance and a point.
     * @param p to check
     * @return distance between the instance and the point
     */
    public double distance(const Vector_3D p) 
    {
        const Vector_3D d = p.subtract(zero);
        const Vector_3D n = Vector_3D(1.0, d, -d.dot_product(direction), direction);
        return n.get_norm();
    }

    /** Compute the shortest distance between the instance and another line.
     * @param line line to check against the instance
     * @return shortest distance between the instance and the line
     */
    public double distance(const Line& line) 
    {

        const Vector_3D normal = Vector_3D.cross_product(direction, line.direction);
        const double n = normal.get_norm();
        if (n < Precision.SAFE_MIN) 
        {
            // lines are parallel
            return distance(line.zero);
        }

        // signed separation of the two parallel planes that contains the lines
        const double offset = line.zero.subtract(zero).dot_product(normal) / n;

        return std::abs(offset);

    }

    /** Compute the point of the instance closest to another line.
     * @param line line to check against the instance
     * @return point of the instance closest to another line
     */
    public Vector_3D closest_point(const Line& line) 
    {

        const double cos = direction.dot_product(line.direction);
        const double n = 1 - cos * cos;
        if (n < Precision.EPSILON) 
        {
            // the lines are parallel
            return zero;
        }

        const Vector_3D delta0 = line.zero.subtract(zero);
        const double& a        = delta0.dot_product(direction);
        const double b        = delta0.dot_product(line.direction);

        return Vector_3D(1, zero, (a - b * cos) / n, direction);

    }

    /** Get the intersection point of the instance and another line.
     * @param line other line
     * @return intersection point of the instance and the other line
     * or NULL if there are no intersection points
     */
    public Vector_3D intersection(const Line& line) 
    {
        const Vector_3D closest = closest_point(line);
        return line.contains(closest) ? closest : NULL;
    }

    /** Build a sub-line covering the whole line.
     * @return a sub-line covering the whole line
     */
    public Sub_Line whole_line() 
    {
        return Sub_Line(this, Intervals_Set(tolerance));
    }

}


