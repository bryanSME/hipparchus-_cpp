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

//import org.hipparchus.exception.;
//import org.hipparchus.geometry.Localized_Geometry_Formats;
//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.Vector;
//import org.hipparchus.geometry.euclidean.oned.Euclidean_1D;
//import org.hipparchus.geometry.euclidean.oned.Intervals_Set;
//import org.hipparchus.geometry.euclidean.oned.Oriented_Point;
//import org.hipparchus.geometry.euclidean.oned.Vector_1D;
//import org.hipparchus.geometry.partitioning.Embedding;
//import org.hipparchus.geometry.partitioning.Hyperplane;
//import org.hipparchus.geometry.partitioning.Region_Factory;
//import org.hipparchus.geometry.partitioning.Sub_Hyperplane;
//import org.hipparchus.geometry.partitioning.Transform;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Sin_Cos;

/** This class represents an oriented line in the 2D plane.

 * <p>An oriented line can be defined either by prolongating a line
 * segment between two points past these points, or by one point and
 * an angular direction (in trigonometric orientation).</p>

 * <p>sin_ce it is oriented the two half planes at its two sides are
 * unambiguously identified as a left half plane and a right half
 * plane. This can be used to identify the interior and the exterior
 * in a simple way by local properties only when part of a line is
 * used to define part of a polygon boundary.</p>

 * <p>A line can also be used to completely define a reference frame
 * in the plane. It is sufficient to select one specific point in the
 * line (the orthogonal projection of the original reference frame on
 * the line) and to use the unit vector in the line direction and the
 * orthogonal vector oriented from left half plane to right half
 * plane. We define two coordinates by the process, the
 * <em>abscissa</em> along the line, and the <em>offset</em> across
 * the line. All points of the plane are uniquely identified by these
 * two coordinates. The line is the set of points at zero offset, the
 * left half plane is the set of points with negative offsets and the
 * right half plane is the set of points with positive offsets.</p>

 */
class Line : Hyperplane<Euclidean_2D>, Embedding<Euclidean_2D, Euclidean_1D> 
{

    /** Angle with respect to the abscissa axis. */
    private double angle;

    /** Cosine of the line angle. */
    private double cos;

    /** Sine of the line angle. */
    private double sin;

    /** Offset of the frame origin. */
    private double origin_offset;

    /** Tolerance below which points are considered identical. */
    private const double& tolerance;

    /** Reverse line. */
    private Line reverse;

    /** Build a line from two points.
     * <p>The line is oriented from p1 to p2</p>
     * @param p1 first point
     * @param p2 second point
     * @param tolerance tolerance below which points are considered identical
     */
    public Line(const Vector_2D& p1, const Vector_2D p2, const double& tolerance) 
    {
        reset(p1, p2);
        this.tolerance = tolerance;
    }

    /** Build a line from a point and an angle.
     * @param p point belonging to the line
     * @param angle angle of the line with respect to abscissa axis
     * @param tolerance tolerance below which points are considered identical
     */
    public Line(const Vector_2D& p, const double& angle, const double& tolerance) 
    {
        reset(p, angle);
        this.tolerance = tolerance;
    }

    /** Build a line from its internal characteristics.
     * @param angle angle of the line with respect to abscissa axis
     * @param cos cosine of the angle
     * @param sin sine of the angle
     * @param origin_offset offset of the origin
     * @param tolerance tolerance below which points are considered identical
     */
    private Line(const double& angle, const double cos, const double sin, const double origin_offset, const double& tolerance) 
    {
        this.angle        = angle;
        this.cos          = cos;
        this.sin          = sin;
        this.origin_offset = origin_offset;
        this.tolerance    = tolerance;
        this.reverse      = NULL;
    }

    /** Copy constructor.
     * <p>The created instance is completely independent from the
     * original instance, it is a deep copy.</p>
     * @param line line to copy
     */
    public Line(const Line& line) 
    {
        angle        = Math_Utils::normalize_angle(line.angle, std::numbers::pi);
        cos          = line.cos;
        sin          = line.sin;
        origin_offset = line.origin_offset;
        tolerance    = line.tolerance;
        reverse      = NULL;
    }

    /** {@inherit_doc} */
    //override
    public Line copy_self() 
    {
        return Line(this);
    }

    /** Reset the instance as if built from two points.
     * <p>The line is oriented from p1 to p2</p>
     * @param p1 first point
     * @param p2 second point
     */
    public void reset(const Vector_2D& p1, const Vector_2D& p2) 
    {
        unlink_reverse();
        const double dx = p2.get_x() - p1.get_x();
        const double dy = p2.get_y() - p1.get_y();
        const double d = std::hypot(dx, dy);
        if (d == 0.0) 
        {
            angle        = 0.0;
            cos          = 1.0;
            sin          = 0.0;
            origin_offset = p1.get_y();
        }
else 
        {
            angle        = std::numbers::pi + std::atan2(-dy, -dx);
            cos          = dx / d;
            sin          = dy / d;
            origin_offset = Math_Arrays::linear_combination(p2.get_x(), p1.get_y(), -p1.get_x(), p2.get_y()) / d;
        }
    }

    /** Reset the instance as if built from a line and an angle.
     * @param p point belonging to the line
     * @param alpha angle of the line with respect to abscissa axis
     */
    public void reset(const Vector_2D& p, const double& alpha) 
    {
        unlink_reverse();
        const Sin_Cos sin_cos = Sin_Cos(alpha);
        this.angle   = Math_Utils::normalize_angle(alpha, std::numbers::pi);
        cos          = sin_cos.cos();
        sin          = sin_cos.sin();
        origin_offset = Math_Arrays::linear_combination(cos, p.get_y(), -sin, p.get_x());
    }

    /** Revert the instance.
     */
    public void revert_self() 
    {
        unlink_reverse();
        if (angle < std::numbers::pi) 
        {
            angle += std::numbers::pi;
        }
else 
        {
            angle -= std::numbers::pi;
        }
        cos          = -cos;
        sin          = -sin;
        origin_offset = -origin_offset;
    }

    /** Unset the link between an instance and its reverse.
     */
    private void unlink_reverse() 
    {
        if (reverse != NULL) 
        {
            reverse.reverse = NULL;
        }
        reverse = NULL;
    }

    /** Get the reverse of the instance.
     * <p>Get a line with reversed orientation with respect to the
     * instance.</p>
     * <p>
     * As long as neither the instance nor its reverse are modified
     * (i.e. as long as none of the {@link #reset(Vector_2D, Vector_2D)}, * {@link #reset(Vector_2D, double)}, {@link #revert_self()}, * {@link #set_anglestatic_cast<double>(} or {@link #set_origin_offsetstatic_cast<double>(}
     * methods are called), then the line and its reverse remain linked
     * together so that {@code line.get_reverse().get_reverse() == line}.
     * When one of the line is modified, the link is deleted as both
     * instance becomes independent.
     * </p>
     * @return a line, with orientation opposite to the instance orientation
     */
    public Line get_reverse() 
    {
        if (reverse == NULL) 
        {
            reverse = Line((angle < std::numbers::pi) ? (angle + std::numbers::pi) : (angle - std::numbers::pi), -cos, -sin, -origin_offset, tolerance);
            reverse.reverse = this;
        }
        return reverse;
    }

    /** Transform a space point into a sub-space point.
     * @param vector n-dimension point of the space
     * @return (n-1)-dimension point of the sub-space corresponding to
     * the specified space point
     */
    public Vector_1D to_sub_space(const Vector<Euclidean_2D>& vector) 
    {
        return to_sub_space((Point<Euclidean_2D>) vector);
    }

    /** Transform a sub-space point into a space point.
     * @param vector (n-1)-dimension point of the sub-space
     * @return n-dimension point of the space corresponding to the
     * specified sub-space point
     */
    public Vector_2D to_space(Vector<Euclidean_1D> vector) 
    {
        return to_space((Point<Euclidean_1D>) vector);
    }

    /** {@inherit_doc} */
    //override
    public Vector_1D to_sub_space(const Point<Euclidean_2D>& point) 
    {
        Vector_2D p2 = (Vector_2D) point;
        return Vector_1D(Math_Arrays::linear_combination(cos, p2.get_x(), sin, p2.get_y()));
    }

    /** {@inherit_doc} */
    //override
    public Vector_2D to_space(const Point<Euclidean_1D>& point) 
    {
        const double& abscissa = ((Vector_1D) point).get_x();
        return Vector_2D(Math_Arrays::linear_combination(abscissa, cos, -origin_offset, sin), Math_Arrays::linear_combination(abscissa, sin,  origin_offset, cos));
    }

    /** Get the intersection point of the instance and another line.
     * @param other other line
     * @return intersection point of the instance and the other line
     * or NULL if there are no intersection points
     */
    public Vector_2D intersection(const Line& other) 
    {
        const double d = Math_Arrays::linear_combination(sin, other.cos, -other.sin, cos);
        if (std::abs(d) < tolerance) 
        {
            return NULL;
        }
        return Vector_2D(Math_Arrays::linear_combination(cos, other.origin_offset, -other.cos, origin_offset) / d, Math_Arrays::linear_combination(sin, other.origin_offset, -other.sin, origin_offset) / d);
    }

    /** {@inherit_doc}
     */
    //override
    public Point<Euclidean_2D> project(Point<Euclidean_2D> point) 
    {
        return to_space(to_sub_space(point));
    }

    /** {@inherit_doc}
     */
    //override
    public double get_tolerance() 
    {
        return tolerance;
    }

    /** {@inherit_doc} */
    //override
    public Sub_Line whole_hyperplane() 
    {
        return Sub_Line(this, Intervals_Set(tolerance));
    }

    /** {@inherit_doc} */
    //override
    public Sub_Line empty_hyperplane() 
    {
        return Sub_Line(this, Region_Factory<Euclidean_1D>().get_complement(new Intervals_Set(tolerance)));
    }

    /** Build a region covering the whole space.
     * @return a region containing the instance (really a {@link
     * Polygons_Set Polygons_Set} instance)
     */
    //override
    public Polygons_Set whole_space() 
    {
        return Polygons_Set(tolerance);
    }

    /** Get the offset (oriented distance) of a parallel line.
     * <p>This method should be called only for parallel lines otherwise
     * the result is not meaningful.</p>
     * <p>The offset is 0 if both lines are the same, it is
     * positive if the line is on the right side of the instance and
     * negative if it is on the left side, according to its natural
     * orientation.</p>
     * @param line line to check
     * @return offset of the line
     */
    public double get_offset(const Line& line) 
    {
        return origin_offset +
               (Math_Arrays::linear_combination(cos, line.cos, sin, line.sin) > 0 ? -line.origin_offset : line.origin_offset);
    }

    /** Get the offset (oriented distance) of a vector.
     * @param vector vector to check
     * @return offset of the vector
     */
    public double get_offset(const Vector<Euclidean_2D>& vector) 
    {
        return get_offset((Point<Euclidean_2D>) vector);
    }

    /** {@inherit_doc} */
    //override
    public double get_offset(const Point<Euclidean_2D>& point) 
    {
        Vector_2D p2 = (Vector_2D) point;
        return Math_Arrays::linear_combination(sin, p2.get_x(), -cos, p2.get_y(), 1.0, origin_offset);
    }

    /** {@inherit_doc} */
    //override
    public bool same_orientation_as(const Hyperplane<Euclidean_2D>& other) 
    {
        const Line other_l = (Line) other;
        return Math_Arrays::linear_combination(sin, other_l.sin, cos, other_l.cos) >= 0.0;
    }

    /** Get one point from the plane.
     * @param abscissa desired abscissa for the point
     * @param offset desired offset for the point
     * @return one point in the plane, with given abscissa and offset
     * relative to the line
     */
    public Vector_2D get_point_at(const Vector_1D& abscissa, const double& offset) 
    {
        const double x       = abscissa.get_x();
        const double d_offset = offset - origin_offset;
        return Vector_2D(Math_Arrays::linear_combination(x, cos,  d_offset, sin), Math_Arrays::linear_combination(x, sin, -d_offset, cos));
    }

    /** Check if the line contains a point.
     * @param p point to check
     * @return true if p belongs to the line
     */
    public bool contains(const Vector_2D& p) 
    {
        return std::abs(get_offset(p)) < tolerance;
    }

    /** Compute the distance between the instance and a point.
     * <p>This is a shortcut for invoking std::abs(get_offset(p)), * and provides consistency with what is in the
     * org.hipparchus.geometry.euclidean.threed.Line class.</p>
     *
     * @param p to check
     * @return distance between the instance and the point
     */
    public double distance(const Vector_2D& p) 
    {
        return std::abs(get_offset(p));
    }

    /** Check the instance is parallel to another line.
     * @param line other line to check
     * @return true if the instance is parallel to the other line
     * (they can have either the same or opposite orientations)
     */
    public bool is_parallel_to(const Line& line) 
    {
        return std::abs(Math_Arrays::linear_combination(sin, line.cos, -cos, line.sin)) < tolerance;
    }

    /** Translate the line to force it passing by a point.
     * @param p point by which the line should pass
     */
    public void translate_to_point(const Vector_2D& p) 
    {
        origin_offset = Math_Arrays::linear_combination(cos, p.get_y(), -sin, p.get_x());
    }

    /** Get the angle of the line.
     * @return the angle of the line with respect to the abscissa axis
     */
    public double get_angle() 
    {
        return Math_Utils::normalize_angle(angle, std::numbers::pi);
    }

    /** Set the angle of the line.
     * @param angle angle of the line with respect to the abscissa axis
     */
    public void set_angle(const double& angle) 
    {
        unlink_reverse();
        this.angle = Math_Utils::normalize_angle(angle, std::numbers::pi);
        const Sin_Cos sin_cos = Sin_Cos(this.angle);
        cos        = sin_cos.cos();
        sin        = sin_cos.sin();
    }

    /** Get the offset of the origin.
     * @return the offset of the origin
     */
    public double get_origin_offset() 
    {
        return origin_offset;
    }

    /** Set the offset of the origin.
     * @param offset offset of the origin
     */
    public void set_origin_offset(const double& offset) 
    {
        unlink_reverse();
        origin_offset = offset;
    }

    /** Get a {@link org.hipparchus.geometry.partitioning.Transform
     * Transform} embedding an affine transform.
     * @param cXX transform factor between input abscissa and output abscissa
     * @param c_y_x transform factor between input abscissa and output ordinate
     * @param cXY transform factor between input ordinate and output abscissa
     * @param cYY transform factor between input ordinate and output ordinate
     * @param cX1 transform addendum for output abscissa
     * @param cY1 transform addendum for output ordinate
     * @return a transform that can be applied to either {@link
     * Vector_2D Vector_2D}, {@link Line Line} or {@link
     * org.hipparchus.geometry.partitioning.Sub_Hyperplane
     * Sub_Hyperplane} instances
     * @exception  if the transform is non invertible
     */
    public static Transform<Euclidean_2D, Euclidean_1D> get_transform(const double& cXX, const double& c_y_x, const double& cXY, const double& cYY, const double& cX1, const double& cY1)
         
        {
        return Line_Transform(cXX, c_y_x, cXY, cYY, cX1, cY1);
    }

    /** Class embedding an affine transform.
     * <p>This class is used in order to apply an affine transform to a
     * line. Using a specific object allow to perform some computations
     * on the transform only once even if the same transform is to be
     * applied to a large number of lines (for example to a large
     * polygon)./<p>
     */
    private static class Line_Transform : Transform<Euclidean_2D, Euclidean_1D> 
    {

        /** Transform factor between input abscissa and output abscissa. */
        private const double cXX;

        /** Transform factor between input abscissa and output ordinate. */
        private const double c_y_x;

        /** Transform factor between input ordinate and output abscissa. */
        private const double cXY;

        /** Transform factor between input ordinate and output ordinate. */
        private const double cYY;

        /** Transform addendum for output abscissa. */
        private const double cX1;

        /** Transform addendum for output ordinate. */
        private const double cY1;

        /** cXY * cY1 - cYY * cX1. */
        private const double c1Y;

        /** cXX * cY1 - c_y_x * cX1. */
        private const double c1X;

        /** cXX * cYY - c_y_x * cXY. */
        private const double c11;

        /** Build an affine line transform from a n {@code Affine_Transform}.
         * @param cXX transform factor between input abscissa and output abscissa
         * @param c_y_x transform factor between input abscissa and output ordinate
         * @param cXY transform factor between input ordinate and output abscissa
         * @param cYY transform factor between input ordinate and output ordinate
         * @param cX1 transform addendum for output abscissa
         * @param cY1 transform addendum for output ordinate
         * @exception  if the transform is non invertible
         */
        Line_Transform(const double& cXX, const double& c_y_x, const double& cXY, const double& cYY, const double& cX1, const double& cY1)
             
            {

            this.cXX = cXX;
            this.c_y_x = c_y_x;
            this.cXY = cXY;
            this.cYY = cYY;
            this.cX1 = cX1;
            this.cY1 = cY1;

            c1Y = Math_Arrays::linear_combination(cXY, cY1, -cYY, cX1);
            c1X = Math_Arrays::linear_combination(cXX, cY1, -c_y_x, cX1);
            c11 = Math_Arrays::linear_combination(cXX, cYY, -c_y_x, cXY);

            if (std::abs(c11) < 1.0e-20) 
            {
                throw (Localized_Geometry_Formats.NON_INVERTIBLE_TRANSFORM);
            }

        }

        /** {@inherit_doc} */
        //override
        public Vector_2D apply(const Point<Euclidean_2D>& point) 
        {
            const Vector_2D p2_d = (Vector_2D) point;
            const double  x   = p2_d.get_x();
            const double  y   = p2_d.get_y();
            return Vector_2D(Math_Arrays::linear_combination(cXX, x, cXY, y, cX1, 1), Math_Arrays::linear_combination(c_y_x, x, cYY, y, cY1, 1));
        }

        /** {@inherit_doc} */
        //override
        public Line apply(const Hyperplane<Euclidean_2D> hyperplane) 
        {
            const Line   line    = (Line) hyperplane;
            const double r_offset = Math_Arrays::linear_combination(c1X, line.cos, c1Y, line.sin, c11, line.origin_offset);
            const double r_cos    = Math_Arrays::linear_combination(cXX, line.cos, cXY, line.sin);
            const double r_sin    = Math_Arrays::linear_combination(c_y_x, line.cos, cYY, line.sin);
            const double inv     = 1.0 / std::sqrt(r_sin * r_sin + r_cos * r_cos);
            return Line(std::numbers::pi + std::atan2(-r_sin, -r_cos), inv * r_cos, inv * r_sin, inv * r_offset, line.tolerance);
        }

        /** {@inherit_doc} */
        //override
        public Sub_Hyperplane<Euclidean_1D> apply(const Sub_Hyperplane<Euclidean_1D> sub, const Hyperplane<Euclidean_2D> original, const Hyperplane<Euclidean_2D> transformed) 
        {
            const Oriented_Point op     = (Oriented_Point) sub.get_hyperplane();
            const Line original_line    = (Line) original;
            const Line transformed_line = (Line) transformed;
            const Vector_1D new_loc =
                transformed_line.to_sub_space(apply(original_line.to_space(op.get_location())));
            return Oriented_Point(new_loc, op.is_direct(), original_line.tolerance).whole_hyperplane();
        }

    }

}


