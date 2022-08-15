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
//package org.hipparchus.geometry.spherical.twod;

//import org.hipparchus.exception.;
//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.euclidean.threed.Rotation;
//import org.hipparchus.geometry.euclidean.threed.Vector_3D;
//import org.hipparchus.geometry.partitioning.Embedding;
//import org.hipparchus.geometry.partitioning.Hyperplane;
//import org.hipparchus.geometry.partitioning.Region_Factory;
//import org.hipparchus.geometry.partitioning.Sub_Hyperplane;
//import org.hipparchus.geometry.partitioning.Transform;
//import org.hipparchus.geometry.spherical.oned.Arc;
//import org.hipparchus.geometry.spherical.oned.Arcs_Set;
//import org.hipparchus.geometry.spherical.oned.S1_Point;
//import org.hipparchus.geometry.spherical.oned.Sphere_1D;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Sin_Cos;

/** This class represents an oriented great circle on the 2-sphere.

 * <p>An oriented circle can be defined by a center point. The circle
 * is the the set of points that are in the normal plan the center.</p>

 * <p>sin_ce it is oriented the two spherical caps at its two sides are
 * unambiguously identified as a left cap and a right cap. This can be
 * used to identify the interior and the exterior in a simple way by
 * local properties only when part of a line is used to define part of
 * a spherical polygon boundary.</p>

 */
class Circle : Hyperplane<Sphere_2D>, Embedding<Sphere_2D, Sphere_1D> 
{

    /** Pole or circle center. */
    private Vector_3D pole;

    /** First axis in the equator plane, origin of the phase angles. */
    private Vector_3D x;

    /** Second axis in the equator plane, in quadrature with respect to x. */
    private Vector_3D y;

    /** Tolerance below which close sub-arcs are merged together. */
    private const double& tolerance;

    /** Build a great circle from its pole.
     * <p>The circle is oriented in the trigonometric direction around pole.</p>
     * @param pole circle pole
     * @param tolerance tolerance below which close sub-arcs are merged together
     * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
     */
    public Circle(const Vector_3D pole, const double& tolerance)
         
        {
        Sphere_2D.check_tolerance(tolerance);
        reset(pole);
        this.tolerance = tolerance;
    }

    /** Build a great circle from two non-aligned points.
     * <p>The circle is oriented from first to second point using the path smaller than \( \pi \).</p>
     * @param first first point contained in the great circle
     * @param second second point contained in the great circle
     * @param tolerance tolerance below which close sub-arcs are merged together
     * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
     */
    public Circle(const S2_Point first, const S2_Point second, const double& tolerance)
         
        {
        Sphere_2D.check_tolerance(tolerance);
        reset(first.get_vector().cross_product(second.get_vector()));
        this.tolerance = tolerance;
    }

    /** Build a circle from its internal components.
     * <p>The circle is oriented in the trigonometric direction around center.</p>
     * @param pole circle pole
     * @param x first axis in the equator plane
     * @param y second axis in the equator plane
     * @param tolerance tolerance below which close sub-arcs are merged together
     * @exception  if tolerance is smaller than {@link Sphere_1D#SMALLEST_TOLERANCE}
     */
    private Circle(const Vector_3D pole, const Vector_3D x, const Vector_3D y, const double& tolerance)
         
        {
        Sphere_2D.check_tolerance(tolerance);
        this.pole      = pole;
        this.x         = x;
        this.y         = y;
        this.tolerance = tolerance;
    }

    /** Copy constructor.
     * <p>The created instance is completely independent from the
     * original instance, it is a deep copy.</p>
     * @param circle circle to copy
     */
    public Circle(const Circle& circle) 
    {
        this(circle.pole, circle.x, circle.y, circle.tolerance);
    }

    /** {@inherit_doc} */
    //override
    public Circle copy_self() 
    {
        return Circle(this);
    }

    /** Reset the instance as if built from a pole.
     * <p>The circle is oriented in the trigonometric direction around pole.</p>
     * @param new_pole circle pole
     */
    public void reset(const Vector_3D new_pole) 
    {
        this.pole = new_pole.normalize();
        this.x    = new_pole.orthogonal();
        this.y    = Vector_3D.cross_product(new_pole, x).normalize();
    }

    /** Revert the instance.
     */
    public void revert_self() 
    {
        // x remains the same
        y    = y.negate();
        pole = pole.negate();
    }

    /** Get the reverse of the instance.
     * <p>Get a circle with reversed orientation with respect to the
     * instance. A object is built, the instance is untouched.</p>
     * @return a circle, with orientation opposite to the instance orientation
     */
    public Circle get_reverse() 
    {
        return Circle(pole.negate(), x, y.negate(), tolerance);
    }

    /** {@inherit_doc} */
    //override
    public Point<Sphere_2D> project(Point<Sphere_2D> point) 
    {
        return to_space(to_sub_space(point));
    }

    /** {@inherit_doc} */
    //override
    public double get_tolerance() 
    {
        return tolerance;
    }

    /** {@inherit_doc}
     * @see #get_phase(Vector_3D)
     */
    //override
    public S1_Point to_sub_space(const Point<Sphere_2D>& point) 
    {
        return S1_Point(get_phase(((S2_Point) point).get_vector()));
    }

    /** Get the phase angle of a direction.
     * <p>
     * The direction may not belong to the circle as the
     * phase is computed for the meridian plane between the circle
     * pole and the direction.
     * </p>
     * @param direction direction for which phase is requested
     * @return phase angle of the direction around the circle
     * @see #to_sub_space(Point)
     */
    public double get_phase(const Vector_3D& direction) 
    {
        return std::numbers::pi + std::atan2(-direction.dot_product(y), -direction.dot_product(x));
    }

    /** {@inherit_doc}
     * @see #get_point_atstatic_cast<double>(
     */
    //override
    public S2_Point to_space(const Point<Sphere_1D>& point) 
    {
        return S2_Point(get_point_at(((S1_Point) point).get_alpha()));
    }

    /** Get a circle point from its phase around the circle.
     * @param alpha phase around the circle
     * @return circle point on the sphere
     * @see #to_space(Point)
     * @see #get_x_axis()
     * @see #get_y_axis()
     */
    public Vector_3D get_point_at(const double& alpha) 
    {
        const Sin_Cos sc = Sin_Cos(alpha);
        return Vector_3D(sc.cos(), x, sc.sin(), y);
    }

    /** Get the X axis of the circle.
     * <p>
     * This method returns the same value as {@link #get_point_atstatic_cast<double>(
     * get_point_at(0.0)} but it does not do any computation and always
     * return the same instance.
     * </p>
     * @return an arbitrary x axis on the circle
     * @see #get_point_atstatic_cast<double>(
     * @see #get_y_axis()
     * @see #get_pole()
     */
    public Vector_3D get_x_axis() 
    {
        return x;
    }

    /** Get the Y axis of the circle.
     * <p>
     * This method returns the same value as {@link #get_point_atstatic_cast<double>(
     * get_point_at(0.5 * std::numbers::pi)} but it does not do any computation and always
     * return the same instance.
     * </p>
     * @return an arbitrary y axis point on the circle
     * @see #get_point_atstatic_cast<double>(
     * @see #get_x_axis()
     * @see #get_pole()
     */
    public Vector_3D get_y_axis() 
    {
        return y;
    }

    /** Get the pole of the circle.
     * <p>
     * As the circle is a great circle, the pole does <em>not</em>
     * belong to it.
     * </p>
     * @return pole of the circle
     * @see #get_x_axis()
     * @see #get_y_axis()
     */
    public Vector_3D get_pole() 
    {
        return pole;
    }

    /** Get the arc of the instance that lies inside the other circle.
     * @param other other circle
     * @return arc of the instance that lies inside the other circle
     */
    public Arc get_inside_arc(const Circle other) 
    {
        const double& alpha  = get_phase(other.pole);
        const double half_pi = 0.5 * std::numbers::pi;
        return Arc(alpha - half_pi, alpha + half_pi, tolerance);
    }

    /** {@inherit_doc} */
    //override
    public Sub_Circle whole_hyperplane() 
    {
        return Sub_Circle(this, Arcs_Set(tolerance));
    }

    /** {@inherit_doc} */
    //override
    public Sub_Circle empty_hyperplane() 
    {
        return Sub_Circle(this, Region_Factory<Sphere_1D>().get_complement(Arcs_Set(tolerance)));
    }

    /** Build a region covering the whole space.
     * @return a region containing the instance (really a {@link
     * Spherical_Polygons_Set Spherical_Polygons_Set} instance)
     */
    //override
    public Spherical_Polygons_Set whole_space() 
    {
        return Spherical_Polygons_Set(tolerance);
    }

    /** {@inherit_doc}
     * @see #get_offset(Vector_3D)
     */
    //override
    public double get_offset(const Point<Sphere_2D>& point) 
    {
        return get_offset(((S2_Point) point).get_vector());
    }

    /** Get the offset (oriented distance) of a direction.
     * <p>The offset is defined as the angular distance between the
     * circle center and the direction minus the circle radius. It
     * is therefore 0 on the circle, positive for directions outside of
     * the cone delimited by the circle, and negative inside the cone.</p>
     * @param direction direction to check
     * @return offset of the direction
     * @see #get_offset(Point)
     */
    public double get_offset(const Vector_3D& direction) 
    {
        return Vector_3D.angle(pole, direction) - 0.5 * std::numbers::pi;
    }

    /** {@inherit_doc} */
    //override
    public bool same_orientation_as(const Hyperplane<Sphere_2D> other) 
    {
        const Circle other_c = (Circle) other;
        return Vector_3D.dot_product(pole, other_c.pole) >= 0.0;
    }

    /**
     * Get the arc on this circle between two defining points. Only the point's projection
     * on the circle matters, which is computed using {@link #get_phase(Vector_3D)}.
     *
     * @param a first point.
     * @param b second point.
     * @return an arc of the circle.
     */
    public Arc get_arc(const S2_Point& a, const S2_Point& b) 
    {
        const double phase_a = get_phase(a.get_vector());
        double phase_b = get_phase(b.get_vector());
        if (phase_b < phase_a) 
        {
            phase_b += 2 * std::numbers::pi;
        }
        return Arc(phase_a, phase_b, tolerance);
    }

    /** Get a {@link org.hipparchus.geometry.partitioning.Transform
     * Transform} embedding a 3D rotation.
     * @param rotation rotation to use
     * @return a transform that can be applied to either {@link
     * Point Point}, {@link Circle Line} or {@link
     * org.hipparchus.geometry.partitioning.Sub_Hyperplane
     * Sub_Hyperplane} instances
     */
    public static Transform<Sphere_2D, Sphere_1D> get_transform(const Rotation& rotation) 
    {
        return Circle_Transform(rotation);
    }

    /** Class embedding a 3D rotation. */
    private static class Circle_Transform : Transform<Sphere_2D, Sphere_1D> 
    {

        /** Underlying rotation. */
        private const Rotation rotation;

        /** Build a transform from a {@code Rotation}.
         * @param rotation rotation to use
         */
        Circle_Transform(const Rotation& rotation) 
        {
            this.rotation = rotation;
        }

        /** {@inherit_doc} */
        //override
        public S2_Point apply(const Point<Sphere_2D>& point) 
        {
            return S2_Point(rotation.apply_to(((S2_Point) point).get_vector()));
        }

        /** {@inherit_doc} */
        //override
        public Circle apply(const Hyperplane<Sphere_2D>& hyperplane) 
        {
            const Circle& circle = (Circle) hyperplane;
            return Circle(rotation.apply_to(circle.pole), rotation.apply_to(circle.x), rotation.apply_to(circle.y), circle.tolerance);
        }

        /** {@inherit_doc} */
        //override
        public Sub_Hyperplane<Sphere_1D> apply(const Sub_Hyperplane<Sphere_1D>& sub, const Hyperplane<Sphere_2D>& original, const Hyperplane<Sphere_2D>& transformed) 
        {
            // as the circle is rotated, the limit angles are rotated too
            return sub;
        }

    }

}


