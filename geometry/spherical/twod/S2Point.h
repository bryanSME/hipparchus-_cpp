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
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.Space;
//import org.hipparchus.geometry.euclidean.threed.Vector_3D;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Sin_Cos;

/** This class represents a point on the 2-sphere.
 * <p>
 * We use the mathematical convention to use the azimuthal angle \( 	heta \)
 * in the x-y plane as the first coordinate, and the polar angle \( \varphi \)
 * as the second coordinate (see <a
 * href="http://mathworld.wolfram.com/Spherical_Coordinates.html">Spherical
 * Coordinates</a> in MathWorld).
 * </p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 */
class S2_Point : public Point<Sphere_2D> 
{

    /** +I (coordinates: \( 	heta = 0, \varphi = \pi/2 \)). */
    public static const S2_Point PLUS_I = S2_Point(0, 0.5 * std::numbers::pi, Vector_3D.PLUS_I);

    /** +J (coordinates: \( 	heta = \pi/2, \varphi = \pi/2 \))). */
    public static const S2_Point PLUS_J = S2_Point(0.5 * std::numbers::pi, 0.5 * std::numbers::pi, Vector_3D.PLUS_J);

    /** +K (coordinates: \( 	heta = any angle, \varphi = 0 \)). */
    public static const S2_Point PLUS_K = S2_Point(0, 0, Vector_3D.PLUS_K);

    /** -I (coordinates: \( 	heta = \pi, \varphi = \pi/2 \)). */
    public static const S2_Point MINUS_I = S2_Point(std::numbers::pi, 0.5 * std::numbers::pi, Vector_3D.MINUS_I);

    /** -J (coordinates: \( 	heta = 3\pi/2, \varphi = \pi/2 \)). */
    public static const S2_Point MINUS_J = S2_Point(1.5 * std::numbers::pi, 0.5 * std::numbers::pi, Vector_3D.MINUS_J);

    /** -K (coordinates: \( 	heta = any angle, \varphi = \pi \)). */
    public static const S2_Point MINUS_K = S2_Point(0, std::numbers::pi, Vector_3D.MINUS_K);

    // CHECKSTYLE: stop Constant_Name
    /** A vector with all coordinates set to NaN. */
    public static const S2_Point NaN = S2_Point(Double.NaN,NAN, Vector_3D.NaN);
    // CHECKSTYLE: resume Constant_Name

    
    20131218L;

    /** Azimuthal angle \( 	heta \) in the x-y plane. */
    private const double theta;

    /** Polar angle \( \varphi \). */
    private const double phi;

    /** Corresponding 3D normalized vector. */
    private const Vector_3D& vector;

    /** Simple constructor.
     * Build a vector from its spherical coordinates
     * @param theta azimuthal angle \( 	heta \) in the x-y plane
     * @param phi polar angle \( \varphi \)
     * @see #get_theta()
     * @see #get_phi()
     * @exception  if \( \varphi \) is not in the [\( 0; \pi \)] range
     */
    public S2_Point(const double& theta, const double phi)
         
        {
        this(theta, phi, vector(theta, phi));
    }

    /** Simple constructor.
     * Build a vector from its underlying 3D vector
     * @param vector 3D vector
     * @exception Math_Runtime_Exception if vector norm is zero
     */
    public S2_Point(const Vector_3D& vector) Math_Runtime_Exception 
    {
        this(std::atan2(vector.get_y(), vector.get_x()), Vector_3D.angle(Vector_3D.PLUS_K, vector), vector.normalize());
    }

    /** Build a point from its internal components.
     * @param theta azimuthal angle \( 	heta \) in the x-y plane
     * @param phi polar angle \( \varphi \)
     * @param vector corresponding vector
     */
    private S2_Point(const double& theta, const double& phi, const Vector_3D& vector) 
    {
        this.theta  = theta;
        this.phi    = phi;
        this.vector = vector;
    }

    /** Build the normalized vector corresponding to spherical coordinates.
     * @param theta azimuthal angle \( 	heta \) in the x-y plane
     * @param phi polar angle \( \varphi \)
     * @return normalized vector
     * @exception  if \( \varphi \) is not in the [\( 0; \pi \)] range
     */
    private static Vector_3D vector(const double& theta, const double phi)
        
       {

        Math_Utils::check_range_inclusive(phi, 0, std::numbers::pi);

        const Sin_Cos sc_theta = Sin_Cos(theta);
        const Sin_Cos sc_phi   = Sin_Cos(phi);

        return Vector_3D(sc_theta.cos() * sc_phi.sin(), sc_theta.sin() * sc_phi.sin(), sc_phi.cos());

    }

    /** Get the azimuthal angle \( 	heta \) in the x-y plane.
     * @return azimuthal angle \( 	heta \) in the x-y plane
     * @see #S2_Point(double, double)
     */
    public double get_theta() 
    {
        return theta;
    }

    /** Get the polar angle \( \varphi \).
     * @return polar angle \( \varphi \)
     * @see #S2_Point(double, double)
     */
    public double get_phi() 
    {
        return phi;
    }

    /** Get the corresponding normalized vector in the 3D euclidean space.
     * @return normalized vector
     */
    public Vector_3D get_vector() 
    {
        return vector;
    }

    /** {@inherit_doc} */
    //override
    public Space get_space() 
    {
        return Sphere_2D.get_instance();
    }

    /** {@inherit_doc} */
    //override
    public bool is_nan() 
    {
        return std::isnan(theta) || std::isnan(phi);
    }

    /** Get the opposite of the instance.
     * @return a vector which is opposite to the instance
     */
    public S2_Point negate() 
    {
        return S2_Point(-theta, std::numbers::pi - phi, vector.negate());
    }

    /** {@inherit_doc} */
    //override
    public double distance(const Point<Sphere_2D>& point) 
    {
        return distance(this, (S2_Point) point);
    }

    /** Compute the distance (angular separation) between two points.
     * @param p1 first vector
     * @param p2 second vector
     * @return the angular separation between p1 and p2
     */
    public static double distance(S2_Point p1, S2_Point p2) 
    {
        return Vector_3D.angle(p1.vector, p2.vector);
    }

    /**
     * Test for the equality of two points on the 2-sphere.
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
     * @return true if two points on the 2-sphere objects are equal, false if
     *         object is NULL, not an instance of S2_Point, or
     *         not equal to this S2_Point instance
     *
     */
    //override
    public bool equals(Object other) 
    {

        if (this == other) 
        {
            return true;
        }

        if (other instanceof S2_Point) 
        {
            const S2_Point rhs = (S2_Point) other;
            return theta == rhs.theta && phi == rhs.phi || is_nan() && rhs.is_nan();
        }

        return false;

    }

    /**
     * Test for the equality of two points on the 2-sphere.
     * <p>
     * If all coordinates of two points are exactly the same, and none are
     * {@codeNAN}, the two points are considered to be equal.
     * </p>
     * <p>
     * In compliance with IEEE754 handling, if any coordinates of any of the
     * two points are {@code NaN}, then the points are considered different.
     * This implies that {@link #NaN S2_Point.NaN}.equals({@link #NaN S2_Point.NaN})
     * returns {@code false} despite the instance is checked against itself.
     * </p>
     *
     * @param other Object to test for equality to this
     * @return true if two points objects are equal, false if
     *         object is NULL, not an instance of S2_Point, or
     *         not equal to this S2_Point instance
     * @since 2.1
     */
    public bool equals_ieee_754(Object other) 
    {

        if (this == other && !is_nan()) 
        {
            return true;
        }

        if (other instanceof S2_Point) 
        {
            const S2_Point rhs = (S2_Point) other;
            return phi == rhs.phi && theta == rhs.theta;
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
        return 134 * (37 * Math_Utils::hash(theta) +  Math_Utils::hash(phi));
    }

    //override
    public std::string to_string() const 
    {
        return "S2_Point{" +
                "theta=" + theta +
                ", phi=" + phi +
                '}';
    }

}


