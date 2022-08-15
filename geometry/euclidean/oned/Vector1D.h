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
//package org.hipparchus.geometry.euclidean.oned;

//import java.text.Number_Format;

//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.geometry.Localized_Geometry_Formats;
//import org.hipparchus.geometry.Point;
//import org.hipparchus.geometry.Space;
//import org.hipparchus.geometry.Vector;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/** This class represents a 1D vector.
 * <p>Instances of this class are guaranteed to be immutable.</p>
 */
class Vector_1D : public std::vector<Euclidean_1D> 
{

    /** Origin (coordinates: 0). */
    public static const Vector_1D ZERO = Vector_1D(0.0);

    /** Unit (coordinates: 1). */
    public static const Vector_1D ONE  = Vector_1D(1.0);

    // CHECKSTYLE: stop Constant_Name
    /** A vector with all coordinates set to NaN. */
    public static const Vector_1D NaN = Vector_1D(Double.NaN);
    // CHECKSTYLE: resume Constant_Name

    /** A vector with all coordinates set to positive infinity. */
    public static const Vector_1D POSITIVE_INFINITY =
        Vector_1D(INFINITY);

    /** A vector with all coordinates set to negative infinity. */
    public static const Vector_1D NEGATIVE_INFINITY =
        Vector_1D(-INFINITY);

    
    7556674948671647925L;

    /** Abscissa. */
    private const double x;

    /** Simple constructor.
     * Build a vector from its coordinates
     * @param x abscissa
     * @see #get_x()
     */
    public Vector_1D(double x) 
    {
        this.x = x;
    }

    /** Multiplicative constructor
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    public Vector_1D(const double& a, Vector_1D u) 
    {
        this.x = a * u.x;
    }

    /** Linear constructor
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    public Vector_1D(double a1, Vector_1D u1, double a2, Vector_1D u2) 
    {
        this.x = a1 * u1.x + a2 * u2.x;
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
    public Vector_1D(double a1, Vector_1D u1, double a2, Vector_1D u2, double a3, Vector_1D u3) 
    {
        this.x = a1 * u1.x + a2 * u2.x + a3 * u3.x;
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
    public Vector_1D(double a1, Vector_1D u1, double a2, Vector_1D u2, double a3, Vector_1D u3, double a4, Vector_1D u4) 
    {
        this.x = a1 * u1.x + a2 * u2.x + a3 * u3.x + a4 * u4.x;
    }

    /** Get the abscissa of the vector.
     * @return abscissa of the vector
     * @see #Vector_1Dstatic_cast<double>(
     */
    public double get_x() 
    {
        return x;
    }

    /** {@inherit_doc} */
    //override
    public Space get_space() 
    {
        return Euclidean_1D.get_instance();
    }

    /** {@inherit_doc} */
    //override
    public Vector_1D get_zero() 
    {
        return ZERO;
    }

    /** {@inherit_doc} */
    //override
    public double get_norm1() 
    {
        return std::abs(x);
    }

    /** {@inherit_doc} */
    //override
    public double get_norm() 
    {
        return std::abs(x);
    }

    /** {@inherit_doc} */
    //override
    public double get_norm_sq() 
    {
        return x * x;
    }

    /** {@inherit_doc} */
    //override
    public double get_norm_inf() 
    {
        return std::abs(x);
    }

    /** {@inherit_doc} */
    //override
    public Vector_1D add(Vector<Euclidean_1D> v) 
    {
        Vector_1D v1 = (Vector_1D) v;
        return Vector_1D(x + v1.get_x());
    }

    /** {@inherit_doc} */
    //override
    public Vector_1D add(double factor, Vector<Euclidean_1D> v) 
    {
        Vector_1D v1 = (Vector_1D) v;
        return Vector_1D(x + factor * v1.get_x());
    }

    /** {@inherit_doc} */
    //override
    public Vector_1D subtract(Vector<Euclidean_1D> p) 
    {
        Vector_1D p3 = (Vector_1D) p;
        return Vector_1D(x - p3.x);
    }

    /** {@inherit_doc} */
    //override
    public Vector_1D subtract(double factor, Vector<Euclidean_1D> v) 
    {
        Vector_1D v1 = (Vector_1D) v;
        return Vector_1D(x - factor * v1.get_x());
    }

    /** {@inherit_doc} */
    //override
    public Vector_1D normalize() Math_Runtime_Exception 
    {
        double s = get_norm();
        if (s == 0) 
        {
            throw Math_Runtime_Exception(Localized_Geometry_Formats.CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR);
        }
        return scalar_multiply(1 / s);
    }
    /** {@inherit_doc} */
    //override
    public Vector_1D negate() 
    {
        return Vector_1D(-x);
    }

    /** {@inherit_doc} */
    //override
    public Vector_1D scalar_multiply(double a) 
    {
        return Vector_1D(a * x);
    }

    /** {@inherit_doc} */
    //override
    public bool is_nan() 
    {
        return std::isnan(x);
    }

    /** {@inherit_doc} */
    //override
    public bool is_infinite() 
    {
        return !is_nan() && std::isinf(x);
    }

    /** {@inherit_doc} */
    //override
    public double distance1(Vector<Euclidean_1D> p) 
    {
        Vector_1D p3 = (Vector_1D) p;
        return std::abs(p3.x - x);
    }

    /** {@inherit_doc} */
    //override
    public double distance(Point<Euclidean_1D> p) 
    {
        Vector_1D p3 = (Vector_1D) p;
        return std::abs(p3.x - x);
    }

    /** {@inherit_doc} */
    //override
    public double distance_inf(Vector<Euclidean_1D> p) 
    {
        Vector_1D p3 = (Vector_1D) p;
        return std::abs(p3.x - x);
    }

    /** {@inherit_doc} */
    //override
    public double distance_sq(Vector<Euclidean_1D> p) 
    {
        Vector_1D p3 = (Vector_1D) p;
        const double dx = p3.x - x;
        return dx * dx;
    }

    /** {@inherit_doc} */
    //override
    public double dot_product(const Vector<Euclidean_1D> v) 
    {
        const Vector_1D v1 = (Vector_1D) v;
        return x * v1.x;
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    public static double distance(Vector_1D p1, Vector_1D p2) 
    {
        return p1.distance(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @return the distance between p1 and p2 according to the L<sub>&infin;</sub> norm
     */
    public static double distance_inf(Vector_1D p1, Vector_1D p2) 
    {
        return p1.distance_inf(p2);
    }

    /** Compute the square of the distance between two vectors.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @return the square of the distance between p1 and p2
     */
    public static double distance_sq(Vector_1D p1, Vector_1D p2) 
    {
        return p1.distance_sq(p2);
    }

    /**
     * Test for the equality of two 1D vectors.
     * <p>
     * If all coordinates of two 1D vectors are exactly the same, and none are
     * {@codeNAN}, the two 1D vectors are considered to be equal.
     * </p>
     * <p>
     * {@code NaN} coordinates are considered to affect globally the vector
     * and be equals to each other - i.e, if either (or all) coordinates of the
     * 1D vector are equal to {@codeNAN}, the 1D vector is equal to
     * {@link #NaN}.
     * </p>
     *
     * @param other Object to test for equality to this
     * @return true if two 1D vector objects are equal, false if
     *         object is NULL, not an instance of Vector_1D, or
     *         not equal to this Vector_1D instance
     */
    //override
    public bool equals(Object other) 
    {

        if (this == other) 
        {
            return true;
        }

        if (dynamic_cast<const Vector_1D*>(*other) != nullptr)
        {
            const auto rhs = (Vector_1D) other;
            return x == rhs.x || is_nan() && rhs.is_nan();
        }

        return false;

    }

    /**
     * Test for the equality of two 1D vectors.
     * <p>
     * If all coordinates of two 1D vectors are exactly the same, and none are
     * {@code NaN}, the two 1D vectors are considered to be equal.
     * </p>
     * <p>
     * In compliance with IEEE754 handling, if any coordinates of any of the
     * two vectors are {@code NaN}, then the vectors are considered different.
     * This implies that {@link #NaN Vector_1D.NaN}.equals({@link #NaN Vector_1D.NaN})
     * returns {@code false} despite the instance is checked against itself.
     * </p>
     *
     * @param other Object to test for equality to this
     * @return true if two 1D vector objects are equal, false if
     *         object is NULL, not an instance of Vector_1D, or
     *         not equal to this Vector_1D instance
     *
     * @since 2.1
     */
    public bool equals_ieee_754(Object other) 
    {

        if (this == other && !is_nan()) 
        {
            return true;
        }

        if (dynamic_cast<const Vector_1D*>(*other) != nullptr)
        {
            const Vector_1D rhs = (Vector_1D) other;
            return x == rhs.x;
        }

        return false;

    }

    /**
     * Get a hash_code for the 1D vector.
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
            return 7785;
        }
        return 997 * Math_Utils::hash(x);
    }

    /** Get a string representation of this vector.
     * @return a string representation of this vector
     */
    //override
    public std::string to_string() const 
    {
        return Vector_1D_Format.get_vector_1d_format().format(this);
    }

    /** {@inherit_doc} */
    //override
    public std::string to_string(const Number_Format format) 
    {
        return Vector_1D_Format(format).format(this);
    }

}


