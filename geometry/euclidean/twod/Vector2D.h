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
#include <algorithm>
#include <string>
#include <cmath>
#include <vector>
#include "Vector2DFormat.h"
#include "../../Space.h"
#include "../../../core/util/MathUtils.h"
#include "../../../core/util/MathArrays.h"

/** This class represents a 2D vector.
 * <p>Instances of this class are guaranteed to be immutable.</p>
 */
class Vector_2D : std::vector<Euclidean_2D> 
{
private:
    /** Abscissa. */
    double my_x;

    /** Ordinate. */
    double my_y;

public:
    /** Origin (coordinates: 0, 0). */
    static const Vector_2D ZERO = Vector_2D(0, 0);

    /** First canonical vector (coordinates: 1, 0).
     * @since 1.6
     */
    static const Vector_2D PLUS_I = Vector_2D(1, 0);

    /** Opposite of the first canonical vector (coordinates: -1, 0).
     * @since 1.6
     */
    static const Vector_2D MINUS_I = Vector_2D(-1, 0);

    /** Second canonical vector (coordinates: 0, 1).
     * @since 1.6
     */
    static const Vector_2D PLUS_J = Vector_2D(0, 1);

    /** Opposite of the second canonical vector (coordinates: 0, -1).
     * @since 1.6
     */
    static const Vector_2D MINUS_J = Vector_2D(0, -1);

    // CHECKSTYLE: stop Constant_Name
    /** A vector with all coordinates set to NaN. */
    static const Vector_2D NaN = Vector_2D(NAN,NAN);
    // CHECKSTYLE: resume Constant_Name

    /** A vector with all coordinates set to positive infinity. */
    static const Vector_2D POSITIVE_INFINITY = Vector_2D(INFINITY, INFINITY);

    /** A vector with all coordinates set to negative infinity. */
    static const Vector_2D NEGATIVE_INFINITY = Vector_2D(-INFINITY, -INFINITY);

    /** Simple constructor.
     * Build a vector from its coordinates
     * @param x abscissa
     * @param y ordinate
     * @see #get_x()
     * @see #get_y()
     */
    Vector_2D(const double& x, const double& y) : my_x{ x }, my_y{ y } {};

    /** Simple constructor.
     * Build a vector from its coordinates
     * @param v coordinates array
     * @exception  if array does not have 2 elements
     * @see #to_array()
     */
    Vector_2D(const std::vector<double>& v)
    {
        if (v.size() != 2) 
        {
            throw std::exception("Not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), 2);
        }
        my_x = v[0];
        my_y = v[1];
    }

    /** Multiplicative constructor
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    Vector_2D(const double& a, const Vector_2D& u) 
    {
        my_x = a * u.get_x();
        my_y = a * u.get_y();
    }

    /** Linear constructor
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    Vector_2D(const double& a1, const Vector_2D& u1, const double& a2, const Vector_2D& u2) 
    {
        my_x = a1 * u1.get_x() + a2 * u2.get_x();
        my_y = a1 * u1.get_y() + a2 * u2.get_y();
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
    Vector_2D(const double& a1, const Vector_2D& u1, const double& a2, const Vector_2D& u2, const double& a3, const Vector_2D& u3) 
    {
        my_x = a1 * u1.x + a2 * u2.x + a3 * u3.x;
        my_y = a1 * u1.y + a2 * u2.y + a3 * u3.y;
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
    Vector_2D(const double& a1, const Vector_2D& u1, const double& a2, const Vector_2D& u2, const double& a3, const Vector_2D& u3, const double& a4, Vector_2D u4) 
    {
        my_x = a1 * u1.get_x() + a2 * u2.get_x() + a3 * u3.get_x() + a4 * u4.get_x();
        my_y = a1 * u1.get_y() + a2 * u2.get_y() + a3 * u3.get_y() + a4 * u4.get_y();
    }

    /** Get the abscissa of the vector.
     * @return abscissa of the vector
     * @see #Vector_2D(double, double)
     */
    double get_x() const
    {
        return my_x;
    }

    /** Get the ordinate of the vector.
     * @return ordinate of the vector
     * @see #Vector_2D(double, double)
     */
    double get_y() const 
    {
        return my_y;
    }

    /** Get the vector coordinates as a dimension 2 array.
     * @return vector coordinates
     * @see #Vector_2D(std::vector<double>)
     */
    std::vector<double> to_array() const
    {
        return std::vector<double> { my_x, my_y };
    }

    /** {@inherit_doc} */
    //override
    Space get_space() 
    {
        return Euclidean_2D.get_instance();
    }

    /** {@inherit_doc} */
    //override
    Vector_2D get_zero() const
    {
        return ZERO;
    }

    /** {@inherit_doc} */
    //override
    double get_norm1() const
    {
        return std::abs(my_x) + std::abs(my_y);
    }

    /** {@inherit_doc} */
    //override
    double get_norm() const
    {
        return std::sqrt(my_x * my_x + my_y * my_y);
    }

    /** {@inherit_doc} */
    //override
    double get_norm_sq() const
    {
        return my_x * my_x + my_y * my_y;
    }

    /** {@inherit_doc} */
    //override
    double get_norm_inf() 
    {
        return std::max(std::abs(my_x), std::abs(my_y));
    }

    /** {@inherit_doc} */
    //override
    Vector_2D add(const std::vector<Euclidean_2D>& v) 
    {
        Vector_2D v2 = (Vector_2D) v;
        return Vector_2D(my_x + v2.get_x(), my_y + v2.get_y());
    }

    /** {@inherit_doc} */
    //override
    Vector_2D add(double factor, std::vector<Euclidean_2D> v) 
    {
        Vector_2D v2 = (Vector_2D) v;
        return Vector_2D(my_x + factor * v2.get_x(), my_y + factor * v2.get_y());
    }

    /** {@inherit_doc} */
    //override
    Vector_2D subtract(std::vector<Euclidean_2D> p) 
    {
        Vector_2D p3 = (Vector_2D) p;
        return Vector_2D(my_x - p3.get_x(), my_y - p3.get_y());
    }

    /** {@inherit_doc} */
    //override
    Vector_2D subtract(double factor, std::vector<Euclidean_2D> v) 
    {
        Vector_2D v2 = (Vector_2D) v;
        return Vector_2D(my_x - factor * v2.get_x(), my_y - factor * v2.get_y());
    }

    /** {@inherit_doc} */
    //override
    Vector_2D normalize() 
    {
        const double s = get_norm();
        if (s == 0) 
        {
            throw std::exception("Not implemented");
            //throw Math_Runtime_Exception(Localized_Geometry_Formats.CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR);
        }
        return scalar_multiply(1 / s);
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
    static double angle(const Vector_2D& v1, const Vector_2D& v2) 
    {

        double norm_product = v1.get_norm() * v2.get_norm();
        if (norm_product == 0) 
        {
            throw std::exception("not implemented");
            //throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }

        double dot = v1.dot_product(v2);
        double threshold = norm_product * 0.9999;
        if (std::abs(dot) > threshold) 
        {
            // the vectors are almost aligned, compute using the sine
            const double n = std::abs(Math_Arrays::linear_combination(v1.x, v2.y, -v1.y, v2.x));
            if (dot >= 0) 
            {
                return std::asin(n / norm_product);
            }
            return std::numbers::pi - std::asin(n / norm_product);
        }

        // the vectors are sufficiently separated to use the cosine
        return std::acos(dot / norm_product);

    }

    /** {@inherit_doc} */
    //override
    Vector_2D negate() 
    {
        return Vector_2D(-x, -y);
    }

    /** {@inherit_doc} */
    //override
    Vector_2D scalar_multiply(double a) 
    {
        return Vector_2D(a * x, a * y);
    }

    /** {@inherit_doc} */
    //override
    bool is_nan() 
    {
        return std::isnan(x) || std::isnan(y);
    }

    /** {@inherit_doc} */
    //override
    bool is_infinite() 
    {
        return !is_nan() && (std::isinf(x) || std::isinf(y));
    }

    /** {@inherit_doc} */
    //override
    double distance1(Vector<Euclidean_2D> p) 
    {
        Vector_2D p3 = (Vector_2D) p;
        const double dx = std::abs(p3.x - x);
        const double dy = std::abs(p3.y - y);
        return dx + dy;
    }

    /** {@inherit_doc} */
    //override
    double distance(Point<Euclidean_2D> p) 
    {
        Vector_2D p3 = (Vector_2D) p;
        const double dx = p3.x - x;
        const double dy = p3.y - y;
        return std::sqrt(dx * dx + dy * dy);
    }

    /** {@inherit_doc} */
    //override
    double distance_inf(Vector<Euclidean_2D> p) 
    {
        Vector_2D p3 = (Vector_2D) p;
        const double dx = std::abs(p3.x - x);
        const double dy = std::abs(p3.y - y);
        return std::max(dx, dy);
    }

    /** {@inherit_doc} */
    //override
    double distance_sq(std::vector<Euclidean_2D> p) 
    {
        Vector_2D p3 = (Vector_2D) p;
        const double dx = p3.x - x;
        const double dy = p3.y - y;
        return dx * dx + dy * dy;
    }

    /** {@inherit_doc} */
    //override
    double dot_product(const std::vector<Euclidean_2D>& v) 
    {
        const Vector_2D v2 = (Vector_2D) v;
        return Math_Arrays::linear_combination(my_x, v2.get_x(), my_y, v2.get_y());
    }

    /**
     * Compute the cross-product of the instance and the given points.
     * <p>
     * The cross product can be used to determine the location of a point
     * with regard to the line formed by (p1, p2) and is calculated as:
     * \[
     *    P = (x_2 - x_1)(y_3 - y_1) - (y_2 - y_1)(x_3 - x_1)
     * \]
     * with \(p3 = (x_3, y_3)\) being this instance.
     * <p>
     * If the result is 0, the points are collinear, i.e. lie on a single straight line L;
     * if it is positive, this point lies to the left, otherwise to the right of the line
     * formed by (p1, p2).
     *
     * @param p1 first point of the line
     * @param p2 second point of the line
     * @return the cross-product
     *
     * @see <a href="http://en.wikipedia.org/wiki/Cross_product">Cross product (Wikipedia)</a>
     */
    double cross_product(const Vector_2D& p1, const Vector_2D& p2) 
    {
        const double x1 = p2.get_x() - p1.get_x();
        const double y1 = get_y() - p1.get_y();
        const double x2 = get_x() - p1.get_x();
        const double y2 = p2.get_y() - p1.get_y();
        return Math_Arrays::linear_combination(x1, y1, -x2, y2);
    }

    /** Compute the distance between two vectors according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm1()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @return the distance between p1 and p2 according to the L<sub>1</sub> norm
     * @since 1.6
     */
    static double distance1(const Vector_2D& p1, const Vector_2D& p2) 
    {
        return p1.distance1(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    static double distance(const Vector_2D& p1, const Vector_2D& p2) 
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
    static double distance_inf(const Vector_2D& p1, const Vector_2D& p2) 
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
    static double distance_sq(const Vector_2D& p1, const Vector_2D& p2) 
    {
        return p1.distance_sq(p2);
    }

    /** Compute the orientation of a triplet of points.
     * @param p first vector of the triplet
     * @param q second vector of the triplet
     * @param r third vector of the triplet
     * @return a positive value if (p, q, r) defines a counterclockwise oriented
     * triangle, a negative value if (p, q, r) defines a clockwise oriented
     * triangle, and 0 if (p, q, r) are collinear or some points are equal
     * @since 1.2
     */
    static double orientation(const Vector_2D& p, const Vector_2D q, const Vector_2D& r) 
    {
        return Math_Arrays::linear_combination(std::vector<double> 
        {
            p.get_x(), -p.get_x(), q.get_x(), -q.get_x(), r.get_x(), -r.get_x()
        }, std::vector<double> 
        {
            q.get_y(),  r.get_y(), r.get_y(),  p.get_y(), p.get_y(),  q.get_y()
        });
    }

    /**
     * Test for the equality of two 2D vectors.
     * <p>
     * If all coordinates of two 2D vectors are exactly the same, and none are
     * {@codeNAN}, the two 2D vectors are considered to be equal.
     * </p>
     * <p>
     * {@code NaN} coordinates are considered to affect globally the vector
     * and be equals to each other - i.e, if either (or all) coordinates of the
     * 2D vector are equal to {@codeNAN}, the 2D vector is equal to
     * {@link #NaN}.
     * </p>
     *
     * @param other Object to test for equality to this
     * @return true if two 2D vector objects are equal, false if
     *         object is NULL, not an instance of Vector_2D, or
     *         not equal to this Vector_2D instance
     */
    //override
    bool equals(const Vector_2D& other)
    {

        if (*this == other) 
        {
            return true;
        }

        if (other instanceof Vector_2D) 
        {
            return x == other.get_x() && my_y == other.get_y() || is_nan() && other.is_nan();
        }

        return false;

    }

    /**
     * Test for the equality of two 2D vectors.
     * <p>
     * If all coordinates of two 2D vectors are exactly the same, and none are
     * {@code NaN}, the two 2D vectors are considered to be equal.
     * </p>
     * <p>
     * In compliance with IEEE754 handling, if any coordinates of any of the
     * two vectors are {@code NaN}, then the vectors are considered different.
     * This implies that {@link #NaN Vector_2D.NaN}.equals({@link #NaN Vector_2D.NaN})
     * returns {@code false} despite the instance is checked against itself.
     * </p>
     *
     * @param other Object to test for equality to this
     * @return true if two 2D vector objects are equal, false if
     *         object is NULL, not an instance of Vector_2D, or
     *         not equal to this Vector_2D instance
     * @since 2.1
     */
    bool equals_ieee_754(const Vector_2D& other)
    {
        if (*this == other && !is_nan()) 
        {
            return true;
        }

        if (other instanceof Vector_2D) 
        {
            return my_x == other.get_x() && my_y == other.get_y();
        }
        return false;
    }

    /**
     * Get a hash_code for the 2D vector.
     * <p>
     * All NaN values have the same hash code.</p>
     *
     * @return a hash code value for this object
     */
    //override
    int hash_code() 
    {
        if (is_nan()) 
        {
            return 542;
        }
        return 122 * (76 * Math_Utils::hash(my_x) +  Math_Utils::hash(my_y));
    }

    /** Get a string representation of this vector.
     * @return a string representation of this vector
     */
    //override
    std::string to_string() const 
    {
        return Vector_2D_Format.get_vector_2d_format().format(*this);
    }

    /** {@inherit_doc} */
    //override
    std::string to_string(const Number_Format& format) 
    {
        return Vector_2D_Format(format).format(*this);
    }
};