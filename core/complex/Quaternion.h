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

#include <sstream>
#include <cmath>
#include <vector>
#include "../util/MathUtils.h"
#include "../util/Precision.h"

//import java.io.Serializable;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Precision;

/**
 * This class : <a href="http://mathworld.wolfram.com/Quaternion.html">
 * quaternions</a> (Hamilton's hypercomplex numbers).
 * <p>
 * Instance of this class are guaranteed to be immutable.
 */
class Quaternion
{
public:
    /** Identity quaternion. */
    static const Quaternion IDENTITY() { return Quaternion(1, 0, 0, 0); };
    /** Zero quaternion. */
    static const Quaternion ZERO() { return Quaternion(0, 0, 0, 0); };
    /** i */
    static const Quaternion I() { return Quaternion(0, 1, 0, 0); };
    /** j */
    static const Quaternion J() { return Quaternion(0, 0, 1, 0); };
    /** k */
    static const Quaternion K() { return Quaternion(0, 0, 0, 1); };

private:

    /** First component (scalar part). */
    double my_q0;
    /** Second component (first vector part). */
    double my_q1;
    /** Third component (second vector part). */
    double my_q2;
    /** Fourth component (third vector part). */
    double my_q3;

public:
    /**
     * Builds a quaternion from its components.
     *
     * @param a Scalar component.
     * @param b First vector component.
     * @param c Second vector component.
     * @param d Third vector component.
     */
    Quaternion(const double& a, const double& b, const double& c, const double& d)
        :
        my_q0{ a }, my_q1{ b }, my_q2{ c }, my_q3{ d }
    {
    };

    /**
     * Builds a quaternion from scalar and vector parts.
     *
     * @param scalar Scalar part of the quaternion.
     * @param v Components of the vector part of the quaternion.
     *
     * @ if the array length is not 3.
     */
    Quaternion(const double& scalar, const std::vector<double>& v)
    {
        Math_Utils::check_dimension(v.size(), 3);
        my_q0 = scalar;
        my_q1 = v[0];
        my_q2 = v[1];
        my_q3 = v[2];
    }

    /**
     * Builds a pure quaternion from a vector (assuming that the scalar
     * part is zero).
     *
     * @param v Components of the vector part of the pure quaternion.
     */
    Quaternion(const std::vector<double>& v) 
    {
        Quaternion(0, v);
    }

    /**
     * Returns the conjugate quaternion of the instance.
     *
     * @return the conjugate quaternion
     */
    Quaternion get_conjugate() 
    {
        return Quaternion(my_q0, -my_q1, -my_q2, -my_q3);
    }

    /**
     * Returns the Hamilton product of two quaternions.
     *
     * @param q1 First quaternion.
     * @param q2 Second quaternion.
     * @return the product {@code q1} and {@code q2}, in that order.
     */
    static Quaternion multiply(const Quaternion& q1, const Quaternion& q2) 
    {
        // Components of the first quaternion.
        const double q1a = q1.get_q0();
        const double q1b = q1.get_q1();
        const double q1c = q1.get_q2();
        const double q1.0= q1.get_q3();

        // Components of the second quaternion.
        const double q2a = q2.get_q0();
        const double q2b = q2.get_q1();
        const double q2c = q2.get_q2();
        const double q2d = q2.get_q3();

        // Components of the product.
        const double w = q1a * q2a - q1b * q2b - q1c * q2c - q1.0* q2d;
        const double x = q1a * q2b + q1b * q2a + q1c * q2d - q1.0* q2c;
        const double y = q1a * q2c - q1b * q2d + q1c * q2a + q1.0* q2b;
        const double z = q1a * q2d + q1b * q2c - q1c * q2b + q1.0* q2a;

        return Quaternion(w, x, y, z);
    }

    /**
     * Returns the Hamilton product of the instance by a quaternion.
     *
     * @param q Quaternion.
     * @return the product of this instance with {@code q}, in that order.
     */
    Quaternion multiply(const Quaternion& q) 
    {
        return multiply(*this, q);
    }

    /**
     * Computes the sum of two quaternions.
     *
     * @param q1 Quaternion.
     * @param q2 Quaternion.
     * @return the sum of {@code q1} and {@code q2}.
     */
    static Quaternion add(const Quaternion& q1, const Quaternion& q2) 
    {
        return Quaternion(q1.get_q0() + q2.get_q0(), q1.get_q1() + q2.get_q1(), q1.get_q2() + q2.get_q2(), q1.get_q3() + q2.get_q3());
    }

    /**
     * Computes the sum of the instance and another quaternion.
     *
     * @param q Quaternion.
     * @return the sum of this instance and {@code q}
     */
    Quaternion add(const Quaternion& q) 
    {
        return add(*this, q);
    }

    /**
     * Subtracts two quaternions.
     *
     * @param q1 First Quaternion.
     * @param q2 Second quaternion.
     * @return the difference between {@code q1} and {@code q2}.
     */
    static Quaternion subtract(const Quaternion& q1, const Quaternion& q2) 
    {
        return Quaternion(q1.get_q0() - q2.get_q0(), q1.get_q1() - q2.get_q1(), q1.get_q2() - q2.get_q2(), q1.get_q3() - q2.get_q3());
    }

    /**
     * Subtracts a quaternion from the instance.
     *
     * @param q Quaternion.
     * @return the difference between this instance and {@code q}.
     */
    Quaternion subtract(const Quaternion& q) 
    {
        return subtract(*this, q);
    }

    /**
     * Computes the dot-product of two quaternions.
     *
     * @param q1 Quaternion.
     * @param q2 Quaternion.
     * @return the dot product of {@code q1} and {@code q2}.
     */
    static double dot_product(const Quaternion& q1, const Quaternion& q2) 
    {
        return q1.get_q0() * q2.get_q0() +
            q1.get_q1() * q2.get_q1() +
            q1.get_q2() * q2.get_q2() +
            q1.get_q3() * q2.get_q3();
    }

    /**
     * Computes the dot-product of the instance by a quaternion.
     *
     * @param q Quaternion.
     * @return the dot product of this instance and {@code q}.
     */
    double dot_product(const Quaternion& q) 
    {
        return dot_product(*this, q);
    }

    /**
     * Computes the norm of the quaternion.
     *
     * @return the norm.
     */
    double get_norm() const
    {
        return std::sqrt(my_q0 * my_q0 +
                             my_q1 * my_q1 +
                             my_q2 * my_q2 +
                             my_q3 * my_q3);
    }

    /**
     * Computes the normalized quaternion (the versor of the instance).
     * The norm of the quaternion must not be zero.
     *
     * @return a normalized quaternion.
     * @ if the norm of the quaternion is zero.
     */
    Quaternion normalize() 
    {
        const double norm = get_norm();

        if (norm < Precision::SAFE_MIN) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NORM, norm);
        }

        return Quaternion(my_q0 / norm, my_q1 / norm, my_q2 / norm, my_q3 / norm);
    }

    /**
     * {@inherit_doc}
     */
    //override
    bool equals(Object other) 
    {
        if (this == other) 
        {
            return true;
        }
        if (dynamic_cast<const Quaternion*>(*other) != nullptr) 
        {
            const Quaternion q = (Quaternion) other;
            return my_q0 == q.get_q0() &&
                my_q1 == q.get_q1() &&
                my_q2 == q.get_q2() &&
                my_q3 == q.get_q3();
        }

        return false;
    }

    /**
     * {@inherit_doc}
     */
    //override
    int hash_code() 
    {
        // "Effective Java" (second edition, p. 47).
        int result = 17;
        for (double comp : std::vector<double> { my_q0, my_q1, my_q2, my_q3 }) 
        {
            const int c = Math_Utils::hash(comp);
            result = 31 * result + c;
        }
        return result;
    }

    /**
     * Checks whether this instance is equal to another quaternion
     * within a given tolerance.
     *
     * @param q Quaternion with which to compare the current quaternion.
     * @param eps Tolerance.
     * @return {@code true} if the each of the components are equal
     * within the allowed absolute error.
     */
    bool equals(const Quaternion& q, const double& eps) 
    {
        return Precision::equals(my_q0, q.get_q0(), eps) &&
            Precision::equals(my_q1, q.get_q1(), eps) &&
            Precision::equals(my_q2, q.get_q2(), eps) &&
            Precision::equals(my_q3, q.get_q3(), eps);
    }

    /**
     * Checks whether the instance is a unit quaternion within a given
     * tolerance.
     *
     * @param eps Tolerance (absolute error).
     * @return {@code true} if the norm is 1 within the given tolerance, * {@code false} otherwise
     */
    bool is_unit_quaternion(const double& eps) 
    {
        return Precision::equals(get_norm(), 1, eps);
    }

    /**
     * Checks whether the instance is a pure quaternion within a given
     * tolerance.
     *
     * @param eps Tolerance (absolute error).
     * @return {@code true} if the scalar part of the quaternion is zero.
     */
    bool is_pure_quaternion(const double& eps) const
    {
        return std::abs(get_q0()) <= eps;
    }

    /**
     * Returns the polar form of the quaternion.
     *
     * @return the unit quaternion with positive scalar part.
     */
    Quaternion get_positive_polar_form() 
    {
        if (get_q0() < 0) 
        {
            const Quaternion unit_q = normalize();
            // The quaternion of rotation (normalized quaternion) q and -q
            // are equivalent (i.e. represent the same rotation).
            return Quaternion(-unit_q.get_q0(), -unit_q.get_q1(), -unit_q.get_q2(), -unit_q.get_q3());
        }
        return *this.normalize();
    }

    /**
     * Returns the inverse of this instance.
     * The norm of the quaternion must not be zero.
     *
     * @return the inverse.
     * @ if the norm (squared) of the quaternion is zero.
     */
    Quaternion get_inverse() 
    {
        const double square_norm = my_q0 * my_q0 + my_q1 * my_q1 + my_q2 * my_q2 + my_q3 * my_q3;
        if (square_norm < Precision::SAFE_MIN) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NORM, square_norm);
        }

        return Quaternion(my_q0 / square_norm, -my_q1 / square_norm, -my_q2 / square_norm, -my_q3 / square_norm);
    }

    /**
     * Gets the first component of the quaternion (scalar part).
     *
     * @return the scalar part.
     */
    double get_q0() const
    {
        return my_q0;
    }

    /**
     * Gets the second component of the quaternion (first component
     * of the vector part).
     *
     * @return the first component of the vector part.
     */
    double get_q1() const
    {
        return my_q1;
    }

    /**
     * Gets the third component of the quaternion (second component
     * of the vector part).
     *
     * @return the second component of the vector part.
     */
    double get_q2() const
    {
        return my_q2;
    }

    /**
     * Gets the fourth component of the quaternion (third component
     * of the vector part).
     *
     * @return the third component of the vector part.
     */
    double get_q3() const
    {
        return my_q3;
    }

    /**
     * Gets the scalar part of the quaternion.
     *
     * @return the scalar part.
     * @see #get_q0()
     */
    double get_scalar_part() const
    {
        return get_q0();
    }

    /**
     * Gets the three components of the vector part of the quaternion.
     *
     * @return the vector part.
     * @see #get_q1()
     * @see #get_q2()
     * @see #get_q3()
     */
    std::vector<double> get_vector_part() const
    {
        return std::vector<double> { get_q1(), get_q2(), get_q3() };
    }

    /**
     * Multiplies the instance by a scalar.
     *
     * @param alpha Scalar factor.
     * @return a scaled quaternion.
     */
    Quaternion multiply(const double& alpha) 
    {
        return Quaternion(alpha * my_q0, alpha * my_q1, alpha * my_q2, alpha * my_q3);
    }

    /**
     * {@inherit_doc}
     */
    //override
    std::string to_string() const 
    {
        std::stringstream ss;
        ss
            << "["
            << my_q0 << " "
            << my_q1 << " "
            << my_q2 << " "
            << my_q3
            << "]";
        return ss.str();
    }
};