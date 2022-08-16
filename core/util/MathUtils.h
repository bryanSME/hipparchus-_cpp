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

//package org.hipparchus.util;

//import java.util.Arrays;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field_Element;
//import org.hipparchus.exception.Localizable;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.;
#include <numbers>
#include <vector>
#include <cmath>
#include "../CalculusFieldElement.hpp"
#include "../exception/MathRuntimeException.h"
#include "../exception/LocalizedCoreFormats.h"
#include <type_traits>
#include "../CalculusFieldElement.hpp"

/**
 * Miscellaneous utility functions.
 *
 * @see Arithmetic_Utils
 * @see Precision
 * @see Math_Arrays
 */
namespace Math_Utils 
{
    /** \(2\pi\) */
    inline static constexpr double TWO_PI{ 2 * std::numbers::pi };

    /** \(\pi^2\) */
    inline static constexpr double PI_SQUARED{ std::numbers::pi * std::numbers::pi };

    /** \(\pi/2\). */
    inline static constexpr double SEMI_PI{ 0.5 * std::numbers::pi };

    /**
     * Returns an integer hash code representing the given double value.
     *
     * @param value the value to be hashed
     * @return the hash code
     */
    inline static int hash(double value)
    {
        return Double.hash_code(value);
    }

    /**
     * Returns {@code true} if the values are equal according to semantics of
     * {@link Double#equals(Object)}.
     *
     * @param x Value
     * @param y Value
     * @return {@code Double(x).equals(new Double(y))}
     */
    inline static bool equals(const double& x, double y)
    {
        return Double(x).equals(new Double(y));
    }

    /**
     * Returns an integer hash code representing the given double array.
     *
     * @param value the value to be hashed (may be NULL)
     * @return the hash code
     */
    inline static int hash(std::vector<double> value)
    {
        return Arrays.hash_code(value);
    }

    /**
     * Normalize an angle in a 2&pi; wide interval around a center value.
     * <p>This method has three main uses:</p>
     * <ul>
     *   <li>normalize an angle between 0 and 2&pi;:<br/>
     *       {@code a = Math_Utils::normalize_angle(a, std::numbers::pi);}</li>
     *   <li>normalize an angle between -&pi; and +&pi;<br/>
     *       {@code a = Math_Utils::normalize_angle(a, 0.0);}</li>
     *   <li>compute the angle between two defining angular positions:<br>
     *       {@code angle = Math_Utils::normalize_angle(end, start) - start;}</li>
     * </ul>
     * <p>Note that due to numerical accuracy and since &pi; cannot be represented
     * exactly, the result interval is <em>closed</em>, it cannot be half-closed
     * as would be more satisfactory in a purely mathematical view.</p>
     * @param a angle to normalize
     * @param center center of the desired 2&pi; interval for the result
     * @return a-2k&pi; with integer k and center-&pi; &lt;= a-2k&pi; &lt;= center+&pi;
     */
    inline static double normalize_angle(const double& a, const double& center)
     {
         return a - TWO_PI * std::floor((a + std::numbers::pi - center) / TWO_PI);
     }

     /**
      * Normalize an angle in a 2&pi; wide interval around a center value.
      * <p>This method has three main uses:</p>
      * <ul>
      *   <li>normalize an angle between 0 and 2&pi;:<br/>
      *       {@code a = Math_Utils::normalize_angle(a, std::numbers::pi);}</li>
      *   <li>normalize an angle between -&pi; and +&pi;<br/>
      *       {@code a = Math_Utils::normalize_angle(a, zero);}</li>
      *   <li>compute the angle between two defining angular positions:<br>
      *       {@code angle = Math_Utils::normalize_angle(end, start).subtract(start);}</li>
      * </ul>
      * <p>Note that due to numerical accuracy and since &pi; cannot be represented
      * exactly, the result interval is <em>closed</em>, it cannot be half-closed
      * as would be more satisfactory in a purely mathematical view.</p>
      * @param <T> the type of the field elements
      * @param a angle to normalize
      * @param center center of the desired 2&pi; interval for the result
      * @return a-2k&pi; with integer k and center-&pi; &lt;= a-2k&pi; &lt;= center+&pi;
      */
     template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
     inline static T normalize_angle(T a, T center)
     {
         return a.subtract(std::floor(a.add(std::numbers::pi).subtract(center).divide(TWO_PI)).multiply(TWO_PI));
     }

     /** Find the maximum of two field elements.
      * @param <T> the type of the field elements
      * @param e1 first element
      * @param e2 second element
      * @return max(a1, e2)
      */
     template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
     inline static T max(const T& e1, const T& e2)
     {
         return e1.subtract(e2).get_real() >= 0 ? e1 : e2;
     }

     /** Find the minimum of two field elements.
      * @param <T> the type of the field elements
      * @param e1 first element
      * @param e2 second element
      * @return min(a1, e2)
      */
     template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
     inline static T min(const T& e1, const T& e2)
     {
         return e1.subtract(e2).get_real() >= 0 ? e2 : e1;
     }

    /**
     * <p>Reduce {@code |a - offset|} to the primary interval
     * {@code [0, |period|)}.</p>
     *
     * <p>Specifically, the value returned is <br/>
     * {@code a - |period| * floor((a - offset) / |period|) - offset}.</p>
     *
     * <p>If any of the parameters are {@code NaN} or infinite, the result is
     * {@code NaN}.</p>
     *
     * @param a Value to reduce.
     * @param period Period.
     * @param offset Value that will be mapped to {@code 0}.
     * @return the value, within the interval {@code [0 |period|)}, * that corresponds to {@code a}.
     */
    inline static double reduce(const double& a, const double& period, const double& offset)
    {
        const double p = std::abs(period);
        return a - p * std::floor((a - offset) / p) - offset;
    }

    /**
     * Returns the first argument with the sign of the second argument.
     *
     * @param magnitude Magnitude of the returned value.
     * @param sign Sign of the returned value.
     * @return a value with magnitude equal to {@code magnitude} and with the
     * same sign as the {@code sign} argument.
     * @Math_Runtime_Exception if {@code magnitude == Byte.MIN_VALUE}
     * and {@code sign >= 0}.
     */
    inline static std::byte copy_sign(const std::byte& magnitude, const std::byte& sign) 
    {
        if ((magnitude >= 0 && sign >= 0) || (magnitude < 0 && sign < 0)) // Sign is OK.
        {
            return magnitude;
        }
        if (sign >= 0 && magnitude == std::numeric_limits<std::byte>::min()) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW);
        }
        return -magnitude; // Flip sign.
    }

    /**
     * Returns the first argument with the sign of the second argument.
     *
     * @param magnitude Magnitude of the returned value.
     * @param sign Sign of the returned value.
     * @return a value with magnitude equal to {@code magnitude} and with the
     * same sign as the {@code sign} argument.
     * @Math_Runtime_Exception if {@code magnitude == Short.MIN_VALUE}
     * and {@code sign >= 0}.
     */
    inline static short copy_sign(const short& magnitude, const short& sign)
    {
        if ((magnitude >= 0 && sign >= 0) ||
            (magnitude < 0 && sign < 0)) { // Sign is OK.
            return magnitude;
        }
        if (sign >= 0 && magnitude == std::numeric_limits<short>::min()) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW);
        }
        return -magnitude; // Flip sign.
    }

    /**
     * Returns the first argument with the sign of the second argument.
     *
     * @param magnitude Magnitude of the returned value.
     * @param sign Sign of the returned value.
     * @return a value with magnitude equal to {@code magnitude} and with the
     * same sign as the {@code sign} argument.
     * @Math_Runtime_Exception if {@code magnitude == std::numeric_limits<int>::min()}
     * and {@code sign >= 0}.
     */
    inline static int copy_sign(const int& magnitude, const int& sign)
    {
        if ((magnitude >= 0 && sign >= 0) ||
            (magnitude < 0 && sign < 0)) { // Sign is OK.
            return magnitude;
        }
        if (sign >= 0 && magnitude == std::numeric_limits<short>::min())
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW);
        }
        return -magnitude; // Flip sign.
    }

    /**
     * Returns the first argument with the sign of the second argument.
     *
     * @param magnitude Magnitude of the returned value.
     * @param sign Sign of the returned value.
     * @return a value with magnitude equal to {@code magnitude} and with the
     * same sign as the {@code sign} argument.
     * @Math_Runtime_Exception if {@code magnitude == long.MIN_VALUE}
     * and {@code sign >= 0}.
     */
    inline static long copy_sign(const long& magnitude, const long& sign)
    {
        if ((magnitude >= 0 && sign >= 0) ||
            (magnitude < 0 && sign < 0)) { // Sign is OK.
            return magnitude;
        }
        if (sign >= 0 && magnitude == std::numeric_limits<short>::min())
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW);
        }

        return -magnitude; // Flip sign.
    }
    /**
     * Check that the argument is a real number.
     *
     * @param x Argument.
     * @ if {@code x} is not a
     * finite real number.
     */
    inline static void check_finite(const double& x)
    {
        if (std::isfinite(x) || std::isnan(x)) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_FINITE_NUMBER, x);
        }
    }

    /**
     * Check that all the elements are real numbers.
     *
     * @param val Arguments.
     * @ if any values of the array is not a
     * finite real number.
     */
    inline static void check_finite(const std::vector<double>& val)
    {
        for (const auto& ele : val)
        {
            if (std::isinfinite(ele) || std::isnan(ele))
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_FINITE_NUMBER, x);
            }
        }
    }

    /**
     * Checks that an object is not NULL.
     *
     * @param o Object to be checked.
     * @param pattern Message pattern.
     * @param args Arguments to replace the placeholders in {@code pattern}.
     * @ if {@code o} is {@code NULL}.
     */
    //inline static void check_not_null(Object o, Localizable pattern, Object ... args)
    //{
    //    if (o == NULL) 
    //    {
    //        throw (pattern, args);
    //    }
    //}

    /**
     * Checks that an object is not NULL.
     *
     * @param o Object to be checked.
     * @ if {@code o} is {@code NULL}.
     */
    //inline static void check_not_null(Object o)
    //{
    //    if (o == NULL) 
    //    {
    //        throw (hipparchus::exception::Localized_Core_Formats_Type::NULL_NOT_ALLOWED);
    //    }
    //}

    /**
     * Checks that the given value is strictly within the range [lo, hi].
     *
     * @param value value to be checked.
     * @param lo the lower bound (inclusive).
     * @param hi the upper bound (inclusive).
     * @ if {@code value} is strictly outside [lo, hi].
     */
    inline static void check_range_inclusive(const long& value, const long& lo, const long& hi)
    {
        if (value < lo || value > hi) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, value, lo, hi);
        }
    }

    /**
     * Checks that the given value is strictly within the range [lo, hi].
     *
     * @param value value to be checked.
     * @param lo the lower bound (inclusive).
     * @param hi the upper bound (inclusive).
     * @ if {@code value} is strictly outside [lo, hi].
     */
    inline static void check_range_inclusive(const double& value, const double& lo, const double& hi)
    {
        if (value < lo || value > hi) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, value, lo, hi);
        }
    }

    /**
     * Checks that the given dimensions match.
     *
     * @param dimension the first dimension.
     * @param other_dimension the second dimension.
     * @ if length != other_length.
     */
    inline static void check_dimension(const int& dimension, const int& other_dimension)
    {
        if (dimension != other_dimension) 
        {
            throw (Localized_Core_Formats_Types::DIMENSIONS_MISMATCH, dimension, other_dimension);
        }
    }

    /**
     * Result class for {@link Math_Utils#two_sum(double, double)} containing the
     * sum and the residual error in the sum.
     */
    class Sum_And_Residual
    {
    private:
        /** Sum. */
        const double my_sum;
        /** Residual error in the sum. */
        const double my_residual;

    public:
        /**
         * Constructs a {@link Sum_And_Residual} instance.
         * @param sum sum
         * @param residual residual error in the sum
         */
        Sum_And_Residual(const double& sum, const double& residual) : my_sum{ sum }, my_residual{ residual }{}

        /**
         * Returns the sum.
         * @return sum
         */
        double get_sum() const
        {
            return my_sum;
        }

        /**
         * Returns the residual error in the sum.
         * @return residual error in the sum
         */
        double get_residual() const
        {
            return my_residual;
        }
    };

    /**
     * Sums {@code a} and {@code b} using Møller's 2Sum algorithm.
     * <p>
     * References:
     * <ul>
     * <li>Møller, Ole. "Quasi double-precision in floating point addition." BIT
     * 5, 37–50 (1965).</li>
     * <li>Shewchuk, Richard J. "Adaptive Precision Floating-Point Arithmetic
     * and Fast Robust Geometric Predicates." Discrete & Computational Geometry
     * 18, 305–363 (1997).</li>
     * <li><a href=
     * "https://en.wikipedia.org/wiki/2Sum">https://en.wikipedia.org/wiki/2Sum</a></li>
     * </ul>
     * @param a first summand
     * @param b second summand
     * @return sum and residual error in the sum
     */
    inline static Sum_And_Residual two_sum(const double& a, const double& b)
    {
        const double s = a + b;
        const double a_prime = s - b;
        const double b_prime = s - a_prime;
        const double delta_a = a - a_prime;
        const double delta_b = b - b_prime;
        const double t = delta_a + delta_b;
        return Math_Utils::Sum_And_Residual(s, t);
        //return Sum_And_Residual(s, t);
    }

    /**
     * Sums {@code a} and {@code b} using Møller's 2Sum algorithm.
     * <p>
     * References:
     * <ul>
     * <li>Møller, Ole. "Quasi double-precision in floating point addition." BIT
     * 5, 37–50 (1965).</li>
     * <li>Shewchuk, Richard J. "Adaptive Precision Floating-Point Arithmetic
     * and Fast Robust Geometric Predicates." Discrete & Computational Geometry
     * 18, 305–363 (1997).</li>
     * <li><a href=
     * "https://en.wikipedia.org/wiki/2Sum">https://en.wikipedia.org/wiki/2Sum</a></li>
     * </ul>
     * @param <T> field element type
     * @param a first summand
     * @param b second summand
     * @return sum and residual error in the sum
     */
    /*inline static <T extends Field_Element<T>> Field_Sum_And_Residual<T> two_sum(const T& a, const T& b)
    {
        const T s = a.add(b);
        const T a_prime = s.subtract(b);
        const T& b_prime = s.subtract(a_prime);
        const T delta_a = a.subtract(a_prime);
        const T delta_b = b.subtract(b_prime);
        const T t = delta_a.add(delta_b);
        return Field_Sum_And_Residual<>(s, t);
    }*/



    ///**
    // * Result class for
    // * {@link Math_Utils#two_sum(Field_Element, Field_Element)} containing
    // * the sum and the residual error in the sum.
    // * @param <T> field element type
    // */
    //template<typename T, typename std::enable_if<std::is_base_of<Field_Element, T>::value>::type* = nullptr>
    //class Field_Sum_And_Residual
    //{
    //private:
    //    /** Sum. */
    //    const T my_sum;
    //    /** Residual error in the sum. */
    //    const T my_residual;

    //    /**
    //     * Constructs a {@link Field_Sum_And_Residual} instance.
    //     * @param sum sum
    //     * @param residual residual error in the sum
    //     */
    //    Field_Sum_And_Residual(const T& sum, const T& residual) : my_sum{ sum }, my_residual{ residual }{}

    //public:
    //    /**
    //     * Returns the sum.
    //     * @return sum
    //     */
    //    T get_sum() const
    //    {
    //        return my_sum;
    //    }

    //    /**
    //     * Returns the residual error in the sum.
    //     * @return residual error in the sum
    //     */
    //    T get_residual() const
    //    {
    //        return my_residual;
    //    }
    //};
};