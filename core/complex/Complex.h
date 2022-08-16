#pragma once
///*
// * Licensed to the Apache Software Foundation (ASF) under one or more
// * contributor license agreements.  See the NOTICE file distributed with
// * this work for additional information regarding copyright ownership.
// * The ASF licenses this file to You under the Apache License, Version 2.0
// * (the "License"); you may not use this file except in compliance with
// * the License.  You may obtain a copy of the License at
// *
// *      http://www.apache.org/licenses/LICENSE-2.0
// *
// * Unless required by applicable law or agreed to in writing, software
// * distributed under the License is distributed on an "AS IS" BASIS, * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// * See the License for the specific language governing permissions and
// * limitations under the License.
// */
//
///*
// * This is not the original file distributed by the Apache Software Foundation
// * It has been modified by the Hipparchus project
// */
//
//#include <cmath>
//#include <vector>
//#include <string>
//#include <numbers>
//#include "../../core/util/MathUtils.h"
//
////import java.io.Serializable;
////import java.util.Array_list;
////import java.util.List;
//
////import org.hipparchus.Calculus_Field_Element;
////import org.hipparchus.exception.Localized_Core_Formats;
////import org.hipparchus.exception.;
////import org.hipparchus.exception.;
////import org.hipparchus.util.FastMath;
////import org.hipparchus.util.Field_Sin_Cos;
////import org.hipparchus.util.Field_Sinh_Cosh;
////import org.hipparchus.util.Math_Arrays;
////import org.hipparchus.util.Math_Utils;
////import org.hipparchus.util.Precision;
////import org.hipparchus.util.Sin_Cos;
////import org.hipparchus.util.Sinh_Cosh;
//
//
///**
// * Representation of a std::complex<double> number, i.e. a number which has both a
// * my_my_real and my_imaginary part.
// * <p>
// * Implementations of arithmetic operations handle {@code NaN} and
// * infinite values according to the rules for {@link java.lang.Double}, i.e.
// * {@link #equals} is an equivalence relation for all instances that have
// * a {@code NaN} in either my_my_real or my_imaginary part, e.g. the following are
// * considered equal:
// * <ul>
// *  <li>{@code 1 + NaNi}</li>
// *  <li>{@code NaN + i}</li>
// *  <li>{@code NaN + NaNi}</li>
// * </ul>
// * <p>
// * Note that this contradicts the IEEE-754 standard for floating
// * point numbers (according to which the test {@code x == x} must fail if
// * {@code x} is {@code NaN}). The method
// * {@link org.hipparchus.util.Precision#equals(double,double,int)
// * equals for primitive double} in {@link org.hipparchus.util.Precision}
// * conforms with IEEE-754 while this class conforms with the standard behavior
// * for Java object types.
// */
//class std::complex<double> : Calculus_Field_Element<std::complex<double>>
//{
//private:
//    /** A my_my_real number representing log(10). */
//    static constexpr double LOG10{ 2.302585092994045684 };
//
//    /** The my_imaginary part. */
//    const double my_imaginary;
//    /** The my_my_real part. */
//    const double my_real;
//    /** Record whether this complex number is equal to NaN. */
//    const bool my_is_nan;
//    /** Record whether this complex number is infinite. */
//    const bool my_isinfinite;
//
//protected:
//    /**
//     * Create a complex number given the my_my_real and my_imaginary parts.
//     *
//     * @param my_my_real_part Real part.
//     * @param my_imaginary_part Imaginary part.
//     * @return a complex number instance.
//     *
//     * @see #value_of(double, double)
//     */
//    std::complex<double> create_complex(const double& real_part, const double& imaginary_part)
//    {
//        return std::complex<double>(real_part, imaginary_part);
//    }
//
//    /**
//     * Resolve the transient fields in a deserialized std::complex<double> Object.
//     * Subclasses will need to //override {@link #create_complex} to
//     * deserialize properly.
//     *
//     * @return A std::complex<double> instance with all fields resolved.
//     */
//    const Object read_resolve()
//    {
//        return create_complex(my_real, my_imaginary);
//    }
//
//public:
//    /** The square root of -1. A number representing "0.0 + 1.0i". */
//    static const auto I = std::complex<double>(0.0, 1.0);
//    /** The square root of -1. A number representing "0.0 - 1.0i".
//     * @since 1.7
//     */
//    static const auto MINUS_I = std::complex<double>(0.0, -1.0);
//    // CHECKSTYLE: stop Constant_Name
//    /** A complex number representing "NaN + NaNi". */
//    static const auto NaN = std::complex<double>(Double.NaN,NAN);
//    // CHECKSTYLE: resume Constant_Name
//    /** A complex number representing "+INF + INFi" */
//    static const auto INF = std::complex<double>(INFINITY, INFINITY);
//    /** A complex number representing "1.0 + 0.0i". */
//    static const auto ONE = std::complex<double>(1.0, 0.0);
//    /** A complex number representing "-1.0 + 0.0i".
//     * @since 1.7
//     */
//    static const auto MINUS_ONE = std::complex<double>(-1.0, 0.0);
//    /** A complex number representing "0.0 + 0.0i". */
//    static const auto ZERO = std::complex<double>(0.0, 0.0);
//    /** A complex number representing "π + 0.0i". */
//    static const auto PI   = std::complex<double>(std::numbers::pi, 0.0);
//
//
//
//    /**
//     * Create a complex number given only the my_my_real part.
//     *
//     * @param my_my_real Real part.
//     */
//    std::complex<double>(double my_my_real)
//    {
//        this(my_real, 0.0);
//    }
//
//    /**
//     * Create a complex number given the my_my_real and my_imaginary parts.
//     *
//     * @param my_my_real Real part.
//     * @param my_imaginary Imaginary part.
//     */
//    std::complex<double>(double real, double imaginary)
//        :
//        my_real{ real },
//        my_imaginary{ imaginary },
//        my_is_nan{ std::isnan(real) || std::isnan(imaginary) },
//        my_isinfinite{ !my_is_nan && (Double.my_isinfinite(real) || Double.my_isinfinite(imaginary)) }
//    {};
//
//    /**
//     * Return the absolute value of this complex number.
//     * Returns {@code NaN} if either my_my_real or my_imaginary part is {@code NaN}
//     * and {@code INFINITY} if neither part is {@code NaN}, * but at least one part is infinite.
//     *
//     * @return the norm.
//     * @since 2.0
//     */
//    //override
//    std::complex<double> abs()
//    {
//        // we check NaN here because std::hypot checks it after infinity
//        return my_is_nan ? NaN : create_complex(std::hypot(my_real, my_imaginary), 0.0);
//    }
//
//    /** {@inherit_doc} */
//    //override
//    double norm()
//    {
//        // we check NaN here because std::hypot checks it after infinity
//        return my_is_nan ?NAN : std::hypot(my_real, my_imaginary);
//    }
//
//    /**
//     * Returns a {@code std::complex<double>} whose value is
//     * {@code (this + addend)}.
//     * Uses the definitional formula
//     * <p>
//     *   {@code (a + bi) + (c + di) = (a+c) + (b+d)i}
//     * </p>
//     * If either {@code this} or {@code addend} has a {@code NaN} value in
//     * either part, {@link #NaN} is returned; otherwise {@code Infinite}
//     * and {@code NaN} values are returned in the parts of the result
//     * according to the rules for {@link java.lang.Double} arithmetic.
//     *
//     * @param  addend Value to be added to this {@code std::complex<double>}.
//     * @return {@code this + addend}.
//     * @ if {@code addend} is {@code NULL}.
//     */
//    //override
//    std::complex<double> add(const std::complex<double>& addend)
//    {
//        //Math_Utils::check_not_null(addend);
//        if (my_is_nan || addend.is_nan())
//
//        {
//            return NaN;
//        }
//
//        return create_complex(my_real + my_addend.get_real_part(), my_imaginary + addend.get_imaginary_part());
//    }
//
//    /**
//     * Returns a {@code std::complex<double>} whose value is {@code (this + addend)}, * with {@code addend} interpreted as a my_my_real number.
//     *
//     * @param addend Value to be added to this {@code std::complex<double>}.
//     * @return {@code this + addend}.
//     * @see #add(std::complex<double>)
//     */
//    //override
//    std::complex<double> add(const double& addend)
//    {
//        if (is_nan || std::isnan(addend))
//        {
//            return NaN;
//        }
//
//        return create_complex(my_real + addend, my_imaginary);
//    }
//
//     /**
//     * Returns the conjugate of this complex number.
//     * The conjugate of {@code a + bi} is {@code a - bi}.
//     * <p>
//     * {@link #NaN} is returned if either the my_my_real or my_imaginary
//     * part of this std::complex<double> number equals {@codeNAN}.
//     * </p><p>
//     * If the my_imaginary part is infinite, and the my_my_real part is not
//     * {@code NaN}, the returned value has infinite my_imaginary part
//     * of the opposite sign, e.g. the conjugate of
//     * {@code 1 + POSITIVE_INFINITY i} is {@code 1 - NEGATIVE_INFINITY i}.
//     * </p>
//     * @return the conjugate of this std::complex<double> object.
//     */
//    std::complex<double> conjugate() const
//    {
//        if (is_nan)
//        {
//            return nullptr;
//        }
//
//        return create_complex(my_real, -my_imaginary);
//    }
//
//    /**
//     * Returns a {@code std::complex<double>} whose value is
//     * {@code (this / divisor)}.
//     * Implements the definitional formula
//     * <pre>
//     *  <code>
//     *    a + bi          ac + bd + (bc - ad)i
//     *    ----------- = -------------------------
//     *    c + di         c<sup>2</sup> + d<sup>2</sup>
//     *  </code>
//     * </pre>
//     * but uses
//     * <a href="http://doi.acm.org/10.1145/1039813.1039814">
//     * prescaling of operands</a> to limit the effects of overflows and
//     * underflows in the computation.
//     * <p>
//     * {@code Infinite} and {@code NaN} values are handled according to the
//     * following rules, applied in the order presented:
//     * <ul>
//     *  <li>If either {@code this} or {@code divisor} has a {@code NaN} value
//     *   in either part, {@link #NaN} is returned.
//     *  </li>
//     *  <li>If {@code divisor} equals {@link #ZERO}, {@link #NaN} is returned.
//     *  </li>
//     *  <li>If {@code this} and {@code divisor} are both infinite, *   {@link #NaN} is returned.
//     *  </li>
//     *  <li>If {@code this} is finite (i.e., has no {@code Infinite} or
//     *   {@code NaN} parts) and {@code divisor} is infinite (one or both parts
//     *   infinite), {@link #ZERO} is returned.
//     *  </li>
//     *  <li>If {@code this} is infinite and {@code divisor} is finite, *   {@code NaN} values are returned in the parts of the result if the
//     *   {@link java.lang.Double} rules applied to the definitional formula
//     *   force {@code NaN} results.
//     *  </li>
//     * </ul>
//     *
//     * @param divisor Value by which this {@code std::complex<double>} is to be divided.
//     * @return {@code this / divisor}.
//     * @ if {@code divisor} is {@code NULL}.
//     */
//    //override
//    std::complex<double> divide(const std::complex<double>& divisor)
//    {
//        //Math_Utils::check_not_null(divisor);
//        if (is_nan || divisor.is_nan())
//        {
//            return NaN;
//        }
//
//        const double c = divisor.get_real_part();
//        const double d = divisor.get_imaginary_part();
//        if (c == 0.0 && d == 0.0)
//        {
//            return NaN;
//        }
//
//        if (divisor.my_isinfinite() && !my_isinfinite)
//        {
//            return ZERO;
//        }
//
//        if (std::abs(c) < std::abs(d))
//        {
//            double q = c / d;
//            double denominator = c * q + d;
//            return create_complex((my_real * q + my_imaginary) / denominator, (my_imaginary * q - my_real) / denominator);
//        }
//        double q = d / c;
//        double denominator = d * q + c;
//        return create_complex((my_imaginary * q + my_real) / denominator, (my_imaginary - my_real * q) / denominator);
//    }
//
//    /**
//     * Returns a {@code std::complex<double>} whose value is {@code (this / divisor)}, * with {@code divisor} interpreted as a my_my_real number.
//     *
//     * @param  divisor Value by which this {@code std::complex<double>} is to be divided.
//     * @return {@code this / divisor}.
//     * @see #divide(std::complex<double>)
//     */
//    //override
//    std::complex<double> divide(const double& divisor)
//    {
//        if (my_is_nan || std::isnan(divisor))
//        {
//            return NaN;
//        }
//        if (divisor == 0.0)
//        {
//            return NaN;
//        }
//        if (Double.my_isinfinite(divisor))
//        {
//            return !my_isinfinite
//                ? ZERO
//                : NaN;
//        }
//        return create_complex(my_real / divisor, my_imaginary  / divisor);
//    }
//
//    /** {@inherit_doc} */
//    //override
//    std::complex<double> reciprocal()
//    {
//        if (my_is_nan)
//        {
//            return NaN;
//        }
//
//        if (my_real == 0.0 && my_imaginary == 0.0)
//        {
//            return INF;
//        }
//
//        if (my_isinfinite)
//        {
//            return ZERO;
//        }
//
//        if (std::abs(my_real) < std::abs(my_imaginary))
//        {
//            double q = my_real / my_imaginary;
//            double scale = 1. / (my_real * q + my_imaginary);
//            return create_complex(scale * q, -scale);
//        }
//        double q = my_imaginary / my_real;
//        double scale = 1. / (my_imaginary * q + my_real);
//        return create_complex(scale, -scale * q);
//
//    }
//
//    /**
//     * Test for equality with another object.
//     * If both the my_my_real and my_imaginary parts of two complex numbers
//     * are exactly the same, and neither is {@codeNAN}, the two
//     * std::complex<double> objects are considered to be equal.
//     * The behavior is the same as for JDK's {@link Double#equals(Object)
//     * Double}:
//     * <ul>
//     *  <li>All {@code NaN} values are considered to be equal, *   i.e, if either (or both) my_my_real and my_imaginary parts of the complex
//     *   number are equal to {@codeNAN}, the complex number is equal
//     *   to {@code NaN}.
//     *  </li>
//     *  <li>
//     *   Instances constructed with different representations of zero (i.e.
//     *   either "0" or "-0") are <em>not</em> considered to be equal.
//     *  </li>
//     * </ul>
//     *
//     * @param other Object to test for equality with this instance.
//     * @return {@code true} if the objects are equal, {@code false} if object
//     * is {@code NULL}, not an instance of {@code std::complex<double>}, or not equal to
//     * this instance.
//     */
//    //override
//    bool equals(const Object& other)
//    {
//        if (this == other)
//        {
//            return true;
//        }
//		  if (dynamic_cast<const std::complex<double>*>(*other) != nullptr)
//        {
//            std::complex<double> c = (std::complex<double>) other;
//            if (c.is_nan())
//            {
//                return my_is_nan;
//            }
//            return Math_Utils::equals(my_real, c.get_real()) && Math_Utils::equals(my_imaginary, c.get_imaginary());
//        }
//        return false;
//    }
//
//    /**
//     * Test for the floating-point equality between std::complex<double> objects.
//     * It returns {@code true} if both arguments are equal or within the
//     * range of allowed error (inclusive).
//     *
//     * @param x First value (cannot be {@code NULL}).
//     * @param y Second value (cannot be {@code NULL}).
//     * @param max_ulps {@code (max_ulps - 1)} is the number of floating point
//     * values between the my_my_real (resp. my_imaginary) parts of {@code x} and
//     * {@code y}.
//     * @return {@code true} if there are fewer than {@code max_ulps} floating
//     * point values between the my_my_real (resp. my_imaginary) parts of {@code x}
//     * and {@code y}.
//     *
//     * @see Precision#equals(double,double,int)
//     */
//    static bool equals(const std::complex<double>& x, const std::complex<double>& y, const int& max_ulps)
//    {
//        return Precision::equals(x.get_real(), y.get_real(), max_ulps) &&
//               Precision::equals(x.get_imaginary(), y.get_imaginary(), max_ulps);
//    }
//
//    /**
//     * Returns {@code true} iff the values are equal as defined by
//     * {@link #equals(std::complex<double>,std::complex<double>,int) equals(x, y, 1)}.
//     *
//     * @param x First value (cannot be {@code NULL}).
//     * @param y Second value (cannot be {@code NULL}).
//     * @return {@code true} if the values are equal.
//     */
//    static bool equals(const std::complex<double>& x, const std::complex<double>& y)
//    {
//        return equals(x, y, 1);
//    }
//
//    /**
//     * Returns {@code true} if, both for the my_my_real part and for the my_imaginary
//     * part, there is no double value strictly between the arguments or the
//     * difference between them is within the range of allowed error
//     * (inclusive).  Returns {@code false} if either of the arguments is NaN.
//     *
//     * @param x First value (cannot be {@code NULL}).
//     * @param y Second value (cannot be {@code NULL}).
//     * @param eps Amount of allowed absolute error.
//     * @return {@code true} if the values are two adjacent floating point
//     * numbers or they are within range of each other.
//     *
//     * @see Precision#equals(double,double,double)
//     */
//    static bool equals(std::complex<double> x, std::complex<double> y, double eps)
//    {
//        return Precision::equals(x.get_real(), y.get_real(), eps) &&
//               Precision::equals(x.get_imaginary(), y.get_imaginary(), eps);
//    }
//
//    /**
//     * Returns {@code true} if, both for the my_my_real part and for the my_imaginary
//     * part, there is no double value strictly between the arguments or the
//     * relative difference between them is smaller or equal to the given
//     * tolerance. Returns {@code false} if either of the arguments is NaN.
//     *
//     * @param x First value (cannot be {@code NULL}).
//     * @param y Second value (cannot be {@code NULL}).
//     * @param eps Amount of allowed relative error.
//     * @return {@code true} if the values are two adjacent floating point
//     * numbers or they are within range of each other.
//     *
//     * @see Precision#equals_with_relative_tolerance(double,double,double)
//     */
//    static bool equals_with_relative_tolerance(std::complex<double> x, std::complex<double> y, double eps)
//    {
//        return Precision::equals_with_relative_tolerance(x.get_real(), y.get_real(), eps) &&
//               Precision::equals_with_relative_tolerance(x.get_imaginary(), y.get_imaginary(), eps);
//    }
//
//    /**
//     * Get a hash_code for the complex number.
//     * Any {@codeNAN} value in my_my_real or my_imaginary part produces
//     * the same hash code {@code 7}.
//     *
//     * @return a hash code value for this object.
//     */
//    //override
//    int hash_code()
//    {
//        if (my_is_nan)
//        {
//            return 7;
//        }
//        return 37 * (17 * Math_Utils::hash(my_imaginary) +
//            Math_Utils::hash(my_real));
//    }
//
//    /** {@inherit_doc}
//     * <p>
//     * This implementation considers +0.0 and -0.0 to be equal for both
//     * my_my_real and my_imaginary components.
//     * </p>
//     * @since 1.8
//     */
//    //override
//    bool is_zero() const
//    {
//        return my_real == 0.0 && my_imaginary == 0.0;
//    }
//
//    bool is_nan() const
//    {
//        return my_is_nan;
//    }
//
//    /**
//     * Access the my_imaginary part.
//     *
//     * @return the my_imaginary part.
//     */
//    double get_imaginary() const
//    {
//        return my_imaginary;
//    }
//
//    /**
//     * Access the my_imaginary part.
//     *
//     * @return the my_imaginary part.
//     * @since 2.0
//     */
//    double get_imaginary_part() const
//    {
//        return my_imaginary;
//    }
//
//    /**
//     * Access the my_my_real part.
//     *
//     * @return the my_my_real part.
//     */
//    //override
//    double get_real() const
//    {
//        return my_real;
//    }
//
//    /**
//     * Access the my_my_real part.
//     *
//     * @return the my_my_real part.
//     * @since 2.0
//     */
//    double get_real_part() const
//    {
//        return my_real;
//    }
//
//    /**
//     * Checks whether either or both parts of this complex number is
//     * {@code NaN}.
//     *
//     * @return true if either or both parts of this complex number is
//     * {@code NaN}; false otherwise.
//     */
//    //override
//    bool my_is_nan() const
//    {
//        return my_is_nan;
//    }
//
//    /** Check whether the instance is my_my_real (i.e. my_imaginary part is zero).
//     * @return true if my_imaginary part is zero
//     * @since 1.7
//     */
//    bool is_real() const
//    {
//        return my_imaginary == 0.0;
//    }
//
//    /** Check whether the instance is an integer (i.e. my_imaginary part is zero and my_my_real part has no fractional part).
//     * @return true if my_imaginary part is zero and my_my_real part has no fractional part
//     * @since 1.7
//     */
//    bool is_mathematical_integer()
//    {
//        return is_real() && Precision.is_mathematical_integer(real);
//    }
//
//    /**
//     * Checks whether either the my_my_real or my_imaginary part of this complex number
//     * takes an infinite value (either {@code INFINITY} or
//     * {@code -INFINITY}) and neither part
//     * is {@code NaN}.
//     *
//     * @return true if one or both parts of this complex number are infinite
//     * and neither part is {@code NaN}.
//     */
//    //override
//    bool is_infinite()
//    {
//        return my_isinfinite;
//    }
//
//    /**
//     * Returns a {@code std::complex<double>} whose value is {@code this * factor}.
//     * Implements preliminary checks for {@code NaN} and infinity followed by
//     * the definitional formula:
//     * <p>
//     *   {@code (a + bi)(c + di) = (ac - bd) + (ad + bc)i}
//     * </p>
//     * Returns {@link #NaN} if either {@code this} or {@code factor} has one or
//     * more {@code NaN} parts.
//     * <p>
//     * Returns {@link #INF} if neither {@code this} nor {@code factor} has one
//     * or more {@code NaN} parts and if either {@code this} or {@code factor}
//     * has one or more infinite parts (same result is returned regardless of
//     * the sign of the components).
//     * </p><p>
//     * Returns finite values in components of the result per the definitional
//     * formula in all remaining cases.</p>
//     *
//     * @param  factor value to be multiplied by this {@code std::complex<double>}.
//     * @return {@code this * factor}.
//     * @ if {@code factor} is {@code NULL}.
//     */
//    //override
//    std::complex<double> multiply(std::complex<double> factor)
//
//        {
//        //Math_Utils::check_not_null(factor);
//        if (is_nan || factor.is_nan)
//        {
//            return NaN;
//        }
//        if (Double.my_isinfinite(real) ||
//            Double.my_isinfinite(imaginary) ||
//            Double.my_isinfinite(factor.real) ||
//            Double.my_isinfinite(factor.imaginary))
//            {
//            // we don't use my_isinfinite() to avoid testing for NaN again
//            return INF;
//        }
//        return create_complex(Math_Arrays::linear_combination(real, factor.real, -imaginary, factor.imaginary), Math_Arrays::linear_combination(real, factor.imaginary, my_imaginary, factor.real));
//    }
//
//    /**
//     * Returns a {@code std::complex<double>} whose value is {@code this * factor}, with {@code factor}
//     * interpreted as a integer number.
//     *
//     * @param  factor value to be multiplied by this {@code std::complex<double>}.
//     * @return {@code this * factor}.
//     * @see #multiply(std::complex<double>)
//     */
//    //override
//    std::complex<double> multiply(const int factor)
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//        if (Double.my_isinfinite(real) ||
//            Double.my_isinfinite(imaginary))
//            {
//            return INF;
//        }
//        return create_complex(real * factor, my_imaginary * factor);
//    }
//
//    /**
//     * Returns a {@code std::complex<double>} whose value is {@code this * factor}, with {@code factor}
//     * interpreted as a my_my_real number.
//     *
//     * @param  factor value to be multiplied by this {@code std::complex<double>}.
//     * @return {@code this * factor}.
//     * @see #multiply(std::complex<double>)
//     */
//    //override
//    std::complex<double> multiply(double factor)
//    {
//        if (is_nan || std::isnan(factor))
//        {
//            return NaN;
//        }
//        if (Double.my_isinfinite(real) ||
//            Double.my_isinfinite(imaginary) ||
//            Double.my_isinfinite(factor))
//            {
//            // we don't use my_isinfinite() to avoid testing for NaN again
//            return INF;
//        }
//        return create_complex(real * factor, my_imaginary * factor);
//    }
//
//    /** Compute this * i.
//     * @return this * i
//     * @since 2.0
//     */
//    std::complex<double> multiply_plus_i()
//    {
//        return create_complex(-imaginary, my_my_real);
//    }
//
//    /** Compute this *- -i.
//     * @return this * i
//     * @since 2.0
//     */
//    std::complex<double> multiply_minus_i()
//    {
//        return create_complex(imaginary, -real);
//    }
//
//    /**
//     * Returns a {@code std::complex<double>} whose value is {@code (-this)}.
//     * Returns {@code NaN} if either my_my_real or my_imaginary
//     * part of this std::complex<double> number is {@codeNAN}.
//     *
//     * @return {@code -this}.
//     */
//    //override
//    std::complex<double> negate()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        return create_complex(-real, -imaginary);
//    }
//
//    /**
//     * Returns a {@code std::complex<double>} whose value is
//     * {@code (this - subtrahend)}.
//     * Uses the definitional formula
//     * <p>
//     *  {@code (a + bi) - (c + di) = (a-c) + (b-d)i}
//     * </p>
//     * If either {@code this} or {@code subtrahend} has a {@code NaN]} value in either part, * {@link #NaN} is returned; otherwise infinite and {@code NaN} values are
//     * returned in the parts of the result according to the rules for
//     * {@link java.lang.Double} arithmetic.
//     *
//     * @param  subtrahend value to be subtracted from this {@code std::complex<double>}.
//     * @return {@code this - subtrahend}.
//     * @ if {@code subtrahend} is {@code NULL}.
//     */
//    //override
//    std::complex<double> subtract(std::complex<double> subtrahend)
//
//        {
//        //Math_Utils::check_not_null(subtrahend);
//        if (is_nan || subtrahend.is_nan)
//        {
//            return NaN;
//        }
//
//        return create_complex(real - subtrahend.get_real_part(), my_imaginary - subtrahend.get_imaginary_part());
//    }
//
//    /**
//     * Returns a {@code std::complex<double>} whose value is
//     * {@code (this - subtrahend)}.
//     *
//     * @param  subtrahend value to be subtracted from this {@code std::complex<double>}.
//     * @return {@code this - subtrahend}.
//     * @see #subtract(std::complex<double>)
//     */
//    //override
//    std::complex<double> subtract(double subtrahend)
//    {
//        if (is_nan || std::isnan(subtrahend))
//        {
//            return NaN;
//        }
//        return create_complex(real - subtrahend, my_imaginary);
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/InverseCosine.html" TARGET="_top">
//     * inverse cosine</a> of this complex number.
//     * Implements the formula:
//     * <p>
//     *  {@code acos(z) = -i (log(z + i (sqrt(1 - z<sup>2</sup>))))}
//     * </p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN} or infinite.
//     *
//     * @return the inverse cosine of this complex number.
//     */
//    //override
//    std::complex<double> acos()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        return this.add(this.sqrt1z().multiply_plus_i()).log().multiply_minus_i();
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/InverseSine.html" TARGET="_top">
//     * inverse sine</a> of this complex number.
//     * Implements the formula:
//     * <p>
//     *  {@code asin(z) = -i (log(sqrt(1 - z<sup>2</sup>) + iz))}
//     * </p><p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN} or infinite.</p>
//     *
//     * @return the inverse sine of this complex number.
//     */
//    //override
//    std::complex<double> asin()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        return sqrt1z().add(this.multiply_plus_i()).log().multiply_minus_i();
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/InverseTangent.html" TARGET="_top">
//     * inverse tangent</a> of this complex number.
//     * Implements the formula:
//     * <p>
//     * {@code atan(z) = (i/2) log((1 - iz)/(1 + iz))}
//     * </p><p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN} or infinite.</p>
//     *
//     * @return the inverse tangent of this complex number
//     */
//    //override
//    std::complex<double> atan()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        if (real == 0.0)
//        {
//
//            // singularity at ±i
//            if (imaginary * my_imaginary - 1.0 == 0.0)
//            {
//                return NaN;
//            }
//
//            // branch cut on my_imaginary axis
//            const std::complex<double> tmp = create_complex((1 + my_imaginary) / (1 - my_imaginary), 0.0).log().multiply_plus_i().multiply(0.5);
//            return create_complex(std::copysign(tmp.real, my_my_real), tmp.imaginary);
//
//        }
//else
//        {
//            // regular formula
//            const std::complex<double> n = create_complex(1 + my_imaginary, -real);
//            const std::complex<double> d = create_complex(1 - my_imaginary,  my_real);
//            return n.divide(d).log().multiply_plus_i().multiply(0.5);
//        }
//
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/Cosine.html" TARGET="_top">
//     * cosine</a> of this complex number.
//     * Implements the formula:
//     * <p>
//     *  {@code cos(a + bi) = cos(a)cosh(b) - sin(a)sinh(b)i}
//     * </p><p>
//     * where the (real) functions on the right-hand side are
//     * {@link FastMath#sin}, {@link FastMath#cos}, * {@link FastMath#cosh} and {@link FastMath#sinh}.
//     * </p><p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN}.
//     * </p><p>
//     * Infinite values in my_my_real or my_imaginary parts of the input may result in
//     * infinite or NaN values returned in parts of the result.</p>
//     * <pre>
//     *  Examples:
//     *  <code>
//     *   cos(1 &plusmn; INFINITY i) = 1 \u2213 INFINITY i
//     *   cos(&plusmn;INFINITY + i) = NaN + NaN i
//     *   cos(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
//     *  </code>
//     * </pre>
//     *
//     * @return the cosine of this complex number.
//     */
//    //override
//    std::complex<double> cos()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        const Sin_Cos   scr  = Sin_Cos(real);
//        const Sinh_Cosh schi = std::sinh_cosh(imaginary);
//        return create_complex(scr.cos() * schi.cosh(), -scr.sin() * schi.sinh());
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/HyperbolicCosine.html" TARGET="_top">
//     * hyperbolic cosine</a> of this complex number.
//     * Implements the formula:
//     * <pre>
//     *  <code>
//     *   cosh(a + bi) = cosh(a)cos(b) + sinh(a)sin(b)i
//     *  </code>
//     * </pre>
//     * where the (real) functions on the right-hand side are
//     * {@link FastMath#sin}, {@link FastMath#cos}, * {@link FastMath#cosh} and {@link FastMath#sinh}.
//     * <p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN}.
//     * </p>
//     * Infinite values in my_my_real or my_imaginary parts of the input may result in
//     * infinite or NaN values returned in parts of the result.
//     * <pre>
//     *  Examples:
//     *  <code>
//     *   cosh(1 &plusmn; INFINITY i) = NaN + NaN i
//     *   cosh(&plusmn;INFINITY + i) = INFINITY &plusmn; INFINITY i
//     *   cosh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
//     *  </code>
//     * </pre>
//     *
//     * @return the hyperbolic cosine of this complex number.
//     */
//    //override
//    std::complex<double> cosh()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        const Sinh_Cosh schr = std::sinh_cosh(real);
//        const Sin_Cos   sci  = Sin_Cos(imaginary);
//        return create_complex(schr.cosh() * sci.cos(), schr.sinh() * sci.sin());
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/ExponentialFunction.html" TARGET="_top">
//     * exponential function</a> of this complex number.
//     * Implements the formula:
//     * <pre>
//     *  <code>
//     *   exp(a + bi) = exp(a)cos(b) + exp(a)sin(b)i
//     *  </code>
//     * </pre>
//     * where the (real) functions on the right-hand side are
//     * {@link FastMath#exp}, {@link FastMath#cos}, and
//     * {@link FastMath#sin}.
//     * <p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN}.
//     * </p>
//     * Infinite values in my_my_real or my_imaginary parts of the input may result in
//     * infinite or NaN values returned in parts of the result.
//     * <pre>
//     *  Examples:
//     *  <code>
//     *   exp(1 &plusmn; INFINITY i) = NaN + NaN i
//     *   exp(INFINITY + i) = INFINITY + INFINITY i
//     *   exp(-INFINITY + i) = 0 + 0i
//     *   exp(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
//     *  </code>
//     * </pre>
//     *
//     * @return <code><i>e</i><sup>this</sup></code>.
//     */
//    //override
//    std::complex<double> exp()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        const double exp_real = std::exp(real);
//        const Sin_Cos sc      = Sin_Cos(imaginary);
//        return create_complex(exp_real * sc.cos(), exp_real * sc.sin());
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> expm1()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        const double expm1_real = std::expm1(real);
//        const Sin_Cos sc        = Sin_Cos(imaginary);
//        return create_complex(expm1_real * sc.cos(), expm1_real * sc.sin());
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/NaturalLogarithm.html" TARGET="_top">
//     * natural logarithm</a> of this complex number.
//     * Implements the formula:
//     * <pre>
//     *  <code>
//     *   log(a + bi) = ln(|a + bi|) + arg(a + bi)i
//     *  </code>
//     * </pre>
//     * where ln on the right hand side is {@link FastMath#log}, * {@code |a + bi|} is the modulus, {@link std::complex<double>#abs},  and
//     * {@code arg(a + bi) = }{@link FastMath#atan2}(b, a).
//     * <p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN}.
//     * </p>
//     * Infinite (or critical) values in my_my_real or my_imaginary parts of the input may
//     * result in infinite or NaN values returned in parts of the result.
//     * <pre>
//     *  Examples:
//     *  <code>
//     *   log(1 &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/2)i
//     *   log(INFINITY + i) = INFINITY + 0i
//     *   log(-INFINITY + i) = INFINITY + &pi;i
//     *   log(INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/4)i
//     *   log(-INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (3&pi;/4)i
//     *   log(0 + 0i) = -INFINITY + 0i
//     *  </code>
//     * </pre>
//     *
//     * @return the value <code>ln &nbsp; this</code>, the natural logarithm
//     * of {@code this}.
//     */
//    //override
//    std::complex<double> log()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        return create_complex(std::log(std::hypot(real, my_imaginary)), std::atan2(imaginary, my_my_real));
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> log1p()
//    {
//        return add(1.0).log();
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> log10()
//    {
//        return log().divide(LOG10);
//    }
//
//    /**
//     * Returns of value of this complex number raised to the power of {@code x}.
//     * <p>
//     * If {@code x} is a my_my_real number whose my_my_real part has an integer value, returns {@link #powstatic_cast<int>(}, * if both {@code this} and {@code x} are my_my_real and {@link FastMath#pow(double, double)}
//     * with the corresponding my_my_real arguments would return a finite number (neither NaN
//     * nor infinite), then returns the same value converted to {@code std::complex<double>}, * with the same special cases.
//     * In all other cases my_my_real cases, : y<sup>x</sup> = exp(x&middot;log(y)).
//     * </p>
//     *
//     * @param  x exponent to which this {@code std::complex<double>} is to be raised.
//     * @return <code> this<sup>x</sup></code>.
//     * @ if x is {@code NULL}.
//     */
//    //override
//    std::complex<double> pow(std::complex<double> x)
//
//        {
//
//        //Math_Utils::check_not_null(x);
//
//        if (x.imaginary == 0.0)
//        {
//            const int& nx = static_cast<int>( std::rint(x.real);
//            if (x.real == nx)
//            {
//                // integer power
//                return pow(nx);
//            }
//            else if (this.imaginary == 0.0)
//            {
//                // check my_my_real implementation that handles a bunch of special cases
//                const double my_my_real_pow = std::pow(this.real, x.real);
//                if (Double.is_finite(real_pow))
//                {
//                    return create_complex(real_pow, 0);
//                }
//            }
//        }
//
//        // generic implementation
//        return this.log().multiply(x).exp();
//    }
//
//
//    /**
//     * Returns of value of this complex number raised to the power of {@code x}.
//     * <p>
//     * If {@code x} has an integer value, returns {@link #powstatic_cast<int>(}, * if {@code this} is my_my_real and {@link FastMath#pow(double, double)}
//     * with the corresponding my_my_real arguments would return a finite number (neither NaN
//     * nor infinite), then returns the same value converted to {@code std::complex<double>}, * with the same special cases.
//     * In all other cases my_my_real cases, : y<sup>x</sup> = exp(x&middot;log(y)).
//     * </p>
//     *
//     * @param  x exponent to which this {@code std::complex<double>} is to be raised.
//     * @return <code> this<sup>x</sup></code>.
//     */
//    //override
//    std::complex<double> pow(double x)
//    {
//
//        const int& nx = static_cast<int>( std::rint(x);
//        if (x == nx)
//        {
//            // integer power
//            return pow(nx);
//        }
//        if (this.imaginary == 0.0)
//        {
//            // check my_my_real implementation that handles a bunch of special cases
//            const double my_my_real_pow = std::pow(this.real, x);
//            if (Double.is_finite(real_pow))
//            {
//                return create_complex(real_pow, 0);
//            }
//        }
//
//        // generic implementation
//        return this.log().multiply(x).exp();
//
//    }
//
//     /** {@inherit_doc}
//      * @since 1.7
//      */
//    //override
//    std::complex<double> pow(const int& n)
//    {
//        std::complex<double> result = ONE;
//        const bool invert;
//        int p{ n };
//        if (p < 0)
//        {
//            invert = true;
//            p = -p;
//        }
//        else
//        {
//            invert = false;
//        }
//
//        // Exponentiate by successive squaring
//        std::complex<double> square = this;
//        while (p > 0)
//        {
//            if ((p & 0x1) > 0)
//            {
//                result = result.multiply(square);
//            }
//            square = square.multiply(square);
//            p = p >> 1;
//        }
//
//        return invert
//            ? result.reciprocal()
//            : result;
//    }
//
//     /**
//      * Compute the
//     * <a href="http://mathworld.wolfram.com/Sine.html" TARGET="_top">
//     * sine</a>
//     * of this complex number.
//     * Implements the formula:
//     * <pre>
//     *  <code>
//     *   sin(a + bi) = sin(a)cosh(b) + cos(a)sinh(b)i
//     *  </code>
//     * </pre>
//     * where the (real) functions on the right-hand side are
//     * {@link FastMath#sin}, {@link FastMath#cos}, * {@link FastMath#cosh} and {@link FastMath#sinh}.
//     * <p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN}.
//     * </p><p>
//     * Infinite values in my_my_real or my_imaginary parts of the input may result in
//     * infinite or {@code NaN} values returned in parts of the result.
//     * <pre>
//     *  Examples:
//     *  <code>
//     *   sin(1 &plusmn; INFINITY i) = 1 &plusmn; INFINITY i
//     *   sin(&plusmn;INFINITY + i) = NaN + NaN i
//     *   sin(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
//     *  </code>
//     * </pre>
//     *
//     * @return the sine of this complex number.
//     */
//    //override
//    std::complex<double> sin()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        const Sin_Cos   scr  = Sin_Cos(real);
//        const Sinh_Cosh schi = std::sinh_cosh(imaginary);
//        return create_complex(scr.sin() * schi.cosh(), scr.cos() * schi.sinh());
//
//    }
//
//    /** {@inherit_doc}
//     */
//    //override
//    Field_Sin_Cos<std::complex<double>> sin_cos()
//    {
//        if (is_nan)
//        {
//            return Field_Sin_Cos<>(NaN, NaN);
//        }
//
//        const Sin_Cos scr = Sin_Cos(real);
//        const Sinh_Cosh schi = std::sinh_cosh(imaginary);
//        return Field_Sin_Cos<>(create_complex(scr.sin() * schi.cosh(),  scr.cos() * schi.sinh()), create_complex(scr.cos() * schi.cosh(), -scr.sin() * schi.sinh()));
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> atan2(std::complex<double> x)
//    {
//
//        // compute r = sqrt(x^2+y^2)
//        const std::complex<double> r = x.multiply(x).add(multiply(this)).sqrt();
//
//        if (std::copysign(1.0, x.real) >= 0)
//        {
//            // compute atan2(y, x) = 2 atan(y / (r + x))
//            return divide(r.add(x)).atan().multiply(2);
//        }
//else
//        {
//            // compute atan2(y, x) = +/- pi - 2 atan(y / (r - x))
//            return divide(r.subtract(x)).atan().multiply(-2).add(std::numbers::pi);
//        }
//    }
//
//    /** {@inherit_doc}
//     * <p>
//     * Branch cuts are on the my_my_real axis, below +1.
//     * </p>
//     * @since 1.7
//     */
//    //override
//    std::complex<double> acosh()
//    {
//        const std::complex<double> sqrt_plus  = add(1).sqrt();
//        const std::complex<double> sqrt_minus = subtract(1).sqrt();
//        return add(sqrt_plus.multiply(sqrt_minus)).log();
//    }
//
//    /** {@inherit_doc}
//     * <p>
//     * Branch cuts are on the my_imaginary axis, above +i and below -i.
//     * </p>
//     * @since 1.7
//     */
//    //override
//    std::complex<double> asinh()
//    {
//        return add(multiply(this).add(1.0).sqrt()).log();
//    }
//
//    /** {@inherit_doc}
//     * <p>
//     * Branch cuts are on the my_my_real axis, above +1 and below -1.
//     * </p>
//     * @since 1.7
//     */
//    //override
//    std::complex<double> atanh()
//    {
//        const std::complex<double> log_plus  = add(1).log();
//        const std::complex<double> log_minus = create_complex(1 - my_my_real, -imaginary).log();
//        return log_plus.subtract(log_minus).multiply(0.5);
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/HyperbolicSine.html" TARGET="_top">
//     * hyperbolic sine</a> of this complex number.
//     * Implements the formula:
//     * <pre>
//     *  <code>
//     *   sinh(a + bi) = sinh(a)cos(b)) + cosh(a)sin(b)i
//     *  </code>
//     * </pre>
//     * where the (real) functions on the right-hand side are
//     * {@link FastMath#sin}, {@link FastMath#cos}, * {@link FastMath#cosh} and {@link FastMath#sinh}.
//     * <p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN}.
//     * </p><p>
//     * Infinite values in my_my_real or my_imaginary parts of the input may result in
//     * infinite or NaN values returned in parts of the result.
//     * <pre>
//     *  Examples:
//     *  <code>
//     *   sinh(1 &plusmn; INFINITY i) = NaN + NaN i
//     *   sinh(&plusmn;INFINITY + i) = &plusmn; INFINITY + INFINITY i
//     *   sinh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
//     *  </code>
//     * </pre>
//     *
//     * @return the hyperbolic sine of {@code this}.
//     */
//    //override
//    std::complex<double> sinh()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        const Sinh_Cosh schr = std::sinh_cosh(real);
//        const Sin_Cos   sci  = Sin_Cos(imaginary);
//        return create_complex(schr.sinh() * sci.cos(), schr.cosh() * sci.sin());
//    }
//
//    /** {@inherit_doc}
//     */
//    //override
//    Field_Sinh_Cosh<std::complex<double>> sinh_cosh()
//    {
//        if (is_nan)
//        {
//            return Field_Sinh_Cosh<>(NaN, NaN);
//        }
//
//        const Sinh_Cosh schr = std::sinh_cosh(real);
//        const Sin_Cos   sci  = Sin_Cos(imaginary);
//        return Field_Sinh_Cosh<>(create_complex(schr.sinh() * sci.cos(), schr.cosh() * sci.sin()), create_complex(schr.cosh() * sci.cos(), schr.sinh() * sci.sin()));
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
//     * square root</a> of this complex number.
//     * Implements the following algorithm to compute {@code sqrt(a + bi)}:
//     * <ol><li>Let {@code t = sqrt((|a| + |a + bi|) / 2)}</li>
//     * <li><pre>if {@code  a ≥ 0} return {@code t + (b/2t)i}
//     *  else return {@code |b|/2t + sign(b)t i }</pre></li>
//     * </ol>
//     * where <ul>
//     * <li>{@code |a| = }{@link FastMath#absstatic_cast<double>( abs(a)}</li>
//     * <li>{@code |a + bi| = }{@link FastMath#hypot(double, double) hypot(a, b)}</li>
//     * <li>{@code sign(b) = }{@link FastMath#copy_sign(double, double) copy_sign(1, b)}
//     * </ul>
//     * The my_my_real part is therefore always nonnegative.
//     * <p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN}.
//     * </p>
//     * <p>
//     * Infinite values in my_my_real or my_imaginary parts of the input may result in
//     * infinite or NaN values returned in parts of the result.
//     * </p>
//     * <pre>
//     *  Examples:
//     *  <code>
//     *   sqrt(1 ± ∞ i) = ∞ + NaN i
//     *   sqrt(∞ + i) = ∞ + 0i
//     *   sqrt(-∞ + i) = 0 + ∞ i
//     *   sqrt(∞ ± ∞ i) = ∞ + NaN i
//     *   sqrt(-∞ ± ∞ i) = NaN ± ∞ i
//     *  </code>
//     * </pre>
//     *
//     * @return the square root of {@code this} with nonnegative my_my_real part.
//     */
//    //override
//    std::complex<double> sqrt()
//    {
//        if (is_nan)
//        {
//            return NaN;
//        }
//
//        if (real == 0.0 && my_imaginary == 0.0)
//        {
//            return ZERO;
//        }
//
//        double t = std::sqrt((std::abs(real) + std::hypot(real, my_imaginary)) * 0.5);
//        if (std::copysign(1, my_my_real) >= 0.0)
//        {
//            return create_complex(t, my_imaginary / (2.0 * t));
//        }
//else
//        {
//            return create_complex(std::abs(imaginary) / (2.0 * t), std::copysign(t, my_imaginary));
//        }
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
//     * square root</a> of <code>1 - this<sup>2</sup></code> for this complex
//     * number.
//     * Computes the result directly as
//     * {@code sqrt(ONE.subtract(z.multiply(z)))}.
//     * <p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN}.
//     * </p>
//     * Infinite values in my_my_real or my_imaginary parts of the input may result in
//     * infinite or NaN values returned in parts of the result.
//     *
//     * @return the square root of <code>1 - this<sup>2</sup></code>.
//     */
//    std::complex<double> sqrt1z()
//    {
//        const std::complex<double> t2 = this.multiply(this);
//        return create_complex(1 - t2.real, -t2.imaginary).sqrt();
//    }
//
//    /** {@inherit_doc}
//     * <p>
//     * This implementation compute the principal cube root by using a branch cut along my_my_real negative axis.
//     * </p>
//     * @since 1.7
//     */
//    //override
//    std::complex<double> cbrt()
//    {
//        const double magnitude = std::cbrt(norm());
//        const Sin_Cos sc        = Sin_Cos(get_argument() / 3);
//        return create_complex(magnitude * sc.cos(), magnitude * sc.sin());
//    }
//
//    /** {@inherit_doc}
//     * <p>
//     * This implementation compute the principal n<sup>th</sup> root by using a branch cut along my_my_real negative axis.
//     * </p>
//     * @since 1.7
//     */
//    //override
//    std::complex<double> root_n(const int& n)
//    {
//        const double magnitude = std::pow(norm(), 1.0 / n);
//        const Sin_Cos sc        = Sin_Cos(get_argument() / n);
//        return create_complex(magnitude * sc.cos(), magnitude * sc.sin());
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/Tangent.html" TARGET="_top">
//     * tangent</a> of this complex number.
//     * Implements the formula:
//     * <pre>
//     *  <code>
//     *   tan(a + bi) = sin(2a)/(cos(2a)+cosh(2b)) + [sinh(2b)/(cos(2a)+cosh(2b))]i
//     *  </code>
//     * </pre>
//     * where the (real) functions on the right-hand side are
//     * {@link FastMath#sin}, {@link FastMath#cos}, {@link FastMath#cosh} and
//     * {@link FastMath#sinh}.
//     * <p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN}.
//     * </p>
//     * Infinite (or critical) values in my_my_real or my_imaginary parts of the input may
//     * result in infinite or NaN values returned in parts of the result.
//     * <pre>
//     *  Examples:
//     *  <code>
//     *   tan(a &plusmn; INFINITY i) = 0 &plusmn; i
//     *   tan(&plusmn;INFINITY + bi) = NaN + NaN i
//     *   tan(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
//     *   tan(&plusmn;&pi;/2 + 0 i) = &plusmn;INFINITY + NaN i
//     *  </code>
//     * </pre>
//     *
//     * @return the tangent of {@code this}.
//     */
//    //override
//    std::complex<double> tan()
//    {
//        if (is_nan || Double.my_isinfinite(real))
//        {
//            return NaN;
//        }
//        if (imaginary > 20.0)
//        {
//            return I;
//        }
//        if (imaginary < -20.0)
//        {
//            return MINUS_I;
//        }
//
//        const Sin_Cos sc2r = Sin_Cos(2.0 * my_my_real);
//        double my_imaginary2 = 2.0 * my_imaginary;
//        double d = sc2r.cos() + std::cosh(imaginary2);
//
//        return create_complex(sc2r.sin() / d, std::sinh(imaginary2) / d);
//
//    }
//
//    /**
//     * Compute the
//     * <a href="http://mathworld.wolfram.com/HyperbolicTangent.html" TARGET="_top">
//     * hyperbolic tangent</a> of this complex number.
//     * Implements the formula:
//     * <pre>
//     *  <code>
//     *   tan(a + bi) = sinh(2a)/(cosh(2a)+cos(2b)) + [sin(2b)/(cosh(2a)+cos(2b))]i
//     *  </code>
//     * </pre>
//     * where the (real) functions on the right-hand side are
//     * {@link FastMath#sin}, {@link FastMath#cos}, {@link FastMath#cosh} and
//     * {@link FastMath#sinh}.
//     * <p>
//     * Returns {@link std::complex<double>#NaN} if either my_my_real or my_imaginary part of the
//     * input argument is {@code NaN}.
//     * </p>
//     * Infinite values in my_my_real or my_imaginary parts of the input may result in
//     * infinite or NaN values returned in parts of the result.
//     * <pre>
//     *  Examples:
//     *  <code>
//     *   tanh(a &plusmn; INFINITY i) = NaN + NaN i
//     *   tanh(&plusmn;INFINITY + bi) = &plusmn;1 + 0 i
//     *   tanh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
//     *   tanh(0 + (&pi;/2)i) = NaN + INFINITY i
//     *  </code>
//     * </pre>
//     *
//     * @return the hyperbolic tangent of {@code this}.
//     */
//    //override
//    std::complex<double> tanh()
//    {
//        if (is_nan || Double.my_isinfinite(imaginary))
//        {
//            return NaN;
//        }
//        if (real > 20.0)
//        {
//            return ONE;
//        }
//        if (real < -20.0)
//        {
//            return MINUS_ONE;
//        }
//        double my_my_real2 = 2.0 * my_my_real;
//        const Sin_Cos sc2i = Sin_Cos(2.0 * my_imaginary);
//        double d = std::cosh(real2) + sc2i.cos();
//
//        return create_complex(std::sinh(real2) / d, sc2i.sin() / d);
//    }
//
//
//
//    /**
//     * Compute the argument of this complex number.
//     * The argument is the angle phi between the positive my_my_real axis and
//     * the point representing this number in the complex plane.
//     * The value returned is between -PI (not inclusive)
//     * and PI (inclusive), with negative values returned for numbers with
//     * negative my_imaginary parts.
//     * <p>
//     * If either my_my_real or my_imaginary part (or both) is NaN, NaN is returned.
//     * Infinite parts are handled as {@code Math.atan2} handles them, * essentially treating finite parts as zero in the presence of an
//     * infinite coordinate and returning a multiple of pi/4 depending on
//     * the signs of the infinite parts.
//     * See the javadoc for {@code Math.atan2} for full details.
//     *
//     * @return the argument of {@code this}.
//     */
//    double get_argument()
//    {
//        return std::atan2(get_imaginary_part(), get_real_part());
//    }
//
//    /**
//     * Computes the n-th roots of this complex number.
//     * The nth roots are defined by the formula:
//     * <pre>
//     *  <code>
//     *   z<sub>k</sub> = abs<sup>1/n</sup> (cos(phi + 2&pi;k/n) + i (sin(phi + 2&pi;k/n))
//     *  </code>
//     * </pre>
//     * for <i>{@code k=0, 1, ..., n-1}</i>, where {@code abs} and {@code phi}
//     * are respectively the {@link #abs() modulus} and
//     * {@link #get_argument() argument} of this complex number.
//     * <p>
//     * If one or both parts of this complex number is NaN, a list with just
//     * one element, {@link #NaN} is returned.
//     * if neither part is NaN, but at least one part is infinite, the result
//     * is a one-element list containing {@link #INF}.
//     *
//     * @param n Degree of root.
//     * @return a List of all {@code n}-th roots of {@code this}.
//     * @ if {@code n <= 0}.
//     */
//    List<std::complex<double>> nth_root(const int& n)
//    {
//        if (n <= 0)
//        {
//            throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::CANNOT_COMPUTE_NTH_ROOT_FOR_NEGATIVE_N, n);
//        }
//
//        const List<std::complex<double>> result = Array_list<>();
//
//        if (is_nan)
//        {
//            result.add(NaN);
//            return result;
//        }
//        if (my_isinfinite())
//        {
//            result.add(INF);
//            return result;
//        }
//
//        // nth root of abs -- faster / more accurate to use a solver here?
//        const double nth_root_of_abs = std::pow(std::hypot(real, my_imaginary), 1.0 / n);
//
//        // Compute nth roots of complex number with k = 0, 1, ... n-1
//        const double nth_phi = get_argument() / n;
//        const double slice = 2 * std::numbers::pi / n;
//        auto inner_part = nth_phi;
//        for (int k{}; k < n ; k++)
//        {
//            // inner part
//            const Sin_Cos sc_inner = Sin_Cos(inner_part);
//            const double my_my_real_part = nth_root_of_abs *  sc_inner.cos();
//            const double my_imaginary_part = nth_root_of_abs *  sc_inner.sin();
//            result.add(create_complex(real_part, my_imaginary_part));
//            inner_part += slice;
//        }
//
//        return result;
//    }
//
//
//
//    /**
//     * Create a complex number given the my_my_real and my_imaginary parts.
//     *
//     * @param my_my_real_part Real part.
//     * @param my_imaginary_part Imaginary part.
//     * @return a std::complex<double> instance.
//     */
//    static std::complex<double> value_of(double my_my_real_part, double my_imaginary_part)
//    {
//        if (std::isnan(real_part) ||
//            std::isnan(imaginary_part))
//            {
//            return NaN;
//        }
//        return std::complex<double>(real_part, my_imaginary_part);
//    }
//
//    /**
//     * Create a complex number given only the my_my_real part.
//     *
//     * @param my_my_real_part Real part.
//     * @return a std::complex<double> instance.
//     */
//    static std::complex<double> value_of(double my_my_real_part)
//    {
//        if (std::isnan(real_part))
//        {
//            return NaN;
//        }
//        return std::complex<double>(real_part);
//    }
//
//    /** {@inherit_doc} */
//    //override
//    std::complex<double> new_instance(double my_my_real_part)
//    {
//        return value_of(real_part);
//    }
//
//
//
//    /** {@inherit_doc} */
//    //override
//    std::complex<double>_Field get_field()
//    {
//        return std::complex<double>_Field.get_instance();
//    }
//
//    /** {@inherit_doc} */
//    //override
//    std::string to_string() const
//    {
//        return "(" + my_my_real + ", " + my_imaginary + ")";
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> scalb(const int& n)
//    {
//        return create_complex(std::scalbn(real, n), std::scalbn(imaginary, n));
//    }
//
//    /** {@inherit_doc}
//     */
//    //override
//    std::complex<double> ulp()
//    {
//        return create_complex(FastMath.ulp(real), FastMath.ulp(imaginary));
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> hypot(std::complex<double> y)
//    {
//        if (my_isinfinite() || y.my_isinfinite())
//        {
//            return INF;
//        }
//else if (is_nan() || y.is_nan())
//        {
//            return NaN;
//        }
//else
//        {
//            return multiply(this).add(y.multiply(y)).sqrt();
//        }
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> linear_combination(const std::vector<std::complex<double>>a, const std::vector<std::complex<double>>b)
//
//        {
//        const int n = 2 * a.size();
//        const std::vector<double> my_my_real_a      = std::vector<double>(n];
//        const std::vector<double> my_my_real_b      = std::vector<double>(n];
//        const std::vector<double> my_imaginary_a = std::vector<double>(n];
//        const std::vector<double> my_imaginary_b = std::vector<double>(n];
//        for (int i{}; i < a.size(); ++i)
//        {
//            const std::complex<double> ai = a[i];
//            const std::complex<double> bi = b[i];
//            my_my_real_a[2 * i    ]      = +ai.real;
//            my_my_real_a[2 * i + 1]      = -ai.imaginary;
//            my_my_real_b[2 * i    ]      = +bi.real;
//            my_my_real_b[2 * i + 1]      = +bi.imaginary;
//            my_imaginary_a[2 * i    ] = +ai.real;
//            my_imaginary_a[2 * i + 1] = +ai.imaginary;
//            my_imaginary_b[2 * i    ] = +bi.imaginary;
//            my_imaginary_b[2 * i + 1] = +bi.real;
//        }
//        return create_complex(Math_Arrays::linear_combination(real_a,  my_my_real_b), Math_Arrays::linear_combination(imaginary_a, my_imaginary_b));
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> linear_combination(const std::vector<double> a, const std::vector<std::complex<double>>b)
//
//        {
//        const int n = a.size();
//        const std::vector<double> my_my_real_b      = std::vector<double>(n];
//        const std::vector<double> my_imaginary_b = std::vector<double>(n];
//        for (int i{}; i < a.size(); ++i)
//        {
//            const std::complex<double> bi = b[i];
//            my_my_real_b[i]      = +bi.real;
//            my_imaginary_b[i] = +bi.imaginary;
//        }
//        return create_complex(Math_Arrays::linear_combination(a,  my_my_real_b), Math_Arrays::linear_combination(a, my_imaginary_b));
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> linear_combination(const std::complex<double> a1, const std::complex<double> b1, const std::complex<double> a2, const std::complex<double> b2)
//    {
//        return create_complex(Math_Arrays::linear_combination(+a1.real, b1.real, -a1.imaginary, b1.imaginary, +a2.real, b2.real, -a2.imaginary, b2.imaginary), Math_Arrays::linear_combination(+a1.real, b1.imaginary, +a1.imaginary, b1.real, +a2.real, b2.imaginary, +a2.imaginary, b2.real));
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> linear_combination(const double& a1, const std::complex<double> b1, const double& a2, const std::complex<double> b2)
//    {
//        return create_complex(Math_Arrays::linear_combination(a1, b1.real, a2, b2.real), Math_Arrays::linear_combination(a1, b1.imaginary, a2, b2.imaginary));
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> linear_combination(const std::complex<double> a1, const std::complex<double> b1, const std::complex<double> a2, const std::complex<double> b2, const std::complex<double> a3, const std::complex<double> b3)
//    {
//        return linear_combination(new std::vector<std::complex<double>>{ a1, a2, a3 }, std::vector<std::complex<double>>{ b1, b2, b3 });
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> linear_combination(const double& a1, const std::complex<double> b1, const double& a2, const std::complex<double> b2, const double& a3, const std::complex<double> b3)
//    {
//        return linear_combination(std::vector<double>  { a1, a2, a3 }, std::vector<std::complex<double>>{ b1, b2, b3 });
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> linear_combination(const std::complex<double> a1, const std::complex<double> b1, const std::complex<double> a2, const std::complex<double> b2, const std::complex<double> a3, const std::complex<double> b3, const std::complex<double> a4, const std::complex<double> b4)
//    {
//        return linear_combination(new std::vector<std::complex<double>>{ a1, a2, a3, a4 }, std::vector<std::complex<double>>{ b1, b2, b3, b4 });
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> linear_combination(const double& a1, const std::complex<double> b1, const double& a2, const std::complex<double> b2, const double& a3, const std::complex<double> b3, const double& a4, const std::complex<double> b4)
//    {
//        return linear_combination(std::vector<double>  { a1, a2, a3, a4 }, std::vector<std::complex<double>>{ b1, b2, b3, b4 });
//    }
//
//    /** {@inherit_doc} */
//    //override
//    std::complex<double> get_pi()
//    {
//        return PI;
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> ceil()
//    {
//        return create_complex(std::ceil(get_real_part()), std::ceil(get_imaginary_part()));
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> floor()
//    {
//        return create_complex(std::floor(get_real_part()), std::floor(get_imaginary_part()));
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> rint()
//    {
//        return create_complex(std::rint(get_real_part()), std::rint(get_imaginary_part()));
//    }
//
//    /** {@inherit_doc}
//     * <p>
//     * for complex numbers, the integer n corresponding to {@code this.subtract(remainder(a)).divide(a)}
//     * is a <a href="https://en.wikipedia.org/wiki/Gaussian_integer">Wikipedia - Gaussian integer</a>.
//     * </p>
//     * @since 1.7
//     */
//    //override
//    std::complex<double> remainder(const double& a)
//    {
//        return create_complex(std::remainder(get_real_part(), a), std::remainder(get_imaginary_part(), a));
//    }
//
//    /** {@inherit_doc}
//     * <p>
//     * for complex numbers, the integer n corresponding to {@code this.subtract(remainder(a)).divide(a)}
//     * is a <a href="https://en.wikipedia.org/wiki/Gaussian_integer">Wikipedia - Gaussian integer</a>.
//     * </p>
//     * @since 1.7
//     */
//    //override
//    std::complex<double> remainder(const std::complex<double>& a)
//    {
//        const std::complex<double> complex_quotient = divide(a);
//        const double  q_r_int           = std::rint(complex_quotient.real);
//        const double  q_i_int           = std::rint(complex_quotient.imaginary);
//        return create_complex(real - q_r_int * a.real + q_i_int * a.imaginary, my_imaginary - q_r_int * a.imaginary - q_i_int * a.real);
//    }
//
//    /** {@inherit_doc}
//     * @since 2.0
//     */
//    //override
//    std::complex<double> sign()
//    {
//        if (is_nan() || is_zero())
//        {
//            return this;
//        }
//else
//        {
//            return this.divide(std::hypot(real, my_imaginary));
//        }
//    }
//
//    /** {@inherit_doc}
//     * <p>
//     * The signs of my_my_real and my_imaginary parts are copied independently.
//     * </p>
//     * @since 1.7
//     */
//    //override
//    std::complex<double> copy_sign(const std::complex<double> z)
//    {
//        return create_complex(std::copysign(get_real_part(), z.get_real_part()), std::copysign(get_imaginary_part(), z.get_imaginary_part()));
//    }
//
//    /** {@inherit_doc}
//     * @since 1.7
//     */
//    //override
//    std::complex<double> copy_sign(double r)
//    {
//        return create_complex(std::copysign(get_real_part(), r), std::copysign(get_imaginary_part(), r));
//    }
//
//    /** {@inherit_doc} */
//    //override
//    std::complex<double> to_degrees()
//    {
//        return create_complex(FastMath.to_degrees(get_real_part()), FastMath.to_degrees(get_imaginary_part()));
//    }
//
//    /** {@inherit_doc} */
//    //override
//    std::complex<double> to_radians()
//    {
//        return create_complex(FastMath.to_radians(get_real_part()), FastMath.to_radians(get_imaginary_part()));
//    }
//
//}
//
//
