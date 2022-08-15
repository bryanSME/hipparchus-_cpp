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

//import java.math.BigDecimal;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.;

/**
 * Utilities for comparing numbers.
 */
class Precision 
{
private:
    /** Exponent offset in IEEE754 representation. */
    static constexpr long EXPONENT_OFFSET{ 1023l };

    /** Offset to order signed double numbers lexicographically. */
    static constexpr long SGN_MASK{ 0x8000000000000000L };
    /** Offset to order signed double numbers lexicographically. */
    static constexpr int SGN_MASK_FLOAT{ 0x80000000 };
    /** Positive zero. */
    static constexpr double POSITIVE_ZERO{};
    /** Positive zero bits. */
    static const long POSITIVE_ZERO_DOUBLE_BITS;// = Double.double_to_raw_long_bits(+0.0);
    /** Negative zero bits. */
    static const long NEGATIVE_ZERO_DOUBLE_BITS;// = Double.double_to_raw_long_bits(-0.0);
    /** Positive zero bits. */
    static const int POSITIVE_ZERO_FLOAT_BITS;// = Float.float_to_raw_int_bits(+0.0f);
    /** Negative zero bits. */
    static const int NEGATIVE_ZERO_FLOAT_BITS;// = Float.float_to_raw_int_bits(-0.0f);
    /** Mask used to extract exponent from double bits. */
    static constexpr long MASK_DOUBLE_EXPONENT{ 0x7ff0000000000000L };
    /** Mask used to extract mantissa from double bits. */
    static constexpr long MASK_DOUBLE_MANTISSA{ 0x000fffffffffffffL };
    /** Mask used to add implicit high order bit for normalized double. */
    static constexpr long IMPLICIT_DOUBLE_HIGH_BIT{ 0x0010000000000000L };
    /** Mask used to extract exponent from float bits. */
    static constexpr int MASK_FLOAT_EXPONENT{ 0x7f800000 };
    /** Mask used to extract mantissa from float bits. */
    static constexpr int MASK_FLOAT_MANTISSA{ 0x007fffff };
    /** Mask used to add implicit high order bit for normalized float. */
    static constexpr int IMPLICIT_FLOAT_HIGH_BIT{ 0x00800000 };

    /**
     * Private constructor.
     */
    Precision() = default;

    /**
     * Rounds the given non-negative value to the "nearest" integer. Nearest is
     * determined by the rounding method specified. Rounding methods are defined
     * in {@link BigDecimal}.
     *
     * @param unscaled Value to round.
     * @param sign Sign of the original, scaled value.
     * @param rounding_method Rounding method, as defined in {@link BigDecimal}.
     * @return the rounded value.
     * @Math_Runtime_Exception if an exact operation is required but result is not exact
     * @ if {@code rounding_method} is not a valid rounding method.
     */
    static double round_unscaled(const double& unscaled, const double& sign, const int& rounding_method)
    {
    switch (rounding_method)
    {
    case BigDecimal.ROUND_CEILING:
        if (sign == -1)
        {
            unscaled = std::floor(FastMath.next_after(unscaled, -INFINITY));
        }
        else
        {
            unscaled = std::ceil(FastMath.next_after(unscaled, INFINITY));
        }
        break;
    case BigDecimal.ROUND_DOWN:
        unscaled = std::floor(FastMath.next_after(unscaled, -INFINITY));
        break;
    case BigDecimal.ROUND_FLOOR:
        if (sign == -1)
        {
            unscaled = std::ceil(FastMath.next_after(unscaled, INFINITY));
        }
        else
        {
            unscaled = std::floor(FastMath.next_after(unscaled, -INFINITY));
        }
        break;
    case BigDecimal.ROUND_HALF_DOWN:
    {
        unscaled = FastMath.next_after(unscaled, -INFINITY);
        double fraction = unscaled - std::floor(unscaled);
        if (fraction > 0.5)
        {
            unscaled = std::ceil(unscaled);
        }
        else
        {
            unscaled = std::floor(unscaled);
        }
        break;
    }
    case BigDecimal.ROUND_HALF_EVEN:
    {
        double fraction = unscaled - std::floor(unscaled);
        if (fraction > 0.5)
        {
            unscaled = std::ceil(unscaled);
        }
        else if (fraction < 0.5)
        {
            unscaled = std::floor(unscaled);
        }
        else
        {
            // The following equality test is intentional and needed for rounding purposes
            if (std::floor(unscaled) / 2.0 == std::floor(std::floor(unscaled) / 2.0)) { // even
                unscaled = std::floor(unscaled);
            }
            else { // odd
                unscaled = std::ceil(unscaled);
            }
        }
        break;
    }
    case BigDecimal.ROUND_HALF_UP:
    {
        unscaled = FastMath.next_after(unscaled, INFINITY);
        double fraction = unscaled - std::floor(unscaled);
        if (fraction >= 0.5)
        {
            unscaled = std::ceil(unscaled);
        }
        else
        {
            unscaled = std::floor(unscaled);
        }
        break;
    }
    case BigDecimal.ROUND_UNNECESSARY:
        if (unscaled != std::floor(unscaled))
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ARITHMETIC_EXCEPTION);
        }
        break;
    case BigDecimal.ROUND_UP:
        // do not round if the discarded fraction is equal to zero
        if (unscaled != std::floor(unscaled))
        {
            unscaled = std::ceil(FastMath.next_after(unscaled, INFINITY));
        }
        break;
    default:
        throw (hipparchus::exception::Localized_Core_Formats_Type::INVALID_ROUNDING_METHOD, rounding_method, "ROUND_CEILING", BigDecimal.ROUND_CEILING, "ROUND_DOWN", BigDecimal.ROUND_DOWN, "ROUND_FLOOR", BigDecimal.ROUND_FLOOR, "ROUND_HALF_DOWN", BigDecimal.ROUND_HALF_DOWN, "ROUND_HALF_EVEN", BigDecimal.ROUND_HALF_EVEN, "ROUND_HALF_UP", BigDecimal.ROUND_HALF_UP, "ROUND_UNNECESSARY", BigDecimal.ROUND_UNNECESSARY, "ROUND_UP", BigDecimal.ROUND_UP);
    }
    return unscaled;
    }

public:

    /**
     * Largest double-precision floating-point number such that
     * {@code 1 + EPSILON} is numerically equal to 1. This value is an upper
     * bound on the relative error due to rounding real numbers to double
     * precision floating-point numbers.
     * <p>
     * In IEEE 754 arithmetic, this is 2<sup>-53</sup>.
     *
     * @see <a href="http://en.wikipedia.org/wiki/Machine_epsilon">Machine epsilon</a>
     */
    static const double EPSILON;

    /**
     * Safe minimum, such that {@code 1 / SAFE_MIN} does not overflow.
     * <br/>
     * In IEEE 754 arithmetic, this is also the smallest normalized
     * number 2<sup>-1022</sup>.
     */
    static const double SAFE_MIN;

    static 
    {
        /*
         *  This was previously expressed as = 0x1.0p-53;
         *  However, OpenJDK (Sparc Solaris) cannot handle such small
         *  constants: MATH-721
         */
        EPSILON = Double.long_bits_to_double((EXPONENT_OFFSET - 53l) << 52);

        /*
         * This was previously expressed as = 0x1.0p-1022;
         * However, OpenJDK (Sparc Solaris) cannot handle such small
         * constants: MATH-721
         */
        SAFE_MIN = Double.long_bits_to_double((EXPONENT_OFFSET - 1022l) << 52);
    }



    /**
     * Compares two numbers given some amount of allowed error.
     *
     * @param x the first number
     * @param y the second number
     * @param eps the amount of error to allow when checking for equality
     * @return <ul><li>0 if  {@link #equals(double, double, double) equals(x, y, eps)}</li>
     *       <li>&lt; 0 if !{@link #equals(double, double, double) equals(x, y, eps)} &amp;&amp; x &lt; y</li>
     *       <li>&gt; 0 if !{@link #equals(double, double, double) equals(x, y, eps)} &amp;&amp; x &gt; y or
     *       either argument is NaN</li></ul>
     */
    static int compare_to(const double& x, const double& y, const double& eps) 
    {
        if (equals(x, y, eps)) 
        {
            return 0;
        }
        if (x < y) 
        {
            return -1;
        }
        return 1;
    }

    /**
     * Compares two numbers given some amount of allowed error.
     * Two float numbers are considered equal if there are {@code (max_ulps - 1)}
     * (or fewer) floating point numbers between them, i.e. two adjacent floating
     * point numbers are considered equal.
     * Adapted from <a
     * href="http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
     * Bruce Dawson</a>. Returns {@code false} if either of the arguments is NaN.
     *
     * @param x first value
     * @param y second value
     * @param max_ulps {@code (max_ulps - 1)} is the number of floating point
     * values between {@code x} and {@code y}.
     * @return <ul><li>0 if  {@link #equals(double, double, int) equals(x, y, max_ulps)}</li>
     *       <li>&lt; 0 if !{@link #equals(double, double, int) equals(x, y, max_ulps)} &amp;&amp; x &lt; y</li>
     *       <li>&gt; 0 if !{@link #equals(double, double, int) equals(x, y, max_ulps)} &amp;&amp; x &gt; y
     *       or either argument is NaN</li></ul>
     */
    static int compare_to(const double& x, const double& y, const int max_ulps) 
    {
        if (equals(x, y, max_ulps)) 
        {
            return 0;
        }
        if (x < y) 
        {
            return -1;
        }
        return 1;
    }

    /**
     * Returns true iff they are equal as defined by
     * {@link #equals(float,float,int) equals(x, y, 1)}.
     *
     * @param x first value
     * @param y second value
     * @return {@code true} if the values are equal.
     */
    static bool equals(const float& x, const float& y) 
    {
        return equals(x, y, 1);
    }

    /**
     * Returns true if both arguments are NaN or they are
     * equal as defined by {@link #equals(float,float) equals(x, y, 1)}.
     *
     * @param x first value
     * @param y second value
     * @return {@code true} if the values are equal or both are NaN.
     */
    static bool equals_including_nan(float x, float y) 
    {
        return (x != x || y != y) ? !(x != x ^ y != y) : equals(x, y, 1);
    }

    /**
     * Returns true if the arguments are equal or within the range of allowed
     * error (inclusive).  Returns {@code false} if either of the arguments
     * is NaN.
     *
     * @param x first value
     * @param y second value
     * @param eps the amount of absolute error to allow.
     * @return {@code true} if the values are equal or within range of each other.
     */
    static bool equals(float x, float y, float eps) 
    {
        return equals(x, y, 1) || std::abs(y - x) <= eps;
    }

    /**
     * Returns true if the arguments are both NaN, are equal, or are within the range
     * of allowed error (inclusive).
     *
     * @param x first value
     * @param y second value
     * @param eps the amount of absolute error to allow.
     * @return {@code true} if the values are equal or within range of each other, * or both are NaN.
     */
    static bool equals_including_nan(float x, float y, float eps) 
    {
        return equals_including_nan(x, y) || (std::abs(y - x) <= eps);
    }

    /**
     * Returns true if the arguments are equal or within the range of allowed
     * error (inclusive).
     * Two float numbers are considered equal if there are {@code (max_ulps - 1)}
     * (or fewer) floating point numbers between them, i.e. two adjacent floating
     * point numbers are considered equal.
     * Adapted from <a
     * href="http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
     * Bruce Dawson</a>.  Returns {@code false} if either of the arguments is NaN.
     *
     * @param x first value
     * @param y second value
     * @param max_ulps {@code (max_ulps - 1)} is the number of floating point
     * values between {@code x} and {@code y}.
     * @return {@code true} if there are fewer than {@code max_ulps} floating
     * point values between {@code x} and {@code y}.
     */
    static bool equals(const float& x, const float& y, const int& max_ulps) 
    {

        const int x_int = Float.float_to_raw_int_bits(x);
        const int y_int = Float.float_to_raw_int_bits(y);

        const bool is_equal;
        if (((x_int ^ y_int) & SGN_MASK_FLOAT) == 0) 
        {
            // number have same sign, there is no risk of overflow
            is_equal = std::abs(x_int - y_int) <= max_ulps;
        }
        else 
        {
            // number have opposite signs, take care of overflow
            const int delta_plus;
            const int delta_minus;
            if (x_int < y_int) 
            {
                delta_plus  = y_int - POSITIVE_ZERO_FLOAT_BITS;
                delta_minus = x_int - NEGATIVE_ZERO_FLOAT_BITS;
            }
            else 
            {
                delta_plus  = x_int - POSITIVE_ZERO_FLOAT_BITS;
                delta_minus = y_int - NEGATIVE_ZERO_FLOAT_BITS;
            }

            if (delta_plus > max_ulps) 
            {
                is_equal = false;
            }
            else 
            {
                is_equal = delta_minus <= (max_ulps - delta_plus);
            }

        }

        return is_equal && !Float.is_nan(x) && !Float.is_nan(y);

    }

    /**
     * Returns true if the arguments are both NaN or if they are equal as defined
     * by {@link #equals(float,float,int) equals(x, y, max_ulps)}.
     *
     * @param x first value
     * @param y second value
     * @param max_ulps {@code (max_ulps - 1)} is the number of floating point
     * values between {@code x} and {@code y}.
     * @return {@code true} if both arguments are NaN or if there are less than
     * {@code max_ulps} floating point values between {@code x} and {@code y}.
     */
    static bool equals_including_nan(const float& x, const float& y, const int& max_ulps) 
    {
        return (x != x || y != y) ? !(x != x ^ y != y) : equals(x, y, max_ulps);
    }

    /**
     * Returns true iff they are equal as defined by
     * {@link #equals(double,double,int) equals(x, y, 1)}.
     *
     * @param x first value
     * @param y second value
     * @return {@code true} if the values are equal.
     */
    static bool equals(const double& x, const double& y) 
    {
        return equals(x, y, 1);
    }

    /**
     * Returns true if the arguments are both NaN or they are
     * equal as defined by {@link #equals(double,double) equals(x, y, 1)}.
     *
     * @param x first value
     * @param y second value
     * @return {@code true} if the values are equal or both are NaN.
     */
    static bool equals_including_nan(const double& x, const double& y) 
    {
        return (x != x || y != y) ? !(x != x ^ y != y) : equals(x, y, 1);
    }

    /**
     * Returns {@code true} if there is no double value strictly between the
     * arguments or the difference between them is within the range of allowed
     * error (inclusive). Returns {@code false} if either of the arguments
     * is NaN.
     *
     * @param x First value.
     * @param y Second value.
     * @param eps Amount of allowed absolute error.
     * @return {@code true} if the values are two adjacent floating point
     * numbers or they are within range of each other.
     */
    static bool equals(const double& x, const double& y, const double& eps) 
    {
        return equals(x, y, 1) || std::abs(y - x) <= eps;
    }

    /**
     * Returns {@code true} if there is no double value strictly between the
     * arguments or the relative difference between them is less than or equal
     * to the given tolerance. Returns {@code false} if either of the arguments
     * is NaN.
     *
     * @param x First value.
     * @param y Second value.
     * @param eps Amount of allowed relative error.
     * @return {@code true} if the values are two adjacent floating point
     * numbers or they are within range of each other.
     */
    static bool equals_with_relative_tolerance(const double& x, const double& y, const double& eps) 
    {
        if (equals(x, y, 1)) 
        {
            return true;
        }

        const double absolute_max = std::max(std::abs(x), std::abs(y));
        const double relative_difference = std::abs((x - y) / absolute_max);

        return relative_difference <= eps;
    }

    /**
     * Returns true if the arguments are both NaN, are equal or are within the range
     * of allowed error (inclusive).
     *
     * @param x first value
     * @param y second value
     * @param eps the amount of absolute error to allow.
     * @return {@code true} if the values are equal or within range of each other, * or both are NaN.
     */
    static bool equals_including_nan(const double& x, const double& y, const double& eps) 
    {
        return equals_including_nan(x, y) || (std::abs(y - x) <= eps);
    }

    /**
     * Returns true if the arguments are equal or within the range of allowed
     * error (inclusive).
     * <p>
     * Two float numbers are considered equal if there are {@code (max_ulps - 1)}
     * (or fewer) floating point numbers between them, i.e. two adjacent
     * floating point numbers are considered equal.
     * </p>
     * <p>
     * Adapted from <a
     * href="http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
     * Bruce Dawson</a>. Returns {@code false} if either of the arguments is NaN.
     * </p>
     *
     * @param x first value
     * @param y second value
     * @param max_ulps {@code (max_ulps - 1)} is the number of floating point
     * values between {@code x} and {@code y}.
     * @return {@code true} if there are fewer than {@code max_ulps} floating
     * point values between {@code x} and {@code y}.
     */
    static bool equals(const double& x, const double& y, const int& max_ulps) 
    {

        const long x_int = Double.double_to_raw_long_bits(x);
        const long y_int = Double.double_to_raw_long_bits(y);

        const bool is_equal;
        if (((x_int ^ y_int) & SGN_MASK) == 0l) 
        {
            // number have same sign, there is no risk of overflow
            is_equal = std::abs(x_int - y_int) <= max_ulps;
        }
        else 
        {
            // number have opposite signs, take care of overflow
            const long delta_plus;
            const long delta_minus;
            if (x_int < y_int) 
            {
                delta_plus  = y_int - POSITIVE_ZERO_DOUBLE_BITS;
                delta_minus = x_int - NEGATIVE_ZERO_DOUBLE_BITS;
            }
            else 
            {
                delta_plus  = x_int - POSITIVE_ZERO_DOUBLE_BITS;
                delta_minus = y_int - NEGATIVE_ZERO_DOUBLE_BITS;
            }

            if (delta_plus > max_ulps) 
            {
                is_equal = false;
            }
            else 
            {
                is_equal = delta_minus <= (max_ulps - delta_plus);
            }

        }

        return is_equal && !std::isnan(x) && !std::isnan(y);

    }

    /**
     * Returns true if both arguments are NaN or if they are equal as defined
     * by {@link #equals(double,double,int) equals(x, y, max_ulps)}.
     *
     * @param x first value
     * @param y second value
     * @param max_ulps {@code (max_ulps - 1)} is the number of floating point
     * values between {@code x} and {@code y}.
     * @return {@code true} if both arguments are NaN or if there are less than
     * {@code max_ulps} floating point values between {@code x} and {@code y}.
     */
    static bool equals_including_nan(const double& x, const double& y, const int& max_ulps) 
    {
        return (x != x || y != y) ? !(x != x ^ y != y) : equals(x, y, max_ulps);
    }

    /**
     * Rounds the given value to the specified number of decimal places.
     * The value is rounded using the {@link BigDecimal#ROUND_HALF_UP} method.
     *
     * @param x Value to round.
     * @param scale Number of digits to the right of the decimal point.
     * @return the rounded value.
     */
    static double round(const double& x, const int& scale) 
    {
        return round(x, scale, BigDecimal.ROUND_HALF_UP);
    }

    /**
     * Rounds the given value to the specified number of decimal places.
     * The value is rounded using the given method which is any method defined
     * in {@link BigDecimal}.
     * If {@code x} is infinite or {@code NaN}, then the value of {@code x} is
     * returned unchanged, regardless of the other parameters.
     *
     * @param x Value to round.
     * @param scale Number of digits to the right of the decimal point.
     * @param rounding_method Rounding method as defined in {@link BigDecimal}.
     * @return the rounded value.
     * @Arithmetic_Exception if {@code rounding_method == ROUND_UNNECESSARY}
     * and the specified scaling operation would require rounding.
     * @Illegal_Argument_Exception if {@code rounding_method} does not
     * represent a valid rounding mode.
     */
    static double round(const double& x, const int& scale, const int& rounding_method) 
    {
        try 
        {
            const double rounded = (new BigDecimal(std::to_string(x))
                   .set_scale(scale, rounding_method))
                   .double_value();
            // MATH-1089: negative values rounded to zero should result in negative zero
            return rounded == POSITIVE_ZERO ? POSITIVE_ZERO * x : rounded;
        }
        catch (Number_FormatException ex) 
        {
            if (std::isinf(x)) 
            {
                return x;
            }
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    /**
     * Rounds the given value to the specified number of decimal places.
     * The value is rounded using the {@link BigDecimal#ROUND_HALF_UP} method.
     *
     * @param x Value to round.
     * @param scale Number of digits to the right of the decimal point.
     * @return the rounded value.
     */
    static float round(const float& x, const int& scale) 
    {
        return round(x, scale, BigDecimal.ROUND_HALF_UP);
    }

    /**
     * Rounds the given value to the specified number of decimal places.
     * The value is rounded using the given method which is any method defined
     * in {@link BigDecimal}.
     *
     * @param x Value to round.
     * @param scale Number of digits to the right of the decimal point.
     * @param rounding_method Rounding method as defined in {@link BigDecimal}.
     * @return the rounded value.
     * @Math_Runtime_Exception if an exact operation is required but result is not exact
     * @ if {@code rounding_method} is not a valid rounding method.
     */
    static float round(const float& x, const int& scale, const int& rounding_method)
    {
        const float sign = std::copysign(1, x);
        const float factor = (float) std::pow(10.0f, scale) * sign;
        return (float) round_unscaled(x * factor, sign, rounding_method) / factor;
    }

    /** Check is x is a mathematical integer.
     * @param x number to check
     * @return true if x is a mathematical integer
     * @since 1.7
     */
    static bool is_mathematical_integer(const double& x) 
    {
        throw std::exception("not implemented");
        //const long bits = Double.double_to_raw_long_bits(x);
        //const int raw_exp = static_cast<int>( ((bits & MASK_DOUBLE_EXPONENT) >> 52);
        //if (raw_exp == 2047) 
        //{
        //    // NaN or infinite
        //    return false;
        //}

        //// a double that may have a fractional part
        //const long raw_mantissa = bits & MASK_DOUBLE_MANTISSA;
        //const long full_mantissa = raw_exp > 0
        //    ? (IMPLICIT_DOUBLE_HIGH_BIT | raw_mantissa)
        //    : raw_mantissa;
        //const long fractional_mask = (IMPLICIT_DOUBLE_HIGH_BIT | MASK_DOUBLE_MANTISSA) >> std::min(53, std::max(0, raw_exp - 1022));
        //return (full_mantissa & fractional_mask) == 0l;   
    }

    /** Check is x is a mathematical integer.
     * @param x number to check
     * @return true if x is a mathematical integer
     * @since 1.7
     */
    static bool is_mathematical_integer(const float& x) 
    {
        throw std::exception("not implemented");
        //const int bits = Float.float_to_raw_int_bits(x);
        //const int raw_exp = static_cast<int>( ((bits & MASK_FLOAT_EXPONENT) >> 23);
        //if (raw_exp == 255) 
        //{
        //    // NaN or infinite
        //    return false;
        //}

        //// a float that may have a fractional part
        //const int raw_mantissa    = bits & MASK_FLOAT_MANTISSA;
        //const int full_mantissa   = raw_exp > 0 ? (IMPLICIT_FLOAT_HIGH_BIT | raw_mantissa) : raw_mantissa;
        //const int fractional_mask = (IMPLICIT_FLOAT_HIGH_BIT | MASK_FLOAT_MANTISSA) >> std::min(24, std::max(0, raw_exp - 126));
        //return (full_mantissa & fractional_mask) == 0;
    }

    /**
     * Computes a number {@code delta} close to {@code original_delta} with
     * the property that <pre><code>
     *   x + delta - x
     * </code></pre>
     * is exactly machine-representable.
     * This is useful when computing numerical derivatives, in order to reduce
     * roundoff errors.
     *
     * @param x Value.
     * @param original_delta Offset value.
     * @return a number {@code delta} so that {@code x + delta} and {@code x}
     * differ by a representable floating number.
     */
    static double representable_delta(const double& x, const double& original_delta)
    {
        return x + original_delta - x;
    }
};