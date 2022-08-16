#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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


//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
//import org.hipparchus.util.Field_Sinh_Cosh;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Precision;
#include <type_traits>
#include <cmath>
#include "../CalculusFieldElement.hpp"


/**
 * Representation of a std::complex<double> number, i.e. a number which has both a
 * real and imaginary part.
 * <p>
 * Implementations of arithmetic operations handle {@code NaN} and
 * infinite values according to the rules for {@link java.lang.Double}, i.e.
 * {@link #equals} is an equivalence relation for all instances that have
 * a {@code NaN} in either real or imaginary part, e.g. the following are
 * considered equal:
 * <ul>
 *  <li>{@code 1 + NaNi}</li>
 *  <li>{@code NaN + i}</li>
 *  <li>{@code NaN + NaNi}</li>
 * </ul>
 * <p>
 * Note that this contradicts the IEEE-754 standard for floating
 * point numbers (according to which the test {@code x == x} must fail if
 * {@code x} is {@code NaN}). The method
 * {@link org.hipparchus.util.Precision#equals(double,double,int)
 * equals for primitive double} in {@link org.hipparchus.util.Precision}
 * conforms with IEEE-754 while this class conforms with the standard behavior
 * for Java object types.
 * @param <T> the type of the field elements
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<Field_Complex<T>>, T>::value>::type* = nullptr>
class Field_Complex : Calculus_Field_Element<Field_Complex<T>>  
{
private:
    /** A real number representing log(10). */
    static constexpr double LOG10{ 2.302585092994045684 };

    /** The imaginary part. */
    const T my_imaginary;

    /** The real part. */
    const T my_real;

    /** Record whether this complex number is equal to NaN. */
    const bool my_is_nan;

    /** Record whether this complex number is infinite. */
    const bool my_is_infinite;

protected:
    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param real_part Real part.
     * @param imaginary_part Imaginary part.
     * @return a complex number instance.
     *
     * @see #value_of(Calculus_Field_Element, Calculus_Field_Element)
     */
    Field_Complex<T> create_complex(const T real_part, const T imaginary_part)
    {
        return Field_Complex<double><>(real_part, imaginary_part);
    }

public:

    /**
     * Create a complex number given only the real part.
     *
     * @param real Real part.
     */
    Field_Complex<double>(T real) 
    {
        Field_Complex<double>(real, real.get_field().get_zero());
    }

    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     */
    Field_Complex<double>(const T& real, const T& imaginary) : my_real{ real }, my_imaginary{ imaginary }
    {
        my_is_nan = real.get_is_nan() || imaginary.get_is_nan();
        my_is_infinite = !is_nan && (real.get_is_infinite() || imaginary.get_is_infinite());
    }

    /** Get the square root of -1.
     * @param field field the complex components belong to
     * @return number representing "0.0 + 1.0i"
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex<T> get_i(const Field<T>& field) 
    {
        return Field_Complex<double><>(field.get_zero(), field.get_one());
    }

    /** Get the square root of -1.
     * @param field field the complex components belong to
     * @return number representing "0.0 _ 1.0i"
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex<T> get_minus_i(const Field<T> field) 
    {
        return Field_Complex<double><>(field.get_zero(), field.get_one().negate());
    }

    /** Get a complex number representing "NaN + NaNi".
     * @param field field the complex components belong to
     * @return complex number representing "NaN + NaNi"
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex<T> get_nan(const Field<T> field) 
    {

        return Field_Complex<double><>(field.get_zero().add(NAN), field.get_zero().add(NAN));
    }

    /** Get a complex number representing "+INF + INFi".
     * @param field field the complex components belong to
     * @return complex number representing "+INF + INFi"
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex<T> get_inf(const Field<T> field) 
    {
        return Field_Complex<double><>(field.get_zero().add(INFINITY), field.get_zero().add(INFINITY));
    }

    /** Get a complex number representing "1.0 + 0.0i".
     * @param field field the complex components belong to
     * @return complex number representing "1.0 + 0.0i"
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex<T> get_one(const Field<T> field) 
    {
        return Field_Complex<double><>(field.get_one(), field.get_zero());
    }

    /** Get a complex number representing "-1.0 + 0.0i".
     * @param field field the complex components belong to
     * @return complex number representing "-1.0 + 0.0i"
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex<T> get_minus_one(const Field<T> field) 
    {
        return Field_Complex<double><>(field.get_one().negate(), field.get_zero());
    }

    /** Get a complex number representing "0.0 + 0.0i".
     * @param field field the complex components belong to
     * @return complex number representing "0.0 + 0.0i
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex<T> get_zero(const Field<T> field) 
    {
        return Field_Complex<double><>(field.get_zero(), field.get_zero());
    }

    /** Get a complex number representing "\xcf\x80 + 0.0i".
     * @param field field the complex components belong to
     * @return complex number representing "\xcf\x80 + 0.0i
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex<T> get_pi(const Field<T> field) 
    {
        return Field_Complex<double><>(field.get_zero().get_pi(), field.get_zero());
    }

    /**
     * Return the absolute value of this complex number.
     * Returns {@code NaN} if either real or imaginary part is {@code NaN}
     * and {@code INFINITY} if neither part is {@code NaN}, * but at least one part is infinite.
     *
     * @return the absolute value.
     */
    //override
    Field_Complex<T> abs() 
    {
        // we check NaN here because std::hypot checks it after infinity
        return is_nan ? get_nan(get_parts_field()) : create_complex(std::hypot(my_real, my_imaginary), get_parts_field().get_zero());
    }

    T get_real() const
    {
        return my_real;
    }

    T get_imaginary() const
    {
        return my_imaginary;
    }

    /**
     * Returns a {@code std::complex<double>} whose value is
     * {@code (this + addend)}.
     * Uses the definitional formula
     * <p>
     *   {@code (a + bi) + (c + di) = (a+c) + (b+d)i}
     * </p>
     * If either {@code this} or {@code addend} has a {@code NaN} value in
     * either part, {@link #get_nan(Field)} is returned; otherwise {@code Infinite}
     * and {@code NaN} values are returned in the parts of the result
     * according to the rules for {@link java.lang.Double} arithmetic.
     *
     * @param  addend Value to be added to this {@code std::complex<double>}.
     * @return {@code this + addend}.
     * @ if {@code addend} is {@code NULL}.
     */
    //override
    Field_Complex<T> add(Field_Complex<T> addend)  
    {
        //Math_Utils::check_not_null(addend);
        if (is_nan || addend.is_nan) 
        {
            return get_nan(get_parts_field());
        }

        return create_complex(my_real.add(addend.get_real_part()), my_imaginary.add(addend.get_imaginary_part()));
    }

    /**
     * Returns a {@code std::complex<double>} whose value is {@code (this + addend)}, * with {@code addend} interpreted as a real number.
     *
     * @param addend Value to be added to this {@code std::complex<double>}.
     * @return {@code this + addend}.
     * @see #add(Field_Complex<double>)
     */
    Field_Complex<T> add(T addend) 
    {
        if (is_nan || addend.is_nan()) 
        {
            return get_nan(get_parts_field());
        }

        return create_complex(my_real.add(addend), my_imaginary);
    }

    /**
     * Returns a {@code std::complex<double>} whose value is {@code (this + addend)}, * with {@code addend} interpreted as a real number.
     *
     * @param addend Value to be added to this {@code std::complex<double>}.
     * @return {@code this + addend}.
     * @see #add(Field_Complex<double>)
     */
    //override
    Field_Complex<T> add(const double& addend) 
    {
        if (is_nan || std::isnan(addend)) 
        {
            return get_nan(get_parts_field());
        }

        return create_complex(my_real.add(addend), my_imaginary);
    }

    /**
     * Returns the conjugate of this complex number.
     * The conjugate of {@code a + bi} is {@code a - bi}.
     * <p>
     * {@link #get_nan(Field)} is returned if either the real or imaginary
     * part of this std::complex<double> number equals {@codeNAN}.
     * </p><p>
     * If the imaginary part is infinite, and the real part is not
     * {@code NaN}, the returned value has infinite imaginary part
     * of the opposite sign, e.g. the conjugate of
     * {@code 1 + POSITIVE_INFINITY i} is {@code 1 - NEGATIVE_INFINITY i}.
     * </p>
     * @return the conjugate of this std::complex<double> object.
     */
    Field_Complex<T> conjugate() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        return create_complex(real, my_imaginary.negate());
    }

    /**
     * Returns a {@code std::complex<double>} whose value is
     * {@code (this / divisor)}.
     * Implements the definitional formula
     * <pre>
     *  <code>
     *    a + bi          ac + bd + (bc - ad)i
     *    ----------- = -------------------------
     *    c + di         c<sup>2</sup> + d<sup>2</sup>
     *  </code>
     * </pre>
     * but uses
     * <a href="http://doi.acm.org/10.1145/1039813.1039814">
     * prescaling of operands</a> to limit the effects of overflows and
     * underflows in the computation.
     * <p>
     * {@code Infinite} and {@code NaN} values are handled according to the
     * following rules, applied in the order presented:
     * <ul>
     *  <li>If either {@code this} or {@code divisor} has a {@code NaN} value
     *   in either part, {@link #get_nan(Field)} is returned.
     *  </li>
     *  <li>If {@code divisor} equals {@link #get_zero(Field)}, {@link #get_nan(Field)} is returned.
     *  </li>
     *  <li>If {@code this} and {@code divisor} are both infinite, *   {@link #get_nan(Field)} is returned.
     *  </li>
     *  <li>If {@code this} is finite (i.e., has no {@code Infinite} or
     *   {@code NaN} parts) and {@code divisor} is infinite (one or both parts
     *   infinite), {@link #get_zero(Field)} is returned.
     *  </li>
     *  <li>If {@code this} is infinite and {@code divisor} is finite, *   {@code NaN} values are returned in the parts of the result if the
     *   {@link java.lang.Double} rules applied to the definitional formula
     *   force {@code NaN} results.
     *  </li>
     * </ul>
     *
     * @param divisor Value by which this {@code std::complex<double>} is to be divided.
     * @return {@code this / divisor}.
     * @ if {@code divisor} is {@code NULL}.
     */
    //override
    Field_Complex<T> divide(const Field_Complex<T>& divisor)
    {
        //Math_Utils::check_not_null(divisor);
        if (is_nan || divisor.is_nan) 
        {
            return get_nan(get_parts_field());
        }

        const T c = divisor.get_real_part();
        const T d = divisor.get_imaginary_part();
        if (c.is_zero() && d.is_zero()) 
        {
            return get_nan(get_parts_field());
        }

        if (divisor.is_infinite() && !is_infinite()) 
        {
            return get_zero(get_parts_field());
        }

        if (std::abs(c).get_real() < std::abs(d).get_real()) 
        {
            T q = c.divide(d);
            T inv_den = c.multiply(q).add(d).reciprocal();
            return create_complex(my_real.multiply(q).add(my_imaginary).multiply(inv_den), my_imaginary.multiply(q).subtract(my_real).multiply(inv_den));
        }

        T q = d.divide(c);
        T inv_den = d.multiply(q).add(c).reciprocal();
        return create_complex(my_imaginary.multiply(q).add(my_real).multiply(inv_den), my_imaginary.subtract(my_real.multiply(q)).multiply(inv_den));
        
    }

    /**
     * Returns a {@code std::complex<double>} whose value is {@code (this / divisor)}, * with {@code divisor} interpreted as a real number.
     *
     * @param  divisor Value by which this {@code std::complex<double>} is to be divided.
     * @return {@code this / divisor}.
     * @see #divide(Field_Complex<double>)
     */
    Field_Complex<T> divide(const T& divisor) 
    {
        if (is_nan || divisor.is_nan()) 
        {
            return get_nan(get_parts_field());
        }
        if (divisor.is_zero()) 
        {
            return get_nan(get_parts_field());
        }
        if (divisor.get_is_infinite()) 
        {
            return !std::isinfinite()
                ? get_zero(get_parts_field())
                : get_nan(get_parts_field());
        }
        return create_complex(my_real.divide(divisor), my_imaginary.divide(divisor));
    }

    /**
     * Returns a {@code std::complex<double>} whose value is {@code (this / divisor)}, * with {@code divisor} interpreted as a real number.
     *
     * @param  divisor Value by which this {@code std::complex<double>} is to be divided.
     * @return {@code this / divisor}.
     * @see #divide(Field_Complex<double>)
     */
    //override
    Field_Complex<T> divide(const double& divisor) 
    {
        if (my_is_nan || std::isnan(divisor)) 
        {
            return get_nan(get_parts_field());
        }
        if (divisor == 0.0) 
        {
            return get_nan(get_parts_field());
        }
        if (std::isinf(divisor)) 
        {
            return !std::isinfinite()
                ? get_zero(get_parts_field())
                : get_nan(get_parts_field());
        }
        return create_complex(my_real.divide(divisor), my_imaginary.divide(divisor));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> reciprocal() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        if (real.is_zero() && imaginary.is_zero()) 
        {
            return get_inf(get_parts_field());
        }

        if (std::isinfinite) 
        {
            return get_zero(get_parts_field());
        }

        if (std::abs(real).get_real() < std::abs(imaginary).get_real()) 
        {
            T q = real.divide(imaginary);
            T scale = real.multiply(q).add(imaginary).reciprocal();
            return create_complex(scale.multiply(q), scale.negate());
        }
            T q = imaginary.divide(real);
            T scale = imaginary.multiply(q).add(real).reciprocal();
            return create_complex(scale, scale.negate().multiply(q));
    }

    /**
     * Test for equality with another object.
     * If both the real and imaginary parts of two complex numbers
     * are exactly the same, and neither is {@codeNAN}, the two
     * std::complex<double> objects are considered to be equal.
     * The behavior is the same as for JDK's {@link Double#equals(Object)
     * Double}:
     * <ul>
     *  <li>All {@code NaN} values are considered to be equal, *   i.e, if either (or both) real and imaginary parts of the complex
     *   number are equal to {@codeNAN}, the complex number is equal
     *   to {@code NaN}.
     *  </li>
     *  <li>
     *   Instances constructed with different representations of zero (i.e.
     *   either "0" or "-0") are <em>not</em> considered to be equal.
     *  </li>
     * </ul>
     *
     * @param other Object to test for equality with this instance.
     * @return {@code true} if the objects are equal, {@code false} if object
     * is {@code NULL}, not an instance of {@code std::complex<double>}, or not equal to
     * this instance.
     */
    //override
    bool equals(Object other) 
    {
        if (this == other) 
        {
            return true;
        }
        if (dynamic_cast<const Field_Complex<double>*>(*other) != nullptr)
        {
            //@Suppress_Warnings("unchecked")
            Field_Complex<T> c = (Field_Complex<T>) other;
            if (c.is_nan) 
            {
                return my_is_nan;
            }
                return my_real.equals(c.get_real) && imaginary.equals(c.imaginary);
            }
        return false;
    }

    /**
     * Test for the floating-point equality between std::complex<double> objects.
     * It returns {@code true} if both arguments are equal or within the
     * range of allowed error (inclusive).
     *
     * @param x First value (cannot be {@code NULL}).
     * @param y Second value (cannot be {@code NULL}).
     * @param max_ulps {@code (max_ulps - 1)} is the number of floating point
     * values between the real (resp. imaginary) parts of {@code x} and
     * {@code y}.
     * @param <T> the type of the field elements
     * @return {@code true} if there are fewer than {@code max_ulps} floating
     * point values between the real (resp. imaginary) parts of {@code x}
     * and {@code y}.
     *
     * @see Precision#equals(double,double,int)
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static bool equals(Field_Complex<T> x, Field_Complex<T> y, int max_ulps) 
    {
        return Precision::equals(x.real.get_real(), y.real.get_real(), max_ulps) &&
               Precision::equals(x.imaginary.get_real(), y.imaginary.get_real(), max_ulps);
    }

    /**
     * Returns {@code true} iff the values are equal as defined by
     * {@link #equals(Field_Complex<double>,Field_Complex<double>,int) equals(x, y, 1)}.
     *
     * @param x First value (cannot be {@code NULL}).
     * @param y Second value (cannot be {@code NULL}).
     * @param <T> the type of the field elements
     * @return {@code true} if the values are equal.
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static bool equals(Field_Complex<T> x, Field_Complex<T> y) 
    {
        return equals(x, y, 1);
    }

    /**
     * Returns {@code true} if, both for the real part and for the imaginary
     * part, there is no T value strictly between the arguments or the
     * difference between them is within the range of allowed error
     * (inclusive).  Returns {@code false} if either of the arguments is NaN.
     *
     * @param x First value (cannot be {@code NULL}).
     * @param y Second value (cannot be {@code NULL}).
     * @param eps Amount of allowed absolute error.
     * @param <T> the type of the field elements
     * @return {@code true} if the values are two adjacent floating point
     * numbers or they are within range of each other.
     *
     * @see Precision#equals(double,double,double)
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static bool equals(Field_Complex<T> x, Field_Complex<T> y, double eps) 
    {
        return Precision::equals(x.real.get_real(), y.real.get_real(), eps) &&
               Precision::equals(x.imaginary.get_real(), y.imaginary.get_real(), eps);
    }

    /**
     * Returns {@code true} if, both for the real part and for the imaginary
     * part, there is no T value strictly between the arguments or the
     * relative difference between them is smaller or equal to the given
     * tolerance. Returns {@code false} if either of the arguments is NaN.
     *
     * @param x First value (cannot be {@code NULL}).
     * @param y Second value (cannot be {@code NULL}).
     * @param eps Amount of allowed relative error.
     * @param <T> the type of the field elements
     * @return {@code true} if the values are two adjacent floating point
     * numbers or they are within range of each other.
     *
     * @see Precision#equals_with_relative_tolerance(double,double,double)
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static bool equals_with_relative_tolerance(Field_Complex<T> x, Field_Complex<T> y, double eps) 
    {
        return Precision::equals_with_relative_tolerance(x.real.get_real(), y.real.get_real(), eps) &&
               Precision::equals_with_relative_tolerance(x.imaginary.get_real(), y.imaginary.get_real(), eps);
    }

    /**
     * Get a hash_code for the complex number.
     * Any {@codeNAN} value in real or imaginary part produces
     * the same hash code {@code 7}.
     *
     * @return a hash code value for this object.
     */
    //override
    int hash_code() 
    {
        if (my_is_nan) 
        {
            return 7;
        }
        return 37 * (17 * my_imaginary.hash_code() + my_real.hash_code());
    }

    /** {@inherit_doc}
     * <p>
     * This implementation considers +0.0 and -0.0 to be equal for both
     * real and imaginary components.
     * </p>
     */
    //override
    bool is_zero() const
    {
        return real.is_zero() && imaginary.is_zero();
    }

    /**
     * Access the imaginary part.
     *
     * @return the imaginary part.
     */
    T get_imaginary() const
    {
        return my_imaginary;
    }

    /**
     * Access the imaginary part.
     *
     * @return the imaginary part.
     */
    T get_imaginary_part() 
    {
        return my_imaginary;
    }

    /**
     * Access the real part.
     *
     * @return the real part.
     */
    //override
    double get_real() 
    {
        return real.get_real();
    }

    /**
     * Access the real part.
     *
     * @return the real part.
     */
    T get_real_part() 
    {
        return real;
    }

    /**
     * Checks whether either or both parts of this complex number is
     * {@code NaN}.
     *
     * @return true if either or both parts of this complex number is
     * {@code NaN}; false otherwise.
     */
    //override
    bool is_nan() 
    {
        return is_nan;
    }

    /** Check whether the instance is real (i.e. imaginary part is zero).
     * @return true if imaginary part is zero
      */
    bool is_real() 
    {
        return imaginary.is_zero();
    }

    /** Check whether the instance is an integer (i.e. imaginary part is zero and real part has no fractional part).
     * @return true if imaginary part is zero and real part has no fractional part
     */
    bool is_mathematical_integer() 
    {
        return is_real() && Precision.is_mathematical_integer(real.get_real());
    }

    /**
     * Checks whether either the real or imaginary part of this complex number
     * takes an infinite value (either {@code INFINITY} or
     * {@code -INFINITY}) and neither part
     * is {@code NaN}.
     *
     * @return true if one or both parts of this complex number are infinite
     * and neither part is {@code NaN}.
     */
    //override
    bool is_infinite() const
    {
        return my_is_infinite;
    }

    /**
     * Returns a {@code std::complex<double>} whose value is {@code this * factor}.
     * Implements preliminary checks for {@code NaN} and infinity followed by
     * the definitional formula:
     * <p>
     *   {@code (a + bi)(c + di) = (ac - bd) + (ad + bc)i}
     * </p>
     * Returns {@link #get_nan(Field)} if either {@code this} or {@code factor} has one or
     * more {@code NaN} parts.
     * <p>
     * Returns {@link #get_inf(Field)} if neither {@code this} nor {@code factor} has one
     * or more {@code NaN} parts and if either {@code this} or {@code factor}
     * has one or more infinite parts (same result is returned regardless of
     * the sign of the components).
     * </p><p>
     * Returns finite values in components of the result per the definitional
     * formula in all remaining cases.</p>
     *
     * @param  factor value to be multiplied by this {@code std::complex<double>}.
     * @return {@code this * factor}.
     * @ if {@code factor} is {@code NULL}.
     */
    //override
    Field_Complex<T> multiply(Field_Complex<T> factor)
         
        {
        //Math_Utils::check_not_null(factor);
        if (is_nan || factor.is_nan) 
        {
            return get_nan(get_parts_field());
        }
        if (real.std::isinfinite() ||
            imaginary.std::isinfinite() ||
            factor.real.std::isinfinite() ||
            factor.imaginary.std::isinfinite()) 
            {
            // we don't use std::isinfinite() to avoid testing for NaN again
            return get_inf(get_parts_field());
        }
        return create_complex(real.linear_combination(real, factor.real, imaginary.negate(), factor.imaginary), real.linear_combination(real, factor.imaginary, imaginary, factor.real));
    }

    /**
     * Returns a {@code std::complex<double>} whose value is {@code this * factor}, with {@code factor}
     * interpreted as a integer number.
     *
     * @param  factor value to be multiplied by this {@code std::complex<double>}.
     * @return {@code this * factor}.
     * @see #multiply(Field_Complex<double>)
     */
    //override
    Field_Complex<T> multiply(const int factor) 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }
        if (real.std::isinfinite() || imaginary.std::isinfinite()) 
        {
            return get_inf(get_parts_field());
        }
        return create_complex(real.multiply(factor), imaginary.multiply(factor));
    }

    /**
     * Returns a {@code std::complex<double>} whose value is {@code this * factor}, with {@code factor}
     * interpreted as a real number.
     *
     * @param  factor value to be multiplied by this {@code std::complex<double>}.
     * @return {@code this * factor}.
     * @see #multiply(Field_Complex<double>)
     */
    //override
    Field_Complex<T> multiply(double factor) 
    {
        if (is_nan || std::isnan(factor)) 
        {
            return get_nan(get_parts_field());
        }
        if (real.std::isinfinite() ||
            imaginary.std::isinfinite() ||
            std::isinf(factor)) 
            {
            // we don't use std::isinfinite() to avoid testing for NaN again
            return get_inf(get_parts_field());
        }
        return create_complex(real.multiply(factor), imaginary.multiply(factor));
    }

    /**
     * Returns a {@code std::complex<double>} whose value is {@code this * factor}, with {@code factor}
     * interpreted as a real number.
     *
     * @param  factor value to be multiplied by this {@code std::complex<double>}.
     * @return {@code this * factor}.
     * @see #multiply(Field_Complex<double>)
     */
    Field_Complex<T> multiply(T factor) 
    {
        if (is_nan || factor.is_nan()) 
        {
            return get_nan(get_parts_field());
        }
        if (real.std::isinfinite() ||
            imaginary.std::isinfinite() ||
            factor.std::isinfinite()) 
            {
            // we don't use std::isinfinite() to avoid testing for NaN again
            return get_inf(get_parts_field());
        }
        return create_complex(real.multiply(factor), imaginary.multiply(factor));
    }

    /** Compute this * i.
     * @return this * i
     * @since 2.0
     */
    Field_Complex<T> multiply_plus_i() 
    {
        return create_complex(imaginary.negate(), real);
    }

    /** Compute this *- -i.
     * @return this * i
     * @since 2.0
     */
    Field_Complex<T> multiply_minus_i() 
    {
        return create_complex(imaginary, real.negate());
    }

    /**
     * Returns a {@code std::complex<double>} whose value is {@code (-this)}.
     * Returns {@code NaN} if either real or imaginary
     * part of this std::complex<double> number is {@codeNAN}.
     *
     * @return {@code -this}.
     */
    //override
    Field_Complex<T> negate() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        return create_complex(real.negate(), imaginary.negate());
    }

    /**
     * Returns a {@code std::complex<double>} whose value is
     * {@code (this - subtrahend)}.
     * Uses the definitional formula
     * <p>
     *  {@code (a + bi) - (c + di) = (a-c) + (b-d)i}
     * </p>
     * If either {@code this} or {@code subtrahend} has a {@code NaN]} value in either part, * {@link #get_nan(Field)} is returned; otherwise infinite and {@code NaN} values are
     * returned in the parts of the result according to the rules for
     * {@link java.lang.Double} arithmetic.
     *
     * @param  subtrahend value to be subtracted from this {@code std::complex<double>}.
     * @return {@code this - subtrahend}.
     * @ if {@code subtrahend} is {@code NULL}.
     */
    //override
    Field_Complex<T> subtract(Field_Complex<T> subtrahend)
         
        {
        //Math_Utils::check_not_null(subtrahend);
        if (is_nan || subtrahend.is_nan) 
        {
            return get_nan(get_parts_field());
        }

        return create_complex(real.subtract(subtrahend.get_real_part()), imaginary.subtract(subtrahend.get_imaginary_part()));
    }

    /**
     * Returns a {@code std::complex<double>} whose value is
     * {@code (this - subtrahend)}.
     *
     * @param  subtrahend value to be subtracted from this {@code std::complex<double>}.
     * @return {@code this - subtrahend}.
     * @see #subtract(Field_Complex<double>)
     */
    //override
    Field_Complex<T> subtract(double subtrahend) 
    {
        if (is_nan || std::isnan(subtrahend)) 
        {
            return get_nan(get_parts_field());
        }
        return create_complex(real.subtract(subtrahend), imaginary);
    }

    /**
     * Returns a {@code std::complex<double>} whose value is
     * {@code (this - subtrahend)}.
     *
     * @param  subtrahend value to be subtracted from this {@code std::complex<double>}.
     * @return {@code this - subtrahend}.
     * @see #subtract(Field_Complex<double>)
     */
    Field_Complex<T> subtract(T subtrahend) 
    {
        if (is_nan || subtrahend.is_nan()) 
        {
            return get_nan(get_parts_field());
        }
        return create_complex(real.subtract(subtrahend), imaginary);
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/InverseCosine.html" TARGET="_top">
     * inverse cosine</a> of this complex number.
     * Implements the formula:
     * <p>
     *  {@code acos(z) = -i (log(z + i (sqrt(1 - z<sup>2</sup>))))}
     * </p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.
     *
     * @return the inverse cosine of this complex number.
     */
    //override
    Field_Complex<T> acos() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        return this.add(this.sqrt1z().multiply_plus_i()).log().multiply_minus_i();
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/InverseSine.html" TARGET="_top">
     * inverse sine</a> of this complex number.
     * Implements the formula:
     * <p>
     *  {@code asin(z) = -i (log(sqrt(1 - z<sup>2</sup>) + iz))}
     * </p><p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.</p>
     *
     * @return the inverse sine of this complex number.
     */
    //override
    Field_Complex<T> asin() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        return sqrt1z().add(this.multiply_plus_i()).log().multiply_minus_i();
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/InverseTangent.html" TARGET="_top">
     * inverse tangent</a> of this complex number.
     * Implements the formula:
     * <p>
     * {@code atan(z) = (i/2) log((1 - iz)/(1 + iz))}
     * </p><p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.</p>
     *
     * @return the inverse tangent of this complex number
     */
    //override
    Field_Complex<T> atan() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        const T one = get_parts_field().get_one();
        if (real.is_zero()) 
        {

            // singularity at \xc2\xb1i
            if (imaginary.multiply(imaginary).subtract(one).is_zero()) 
            {
                return get_nan(get_parts_field());
            }

            // branch cut on imaginary axis
            const T zero = get_parts_field().get_zero();
            const Field_Complex<T> tmp = create_complex(one.add(imaginary).divide(one.subtract(imaginary)), zero).
                                        log().multiply_plus_i().multiply(0.5);
            return create_complex(std::copysign(tmp.real, real), tmp.imaginary);

        }
else 
        {
            // regular formula
            const Field_Complex<T> n = create_complex(one.add(imaginary), real.negate());
            const Field_Complex<T> d = create_complex(one.subtract(imaginary),  real);
            return n.divide(d).log().multiply_plus_i().multiply(0.5);
        }

    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/Cosine.html" TARGET="_top">
     * cosine</a> of this complex number.
     * Implements the formula:
     * <p>
     *  {@code cos(a + bi) = cos(a)cosh(b) - sin(a)sinh(b)i}
     * </p><p>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos}, * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * </p><p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p><p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.</p>
     * <pre>
     *  Examples:
     *  <code>
     *   cos(1 &plusmn; INFINITY i) = 1 \u2213 INFINITY i
     *   cos(&plusmn;INFINITY + i) = NaN + NaN i
     *   cos(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the cosine of this complex number.
     */
    //override
    Field_Complex<T> cos() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        const Field_Sin_Cos<T>   scr  = Sin_Cos(real);
        const Field_Sinh_Cosh<T> schi = std::sinh_cosh(imaginary);
        return create_complex(scr.cos().multiply(schi.cosh()), scr.sin().negate().multiply(schi.sinh()));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/HyperbolicCosine.html" TARGET="_top">
     * hyperbolic cosine</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   cosh(a + bi) = cosh(a)cos(b) + sinh(a)sin(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos}, * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * <p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   cosh(1 &plusmn; INFINITY i) = NaN + NaN i
     *   cosh(&plusmn;INFINITY + i) = INFINITY &plusmn; INFINITY i
     *   cosh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the hyperbolic cosine of this complex number.
     */
    //override
    Field_Complex<T> cosh() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        const Field_Sinh_Cosh<T> schr = std::sinh_cosh(real);
        const Field_Sin_Cos<T>   sci  = Sin_Cos(imaginary);
        return create_complex(schr.cosh().multiply(sci.cos()), schr.sinh().multiply(sci.sin()));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/ExponentialFunction.html" TARGET="_top">
     * exponential function</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   exp(a + bi) = exp(a)cos(b) + exp(a)sin(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#exp}, {@link FastMath#cos}, and
     * {@link FastMath#sin}.
     * <p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   exp(1 &plusmn; INFINITY i) = NaN + NaN i
     *   exp(INFINITY + i) = INFINITY + INFINITY i
     *   exp(-INFINITY + i) = 0 + 0i
     *   exp(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return <code><i>e</i><sup>this</sup></code>.
     */
    //override
    Field_Complex<T> exp() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        const T              exp_real = std::exp(real);
        const Field_Sin_Cos<T> sc      = Sin_Cos(imaginary);
        return create_complex(exp_real.multiply(sc.cos()), exp_real.multiply(sc.sin()));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> expm1() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        const T              expm1_real = std::expm1(real);
        const Field_Sin_Cos<T> sc        = Sin_Cos(imaginary);
        return create_complex(expm1_real.multiply(sc.cos()), expm1_real.multiply(sc.sin()));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/NaturalLogarithm.html" TARGET="_top">
     * natural logarithm</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   log(a + bi) = ln(|a + bi|) + arg(a + bi)i
     *  </code>
     * </pre>
     * where ln on the right hand side is {@link FastMath#log}, * {@code |a + bi|} is the modulus, {@link #abs},  and
     * {@code arg(a + bi) = }{@link FastMath#atan2}(b, a).
     * <p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite (or critical) values in real or imaginary parts of the input may
     * result in infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   log(1 &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/2)i
     *   log(INFINITY + i) = INFINITY + 0i
     *   log(-INFINITY + i) = INFINITY + &pi;i
     *   log(INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/4)i
     *   log(-INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (3&pi;/4)i
     *   log(0 + 0i) = -INFINITY + 0i
     *  </code>
     * </pre>
     *
     * @return the value <code>ln &nbsp; this</code>, the natural logarithm
     * of {@code this}.
     */
    //override
    Field_Complex<T> log() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        return create_complex(std::log(std::hypot(real, imaginary)), std::atan2(imaginary, real));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> log1p() 
    {
        return add(1.0).log();
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> log10() 
    {
        return log().divide(LOG10);
    }

    /**
     * Returns of value of this complex number raised to the power of {@code x}.
     * <p>
     * If {@code x} is a real number whose real part has an integer value, returns {@link #powstatic_cast<int>(}, * if both {@code this} and {@code x} are real and {@link FastMath#pow(double, double)}
     * with the corresponding real arguments would return a finite number (neither NaN
     * nor infinite), then returns the same value converted to {@code std::complex<double>}, * with the same special cases.
     * In all other cases real cases, : y<sup>x</sup> = exp(x&middot;log(y)).
     * </p>
     *
     * @param  x exponent to which this {@code std::complex<double>} is to be raised.
     * @return <code> this<sup>x</sup></code>.
     * @ if x is {@code NULL}.
     */
    //override
    Field_Complex<T> pow(Field_Complex<T> x)
         
        {

        //Math_Utils::check_not_null(x);

        if (x.imaginary.is_zero()) 
        {
            const int& nx = static_cast<int>( std::rint(x.real.get_real());
            if (x.real.get_real() == nx) 
            {
                // integer power
                return pow(nx);
            }
else if (this.imaginary.is_zero()) 
            {
                // check real implementation that handles a bunch of special cases
                const T real_pow = std::pow(this.real, x.real);
                if (real_pow.is_finite()) 
                {
                    return create_complex(real_pow, get_parts_field().get_zero());
                }
            }
        }

        // generic implementation
        return this.log().multiply(x).exp();

    }


    /**
     * Returns of value of this complex number raised to the power of {@code x}.
     * <p>
     * If {@code x} has an integer value, returns {@link #powstatic_cast<int>(}, * if {@code this} is real and {@link FastMath#pow(double, double)}
     * with the corresponding real arguments would return a finite number (neither NaN
     * nor infinite), then returns the same value converted to {@code std::complex<double>}, * with the same special cases.
     * In all other cases real cases, : y<sup>x</sup> = exp(x&middot;log(y)).
     * </p>
     *
     * @param  x exponent to which this {@code std::complex<double>} is to be raised.
     * @return <code> this<sup>x</sup></code>.
     */
    Field_Complex<T> pow(T x) 
    {

        const int& nx = static_cast<int>( std::rint(x.get_real());
        if (x.get_real() == nx) 
        {
            // integer power
            return pow(nx);
        }
else if (this.imaginary.is_zero()) 
        {
            // check real implementation that handles a bunch of special cases
            const T real_pow = std::pow(this.real, x);
            if (real_pow.is_finite()) 
            {
                return create_complex(real_pow, get_parts_field().get_zero());
            }
        }

        // generic implementation
        return this.log().multiply(x).exp();

    }

    /**
     * Returns of value of this complex number raised to the power of {@code x}.
     * <p>
     * If {@code x} has an integer value, returns {@link #powstatic_cast<int>(}, * if {@code this} is real and {@link FastMath#pow(double, double)}
     * with the corresponding real arguments would return a finite number (neither NaN
     * nor infinite), then returns the same value converted to {@code std::complex<double>}, * with the same special cases.
     * In all other cases real cases, : y<sup>x</sup> = exp(x&middot;log(y)).
     * </p>
     *
     * @param  x exponent to which this {@code std::complex<double>} is to be raised.
     * @return <code> this<sup>x</sup></code>.
     */
    //override
    Field_Complex<T> pow(double x) 
    {

        const int& nx = static_cast<int>( std::rint(x);
        if (x == nx) 
        {
            // integer power
            return pow(nx);
        }
else if (this.imaginary.is_zero()) 
        {
            // check real implementation that handles a bunch of special cases
            const T real_pow = std::pow(this.real, x);
            if (real_pow.is_finite()) 
            {
                return create_complex(real_pow, get_parts_field().get_zero());
            }
        }

        // generic implementation
        return this.log().multiply(x).exp();

    }

     /** {@inherit_doc} */
    //override
    Field_Complex<T> pow(const int& n) 
    {

        Field_Complex<T> result = get_field().get_one();
        const bool invert;
        int p = n;
        if (p < 0) 
        {
            invert = true;
            p = -p;
        }
else 
        {
            invert = false;
        }

        // Exponentiate by successive squaring
        Field_Complex<T> square = this;
        while (p > 0) 
        {
            if ((p & 0x1) > 0) 
            {
                result = result.multiply(square);
            }
            square = square.multiply(square);
            p = p >> 1;
        }

        return invert ? result.reciprocal() : result;

    }

     /**
      * Compute the
     * <a href="http://mathworld.wolfram.com/Sine.html" TARGET="_top">
     * sine</a>
     * of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   sin(a + bi) = sin(a)cosh(b) + cos(a)sinh(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos}, * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * <p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p><p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or {@code NaN} values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   sin(1 &plusmn; INFINITY i) = 1 &plusmn; INFINITY i
     *   sin(&plusmn;INFINITY + i) = NaN + NaN i
     *   sin(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the sine of this complex number.
     */
    //override
    Field_Complex<T> sin() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        const Field_Sin_Cos<T>   scr  = Sin_Cos(real);
        const Field_Sinh_Cosh<T> schi = std::sinh_cosh(imaginary);
        return create_complex(scr.sin().multiply(schi.cosh()), scr.cos().multiply(schi.sinh()));

    }

    /** {@inherit_doc}
     */
    //override
    Field_Sin_Cos<Field_Complex<T>> sin_cos() 
    {
        if (is_nan) 
        {
            return Field_Sin_Cos<>(get_nan(get_parts_field()), get_nan(get_parts_field()));
        }

        const Field_Sin_Cos<T>   scr = Sin_Cos(real);
        const Field_Sinh_Cosh<T> schi = std::sinh_cosh(imaginary);
        return Field_Sin_Cos<>(create_complex(scr.sin().multiply(schi.cosh()), scr.cos().multiply(schi.sinh())), create_complex(scr.cos().multiply(schi.cosh()), scr.sin().negate().multiply(schi.sinh())));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> atan2(Field_Complex<T> x) 
    {

        // compute r = sqrt(x^2+y^2)
        const Field_Complex<T> r = x.multiply(x).add(multiply(this)).sqrt();

        if (x.real.get_real() >= 0) 
        {
            // compute atan2(y, x) = 2 atan(y / (r + x))
            return divide(r.add(x)).atan().multiply(2);
        }
else 
        {
            // compute atan2(y, x) = +/- pi - 2 atan(y / (r - x))
            return divide(r.subtract(x)).atan().multiply(-2).add(x.real.get_pi());
        }
    }

    /** {@inherit_doc}
     * <p>
     * Branch cuts are on the real axis, below +1.
     * </p>
     */
    //override
    Field_Complex<T> acosh() 
    {
        const Field_Complex<T> sqrt_plus  = add(1).sqrt();
        const Field_Complex<T> sqrt_minus = subtract(1).sqrt();
        return add(sqrt_plus.multiply(sqrt_minus)).log();
    }

    /** {@inherit_doc}
     * <p>
     * Branch cuts are on the imaginary axis, above +i and below -i.
     * </p>
     */
    //override
    Field_Complex<T> asinh() 
    {
        return add(multiply(this).add(1.0).sqrt()).log();
    }

    /** {@inherit_doc}
     * <p>
     * Branch cuts are on the real axis, above +1 and below -1.
     * </p>
     */
    //override
    Field_Complex<T> atanh() 
    {
        const Field_Complex<T> log_plus  = add(1).log();
        const Field_Complex<T> log_minus = create_complex(get_parts_field().get_one().subtract(real), imaginary.negate()).log();
        return log_plus.subtract(log_minus).multiply(0.5);
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/HyperbolicSine.html" TARGET="_top">
     * hyperbolic sine</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   sinh(a + bi) = sinh(a)cos(b)) + cosh(a)sin(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos}, * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * <p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p><p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   sinh(1 &plusmn; INFINITY i) = NaN + NaN i
     *   sinh(&plusmn;INFINITY + i) = &plusmn; INFINITY + INFINITY i
     *   sinh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the hyperbolic sine of {@code this}.
     */
    //override
    Field_Complex<T> sinh() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        const Field_Sinh_Cosh<T> schr = std::sinh_cosh(real);
        const Field_Sin_Cos<T>   sci  = Sin_Cos(imaginary);
        return create_complex(schr.sinh().multiply(sci.cos()), schr.cosh().multiply(sci.sin()));
    }

    /** {@inherit_doc}
     */
    //override
    Field_Sinh_Cosh<Field_Complex<T>> sinh_cosh() 
    {
        if (is_nan) 
        {
            return Field_Sinh_Cosh<>(get_nan(get_parts_field()), get_nan(get_parts_field()));
        }

        const Field_Sinh_Cosh<T> schr = std::sinh_cosh(real);
        const Field_Sin_Cos<T>   sci  = Sin_Cos(imaginary);
        return Field_Sinh_Cosh<>(create_complex(schr.sinh().multiply(sci.cos()), schr.cosh().multiply(sci.sin())), create_complex(schr.cosh().multiply(sci.cos()), schr.sinh().multiply(sci.sin())));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
     * square root</a> of this complex number.
     * Implements the following algorithm to compute {@code sqrt(a + bi)}:
     * <ol><li>Let {@code t = sqrt((|a| + |a + bi|) / 2)}</li>
     * <li><pre>if {@code  a \xe2\x89\xa5 0} return {@code t + (b/2t)i}
     *  else return {@code |b|/2t + sign(b)t i }</pre></li>
     * </ol>
     * where <ul>
     * <li>{@code |a| = }{@link FastMath#abs(Calculus_Field_Element) abs(a)}</li>
     * <li>{@code |a + bi| = }{@link FastMath#hypot(Calculus_Field_Element, Calculus_Field_Element) hypot(a, b)}</li>
     * <li>{@code sign(b) = }{@link FastMath#copy_sign(Calculus_Field_Element, Calculus_Field_Element) copy_sign(1, b)}
     * </ul>
     * The real part is therefore always nonnegative.
     * <p>
     * Returns {@link #get_nan(Field) NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * <p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * </p>
     * <pre>
     *  Examples:
     *  <code>
     *   sqrt(1 \xc2\xb1 \xe2\x88\x9e i) = \xe2\x88\x9e + NaN i
     *   sqrt(\xe2\x88\x9e + i) = \xe2\x88\x9e + 0i
     *   sqrt(-\xe2\x88\x9e + i) = 0 + \xe2\x88\x9e i
     *   sqrt(\xe2\x88\x9e \xc2\xb1 \xe2\x88\x9e i) = \xe2\x88\x9e + NaN i
     *   sqrt(-\xe2\x88\x9e \xc2\xb1 \xe2\x88\x9e i) = NaN \xc2\xb1 \xe2\x88\x9e i
     *  </code>
     * </pre>
     *
     * @return the square root of {@code this} with nonnegative real part.
     */
    //override
    Field_Complex<T> sqrt() 
    {
        if (is_nan) 
        {
            return get_nan(get_parts_field());
        }

        if (is_zero()) 
        {
            return get_zero(get_parts_field());
        }

        T t = std::sqrt((std::abs(real).add(std::hypot(real, imaginary))).multiply(0.5));
        if (real.get_real() >= 0.0) 
        {
            return create_complex(t, imaginary.divide(t.multiply(2)));
        }
else 
        {
            return create_complex(std::abs(imaginary).divide(t.multiply(2)), std::copysign(t, imaginary));
        }
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
     * square root</a> of <code>1 - this<sup>2</sup></code> for this complex
     * number.
     * Computes the result directly as
     * {@code sqrt(ONE.subtract(z.multiply(z)))}.
     * <p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     *
     * @return the square root of <code>1 - this<sup>2</sup></code>.
     */
    Field_Complex<T> sqrt1z() 
    {
        const Field_Complex<T> t2 = this.multiply(this);
        return create_complex(get_parts_field().get_one().subtract(t2.real), t2.imaginary.negate()).sqrt();
    }

    /** {@inherit_doc}
     * <p>
     * This implementation compute the principal cube root by using a branch cut along real negative axis.
     * </p>
     */
    //override
    Field_Complex<T> cbrt() 
    {
        const T              magnitude = std::cbrt(abs().get_real_part());
        const Field_Sin_Cos<T> sc        = Sin_Cos(get_argument().divide(3));
        return create_complex(magnitude.multiply(sc.cos()), magnitude.multiply(sc.sin()));
    }

    /** {@inherit_doc}
     * <p>
     * This implementation compute the principal n<sup>th</sup> root by using a branch cut along real negative axis.
     * </p>
     */
    //override
    Field_Complex<T> root_n(const int& n) 
    {
        const T              magnitude = std::pow(abs().get_real_part(), 1.0 / n);
        const Field_Sin_Cos<T> sc        = Sin_Cos(get_argument().divide(n));
        return create_complex(magnitude.multiply(sc.cos()), magnitude.multiply(sc.sin()));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/Tangent.html" TARGET="_top">
     * tangent</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   tan(a + bi) = sin(2a)/(cos(2a)+cosh(2b)) + [sinh(2b)/(cos(2a)+cosh(2b))]i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos}, {@link FastMath#cosh} and
     * {@link FastMath#sinh}.
     * <p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite (or critical) values in real or imaginary parts of the input may
     * result in infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   tan(a &plusmn; INFINITY i) = 0 &plusmn; i
     *   tan(&plusmn;INFINITY + bi) = NaN + NaN i
     *   tan(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *   tan(&plusmn;&pi;/2 + 0 i) = &plusmn;INFINITY + NaN i
     *  </code>
     * </pre>
     *
     * @return the tangent of {@code this}.
     */
    //override
    Field_Complex<T> tan() 
    {
        if (is_nan || real.std::isinfinite()) 
        {
            return get_nan(get_parts_field());
        }
        if (imaginary.get_real() > 20.0) 
        {
            return get_i(get_parts_field());
        }
        if (imaginary.get_real() < -20.0) 
        {
            return get_minus_i(get_parts_field());
        }

        const Field_Sin_Cos<T> sc2r = Sin_Cos(real.multiply(2));
        T imaginary2 = imaginary.multiply(2);
        T d = sc2r.cos().add(std::cosh(imaginary2));

        return create_complex(sc2r.sin().divide(d), std::sinh(imaginary2).divide(d));

    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/HyperbolicTangent.html" TARGET="_top">
     * hyperbolic tangent</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   tan(a + bi) = sinh(2a)/(cosh(2a)+cos(2b)) + [sin(2b)/(cosh(2a)+cos(2b))]i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos}, {@link FastMath#cosh} and
     * {@link FastMath#sinh}.
     * <p>
     * Returns {@link #get_nan(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   tanh(a &plusmn; INFINITY i) = NaN + NaN i
     *   tanh(&plusmn;INFINITY + bi) = &plusmn;1 + 0 i
     *   tanh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *   tanh(0 + (&pi;/2)i) = NaN + INFINITY i
     *  </code>
     * </pre>
     *
     * @return the hyperbolic tangent of {@code this}.
     */
    //override
    Field_Complex<T> tanh() 
    {
        if (is_nan || imaginary.std::isinfinite()) 
        {
            return get_nan(get_parts_field());
        }
        if (real.get_real() > 20.0) 
        {
            return get_one(get_parts_field());
        }
        if (real.get_real() < -20.0) 
        {
            return get_minus_one(get_parts_field());
        }
        T real2 = real.multiply(2);
        const Field_Sin_Cos<T> sc2i = Sin_Cos(imaginary.multiply(2));
        T d = std::cosh(real2).add(sc2i.cos());

        return create_complex(std::sinh(real2).divide(d), sc2i.sin().divide(d));
    }



    /**
     * Compute the argument of this complex number.
     * The argument is the angle phi between the positive real axis and
     * the point representing this number in the complex plane.
     * The value returned is between -PI (not inclusive)
     * and PI (inclusive), with negative values returned for numbers with
     * negative imaginary parts.
     * <p>
     * If either real or imaginary part (or both) is NaN, NaN is returned.
     * Infinite parts are handled as {@code Math.atan2} handles them, * essentially treating finite parts as zero in the presence of an
     * infinite coordinate and returning a multiple of pi/4 depending on
     * the signs of the infinite parts.
     * See the javadoc for {@code Math.atan2} for full details.
     *
     * @return the argument of {@code this}.
     */
    T get_argument() 
    {
        return std::atan2(get_imaginary_part(), get_real_part());
    }

    /**
     * Computes the n-th roots of this complex number.
     * The nth roots are defined by the formula:
     * <pre>
     *  <code>
     *   z<sub>k</sub> = abs<sup>1/n</sup> (cos(phi + 2&pi;k/n) + i (sin(phi + 2&pi;k/n))
     *  </code>
     * </pre>
     * for <i>{@code k=0, 1, ..., n-1}</i>, where {@code abs} and {@code phi}
     * are respectively the {@link #abs() modulus} and
     * {@link #get_argument() argument} of this complex number.
     * <p>
     * If one or both parts of this complex number is NaN, a list with just
     * one element, {@link #get_nan(Field)} is returned.
     * if neither part is NaN, but at least one part is infinite, the result
     * is a one-element list containing {@link #get_inf(Field)}.
     *
     * @param n Degree of root.
     * @return a List of all {@code n}-th roots of {@code this}.
     * @ if {@code n <= 0}.
     */
    List<Field_Complex<T>> nth_root(const int& n)  
    {

        if (n <= 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::CANNOT_COMPUTE_NTH_ROOT_FOR_NEGATIVE_N, n);
        }

        const List<Field_Complex<T>> result = Array_list<>();

        if (is_nan) 
        {
            result.add(get_nan(get_parts_field()));
            return result;
        }
        if (std::isinfinite()) 
        {
            result.add(get_inf(get_parts_field()));
            return result;
        }

        // nth root of abs -- faster / more accurate to use a solver here?
        const T nth_root_of_abs = std::pow(std::hypot(real, imaginary), 1.0 / n);

        // Compute nth roots of complex number with k = 0, 1, ... n-1
        const T nth_phi = get_argument().divide(n);
        const double slice = 2 * std::numbers::pi / n;
        T inner_part = nth_phi;
        for (int k{}; k < n ; k++) 
        {
            // inner part
            const Field_Sin_Cos<T> sc_inner = Sin_Cos(inner_part);
            const T real_part = nth_root_of_abs.multiply(sc_inner.cos());
            const T imaginary_part = nth_root_of_abs.multiply(sc_inner.sin());
            result.add(create_complex(real_part, imaginary_part));
            inner_part = inner_part.add(slice);
        }

        return result;
    }



    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param real_part Real part.
     * @param imaginary_part Imaginary part.
     * @param <T> the type of the field elements
     * @return a std::complex<double> instance.
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex<T>
        value_of(T real_part, T imaginary_part) 
        {
        if (real_part.is_nan() || imaginary_part.is_nan()) 
        {
            return get_nan(real_part.get_field());
        }
        return Field_Complex<double><>(real_part, imaginary_part);
    }

    /**
     * Create a complex number given only the real part.
     *
     * @param real_part Real part.
     * @param <T> the type of the field elements
     * @return a std::complex<double> instance.
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Complex<T>
        value_of(T real_part) 
        {
        if (real_part.is_nan()) 
        {
            return get_nan(real_part.get_field());
        }
        return Field_Complex<double><>(real_part);
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> new_instance(double real_part) 
    {
        return value_of(get_parts_field().get_zero().new_instance(real_part));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<double>_Field<T> get_field() 
    {
        return Field_Complex<double>_Field.get_field(get_parts_field());
    }

    /** Get the {@link Field} the real and imaginary parts belong to.
     * @return {@link Field} the real and imaginary parts belong to
     */
    Field<T> get_parts_field() 
    {
        return real.get_field();
    }

    /** {@inherit_doc} */
    //override
    std::string to_string() const 
    {
        return "(" + real + ", " + imaginary + ")";
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> scalb(const int& n) 
    {
        return create_complex(std::scalbn(real, n), std::scalbn(imaginary, n));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> ulp() 
    {
        return create_complex(FastMath.ulp(real), FastMath.ulp(imaginary));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> hypot(Field_Complex<T> y) 
    {
        if (std::isinfinite() || y.std::isinfinite()) 
        {
            return get_inf(get_parts_field());
        }
else if (is_nan() || y.is_nan()) 
        {
            return get_nan(get_parts_field());
        }
else 
        {
            return multiply(this).add(y.multiply(y)).sqrt();
        }
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> linear_combination(const Field_Complex<T>[] a, const Field_Complex<T>[] b)
         
        {
        const int n = 2 * a.size();
        const std::vector<T> real_a      = Math_Arrays::build_array(get_parts_field(), n);
        const std::vector<T> real_b      = Math_Arrays::build_array(get_parts_field(), n);
        const std::vector<T> imaginary_a = Math_Arrays::build_array(get_parts_field(), n);
        const std::vector<T> imaginary_b = Math_Arrays::build_array(get_parts_field(), n);
        for (int i{}; i < a.size(); ++i)  
        {
            const Field_Complex<T> ai = a[i];
            const Field_Complex<T> bi = b[i];
            real_a[2 * i    ]      = ai.real;
            real_a[2 * i + 1]      = ai.imaginary.negate();
            real_b[2 * i    ]      = bi.real;
            real_b[2 * i + 1]      = bi.imaginary;
            imaginary_a[2 * i    ] = ai.real;
            imaginary_a[2 * i + 1] = ai.imaginary;
            imaginary_b[2 * i    ] = bi.imaginary;
            imaginary_b[2 * i + 1] = bi.real;
        }
        return create_complex(real.linear_combination(real_a,  real_b), real.linear_combination(imaginary_a, imaginary_b));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> linear_combination(const std::vector<double> a, const Field_Complex<T>[] b)
         
        {
        const int n = a.size();
        const std::vector<T> real_b      = Math_Arrays::build_array(get_parts_field(), n);
        const std::vector<T> imaginary_b = Math_Arrays::build_array(get_parts_field(), n);
        for (int i{}; i < a.size(); ++i)  
        {
            const Field_Complex<T> bi = b[i];
            real_b[i]      = bi.real;
            imaginary_b[i] = bi.imaginary;
        }
        return create_complex(real.linear_combination(a,  real_b), real.linear_combination(a, imaginary_b));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> linear_combination(const Field_Complex<T> a1, const Field_Complex<T> b1, const Field_Complex<T> a2, const Field_Complex<T> b2) 
    {
        return create_complex(real.linear_combination(a1.real, b1.real, a1.imaginary.negate(), b1.imaginary, a2.real, b2.real, a2.imaginary.negate(), b2.imaginary), real.linear_combination(a1.real, b1.imaginary, a1.imaginary, b1.real, a2.real, b2.imaginary, a2.imaginary, b2.real));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> linear_combination(const double& a1, const Field_Complex<T> b1, const double& a2, const Field_Complex<T> b2) 
    {
        return create_complex(real.linear_combination(a1, b1.real, a2, b2.real), real.linear_combination(a1, b1.imaginary, a2, b2.imaginary));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> linear_combination(const Field_Complex<T> a1, const Field_Complex<T> b1, const Field_Complex<T> a2, const Field_Complex<T> b2, const Field_Complex<T> a3, const Field_Complex<T> b3) 
    {
        Field_Complex<T>[] a = Math_Arrays::build_array(get_field(), 3);
        a[0] = a1;
        a[1] = a2;
        a[2] = a3;
        Field_Complex<T>[] b = Math_Arrays::build_array(get_field(), 3);
        b[0] = b1;
        b[1] = b2;
        b[2] = b3;
        return linear_combination(a, b);
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> linear_combination(const double& a1, const Field_Complex<T> b1, const double& a2, const Field_Complex<T> b2, const double& a3, const Field_Complex<T> b3) 
    {
        Field_Complex<T>[] b = Math_Arrays::build_array(get_field(), 3);
        b[0] = b1;
        b[1] = b2;
        b[2] = b3;
        return linear_combination(std::vector<double>  { a1, a2, a3 }, b);
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> linear_combination(const Field_Complex<T> a1, const Field_Complex<T> b1, const Field_Complex<T> a2, const Field_Complex<T> b2, const Field_Complex<T> a3, const Field_Complex<T> b3, const Field_Complex<T> a4, const Field_Complex<T> b4) 
    {
        Field_Complex<T>[] a = Math_Arrays::build_array(get_field(), 4);
        a[0] = a1;
        a[1] = a2;
        a[2] = a3;
        a[3] = a4;
        Field_Complex<T>[] b = Math_Arrays::build_array(get_field(), 4);
        b[0] = b1;
        b[1] = b2;
        b[2] = b3;
        b[3] = b4;
        return linear_combination(a, b);
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> linear_combination(const double& a1, const Field_Complex<T> b1, const double& a2, const Field_Complex<T> b2, const double& a3, const Field_Complex<T> b3, const double& a4, const Field_Complex<T> b4) 
    {
        Field_Complex<T>[] b = Math_Arrays::build_array(get_field(), 4);
        b[0] = b1;
        b[1] = b2;
        b[2] = b3;
        b[3] = b4;
        return linear_combination(std::vector<double>  { a1, a2, a3, a4 }, b);
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> ceil() 
    {
        return create_complex(std::ceil(get_real_part()), std::ceil(get_imaginary_part()));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> floor() 
    {
        return create_complex(std::floor(get_real_part()), std::floor(get_imaginary_part()));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> rint() 
    {
        return create_complex(std::rint(get_real_part()), std::rint(get_imaginary_part()));
    }

    /** {@inherit_doc}
     * <p>
     * for complex numbers, the integer n corresponding to {@code this.subtract(remainder(a)).divide(a)}
     * is a <a href="https://en.wikipedia.org/wiki/Gaussian_integer">Wikipedia - Gaussian integer</a>.
     * </p>
     */
    //override
    Field_Complex<T> remainder(const double& a) 
    {
        return create_complex(std::remainder(get_real_part(), a), std::remainder(get_imaginary_part(), a));
    }

    /** {@inherit_doc}
     * <p>
     * for complex numbers, the integer n corresponding to {@code this.subtract(remainder(a)).divide(a)}
     * is a <a href="https://en.wikipedia.org/wiki/Gaussian_integer">Wikipedia - Gaussian integer</a>.
     * </p>
     */
    //override
    Field_Complex<T> remainder(const Field_Complex<T> a) 
    {
        const Field_Complex<T> complex_quotient = divide(a);
        const T  q_r_int           = std::rint(complex_quotient.real);
        const T  q_i_int           = std::rint(complex_quotient.imaginary);
        return create_complex(real.subtract(q_r_int.multiply(a.real)).add(q_i_int.multiply(a.imaginary)), imaginary.subtract(q_r_int.multiply(a.imaginary)).subtract(q_i_int.multiply(a.real)));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> sign() 
    {
        if (is_nan() || is_zero()) 
        {
            return this;
        }
else 
        {
            return this.divide(std::hypot(real, imaginary));
        }
    }

    /** {@inherit_doc}
     * <p>
     * The signs of real and imaginary parts are copied independently.
     * </p>
     */
    //override
    Field_Complex<T> copy_sign(const Field_Complex<T> z) 
    {
        return create_complex(std::copysign(get_real_part(), z.get_real_part()), std::copysign(get_imaginary_part(), z.get_imaginary_part()));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> copy_sign(double r) 
    {
        return create_complex(std::copysign(get_real_part(), r), std::copysign(get_imaginary_part(), r));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> to_degrees() 
    {
        return create_complex(FastMath.to_degrees(get_real_part()), FastMath.to_degrees(get_imaginary_part()));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> to_radians() 
    {
        return create_complex(FastMath.to_radians(get_real_part()), FastMath.to_radians(get_imaginary_part()));
    }

    /** {@inherit_doc} */
    //override
    Field_Complex<T> get_pi() 
    {
        return get_pi(get_parts_field());
    }

}


