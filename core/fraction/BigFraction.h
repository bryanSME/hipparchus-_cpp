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
//package org.hipparchus.fraction;

//import java.io.Serializable;
//import java.math.BigDecimal;
//import java.math.BigInteger;
//import java.util.function.Function;
//import java.util.stream.Stream;

//import org.hipparchus.Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.fraction.Convergents_Iterator.Convergence_Step;
//import org.hipparchus.util.Arithmetic_Utils;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Pair;
//import org.hipparchus.util.Precision;
#include "../util/BigInteger.h"
#include "../"
/**
 * Representation of a rational number without any overflow. This class is
 * immutable.
 *
 */
class Big_Fraction : public Field_Element<Big_Fraction>, public Comparable<Big_Fraction>
{
public:
    /** A fraction representing "2 / 1". */
    static const Big_Fraction TWO = Big_Fraction(2);

    /** A fraction representing "1". */
    static const Big_Fraction ONE = Big_Fraction(1);

    /** A fraction representing "0". */
    static const Big_Fraction ZERO = Big_Fraction(0);

    /** A fraction representing "-1 / 1". */
    static const Big_Fraction MINUS_ONE = Big_Fraction(-1);

    /** A fraction representing "4/5". */
    static const Big_Fraction FOUR_FIFTHS = Big_Fraction(4, 5);

    /** A fraction representing "1/5". */
    static const Big_Fraction ONE_FIFTH = Big_Fraction(1, 5);

    /** A fraction representing "1/2". */
    static const Big_Fraction ONE_HALF = Big_Fraction(1, 2);

    /** A fraction representing "1/4". */
    static const Big_Fraction ONE_QUARTER = Big_Fraction(1, 4);

    /** A fraction representing "1/3". */
    static const Big_Fraction ONE_THIRD = Big_Fraction(1, 3);

    /** A fraction representing "3/5". */
    static const Big_Fraction THREE_FIFTHS = Big_Fraction(3, 5);

    /** A fraction representing "3/4". */
    static const Big_Fraction THREE_QUARTERS = Big_Fraction(3, 4);

    /** A fraction representing "2/5". */
    static const Big_Fraction TWO_FIFTHS = Big_Fraction(2, 5);

    /** A fraction representing "2/4". */
    static const Big_Fraction TWO_QUARTERS = Big_Fraction(2, 4);

    /** A fraction representing "2/3". */
    static const Big_Fraction TWO_THIRDS = Big_Fraction(2, 3);

    
    static const long serial_version_uid = -5630213147331578515L;

    /** <code>BigInteger</code> representation of 100. */
    static const BigInteger ONE_HUNDRED = Bigstatic_cast<int>(100);

    /** Convert a convergence step to the corresponding double fraction. */
    static const Function<Convergence_Step, Big_Fraction> STEP_TO_FRACTION = //
            s -> Big_Fraction(s.get_numerator(), s.get_denominator());

private:
    /** The numerator. */
    const BigInteger my_numerator;

    /** The denominator. */
    const BigInteger my_denominator;

public:
    /**
     * <p>
     * Create a {@link Big_Fraction} equivalent to the passed {@code BigInteger}, ie
     * "num / 1".
     * </p>
     *
     * @param num
     *            the numerator.
     */
    Big_Fraction(const BigInteger& num) 
    {
        this(num, BigInteger::ONE);
    }

    /**
     * Create a {@link Big_Fraction} given the numerator and denominator as
     * {@code BigInteger}. The {@link Big_Fraction} is reduced to lowest terms.
     *
     * @param num the numerator, must not be {@code NULL}.
     * @param den the denominator, must not be {@code NULL}.
     * @ if the denominator is zero.
     * @Null_Argument_Exception if either of the arguments is NULL
     */
    Big_Fraction(BigInteger num, BigInteger den) 
    {
        //Math_Utils::check_not_null(num, hipparchus::exception::Localized_Core_Formats_Type::NUMERATOR);
        //Math_Utils::check_not_null(den, hipparchus::exception::Localized_Core_Formats_Type::DENOMINATOR);
        if (den.signum() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR);
        }
        if (num.signum() == 0) 
        {
            numerator   = BigInteger.ZERO;
            denominator = BigInteger.ONE;
        }
        else 
        {

            // reduce numerator and denominator by greatest common denominator
            const BigInteger gcd = num.gcd(den);
            if (BigInteger.ONE.compare_to(gcd) < 0) 
            {
                num = num.divide(gcd);
                den = den.divide(gcd);
            }

            // move sign to numerator
            if (den.signum() == -1) 
            {
                num = num.negate();
                den = den.negate();
            }

            // store the values in the const fields
            numerator   = num;
            denominator = den;

        }
    }

    /**
     * Create a fraction given the double value.
     * <p>
     * This constructor behaves <em>differently</em> from
     * {@link #Big_Fraction(double, double, int)}. It converts the double value
     * exactly, considering its internal bits representation. This works for all
     * values except NaN and infinities and does not requires any loop or
     * convergence threshold.
     * </p>
     * <p>
     * sin_ce this conversion is exact and since double numbers are sometimes
     * approximated, the fraction created may seem strange in some cases. For example, * calling <code>new Big_Fraction(1.0 / 3.0)</code> does <em>not</em> create
     * the fraction 1/3, but the fraction 6004799503160661 / 18014398509481984
     * because the double number passed to the constructor is not exactly 1/3
     * (this number cannot be stored exactly in IEEE754).
     * </p>
     * @see #Big_Fraction(double, double, int)
     * @param value the double value to convert to a fraction.
     * @exception  if value is NaN or infinite
     */
    Big_Fraction(const double value)  
    {
        if (std::isnan(value)) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NAN_VALUE_CONVERSION);
        }
        if (std::isinf(value)) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INFINITE_VALUE_CONVERSION);
        }

        // compute m and k such that value = m * 2^k
        const long bits     = Double.double_to_long_bits(value);
        const long sign     = bits & 0x8000000000000000L;
        const long exponent = bits & 0x7ff0000000000000L;
        long m              = bits & 0x000fffffffffffffL;
        if (exponent != 0) 
        {
            // this was a normalized number, add the implicit most significant bit
            m |= 0x0010000000000000L;
        }
        if (sign != 0) 
        {
            m = -m;
        }
        int k = (static_cast<int>( (exponent >> 52)) - 1075;
        while (((m & 0x001ffffffffffffeL) != 0) && ((m & 0x1) == 0)) 
        {
            m >>= 1;
            ++k;
        }

        if (k < 0) 
        {
            numerator   = Bigstatic_cast<int>(m);
            denominator = BigInteger.ZERO.flip_bit(-k);
        }
else 
        {
            numerator   = Bigstatic_cast<int>(m).multiply(BigInteger.ZERO.flip_bit(k));
            denominator = BigInteger.ONE;
        }

    }

    /**
     * Create a fraction given the double value and maximum error allowed.
     * <p>
     * References:
     * <ul>
     * <li><a href="http://mathworld.wolfram.com/Continued_Fraction.html">
     * Continued Fraction</a> equations (11) and (22)-(26)</li>
     * </ul>
     * </p>
     *
     * @param value
     *            the double value to convert to a fraction.
     * @param epsilon
     *            maximum error allowed. The resulting fraction is within
     *            <code>epsilon</code> of <code>value</code>, in absolute terms.
     * @param max_iterations
     *            maximum number of convergents.
     * @Math_Illegal_State_Exception
     *             if the continued fraction failed to converge.
     * @see #Big_Fractionstatic_cast<double>(
     */
    Big_Fraction(const double& value, const double epsilon, const int max_iterations)
        Math_Illegal_State_Exception 
        {
        Convergence_Step converged = Convergents_Iterator.convergent(value, max_iterations, s -> 
        {
            const double quotient = s.get_fraction_value();
            return Precision.equals(quotient, value, 1) || std::abs(quotient - value) < epsilon;
        }).get_key();
        if (std::abs(converged.get_fraction_value() - value) < epsilon) 
        {
            this.numerator = Bigstatic_cast<int>(converged.get_numerator());
            this.denominator = Bigstatic_cast<int>(converged.get_denominator());
        }
else 
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FAILED_FRACTION_CONVERSION, value, max_iterations);
        }
    }

    /**
     * Create a fraction given the double value and maximum denominator.
     * <p>
     * References:
     * <ul>
     * <li><a href="http://mathworld.wolfram.com/Continued_Fraction.html">
     * Continued Fraction</a> equations (11) and (22)-(26)</li>
     * </ul>
     * </p>
     *
     * @param value
     *            the double value to convert to a fraction.
     * @param max_denominator
     *            The maximum allowed value for denominator.
     * @Math_Illegal_State_Exception
     *             if the continued fraction failed to converge.
     */
    Big_Fraction(const double& value, const long max_denominator)
        Math_Illegal_State_Exception 
        {
        const int max_iterations = 100;
        Convergence_Step[] last_valid = Convergence_Step[1];
        Convergents_Iterator.convergent(value, max_iterations, s -> 
        {
            if (s.get_denominator() < max_denominator) 
            {
                last_valid[0] = s;
            }
            return Precision.equals(s.get_fraction_value(), value, 1);
        });
        if (last_valid[0] != NULL) 
        {
            this.numerator   = Bigstatic_cast<int>(last_valid[0].get_numerator());
            this.denominator = Bigstatic_cast<int>(last_valid[0].get_denominator());
        }
else 
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FAILED_FRACTION_CONVERSION, value, max_iterations);
        }
    }

    /**
     * <p>
     * Create a {@link Big_Fraction} equivalent to the passed {@code int}, ie
     * "num / 1".
     * </p>
     *
     * @param num
     *            the numerator.
     */
    Big_Fraction(const int& num) 
    {
        this(Bigstatic_cast<int>(num), BigInteger.ONE);
    }

    /**
     * <p>
     * Create a {@link Big_Fraction} given the numerator and denominator as simple
     * {@code int}. The {@link Big_Fraction} is reduced to lowest terms.
     * </p>
     *
     * @param num
     *            the numerator.
     * @param den
     *            the denominator.
     */
    Big_Fraction(const int& num, const int den) 
    {
        this(Bigstatic_cast<int>(num), Bigstatic_cast<int>(den));
    }

    /**
     * <p>
     * Create a {@link Big_Fraction} equivalent to the passed long, ie "num / 1".
     * </p>
     *
     * @param num
     *            the numerator.
     */
    Big_Fraction(const long num) 
    {
        this(Bigstatic_cast<int>(num), BigInteger.ONE);
    }

    /**
     * <p>
     * Create a {@link Big_Fraction} given the numerator and denominator as simple
     * {@code long}. The {@link Big_Fraction} is reduced to lowest terms.
     * </p>
     *
     * @param num
     *            the numerator.
     * @param den
     *            the denominator.
     */
    Big_Fraction(const long num, const long den) 
    {
        this(Bigstatic_cast<int>(num), Bigstatic_cast<int>(den));
    }

    /**
     * A test to determine if a series of fractions has converged.
     */
    //@Functional_Interface
    class Convergence_Test 
    {
        /**
         * Evaluates if the fraction formed by {@code numerator/denominator} satisfies
         * this convergence test.
         *
         * @param numerator   the numerator
         * @param denominator the denominator
         * @return if this convergence test is satisfied
         */
        bool test(const long& numerator, const long& denominator);
    }

    /** Generate a {@link Stream stream} of convergents from a real number.
     * @param value value to approximate
     * @param max_convergents maximum number of convergents.
     * @return stream of {@link Big_Fraction} convergents approximating  {@code value}
     * @since 2.1
     */
    static Stream<Big_Fraction> convergents(const double& value, const int& max_convergents) 
    {
        return Convergents_Iterator.convergents(value, max_convergents).map(STEP_TO_FRACTION);
    }

    /**
     * Returns the last element of the series of convergent-steps to approximate the
     * given value.
     * <p>
     * The series terminates either at the first step that satisfies the given
     * {@code convergence_test} or after at most {@code max_convergents} elements. The
     * returned Pair consists of that terminal {@link Big_Fraction} and a
     * {@link Boolean} that indicates if it satisfies the given convergence tests.
     * If the returned pair's value is {@code false} the element at position
     * {@code max_convergents} was examined but failed to satisfy the
     * {@code convergence_test}. A caller can then decide to accept the result
     * nevertheless or to discard it. This method is usually faster than
     * {@link #convergents(double, int)} if only the terminal element is of
     * interest.
     *
     * @param value           value to approximate
     * @param max_convergents  maximum number of convergents to examine
     * @param convergence_test the test if the series has converged at a step
     * @return the pair of last element of the series of convergents and a bool
     *         indicating if that element satisfies the specified convergent test
     */
    static Pair<Big_Fraction, Boolean> convergent(const double& value, const int& max_convergents, const Convergence_Test& convergence_test) 
    {
        std::pair<Convergence_Step, Boolean> converged = Convergents_Iterator.convergent(value, max_convergents, s -> convergence_test.test(s.get_numerator(), s.get_denominator()));
        return Pair.create(STEP_TO_FRACTION.apply(converged.get_key()), converged.get_value());
    }

    /** {@inherit_doc} */
    //override
    double get_real() 
    {
        return double_value();
    }

    /**
     * <p>
     * Creates a <code>Big_Fraction</code> instance with the 2 parts of a fraction
     * Y/Z.
     * </p>
     *
     * <p>
     * Any negative signs are resolved to be on the numerator.
     * </p>
     *
     * @param numerator
     *            the numerator, for example the three in 'three sevenths'.
     * @param denominator
     *            the denominator, for example the seven in 'three sevenths'.
     * @return a fraction instance, with the numerator and denominator
     *         reduced.
     * @Arithmetic_Exception
     *             if the denominator is <code>zero</code>.
     */
    static Big_Fraction get_reduced_fraction(const int& numerator, const int& denominator) 
    {
        if (numerator == 0) 
        {
            return ZERO; // normalize zero.
        }

        return Big_Fraction(numerator, denominator);
    }

    /**
     * <p>
     * Returns the absolute value of this {@link Big_Fraction}.
     * </p>
     *
     * @return the absolute value as a {@link Big_Fraction}.
     */
    Big_Fraction abs() 
    {
        return (numerator.signum() == 1)
            ? *this
            : negate();
    }

    /** Check if a fraction is an integer.
     * @return true of fraction is an integer
     */
    bool is_integer() 
    {
        return denominator.equals(BigInteger.ONE);
    }

    /** Returns the signum function of this {@link Big_Fraction}.
     * <p>
     * The return value is -1 if the specified value is negative;
     * 0 if the specified value is zero; and 1 if the specified value is positive.
     * </p>
     * @return the signum function of this {@link Big_Fraction}
     * @since 1.7
     */
    int signum() 
    {
        return numerator.signum();
    }

    /**
     * <p>
     * Adds the value of this fraction to the passed {@link BigInteger}, * returning the result in reduced form.
     * </p>
     *
     * @param bg
     *            the {@link BigInteger} to add, must'nt be <code>null</code>.
     * @return a <code>Big_Fraction</code> instance with the resulting values.
     * @Null_Argument_Exception
     *             if the {@link BigInteger} is <code>null</code>.
     */
    Big_Fraction add(const BigInteger bg) Null_Argument_Exception 
    {
        //Math_Utils::check_not_null(bg);

        if (numerator.signum() == 0) 
        {
            return Big_Fraction(bg);
        }
        if (bg.signum() == 0) 
        {
            return this;
        }

        return Big_Fraction(numerator.add(denominator.multiply(bg)), denominator);
    }

    /**
     * <p>
     * Adds the value of this fraction to the passed {@code integer}, returning
     * the result in reduced form.
     * </p>
     *
     * @param i
     *            the {@code integer} to add.
     * @return a <code>Big_Fraction</code> instance with the resulting values.
     */
    Big_Fraction add(const int& i) 
    {
        return add(Bigstatic_cast<int>(i));
    }

    /**
     * <p>
     * Adds the value of this fraction to the passed {@code long}, returning
     * the result in reduced form.
     * </p>
     *
     * @param l
     *            the {@code long} to add.
     * @return a <code>Big_Fraction</code> instance with the resulting values.
     */
    Big_Fraction add(const long l) 
    {
        return add(Bigstatic_cast<int>(l));
    }

    /**
     * <p>
     * Adds the value of this fraction to another, returning the result in
     * reduced form.
     * </p>
     *
     * @param fraction
     *            the {@link Big_Fraction} to add, must not be <code>null</code>.
     * @return a {@link Big_Fraction} instance with the resulting values.
     * @Null_Argument_Exception if the {@link Big_Fraction} is {@code NULL}.
     */
    //override
    Big_Fraction add(const Big_Fraction& fraction) 
    {
        //Math_Utils::check_not_null(fraction, hipparchus::exception::Localized_Core_Formats_Type::FRACTION);
        if (fraction.numerator.signum() == 0) 
        {
            return this;
        }
        if (numerator.signum() == 0) 
        {
            return fraction;
        }

        BigInteger num;
        BigInteger den;
        if (denominator.equals(fraction.denominator)) 
        {
            num = numerator.add(fraction.numerator);
            den = denominator;
        }
else 
        {
            num = (numerator.multiply(fraction.denominator)).add((fraction.numerator).multiply(denominator));
            den = denominator.multiply(fraction.denominator);
        }

        if (num.signum() == 0) 
        {
            return ZERO;
        }

        return Big_Fraction(num, den);

    }

    /**
     * <p>
     * Gets the fraction as a <code>BigDecimal</code>. This calculates the
     * fraction as the numerator divided by denominator.
     * </p>
     *
     * @return the fraction as a <code>BigDecimal</code>.
     * @Arithmetic_Exception
     *             if the exact quotient does not have a terminating decimal
     *             expansion.
     * @see BigDecimal
     */
    BigDecimal big_decimal_value() 
    {
        return BigDecimal(numerator).divide(new BigDecimal(denominator));
    }

    /**
     * <p>
     * Gets the fraction as a <code>BigDecimal</code> following the passed
     * rounding mode. This calculates the fraction as the numerator divided by
     * denominator.
     * </p>
     *
     * @param rounding_mode
     *            rounding mode to apply. see {@link BigDecimal} constants.
     * @return the fraction as a <code>BigDecimal</code>.
     * @Illegal_Argument_Exception
     *             if {@code rounding_mode} does not represent a valid rounding
     *             mode.
     * @see BigDecimal
     */
    BigDecimal big_decimal_value(const int rounding_mode) 
    {
        return BigDecimal(numerator).divide(new BigDecimal(denominator), rounding_mode);
    }

    /**
     * <p>
     * Gets the fraction as a <code>BigDecimal</code> following the passed scale
     * and rounding mode. This calculates the fraction as the numerator divided
     * by denominator.
     * </p>
     *
     * @param scale
     *            scale of the <code>BigDecimal</code> quotient to be returned.
     *            see {@link BigDecimal} for more information.
     * @param rounding_mode
     *            rounding mode to apply. see {@link BigDecimal} constants.
     * @return the fraction as a <code>BigDecimal</code>.
     * @see BigDecimal
     */
    BigDecimal big_decimal_value(const int scale, const int rounding_mode) 
    {
        return BigDecimal(numerator).divide(new BigDecimal(denominator), scale, rounding_mode);
    }

    /**
     * <p>
     * Compares this object to another based on size.
     * </p>
     *
     * @param object
     *            the object to compare to, must not be <code>null</code>.
     * @return -1 if this is less than {@code object}, +1 if this is greater
     *         than {@code object}, 0 if they are equal.
     * @see java.lang.Comparable#compare_to(java.lang.Object)
     */
    //override
    int compare_to(const Big_Fraction object) 
    {
        int lhs_sig_num = numerator.signum();
        int rhs_sig_num = object.numerator.signum();

        if (lhs_sig_num != rhs_sig_num) 
        {
            return (lhs_sig_num > rhs_sig_num) ? 1 : -1;
        }
        if (lhs_sig_num == 0) 
        {
            return 0;
        }

        BigInteger nOd = numerator.multiply(object.denominator);
        BigInteger dOn = denominator.multiply(object.numerator);
        return nOd.compare_to(dOn);
    }

    /**
     * <p>
     * Divide the value of this fraction by the passed {@code BigInteger}, * ie {@code this * 1 / bg}, returning the result in reduced form.
     * </p>
     *
     * @param bg the {@code BigInteger} to divide by, must not be {@code NULL}
     * @return a {@link Big_Fraction} instance with the resulting values
     * @Null_Argument_Exception if the {@code BigInteger} is {@code NULL}
     * @Math_Runtime_Exception if the fraction to divide by is zero
     */
    Big_Fraction divide(const BigInteger bg) 
    {
        //Math_Utils::check_not_null(bg);
        if (bg.signum() == 0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR);
        }
        if (numerator.signum() == 0) 
        {
            return ZERO;
        }
        return Big_Fraction(numerator, denominator.multiply(bg));
    }

    /**
     * <p>
     * Divide the value of this fraction by the passed {@code int}, ie
     * {@code this * 1 / i}, returning the result in reduced form.
     * </p>
     *
     * @param i the {@code int} to divide by
     * @return a {@link Big_Fraction} instance with the resulting values
     * @Math_Runtime_Exception if the fraction to divide by is zero
     */
    Big_Fraction divide(const int& i) 
    {
        return divide(Bigstatic_cast<int>(i));
    }

    /**
     * <p>
     * Divide the value of this fraction by the passed {@code long}, ie
     * {@code this * 1 / l}, returning the result in reduced form.
     * </p>
     *
     * @param l the {@code long} to divide by
     * @return a {@link Big_Fraction} instance with the resulting values
     * @Math_Runtime_Exception if the fraction to divide by is zero
     */
    Big_Fraction divide(const long l) 
    {
        return divide(Bigstatic_cast<int>(l));
    }

    /**
     * <p>
     * Divide the value of this fraction by another, returning the result in
     * reduced form.
     * </p>
     *
     * @param fraction Fraction to divide by, must not be {@code NULL}.
     * @return a {@link Big_Fraction} instance with the resulting values.
     * @Null_Argument_Exception if the {@code fraction} is {@code NULL}.
     * @Math_Runtime_Exception if the fraction to divide by is zero
     */
    //override
    Big_Fraction divide(const Big_Fraction& fraction) 
    {
        //Math_Utils::check_not_null(fraction, hipparchus::exception::Localized_Core_Formats_Type::FRACTION);
        if (fraction.numerator.signum() == 0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR);
        }
        if (numerator.signum() == 0) 
        {
            return ZERO;
        }

        return multiply(fraction.reciprocal());
    }

    /**
     * <p>
     * Gets the fraction as a {@code double}. This calculates the fraction as
     * the numerator divided by denominator.
     * </p>
     *
     * @return the fraction as a {@code double}
     * @see java.lang.Number#double_value()
     */
    //override
    double double_value() 
    {
        double result = numerator.double_value() / denominator.double_value();
        if (std::isnan(result)) 
        {
            // Numerator and/or denominator must be out of range:
            // Calculate how far to shift them to put them in range.
            int shift = std::max(numerator.bit_length(), denominator.bit_length()) - FastMath.get_exponent(Double.MAX_VALUE);
            result = numerator.shift_right(shift).double_value() /
                denominator.shift_right(shift).double_value();
        }
        return result;
    }

    /**
     * <p>
     * Test for the equality of two fractions. If the lowest term numerator and
     * denominators are the same for both fractions, the two fractions are
     * considered to be equal.
     * </p>
     *
     * @param other
     *            fraction to test for equality to this fraction, can be
     *            <code>null</code>.
     * @return true if two fractions are equal, false if object is
     *         <code>null</code>, not an instance of {@link Big_Fraction}, or not
     *         equal to this fraction instance.
     * @see java.lang.Object#equals(java.lang.Object)
     */
    //override
    bool equals(const Object& other) 
    {
        bool ret = false;

        if (this == other) 
        {
            ret = true;
        }
        if (dynamic_cast<const Big_Fraction*>(*other) != nullptr)
        {
            Big_Fraction rhs = (Big_Fraction) other;
            ret = numerator.equals(rhs.numerator) && denominator.equals(rhs.denominator);
        }

        return ret;
    }

    /**
     * <p>
     * Gets the fraction as a {@code float}. This calculates the fraction as
     * the numerator divided by denominator.
     * </p>
     *
     * @return the fraction as a {@code float}.
     * @see java.lang.Number#float_value()
     */
    //override
    float float_value() 
    {
        float result = numerator.float_value() / denominator.float_value();
        if (std::isnan(result)) 
        {
            // Numerator and/or denominator must be out of range:
            // Calculate how far to shift them to put them in range.
            int shift = std::max(numerator.bit_length(), denominator.bit_length()) - FastMath.get_exponent(Float.MAX_VALUE);
            result = numerator.shift_right(shift).float_value() /
                denominator.shift_right(shift).float_value();
        }
        return result;
    }

    /**
     * <p>
     * Access the denominator as a <code>BigInteger</code>.
     * </p>
     *
     * @return the denominator as a <code>BigInteger</code>.
     */
    BigInteger get_denominator() 
    {
        return denominator;
    }

    /**
     * <p>
     * Access the denominator as a {@code int}.
     * </p>
     *
     * @return the denominator as a {@code int}.
     */
    int get_denominator_as_int() 
    {
        return denominator.int_value();
    }

    /**
     * <p>
     * Access the denominator as a {@code long}.
     * </p>
     *
     * @return the denominator as a {@code long}.
     */
    long get_denominator_as_long() 
    {
        return denominator.long_value();
    }

    /**
     * <p>
     * Access the numerator as a <code>BigInteger</code>.
     * </p>
     *
     * @return the numerator as a <code>BigInteger</code>.
     */
    BigInteger get_numerator() 
    {
        return numerator;
    }

    /**
     * <p>
     * Access the numerator as a {@code int}.
     * </p>
     *
     * @return the numerator as a {@code int}.
     */
    int get_numerator_as_int() 
    {
        return numerator.int_value();
    }

    /**
     * <p>
     * Access the numerator as a {@code long}.
     * </p>
     *
     * @return the numerator as a {@code long}.
     */
    long get_numerator_as_long() 
    {
        return numerator.long_value();
    }

    /**
     * <p>
     * Gets a hash_code for the fraction.
     * </p>
     *
     * @return a hash code value for this object.
     * @see java.lang.Object#hash_code()
     */
    //override
    int hash_code() 
    {
        return 37 * (37 * 17 + numerator.hash_code()) + denominator.hash_code();
    }

    /**
     * <p>
     * Gets the fraction as an {@code int}. This returns the whole number part
     * of the fraction.
     * </p>
     *
     * @return the whole number fraction part.
     * @see java.lang.Number#int_value()
     */
    //override
    int int_value() 
    {
        return numerator.divide(denominator).int_value();
    }

    /**
     * <p>
     * Gets the fraction as a {@code long}. This returns the whole number part
     * of the fraction.
     * </p>
     *
     * @return the whole number fraction part.
     * @see java.lang.Number#long_value()
     */
    //override
    long long_value() 
    {
        return numerator.divide(denominator).long_value();
    }

    /**
     * <p>
     * Multiplies the value of this fraction by the passed
     * <code>BigInteger</code>, returning the result in reduced form.
     * </p>
     *
     * @param bg the {@code BigInteger} to multiply by.
     * @return a {@code Big_Fraction} instance with the resulting values.
     * @Null_Argument_Exception if {@code bg} is {@code NULL}.
     */
    Big_Fraction multiply(const BigInteger bg) 
    {
        //Math_Utils::check_not_null(bg);
        if (numerator.signum() == 0 || bg.signum() == 0) 
        {
            return ZERO;
        }
        return Big_Fraction(bg.multiply(numerator), denominator);
    }

    /**
     * <p>
     * Multiply the value of this fraction by the passed {@code int}, returning
     * the result in reduced form.
     * </p>
     *
     * @param i
     *            the {@code int} to multiply by.
     * @return a {@link Big_Fraction} instance with the resulting values.
     */
    //override
    Big_Fraction multiply(const int& i) 
    {
        if (i == 0 || numerator.signum() == 0) 
        {
            return ZERO;
        }

        return multiply(Bigstatic_cast<int>(i));
    }

    /**
     * <p>
     * Multiply the value of this fraction by the passed {@code long}, * returning the result in reduced form.
     * </p>
     *
     * @param l
     *            the {@code long} to multiply by.
     * @return a {@link Big_Fraction} instance with the resulting values.
     */
    Big_Fraction multiply(const long l) 
    {
        if (l == 0 || numerator.signum() == 0) 
        {
            return ZERO;
        }

        return multiply(Bigstatic_cast<int>(l));
    }

    /**
     * <p>
     * Multiplies the value of this fraction by another, returning the result in
     * reduced form.
     * </p>
     *
     * @param fraction Fraction to multiply by, must not be {@code NULL}.
     * @return a {@link Big_Fraction} instance with the resulting values.
     * @Null_Argument_Exception if {@code fraction} is {@code NULL}.
     */
    //override
    Big_Fraction multiply(const Big_Fraction& fraction) 
    {
        //Math_Utils::check_not_null(fraction, hipparchus::exception::Localized_Core_Formats_Type::FRACTION);
        if (numerator.signum() == 0 ||
            fraction.numerator.signum() == 0) 
            {
            return ZERO;
        }
        return Big_Fraction(numerator.multiply(fraction.numerator), denominator.multiply(fraction.denominator));
    }

    /**
     * <p>
     * Return the additive inverse of this fraction, returning the result in
     * reduced form.
     * </p>
     *
     * @return the negation of this fraction.
     */
    //override
    Big_Fraction negate() 
    {
        return Big_Fraction(numerator.negate(), denominator);
    }

    /**
     * <p>
     * Gets the fraction percentage as a {@code double}. This calculates the
     * fraction as the numerator divided by denominator multiplied by 100.
     * </p>
     *
     * @return the fraction percentage as a {@code double}.
     */
    double percentage_value() 
    {
        return multiply(ONE_HUNDRED).double_value();
    }

    /**
     * <p>
     * Returns a {@code Big_Fraction} whose value is
     * {@code (this<sup>exponent</sup>)}, returning the result in reduced form.
     * </p>
     *
     * @param exponent
     *            exponent to which this {@code Big_Fraction} is to be
     *            raised.
     * @return <tt>this<sup>exponent</sup></tt>.
     */
    Big_Fraction pow(const int exponent) 
    {
        if (exponent == 0) 
        {
            return ONE;
        }
        if (numerator.signum() == 0) 
        {
            return this;
        }

        if (exponent < 0) 
        {
            return Big_Fraction(denominator.pow(-exponent), numerator.pow(-exponent));
        }
        return Big_Fraction(numerator.pow(exponent), denominator.pow(exponent));
    }

    /**
     * <p>
     * Returns a <code>Big_Fraction</code> whose value is
     * <tt>(this<sup>exponent</sup>)</tt>, returning the result in reduced form.
     * </p>
     *
     * @param exponent
     *            exponent to which this <code>Big_Fraction</code> is to be raised.
     * @return <tt>this<sup>exponent</sup></tt> as a <code>Big_Fraction</code>.
     */
    Big_Fraction pow(const long exponent) 
    {
        if (exponent == 0) 
        {
            return ONE;
        }
        if (numerator.signum() == 0) 
        {
            return this;
        }

        if (exponent < 0) 
        {
            return Big_Fraction(Arithmetic_Utils.pow(denominator, -exponent), Arithmetic_Utils.pow(numerator,   -exponent));
        }
        return Big_Fraction(Arithmetic_Utils.pow(numerator,   exponent), Arithmetic_Utils.pow(denominator, exponent));
    }

    /**
     * <p>
     * Returns a <code>Big_Fraction</code> whose value is
     * <tt>(this<sup>exponent</sup>)</tt>, returning the result in reduced form.
     * </p>
     *
     * @param exponent
     *            exponent to which this <code>Big_Fraction</code> is to be raised.
     * @return <tt>this<sup>exponent</sup></tt> as a <code>Big_Fraction</code>.
     */
    Big_Fraction pow(const BigInteger exponent) 
    {
        if (exponent.signum() == 0) 
        {
            return ONE;
        }
        if (numerator.signum() == 0) 
        {
            return this;
        }

        if (exponent.signum() == -1) 
        {
            const BigInteger e_neg = exponent.negate();
            return Big_Fraction(Arithmetic_Utils.pow(denominator, e_neg), Arithmetic_Utils.pow(numerator,   e_neg));
        }
        return Big_Fraction(Arithmetic_Utils.pow(numerator,   exponent), Arithmetic_Utils.pow(denominator, exponent));
    }

    /**
     * <p>
     * Returns a <code>double</code> whose value is
     * <tt>(this<sup>exponent</sup>)</tt>, returning the result in reduced form.
     * </p>
     *
     * @param exponent
     *            exponent to which this <code>Big_Fraction</code> is to be raised.
     * @return <tt>this<sup>exponent</sup></tt>.
     */
    double pow(const double exponent) 
    {
        return std::pow(numerator.double_value(),   exponent) /
               std::pow(denominator.double_value(), exponent);
    }

    /**
     * <p>
     * Return the multiplicative inverse of this fraction.
     * </p>
     *
     * @return the reciprocal fraction.
     */
    //override
    Big_Fraction reciprocal() 
    {
        return Big_Fraction(denominator, numerator);
    }

    /**
     * <p>
     * Reduce this <code>Big_Fraction</code> to its lowest terms.
     * </p>
     *
     * @return the reduced <code>Big_Fraction</code>. It doesn't change anything if
     *         the fraction can be reduced.
     */
    Big_Fraction reduce() 
    {
        const BigInteger gcd = numerator.gcd(denominator);

        if (BigInteger.ONE.compare_to(gcd) < 0) 
        {
            return Big_Fraction(numerator.divide(gcd), denominator.divide(gcd));
        }
else 
        {
            return this;
        }
    }

    /**
     * <p>
     * Subtracts the value of an {@link BigInteger} from the value of this
     * {@code Big_Fraction}, returning the result in reduced form.
     * </p>
     *
     * @param bg the {@link BigInteger} to subtract, cannot be {@code NULL}.
     * @return a {@code Big_Fraction} instance with the resulting values.
     * @Null_Argument_Exception if the {@link BigInteger} is {@code NULL}.
     */
    Big_Fraction subtract(const BigInteger bg) 
    {
        //Math_Utils::check_not_null(bg);
        if (bg.signum() == 0) 
        {
            return this;
        }
        if (numerator.signum() == 0) 
        {
            return Big_Fraction(bg.negate());
        }

        return Big_Fraction(numerator.subtract(denominator.multiply(bg)), denominator);
    }

    /**
     * <p>
     * Subtracts the value of an {@code integer} from the value of this
     * {@code Big_Fraction}, returning the result in reduced form.
     * </p>
     *
     * @param i the {@code integer} to subtract.
     * @return a {@code Big_Fraction} instance with the resulting values.
     */
    Big_Fraction subtract(const int& i) 
    {
        return subtract(Bigstatic_cast<int>(i));
    }

    /**
     * <p>
     * Subtracts the value of a {@code long} from the value of this
     * {@code Big_Fraction}, returning the result in reduced form.
     * </p>
     *
     * @param l the {@code long} to subtract.
     * @return a {@code Big_Fraction} instance with the resulting values.
     */
    Big_Fraction subtract(const long& l) 
    {
        return subtract(Bigstatic_cast<int>(l));
    }

    /**
     * <p>
     * Subtracts the value of another fraction from the value of this one, * returning the result in reduced form.
     * </p>
     *
     * @param fraction {@link Big_Fraction} to subtract, must not be {@code NULL}.
     * @return a {@link Big_Fraction} instance with the resulting values
     * @Null_Argument_Exception if the {@code fraction} is {@code NULL}.
     */
    //override
    Big_Fraction subtract(const Big_Fraction& fraction) 
    {
        //Math_Utils::check_not_null(fraction, hipparchus::exception::Localized_Core_Formats_Type::FRACTION);
        if (fraction.numerator.signum() == 0) 
        {
            return this;
        }
        if (numerator.signum() == 0) 
        {
            return fraction.negate();
        }

        BigInteger num;
        BigInteger den;
        if (denominator.equals(fraction.denominator)) 
        {
            num = numerator.subtract(fraction.numerator);
            den = denominator;
        }
else 
        {
            num = (numerator.multiply(fraction.denominator)).subtract((fraction.numerator).multiply(denominator));
            den = denominator.multiply(fraction.denominator);
        }
        return Big_Fraction(num, den);

    }

    /**
     * <p>
     * Returns the <code>std::string</code> representing this fraction, ie
     * "num / dem" or just "num" if the denominator is one.
     * </p>
     *
     * @return a string representation of the fraction.
     * @see java.lang.Object#to_string()
     */
    //override
    std::string to_string() const 
    {
        if (BigInteger.ONE.equals(denominator)) 
        {
            return numerator.to_string();
        }
else if (BigInteger.ZERO.equals(numerator)) 
        {
            return "0";
        }
else 
        {
            return numerator + " / " + denominator;
        }
    }

    /** {@inherit_doc} */
    //override
    Big_Fraction_Field get_field() 
    {
        return Big_Fraction_Field.get_instance();
    }

}


