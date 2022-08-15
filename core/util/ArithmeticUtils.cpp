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

#include <limits>
#include <algorithm>
#include "BigInteger.h"
#include "ArithmeticUtils.h"

int Arithmetic_Utils::add_and_check(const int& x, const int& y)
{
    const long s = static_cast<long>(x) + static_cast<long>(y);

    if (s < std::numeric_limits<int>::min() || s > std::numeric_limits<int>::max())
    {
        throw Math_Runtime_Exception(Localized_Core_Formats.OVERFLOW_IN_ADDITION, x, y);
    }
    return static_cast<int>(s);
}

/**
* Add two long integers, checking for overflow.
*
* @param a an addend
* @param b an addend
* @return the sum {@code a+b}
* @Math_Runtime_Exception if the result can not be represented as an long
*/
long Arithmetic_Utils::add_and_check(const long& a, const long& b)
{
    return add_and_check(a, b, Localized_Core_Formats.OVERFLOW_IN_ADDITION);
}

/**
* Computes the greatest common divisor of the absolute value of two
* numbers, using a modified version of the "binary gcd" method.
* See Knuth 4.5.2 algorithm B.
* The algorithm is due to Josef Stein (1961).
* <br/>
* Special cases:
* <ul>
*  <li>The invocations
*   {@code gcd(std::numeric_limits<int>::min(), std::numeric_limits<int>::min())}, *   {@code gcd(std::numeric_limits<int>::min(), 0)} and
*   {@code gcd(0, std::numeric_limits<int>::min())} throw an
*   {@code Arithmetic_Exception}, because the result would be 2^31, which
*   is too large for an int value.</li>
*  <li>The result of {@code gcd(x, x)}, {@code gcd(0, x)} and
*   {@code gcd(x, 0)} is the absolute value of {@code x}, except
*   for the special cases above.</li>
*  <li>The invocation {@code gcd(0, 0)} is the only one which returns
*   {@code 0}.</li>
* </ul>
*
* @param p Number.
* @param q Number.
* @return the greatest common divisor (never negative).
* @Math_Runtime_Exception if the result cannot be represented as
* a non-negative {@code int} value.
*/
int Arithmetic_Utils::gcd(const int& p, const int& q)
{
    int a = p;
    int b = q;
    if (a == 0 ||
        b == 0) 
        {
        if (a == std::numeric_limits<int>::min() ||
            b == std::numeric_limits<int>::min())
            {
            throw Math_Runtime_Exception(Localized_Core_Formats.GCD_OVERFLOW_32_BITS, p, q);
        }
        return std::abs(a + b);
    }

    long al = a;
    long bl = b;
    bool use_long = false;
    if (a < 0) 
    {
        if(std::numeric_limits<int>::min() == a) 
        {
            use_long = true;
        }
        else 
        {
            a = -a;
        }
        al = -al;
    }
    if (b < 0) 
    {
        if (std::numeric_limits<int>::min() == b) 
        {
            use_long = true;
        }
        else 
        {
            b = -b;
        }
        bl = -bl;
    }
    if (use_long) 
    {
        if(al == bl) 
        {
            throw Math_Runtime_Exception(Localized_Core_Formats.GCD_OVERFLOW_32_BITS, p, q);
        }
        long blbu = bl;
        bl = al;
        al = blbu % al;
        if (al == 0) 
        {
            if (bl > std::numeric_limits<int>::max()) 
            {
                throw Math_Runtime_Exception(Localized_Core_Formats.GCD_OVERFLOW_32_BITS, p, q);
            }
            return static_cast<int>( bl;
        }
        blbu = bl;

        // Now "al" and "bl" fit in an "int".
        b = static_cast<int>( al;
        a = static_cast<int>( (blbu % al);
    }

    return gcd_positive(a, b);
}

/**
* Computes the greatest common divisor of two <em>positive</em> numbers
* (this precondition is <em>not</em> checked and the result is undefined
* if not fulfilled) using the "binary gcd" method which avoids division
* and modulo operations.
* See Knuth 4.5.2 algorithm B.
* The algorithm is due to Josef Stein (1961).
* <p>
* Special cases:
* <ul>
*  <li>The result of {@code gcd(x, x)}, {@code gcd(0, x)} and
*   {@code gcd(x, 0)} is the value of {@code x}.</li>
*  <li>The invocation {@code gcd(0, 0)} is the only one which returns
*   {@code 0}.</li>
* </ul>
*
* @param a Positive number.
* @param b Positive number.
* @return the greatest common divisor.
*/
int Arithmetic_Utils::gcd_positive(const int& a, const int& b)
{
    if (a == 0) 
    {
        return b;
    }
    else if (b == 0) 
    {
        return a;
    }

    // Make "a" and "b" odd, keeping track of common power of 2.
    const int a_twos = Integer.number_of_trailing_zeros(a);
    a >>= a_twos;
    const int b_twos = Integer.number_of_trailing_zeros(b);
    b >>= b_twos;
    const int shift = std::min(a_twos, b_twos);

    // "a" and "b" are positive.
    // If a > b then "gdc(a, b)" is equal to "gcd(a - b, b)".
    // If a < b then "gcd(a, b)" is equal to "gcd(b - a, a)".
    // Hence, in the successive iterations:
    //  "a" becomes the absolute difference of the current values, //  "b" becomes the minimum of the current values.
    while (a != b) 
    {
        const int delta = a - b;
        b = std::min(a, b);
        a = std::abs(delta);

        // Remove any power of 2 in "a" ("b" is guaranteed to be odd).
        a >>= Integer.number_of_trailing_zeros(a);
    }

    // Recover the common power of 2.
    return a << shift;
}

/**
* Gets the greatest common divisor of the absolute value of two numbers, * using the "binary gcd" method which avoids division and modulo
* operations. See Knuth 4.5.2 algorithm B. This algorithm is due to Josef
* Stein (1961).
* <p>
* Special cases:
* <ul>
* <li>The invocations
* {@code gcd(std::numeric_limits<long>::min(), std::numeric_limits<long>::min())}, * {@code gcd(std::numeric_limits<long>::min(), 0L)} and
* {@code gcd(0L, std::numeric_limits<long>::min())} throw an
* {@code Arithmetic_Exception}, because the result would be 2^63, which
* is too large for a long value.</li>
* <li>The result of {@code gcd(x, x)}, {@code gcd(0L, x)} and
* {@code gcd(x, 0L)} is the absolute value of {@code x}, except
* for the special cases above.
* <li>The invocation {@code gcd(0L, 0L)} is the only one which returns
* {@code 0L}.</li>
* </ul>
*
* @param p Number.
* @param q Number.
* @return the greatest common divisor, never negative.
* @Math_Runtime_Exception if the result cannot be represented as
* a non-negative {@code long} value.
*/
long Arithmetic_Utils::gcd(const long& p, const long& q)
{
    long u = p;
    long v = q;
    if ((u == 0) || (v == 0)) 
    {
        if ((u == std::numeric_limits<long>::min()) || (v == std::numeric_limits<long>::min()))
        {
            throw Math_Runtime_Exception(Localized_Core_Formats.GCD_OVERFLOW_64_BITS, p, q);
        }
        return std::abs(u) + std::abs(v);
    }
    // keep u and v negative, as negative integers range down to
    // -2^63, while positive numbers can only be as large as 2^63-1
    // (i.e. we can't necessarily negate a negative number without
    // overflow)
    /* assert u!=0 && v!=0; */
    if (u > 0) 
    {
        u = -u;
    } // make u negative
    if (v > 0) 
    {
        v = -v;
    } // make v negative
    // B1. [Find power of 2]
    int k = 0;
    while ((u & 1) == 0 && (v & 1) == 0 && k < 63) { // while u and v are
                                                        // both even...
        u /= 2;
        v /= 2;
        k++; // cast out twos.
    }
    if (k == 63) 
    {
        throw Math_Runtime_Exception(Localized_Core_Formats.GCD_OVERFLOW_64_BITS, p, q);
    }
    // B2. Initialize: u and v have been divided by 2^k and at least
    // one is odd.
    long t = ((u & 1) == 1) ? v : -(u / 2)/* B3 */;
    // t negative: u was odd, v may be even (t replaces v)
    // t positive: u was even, v is odd (t replaces u)
    do 
    {
        /* assert u<0 && v<0; */
        // B4/B3: cast out twos from t.
        while ((t & 1) == 0) { // while t is even..
            t /= 2; // cast out twos
        }
        // B5 [reset max(u,v)]
        if (t > 0) 
        {
            u = -t;
        }
        else 
        {
            v = t;
        }
        // B6/B3. at this point both u and v should be odd.
        t = (v - u) / 2;
        // |u| larger: t positive (replace u)
        // |v| larger: t negative (replace v)
    } while (t != 0);
    return -u * (1L << k); // gcd is u*2^k
}

/**
* Returns the least common multiple of the absolute value of two numbers, * using the formula {@code lcm(a,b) = (a / gcd(a,b)) * b}.
* <p>
* Special cases:
* <ul>
* <li>The invocations {@code lcm(std::numeric_limits<int>::min(), n)} and
* {@code lcm(n, std::numeric_limits<int>::min())}, where {@code abs(n)} is a
* power of 2, throw an {@code Arithmetic_Exception}, because the result
* would be 2^31, which is too large for an int value.</li>
* <li>The result of {@code lcm(0, x)} and {@code lcm(x, 0)} is
* {@code 0} for any {@code x}.
* </ul>
*
* @param a Number.
* @param b Number.
* @return the least common multiple, never negative.
* @Math_Runtime_Exception if the result cannot be represented as
* a non-negative {@code int} value.
*/
int Arithmetic_Utils::lcm(const int& a, int b)
{
    if (a == 0 || b == 0)
    {
        return 0;
    }
    int lcm = std::abs(Arithmetic_Utils.mul_and_check(a / gcd(a, b), b));
    if (lcm == std::numeric_limits<int>::min()) 
    {
        throw Math_Runtime_Exception(Localized_Core_Formats.LCM_OVERFLOW_32_BITS, a, b);
    }
    return lcm;
}

/**
    * Returns the least common multiple of the absolute value of two numbers, * using the formula {@code lcm(a,b) = (a / gcd(a,b)) * b}.
    * <p>
    * Special cases:
    * <ul>
    * <li>The invocations {@code lcm(std::numeric_limits<long>::min(), n)} and
    * {@code lcm(n, std::numeric_limits<long>::min())}, where {@code abs(n)} is a
    * power of 2, throw an {@code Arithmetic_Exception}, because the result
    * would be 2^63, which is too large for an int value.</li>
    * <li>The result of {@code lcm(0L, x)} and {@code lcm(x, 0L)} is
    * {@code 0L} for any {@code x}.
    * </ul>
    *
    * @param a Number.
    * @param b Number.
    * @return the least common multiple, never negative.
    * @Math_Runtime_Exception if the result cannot be represented
    * as a non-negative {@code long} value.
    */
long lcm(long a, long b) Math_Runtime_Exception 
{
    if (a == 0 || b == 0)
    {
        return 0;
    }
    long lcm = std::abs(Arithmetic_Utils.mul_and_check(a / gcd(a, b), b));
    if (lcm == std::numeric_limits<long>::min())
    {
        throw Math_Runtime_Exception(Localized_Core_Formats.LCM_OVERFLOW_64_BITS, a, b);
    }
    return lcm;
}

/**
    * Multiply two integers, checking for overflow.
    *
    * @param x Factor.
    * @param y Factor.
    * @return the product {@code x * y}.
    * @Math_Runtime_Exception if the result can not be
    * represented as an {@code int}.
    */
int mul_and_check(const int& x, int y) Math_Runtime_Exception 
{
    long m = (static_cast<long>(x) * (static_cast<long>(y);
    if (m < std::numeric_limits<int>::min() || m > std::numeric_limits<int>::max()) 
    {
        throw Math_Runtime_Exception(Localized_Core_Formats.ARITHMETIC_EXCEPTION);
    }
    return static_cast<int>(m;
}

    /**
     * Multiply two long integers, checking for overflow.
     *
     * @param a Factor.
     * @param b Factor.
     * @return the product {@code a * b}.
     * @Math_Runtime_Exception if the result can not be represented
     * as a {@code long}.
     */
long Arithmetic_Utils::mul_and_check(long a, long b)
{
    long ret;
    if (a > b) 
    {
        // use symmetry to reduce boundary cases
        ret = mul_and_check(b, a);
    }
    else 
    {
        if (a < 0) 
        {
            if (b < 0) 
            {
                // check for positive overflow with negative a, negative b
                if (a >= std::numeric_limits<long>::max() / b) 
                {
                    ret = a * b;
                }
                else 
                {
                    throw Math_Runtime_Exception(Localized_Core_Formats.ARITHMETIC_EXCEPTION);
                }
            }
            else if (b > 0) 
            {
                // check for negative overflow with negative a, positive b
                if (std::numeric_limits<long>::min() / b <= a) 
                {
                    ret = a * b;
                }
                else 
                {
                    throw Math_Runtime_Exception(Localized_Core_Formats.ARITHMETIC_EXCEPTION);

                }
            }
            else 
            {
                // assert b == 0
                ret = 0;
            }
        }
        else if (a > 0) 
        {
            // assert a > 0
            // assert b > 0

            // check for positive overflow with positive a, positive b
            if (a <= std::numeric_limits<long>::max() / b) 
            {
                ret = a * b;
            }
            else 
            {
                throw Math_Runtime_Exception(Localized_Core_Formats.ARITHMETIC_EXCEPTION);
            }
        }
        else 
        {
            // assert a == 0
            ret = 0;
        }
    }
    return ret;
}

/**
    * Subtract two integers, checking for overflow.
    *
    * @param x Minuend.
    * @param y Subtrahend.
    * @return the difference {@code x - y}.
    * @Math_Runtime_Exception if the result can not be represented
    * as an {@code int}.
    */
int Arithmetic_Utils::sub_and_check(const int& x, int y)
{
    long s = static_cast<long>(x - static_cast<long>(y;
    if (s < std::numeric_limits<int>::min() || s > std::numeric_limits<int>::max()) 
    {
        throw Math_Runtime_Exception(Localized_Core_Formats.OVERFLOW_IN_SUBTRACTION, x, y);
    }
    return static_cast<int>(s;
}

/**
    * Subtract two long integers, checking for overflow.
    *
    * @param a Value.
    * @param b Value.
    * @return the difference {@code a - b}.
    * @Math_Runtime_Exception if the result can not be represented as a
    * {@code long}.
    */
long Arithmetic_Utils::sub_and_check(const long& a, const long& b)
{
    long ret;
    if (b == std::numeric_limits<long>::min()) 
    {
        if (a < 0) 
        {
            ret = a - b;
        }
        else 
        {
            throw Math_Runtime_Exception(Localized_Core_Formats.OVERFLOW_IN_ADDITION, a, -b);
        }
    }
    else 
    {
        // use additive inverse
        ret = add_and_check(a, -b, Localized_Core_Formats.OVERFLOW_IN_ADDITION);
    }
    return ret;
}

/**
    * Raise an int to an int power.
    *
    * @param k Number to raise.
    * @param e Exponent (must be positive or zero).
    * @return \( k^e \)
    * @ if {@code e < 0}.
    * @Math_Runtime_Exception if the result would overflow.
    */
int Arithmetic_Utils::pow(const int& k, const int& e) 
    {
    if (e < 0) 
    {
        throw (Localized_Core_Formats.EXPONENT, e);
    }

    int exp = e;
    int result = 1;
    int k2p    = k;
    while (true) 
    {
        if ((exp & 0x1) != 0) 
        {
            result = mul_and_check(result, k2p);
        }

        exp >>= 1;
    if (exp == 0) 
    {
        break;
    }

    k2p = mul_and_check(k2p, k2p);
    }

    return result;
}

/**
    * Raise a long to an int power.
    *
    * @param k Number to raise.
    * @param e Exponent (must be positive or zero).
    * @return \( k^e \)
    * @ if {@code e < 0}.
    * @Math_Runtime_Exception if the result would overflow.
    */
long Arithmetic_Utils::pow(const long k, const int& e) 
    {
    if (e < 0) 
    {
        throw (Localized_Core_Formats.EXPONENT, e);
    }

    int exp = e;
    long result = 1;
    long k2p    = k;
    while (true) 
    {
        if ((exp & 0x1) != 0) 
        {
            result = mul_and_check(result, k2p);
        }

        exp >>= 1;
    if (exp == 0) 
    {
        break;
    }

    k2p = mul_and_check(k2p, k2p);
    }

    return result;
}

/**
    * Raise a BigInteger to an int power.
    *
    * @param k Number to raise.
    * @param e Exponent (must be positive or zero).
    * @return k<sup>e</sup>
    * @ if {@code e < 0}.
    */
BigInteger Arithmetic_Utils::pow(const BigInteger& k, const int& e)
{
    if (e < 0) 
    {
        throw (Localized_Core_Formats.EXPONENT, e);
    }

    return k.pow(e);
}

/**
    * Raise a BigInteger to a long power.
    *
    * @param k Number to raise.
    * @param e Exponent (must be positive or zero).
    * @return k<sup>e</sup>
    * @ if {@code e < 0}.
    */
BigInteger Arithmetic_Utils::pow(const BigInteger& k, const long& e)
{
    if (e < 0) 
    {
        throw (Localized_Core_Formats.EXPONENT, e);
    }

    BigInteger result = BigInteger.ONE;
    BigInteger k2p    = k;
    while (e != 0) 
    {
        if ((e & 0x1) != 0) 
        {
            result = result.multiply(k2p);
        }
        k2p = k2p.multiply(k2p);
        e >>= 1;
    }

    return result;

}

/**
    * Raise a BigInteger to a BigInteger power.
    *
    * @param k Number to raise.
    * @param e Exponent (must be positive or zero).
    * @return k<sup>e</sup>
    * @ if {@code e < 0}.
    */
BigInteger Arithmetic_Utils::pow(const BigInteger& k, const BigInteger& e)
{
    if (e.compare_to(BigInteger.ZERO) < 0) 
    {
        throw (Localized_Core_Formats.EXPONENT, e);
    }

    BigInteger result = BigInteger.ONE;
    BigInteger k2p    = k;
    while (!BigInteger.ZERO.equals(e)) 
    {
        if (e.test_bit(0)) 
        {
            result = result.multiply(k2p);
        }
        k2p = k2p.multiply(k2p);
        e = e.shift_right(1);
    }

    return result;
}

/**
    * Add two long integers, checking for overflow.
    *
    * @param a Addend.
    * @param b Addend.
    * @param pattern Pattern to use for any thrown exception.
    * @return the sum {@code a + b}.
    * @Math_Runtime_Exception if the result cannot be represented
    * as a {@code long}.
    */
    long Arithmetic_Utils::add_and_check(const long& a, const long& b, const Localizable& pattern)
    {
        const long result = a + b;
        if (!((a ^ b) < 0 || (a ^ result) >= 0)) 
        {
            throw Math_Runtime_Exception(pattern, a, b);
        }
        return result;
}

/**
    * Returns true if the argument is a power of two.
    *
    * @param n the number to test
    * @return true if the argument is a power of two
    */
bool Arithmetic_Utils::is_power_of_two(const long& n)
{
    return (n > 0) && ((n & (n - 1)) == 0);
}

/**
    * Returns the unsigned remainder from dividing the first argument
    * by the second where each argument and the result is interpreted
    * as an unsigned value.
    * <p>
    * This method does not use the {@code long} datatype.
    *
    * @param dividend the value to be divided
    * @param divisor the value doing the dividing
    * @return the unsigned remainder of the first argument divided by
    * the second argument.
    */
int Arithmetic_Utils::remainder_unsigned(const int& dividend, const int& divisor)
{
    if (divisor >= 0) 
    {
        if (dividend >= 0) 
        {
            return dividend % divisor;
        }
        // The implementation is a Java port of algorithm described in the book
        // "Hacker's Delight" (section "Unsigned short division from signed division").
        int q = ((dividend >>> 1) / divisor) << 1;
        dividend -= q * divisor;
        if (dividend < 0 || dividend >= divisor) 
        {
            dividend -= divisor;
        }
        return dividend;
    }
    return dividend >= 0 || dividend < divisor ? dividend : dividend - divisor;
}

/**
    * Returns the unsigned remainder from dividing the first argument
    * by the second where each argument and the result is interpreted
    * as an unsigned value.
    * <p>
    * This method does not use the {@code BigInteger} datatype.
    *
    * @param dividend the value to be divided
    * @param divisor the value doing the dividing
    * @return the unsigned remainder of the first argument divided by
    * the second argument.
    */
long Arithmetic_Utils::remainder_unsigned(const long& dividend, const long& divisor)
{
    if (divisor >= 0L) 
    {
        if (dividend >= 0L) 
        {
            return dividend % divisor;
        }
        // The implementation is a Java port of algorithm described in the book
        // "Hacker's Delight" (section "Unsigned short division from signed division").
        long q = ((dividend >>> 1) / divisor) << 1;
        dividend -= q * divisor;
        if (dividend < 0L || dividend >= divisor) 
        {
            dividend -= divisor;
        }
        return dividend;
    }
    return dividend >= 0L || dividend < divisor ? dividend : dividend - divisor;
}

/**
    * Returns the unsigned quotient of dividing the first argument by
    * the second where each argument and the result is interpreted as
    * an unsigned value.
    * <p>
    * Note that in two's complement arithmetic, the three other
    * basic arithmetic operations of add, subtract, and multiply are
    * bit-wise identical if the two operands are regarded as both
    * being signed or both being unsigned. Therefore separate {@code
    * add_unsigned}, etc. methods are not provided.
    * <p>
    * This method does not use the {@code long} datatype.
    *
    * @param dividend the value to be divided
    * @param divisor the value doing the dividing
    * @return the unsigned quotient of the first argument divided by
    * the second argument
    */
int Arithmetic_Utils::divide_unsigned(const int& dividend, const int& divisor)
{
    if (divisor >= 0) 
    {
        if (dividend >= 0) 
        {
            return dividend / divisor;
        }
        // The implementation is a Java port of algorithm described in the book
        // "Hacker's Delight" (section "Unsigned short division from signed division").
        int q = ((dividend >>> 1) / divisor) << 1;
        dividend -= q * divisor;
        if (dividend < 0L || dividend >= divisor) 
        {
            q++;
        }
        return q;
    }
    return dividend >= 0 || dividend < divisor ? 0 : 1;
}

/**
    * Returns the unsigned quotient of dividing the first argument by
    * the second where each argument and the result is interpreted as
    * an unsigned value.
    * <p>
    * Note that in two's complement arithmetic, the three other
    * basic arithmetic operations of add, subtract, and multiply are
    * bit-wise identical if the two operands are regarded as both
    * being signed or both being unsigned. Therefore separate {@code
    * add_unsigned}, etc. methods are not provided.
    * <p>
    * This method does not use the {@code BigInteger} datatype.
    *
    * @param dividend the value to be divided
    * @param divisor the value doing the dividing
    * @return the unsigned quotient of the first argument divided by
    * the second argument.
    */
long Arithmetic_Utils::divide_unsigned(const long& dividend, const long& divisor)
{
    if (divisor >= 0) 
    {
        if (dividend >= 0) 
        {
            return dividend / divisor;
        }
        // The implementation is a Java port of algorithm described in the book
        // "Hacker's Delight" (section "Unsigned short division from signed division").
        long q = ((dividend >>> 1) / divisor) << 1;
        dividend -= q * divisor;
        if (dividend < 0L || dividend >= divisor) 
        {
            q++;
        }
        return q;
    }
    return dividend >= 0L || dividend < divisor ? 0L : 1L;
}



