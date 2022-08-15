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

#include <limits>
#include "BigInteger.h"

  /**
   * Some useful, arithmetics related, additions to the built-in functions in
   * {@link Math}.
   */
class Arithmetic_Utils
{
public:
	/**
	 * Add two integers, checking for overflow.
	 *
	 * @param x an addend
	 * @param y an addend
	 * @return the sum {@code x+y}
	 * @Math_Runtime_Exception if the result can not be represented
	 * as an {@code int}.
	 */
	static int add_and_check(const int& x, const int& y);

	/**
	 * Add two long integers, checking for overflow.
	 *
	 * @param a an addend
	 * @param b an addend
	 * @return the sum {@code a+b}
	 * @Math_Runtime_Exception if the result can not be represented as an long
	 */
	static long add_and_check(const long& a, const long& b);

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
	static int gcd(const int& p, const int& q);

	/**
	 * Returns true if the argument is a power of two.
	 *
	 * @param n the number to test
	 * @return true if the argument is a power of two
	 */
	static bool is_power_of_two(const long& n);

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
	static int remainder_unsigned(const int& dividend, const int& divisor);

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
	static long remainder_unsigned(long dividend, long divisor);

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
	static int divide_unsigned(const int& dividend, int divisor);

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
	static long divide_unsigned(long dividend, long divisor);

	/**
	 * Gets the greatest common divisor of the absolute value of two numbers, * using the "binary gcd" method which avoids division and modulo
	 * operations. See Knuth 4.5.2 algorithm B. This algorithm is due to Josef
	 * Stein (1961).
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>The invocations
	 * {@code gcd(long.MIN_VALUE, long.MIN_VALUE)}, * {@code gcd(long.MIN_VALUE, 0L)} and
	 * {@code gcd(0L, long.MIN_VALUE)} throw an
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
	static long gcd(const long p, const long q);

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
	static int lcm(const int& a, int b);

	/**
	 * Returns the least common multiple of the absolute value of two numbers, * using the formula {@code lcm(a,b) = (a / gcd(a,b)) * b}.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>The invocations {@code lcm(long.MIN_VALUE, n)} and
	 * {@code lcm(n, long.MIN_VALUE)}, where {@code abs(n)} is a
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
	static long lcm(long a, long b);

	/**
	 * Multiply two integers, checking for overflow.
	 *
	 * @param x Factor.
	 * @param y Factor.
	 * @return the product {@code x * y}.
	 * @Math_Runtime_Exception if the result can not be
	 * represented as an {@code int}.
	 */
	static int mul_and_check(const int& x, int y);

	/**
	 * Multiply two long integers, checking for overflow.
	 *
	 * @param a Factor.
	 * @param b Factor.
	 * @return the product {@code a * b}.
	 * @Math_Runtime_Exception if the result can not be represented
	 * as a {@code long}.
	 */
	static long mul_and_check(long a, long b);

	/**
	 * Subtract two integers, checking for overflow.
	 *
	 * @param x Minuend.
	 * @param y Subtrahend.
	 * @return the difference {@code x - y}.
	 * @Math_Runtime_Exception if the result can not be represented
	 * as an {@code int}.
	 */
	static int sub_and_check(const int& x, int y);

	/**
	 * Subtract two long integers, checking for overflow.
	 *
	 * @param a Value.
	 * @param b Value.
	 * @return the difference {@code a - b}.
	 * @Math_Runtime_Exception if the result can not be represented as a
	 * {@code long}.
	 */
	static long sub_and_check(long a, long b);

	/**
	 * Raise an int to an int power.
	 *
	 * @param k Number to raise.
	 * @param e Exponent (must be positive or zero).
	 * @return \( k^e \)
	 * @ if {@code e < 0}.
	 * @Math_Runtime_Exception if the result would overflow.
	 */
	static int pow(const int& k, const int e);

	/**
	 * Raise a BigInteger to a long power.
	 *
	 * @param k Number to raise.
	 * @param e Exponent (must be positive or zero).
	 * @return k<sup>e</sup>
	 * @ if {@code e < 0}.
	 */
	static BigInteger pow(const BigInteger& k, const long& e);

	/**
	 * Raise a BigInteger to a BigInteger power.
	 *
	 * @param k Number to raise.
	 * @param e Exponent (must be positive or zero).
	 * @return k<sup>e</sup>
	 * @ if {@code e < 0}.
	 */
	static BigInteger pow(const BigInteger& k, const BigInteger& e);

private:

	/** Private constructor. */
	Arithmetic_Utils()
	{
		//super();
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
	static int gcd_positive(const int& a, const int& b);

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
	//static long add_and_check(const long& a, const long& b, const Localizable& pattern);
};