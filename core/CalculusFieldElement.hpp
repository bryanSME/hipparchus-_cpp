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

#include "FieldElement.h"
#include <vector>

  /**
   * Interface representing a <a href="http://mathworld.wolfram.com/Field.html">field</a>
   * with calculus capabilities (sin, cos, ...).
   * @param <T> the type of the field elements
   * @see Field_Element
   * @since 1.7
   */
template<typename T>
class Calculus_Field_Element : public Field_Element<T>
{
	/** Get the Archimedes constant π.
	 * <p>
	 * Archimedes constant is the ratio of a circle's circumference to its diameter.
	 * </p>
	 * @return Archimedes constant π
	 * @since 2.0
	 */
	virtual T get_pi() = 0;

	/** Create an instance corresponding to a constant real value.
	 * @param value constant real value
	 * @return instance corresponding to a constant real value
	 */
	virtual T new_instance(double value) = 0;

	/** '+' operator.
	 * @param a right hand side parameter of the operator
	 * @return this+a
	 */
	virtual T add(double a) = 0;

	/** '-' operator.
	 * @param a right hand side parameter of the operator
	 * @return this-a
	 */
	virtual T subtract(double a) = 0;

	/** '&times;' operator.
	 * @param a right hand side parameter of the operator
	 * @return this&times;a
	 */
	virtual T multiply(double a) = 0;

	/** '&divide;' operator.
	 * @param a right hand side parameter of the operator
	 * @return this&divide;a
	 */
	virtual T divide(double a) = 0;

	/**
	 * Return the exponent of the instance, removing the bias.
	 * <p>
	 * For double numbers of the form 2<sup>x</sup>, the unbiased
	 * exponent is exactly x.
	 * </p>
	 * @return exponent for the instance, without bias
	 */
	default int get_exponent()
	{
		return FastMath.get_exponent(get_real());
	}

	/**
	 * Multiply the instance by a power of 2.
	 * @param n power of 2
	 * @return this &times; 2<sup>n</sup>
	 */
	virtual T scalb(const int& n) = 0;

	/**
	 * Compute least significant bit (Unit in Last Position) for a number.
	 * @return ulp(this)
	 * @since 2.0
	 */
	virtual T ulp() = 0;

	/**
	 * Returns the hypotenuse of a triangle with sides {@code this} and {@code y}
	 * - sqrt(<i>this</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
	 * avoiding intermediate overflow or underflow.
	 *
	 * <ul>
	 * <li> If either argument is infinite, then the result is positive infinity.</li>
	 * <li> else, if either argument is NaN then the result is NaN.</li>
	 * </ul>
	 *
	 * @param y a value
	 * @return sqrt(<i>this</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
	 * @exception  if number of free parameters or orders are inconsistent
	 */
	virtual T hypot(T y) = 0;

	/** {@inherit_doc} */
	//override
	virtual T reciprocal() = 0;

	/** Square root.
	 * @return square root of the instance
	 */
	virtual T sqrt() = 0;

	/** Cubic root.
	 * @return cubic root of the instance
	 */
	virtual T cbrt() = 0;

	/** N<sup>th</sup> root.
	 * @param n order of the root
	 * @return n<sup>th</sup> root of the instance
	 */
	virtual T root_n(const int& n) = 0;

	/** Power operation.
	 * @param p power to apply
	 * @return this<sup>p</sup>
	 */
	virtual T pow(const double& p) = 0;

	/** Integer power operation.
	 * @param n power to apply
	 * @return this<sup>n</sup>
	 */
	virtual T pow(const int& n) = 0;

	/** Power operation.
	 * @param e exponent
	 * @return this<sup>e</sup>
	 * @exception  if number of free parameters or orders are inconsistent
	 */
	virtual T pow(T e) = 0;

	/** Exponential.
	 * @return exponential of the instance
	 */
	virtual T exp() = 0;

	/** Exponential minus 1.
	 * @return exponential minus one of the instance
	 */
	virtual T expm1() = 0;

	/** Natural logarithm.
	 * @return logarithm of the instance
	 */
	virtual T log() = 0;

	/** Shifted natural logarithm.
	 * @return logarithm of one plus the instance
	 */
	virtual T log1p() = 0;

	/** Base 10 logarithm.
	 * @return base 10 logarithm of the instance
	 */
	virtual T log10() = 0;

	/** Cosine operation.
	 * @return cos(this)
	 */
	virtual T cos() = 0;

	/** Sine operation.
	 * @return sin(this)
	 */
	virtual T sin() = 0;

	/** Combined Sine and Cosine operation.
	 * @return [sin(this), cos(this)]
	 * @since 1.4
	 */
	virtual Field_Sin_Cos<T> sin_cos() = 0;

	/** Tangent operation.
	 * @return tan(this)
	 */
	virtual T tan() = 0;

	/** Arc cosine operation.
	 * @return acos(this)
	 */
	virtual T acos() = 0;

	/** Arc sine operation.
	 * @return asin(this)
	 */
	virtual T asin() = 0;

	/** Arc tangent operation.
	 * @return atan(this)
	 */
	virtual T atan() = 0;

	/** Two arguments arc tangent operation.
	 * <p>
	 * Beware of the order or arguments! As this is based on a
	 * two-arguments functions, in order to be consistent with
	 * arguments order, the instance is the <em>first</em> argument
	 * and the single provided argument is the <em>second</em> argument.
	 * In order to be consistent with programming languages {@code atan2}, * this method computes {@code atan2(this, x)}, i.e. the instance
	 * represents the {@code y} argument and the {@code x} argument is
	 * the one passed as a single argument. This may seem confusing especially
	 * for users of Wolfram alpha, as this site is <em>not</em> consistent
	 * with programming languages {@code atan2} two-arguments arc tangent
	 * and puts {@code x} as its first argument.
	 * </p>
	 * @param x second argument of the arc tangent
	 * @return atan2(this, x)
	 * @exception  if number of free parameters or orders are inconsistent
	 */
	virtual T atan2(T x) = 0;

	/** Hyperbolic cosine operation.
	 * @return cosh(this)
	 */
	virtual T cosh() = 0;

	/** Hyperbolic sine operation.
	 * @return sinh(this)
	 */
	virtual T sinh() = 0;

	/** Combined hyperbolic sine and sosine operation.
	 * @return [sinh(this), cosh(this)]
	 * @since 2.0
	 */
	virtual Field_Sinh_Cosh<T> sinh_cosh() = 0;

	/** Hyperbolic tangent operation.
	 * @return tanh(this)
	 */
	virtual T tanh() = 0;

	/** Inverse hyperbolic cosine operation.
	 * @return acosh(this)
	 */
	virtual T acosh() = 0;

	/** Inverse hyperbolic sine operation.
	 * @return asin(this)
	 */
	virtual T asinh() = 0;

	/** Inverse hyperbolic  tangent operation.
	 * @return atanh(this)
	 */
	virtual T atanh() = 0;

	/** Convert radians to degrees, with error of less than 0.5 ULP
	 *  @return instance converted into degrees
	 */
	virtual T to_degrees() = 0;

	/** Convert degrees to radians, with error of less than 0.5 ULP
	 *  @return instance converted into radians
	 */
	virtual T to_radians() = 0;

	/**
	 * Compute a linear combination.
	 * @param a Factors.
	 * @param b Factors.
	 * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
	 * @ if arrays dimensions don't match
	 */
	virtual T linear_combination(std::vector<T> a, std::vector<T> b) = 0;

	/**
	 * Compute a linear combination.
	 * @param a Factors.
	 * @param b Factors.
	 * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
	 * @ if arrays dimensions don't match
	 */
	virtual T linear_combination(std::vector<double> a, std::vector<T> b) = 0;

	/**
	 * Compute a linear combination.
	 * @param a1 first factor of the first term
	 * @param b1 second factor of the first term
	 * @param a2 first factor of the second term
	 * @param b2 second factor of the second term
	 * @return a<sub>1</sub>&times;b<sub>1</sub> +
	 * a<sub>2</sub>&times;b<sub>2</sub>
	 * @see #linear_combination(Field_Element, Field_Element, Field_Element, Field_Element, Field_Element, Field_Element)
	 * @see #linear_combination(Field_Element, Field_Element, Field_Element, Field_Element, Field_Element, Field_Element, Field_Element, Field_Element)
	 */
	virtual T linear_combination(T a1, T b1, T a2, T b2) = 0;

	/**
	 * Compute a linear combination.
	 * @param a1 first factor of the first term
	 * @param b1 second factor of the first term
	 * @param a2 first factor of the second term
	 * @param b2 second factor of the second term
	 * @return a<sub>1</sub>&times;b<sub>1</sub> +
	 * a<sub>2</sub>&times;b<sub>2</sub>
	 * @see #linear_combination(double, Field_Element, double, Field_Element, double, Field_Element)
	 * @see #linear_combination(double, Field_Element, double, Field_Element, double, Field_Element, double, Field_Element)
	 */
	virtual T linear_combination(double a1, T b1, double a2, T b2) = 0;

	/**
	 * Compute a linear combination.
	 * @param a1 first factor of the first term
	 * @param b1 second factor of the first term
	 * @param a2 first factor of the second term
	 * @param b2 second factor of the second term
	 * @param a3 first factor of the third term
	 * @param b3 second factor of the third term
	 * @return a<sub>1</sub>&times;b<sub>1</sub> +
	 * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub>
	 * @see #linear_combination(Field_Element, Field_Element, Field_Element, Field_Element)
	 * @see #linear_combination(Field_Element, Field_Element, Field_Element, Field_Element, Field_Element, Field_Element, Field_Element, Field_Element)
	 */
	virtual T linear_combination(T a1, T b1, T a2, T b2, T a3, T b3) = 0;

	/**
	 * Compute a linear combination.
	 * @param a1 first factor of the first term
	 * @param b1 second factor of the first term
	 * @param a2 first factor of the second term
	 * @param b2 second factor of the second term
	 * @param a3 first factor of the third term
	 * @param b3 second factor of the third term
	 * @return a<sub>1</sub>&times;b<sub>1</sub> +
	 * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub>
	 * @see #linear_combination(double, Field_Element, double, Field_Element)
	 * @see #linear_combination(double, Field_Element, double, Field_Element, double, Field_Element, double, Field_Element)
	 */
	virtual T linear_combination(double a1, T b1, double a2, T b2, double a3, T b3) = 0;

	/**
	 * Compute a linear combination.
	 * @param a1 first factor of the first term
	 * @param b1 second factor of the first term
	 * @param a2 first factor of the second term
	 * @param b2 second factor of the second term
	 * @param a3 first factor of the third term
	 * @param b3 second factor of the third term
	 * @param a4 first factor of the fourth term
	 * @param b4 second factor of the fourth term
	 * @return a<sub>1</sub>&times;b<sub>1</sub> +
	 * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub> +
	 * a<sub>4</sub>&times;b<sub>4</sub>
	 * @see #linear_combination(Field_Element, Field_Element, Field_Element, Field_Element)
	 * @see #linear_combination(Field_Element, Field_Element, Field_Element, Field_Element, Field_Element, Field_Element)
	 */
	virtual T linear_combination(T a1, T b1, T a2, T b2, T a3, T b3, T a4, T b4) = 0;

	/**
	 * Compute a linear combination.
	 * @param a1 first factor of the first term
	 * @param b1 second factor of the first term
	 * @param a2 first factor of the second term
	 * @param b2 second factor of the second term
	 * @param a3 first factor of the third term
	 * @param b3 second factor of the third term
	 * @param a4 first factor of the fourth term
	 * @param b4 second factor of the fourth term
	 * @return a<sub>1</sub>&times;b<sub>1</sub> +
	 * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub> +
	 * a<sub>4</sub>&times;b<sub>4</sub>
	 * @see #linear_combination(double, Field_Element, double, Field_Element)
	 * @see #linear_combination(double, Field_Element, double, Field_Element, double, Field_Element)
	 */
	virtual T linear_combination(double a1, T b1, double a2, T b2, double a3, T b3, double a4, T b4) = 0;

	/** Get the smallest whole number larger than instance.
	 * @return ceil(this)
	 */
	virtual T ceil() = 0;

	/** Get the largest whole number smaller than instance.
	 * @return floor(this)
	 */
	virtual T floor() = 0;

	/** Get the whole number that is the nearest to the instance, or the even one if x is exactly half way between two integers.
	 * @return a double number r such that r is an integer r - 0.5 &le; this &le; r + 0.5
	 */
	virtual T rint() = 0;

	/** IEEE remainder operator.
	 * @param a right hand side parameter of the operator
	 * @return this - n &times; a where n is the closest integer to this/a
	 */
	virtual T remainder(double a) = 0;

	/** IEEE remainder operator.
	 * @param a right hand side parameter of the operator
	 * @return this - n &times; a where n is the closest integer to this/a
	 */
	virtual T remainder(T a) = 0;

	/** Compute the sign of the instance.
	 * The sign is -1 for negative numbers, +1 for positive numbers and 0 otherwise, * for std::complex<double> number, it is extended on the unit circle (equivalent to z/|z|, * with special handling for 0 and NaN)
	 * @return -1.0, -0.0, +0.0, +1.0 or NaN depending on sign of a
	 */
	T sign() = 0;

	/**
	 * Returns the instance with the sign of the argument.
	 * A NaN {@code sign} argument is treated as positive.
	 *
	 * @param sign the sign for the returned value
	 * @return the instance with the same sign as the {@code sign} argument
	 */
	virtual T copy_sign(T sign) = 0;

	/**
	 * Returns the instance with the sign of the argument.
	 * A NaN {@code sign} argument is treated as positive.
	 *
	 * @param sign the sign for the returned value
	 * @return the instance with the same sign as the {@code sign} argument
	 */
	virtual T copy_sign(double sign) = 0;

	/**
	 * Check if the instance is infinite.
	 * @return true if the instance is infinite
	 */
	bool std::isinfinite()
	{
		return std::isinf(get_real());
	}

	/**
	 * Check if the instance is finite (neither infinite nor NaN).
	 * @return true if the instance is finite (neither infinite nor NaN)
	 * @since 2.0
	 */
	bool is_finite()
	{
		return Double.is_finite(get_real());
	}

	/**
	 * Check if the instance is Not a Number.
	 * @return true if the instance is Not a Number
	 */
	bool is_nan()
	{
		return std::isnan(get_real());
	}

	/** norm.
	 * @return norm(this)
	 * @since 2.0
	 */
	double norm()
	{
		return abs().get_real();
	}

	/** absolute value.
	 * <p>
	 * Just another name for {@link #norm()}
	 * </p>
	 * @return abs(this)
	 */
	virtual T abs() = 0;

	/** Get the closest long to instance real value.
	 * @return closest long to {@link #get_real()}
	 */
	long round()
	{
		return std::round(get_real());
	}
};