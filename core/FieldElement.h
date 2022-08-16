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
  //package org.hipparchus;

  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.exception.;

  /**
   * Interface representing <a href="http://mathworld.wolfram.com/Field.html">field</a> elements.
   * @param <T> the type of the field elements
   * @see Field
   */
template<typename T>
class Field_Element
	//class Field_Element<T extends Field_Element<T>>
{
	/** Get the real value of the number.
	 * @return real value
	 */
	virtual double get_real() = 0;

	/** Compute this + a.
	 * @param a element to add
	 * @return a element representing this + a
	 * @ if {@code a} is {@code NULL}.
	 */
	virtual T add(T a) = 0;

	/** Compute this - a.
	 * @param a element to subtract
	 * @return a element representing this - a
	 * @ if {@code a} is {@code NULL}.
	 */
	virtual T subtract(T a) = 0;

	/**
	 * Returns the additive inverse of {@code this} element.
	 * @return the opposite of {@code this}.
	 */
	virtual T negate() = 0;

	/** Compute n &times; this. Multiplication by an integer number is defined
	 * as the following sum
	 * <center>
	 * n &times; this = &sum;<sub>i=1</sub><sup>n</sup> this.
	 * </center>
	 * @param n Number of times {@code this} must be added to itself.
	 * @return A element representing n &times; this.
	 */
	virtual T multiply(const int& n) = 0;

	/** Compute this &times; a.
	 * @param a element to multiply
	 * @return a element representing this &times; a
	 * @ if {@code a} is {@code NULL}.
	 */
	virtual T multiply(T a) = 0;

	/** Compute this &divide; a.
	 * @param a element to divide by
	 * @return a element representing this &divide; a
	 * @ if {@code a} is {@code NULL}.
	 * @Math_Runtime_Exception if {@code a} is zero
	 */
	virtual T divide(T a) = 0;

	/**
	 * Returns the multiplicative inverse of {@code this} element.
	 * @return the inverse of {@code this}.
	 * @Math_Runtime_Exception if {@code this} is zero
	 */
	virtual T reciprocal() = 0;

	/** Get the {@link Field} to which the instance belongs.
	 * @return {@link Field} to which the instance belongs
	 */
	virtual Field<T> get_field() = 0;

	/** Check if an element is semantically equal to zero.
	 * <p>
	 * The default implementation simply calls {@code equals(get_field().get_zero())}.
	 * However, this may need to be overridden in some cases as due to
	 * compatibility with {@code hash_code()} some classes implements
	 * {@code equals(Object)} in such a way that -0.0 and +0.0 are different, * which may be a problem. It prevents for example identifying a diagonal
	 * element is zero and should be avoided when doing partial pivoting in
	 * LU decomposition.
	 * </p>
	 * @return true if the element is semantically equal to zero
	 * @since 1.8
	 */
	bool is_zero()
	{
		return equals(get_field().get_zero());
	}
};