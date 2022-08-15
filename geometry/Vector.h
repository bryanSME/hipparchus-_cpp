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

#include <type_traits>
#include "Point.hpp"
#include "Space.h"

  /** This interface represents a generic vector in a vectorial space or a point in an affine space.
   * @param <S> Type of the space.
   * @see Space
   * @see Point
   */
template<
	typename S,
	typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Vector //<S extends Space> extends Point<S>
{
	/** Get the NULL vector of the vectorial space or origin point of the affine space.
	 * @return NULL vector of the vectorial space or origin point of the affine space
	 */
	Vector<S> get_zero() = 0;

	/** Get the L<sub>1</sub> norm for the vector.
	 * @return L<sub>1</sub> norm for the vector
	 */
	virtual double get_norm1() = 0;

	/** Get the L<sub>2</sub> norm for the vector.
	 * @return Euclidean norm for the vector
	 */
	virtual double get_norm() = 0;

	/** Get the square of the norm for the vector.
	 * @return square of the Euclidean norm for the vector
	 */
	virtual double get_norm_sq() = 0;

	/** Get the L<sub>&infin;</sub> norm for the vector.
	 * @return L<sub>&infin;</sub> norm for the vector
	 */
	virtual double get_norm_inf() = 0;

	/** Add a vector to the instance.
	 * @param v vector to add
	 * @return a vector
	 */
	virtual Vector<S> add(const Vector<S>& v) = 0;

	/** Add a scaled vector to the instance.
	 * @param factor scale factor to apply to v before adding it
	 * @param v vector to add
	 * @return a vector
	 */
	virtual Vector<S> add(const double& factor, const Vector<S>& v) = 0;

	/** Subtract a vector from the instance.
	 * @param v vector to subtract
	 * @return a vector
	 */
	virtual Vector<S> subtract(const Vector<S>& v) = 0;

	/** Subtract a scaled vector from the instance.
	 * @param factor scale factor to apply to v before subtracting it
	 * @param v vector to subtract
	 * @return a vector
	 */
	virtual Vector<S> subtract(const double& factor, const Vector<S>& v) = 0;

	/** Get the opposite of the instance.
	 * @return a vector which is opposite to the instance
	 */
	virtual Vector<S> negate() = 0;

	/** Get a normalized vector aligned with the instance.
	 * @return a normalized vector
	 * @exception Math_Runtime_Exception if the norm is zero
	 */
	virtual Vector<S> normalize() = 0;

	/** Multiply the instance by a scalar.
	 * @param a scalar
	 * @return a vector
	 */
	virtual Vector<S> scalar_multiply(const double& a) = 0;

	/**
	 * Returns true if any coordinate of this vector is infinite and none are NaN;
	 * false otherwise
	 * @return  true if any coordinate of this vector is infinite and none are NaN;
	 * false otherwise
	 */
	virtual bool std::isinfinite() = 0;

	/** Compute the distance between the instance and another vector according to the L<sub>1</sub> norm.
	 * <p>Calling this method is equivalent to calling:
	 * <code>q.subtract(p).get_norm1()</code> except that no intermediate
	 * vector is built</p>
	 * @param v second vector
	 * @return the distance between the instance and p according to the L<sub>1</sub> norm
	 */
	virtual double distance1(const Vector<S>& v) = 0;

	/** Compute the distance between the instance and another vector according to the L<sub>&infin;</sub> norm.
	 * <p>Calling this method is equivalent to calling:
	 * <code>q.subtract(p).get_norm_inf()</code> except that no intermediate
	 * vector is built</p>
	 * @param v second vector
	 * @return the distance between the instance and p according to the L<sub>&infin;</sub> norm
	 */
	virtual double distance_inf(const Vector<S>& v) = 0;

	/** Compute the square of the distance between the instance and another vector.
	 * <p>Calling this method is equivalent to calling:
	 * <code>q.subtract(p).get_norm_sq()</code> except that no intermediate
	 * vector is built</p>
	 * @param v second vector
	 * @return the square of the distance between the instance and p
	 */
	virtual double distance_sq(const Vector<S>& v) = 0;

	/** Compute the dot-product of the instance and another vector.
	 * @param v second vector
	 * @return the dot product this.v
	 */
	virtual double dot_product(const Vector<S>& v) = 0;

	/** Get a string representation of this vector.
	 * @param format the custom format for components
	 * @return a string representation of this vector
	 */
	virtual std::string to_string(const Number_Format& format) = 0;
};