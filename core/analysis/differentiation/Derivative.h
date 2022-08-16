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
 //package org.hipparchus.analysis.differentiation;

 //import org.hipparchus.Calculus_Field_Element;
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/** Interface representing both the value and the differentials of a function.
 * @param <T> the type of the field elements
 * @since 1.7
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Derivative : public Calculus_Field_Element<T>
{
	/** Get the number of free parameters.
	 * @return number of free parameters
	 */
	virtual int get_free_parameters() = 0;

	/** Get the derivation order.
	 * @return derivation order
	 */
	virtual int get_order() = 0;

	/** Get the value part of the function.
	 * @return value part of the value of the function
	 */
	virtual double get_value() = 0;

	/** Get a partial derivative.
	 * @param orders derivation orders with respect to each variable (if all orders are 0, * the value is returned)
	 * @return partial derivative
	 * @see #get_value()
	 * @exception  if the numbers of variables does not
	 * match the instance
	 * @exception  if sum of derivation orders is larger
	 * than the instance limits
	 */
	virtual double get_partial_derivative(const int& ... orders) = 0;

	/** Compute composition of the instance by a univariate function.
	 * @param f array of value and derivatives of the function at
	 * the current point (i.e. [f({@link #get_value()}), * f'({@link #get_value()}), f''({@link #get_value()})...]).
	 * @return f(this)
	 * @exception  if the number of derivatives
	 * in the array is not equal to {@link #get_order() order} + 1
	 */
	virtual T compose(double... f) = 0;
};