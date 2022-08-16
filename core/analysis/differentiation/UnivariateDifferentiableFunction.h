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
#include "../UnivariateFunction.h"
#include "Derivative.h"
#include <type_traits>

  /** Interface for univariate functions derivatives.
   * <p>This interface represents a simple function which computes
   * both the value and the first derivative of a mathematical function.
   * The derivative is computed with respect to the input variable.</p>
   * @see Univariate_Differentiable_Function
   * @see Univariate_Function_differentiator
   */
class Univariate_Differentiable_Function : public Univariate_Function
{
public:
	/**
	 * Compute the value for the function.
	 * @param x the point for which the function value should be computed
	 * @param <T> the type of the field elements
	 * @return the value
	 * @exception  if {@code x} does not
	 * satisfy the function's constraints (argument out of bound, or unsupported
	 * derivative order for example)
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Derivative<T>, T>::value>::type* = nullptr>
	T value(T x);
};