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

#include "../CalculusFieldUnivariateFunction.h"
#include <type_traits>

 /**
  * Interface for univariate real integration algorithms.
  * @param <T> Type of the field elements.
  * @since 2.0
  */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Univariate_Integrator
{
	/**
	 * Get the relative accuracy.
	 *
	 * @return the accuracy
	 */
	virtual double get_relative_accuracy() = 0;

	/**
	 * Get the absolute accuracy.
	 *
	 * @return the accuracy
	 */
	virtual double get_absolute_accuracy() = 0;

	/**
	 * Get the min limit for the number of iterations.
	 *
	 * @return the actual min limit
	 */
	virtual int get_minimal_iteration_count() = 0;

	/**
	 * Get the upper limit for the number of iterations.
	 *
	 * @return the actual upper limit
	 */
	virtual int get_maximal_iteration_count() = 0;

	/**
	 * Integrate the function in the given interval.
	 *
	 * @param max_eval Maximum number of evaluations.
	 * @param f the integrand function
	 * @param min the lower bound for the interval
	 * @param max the upper bound for the interval
	 * @return the value of integral
	 * @Math_Illegal_State_Exception if the maximum number of function
	 * evaluations is exceeded
	 * @Math_Illegal_State_Exception if the maximum iteration count is exceeded
	 * or the integrator detects convergence problems otherwise
	 * @ if {@code min > max} or the endpoints do not
	 * satisfy the requirements specified by the integrator
	 * @ if {@code f} is {@code NULL}.
	 */
	virtual T integrate(const int& max_eval, Calculus_Field_Univariate_Function<T> f, T min, T max) = 0;

	/**
	 * Get the number of function evaluations of the last run of the integrator.
	 *
	 * @return number of function evaluations
	 */
	virtual int get_evaluations() = 0;

	/**
	 * Get the number of iterations of the last run of the integrator.
	 *
	 * @return number of iterations
	 */
	virtual int get_iterations() = 0;
};