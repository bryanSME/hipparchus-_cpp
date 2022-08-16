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
#include "../UnivariateFunction.h"

  /**
   * Interface for (univariate real) rootfinding algorithms.
   * Implementations will search for only one zero in the given interval.
   *
   * This class is not intended for use outside of the Hipparchus
   * library, regular user should rely on more specific interfaces like
   * {@link Univariate_Solver}, {@link Polynomial_Solver} or {@link
   * Univariate_Differentiable_Solver}.
   * @param <F> Type of function to solve.
   *
   * @see Univariate_Solver
   * @see Polynomial_Solver
   * @see Univariate_Differentiable_Solver
   */
template<typename F, typename std::enable_if<std::is_base_of<Univariate_Function, F>::value>::type* = nullptr>
class Base_Univariate_Solver
{
	/**
	 * Get the maximum number of function evaluations.
	 *
	 * @return the maximum number of function evaluations.
	 */
	virtual int get_max_evaluations() = 0;

	/**
	 * Get the number of evaluations of the objective function.
	 * The number of evaluations corresponds to the last call to the
	 * {@code optimize} method. It is 0 if the method has not been
	 * called yet.
	 *
	 * @return the number of evaluations of the objective function.
	 */
	virtual int get_evaluations() = 0;

	/**
	 * Get the absolute accuracy of the solver.  Solutions returned by the
	 * solver should be accurate to this tolerance, i.e., if &epsilon; is the
	 * absolute accuracy of the solver and {@code v} is a value returned by
	 * one of the {@code solve} methods, then a root of the function should
	 * exist somewhere in the interval ({@code v} - &epsilon;, {@code v} + &epsilon;).
	 *
	 * @return the absolute accuracy.
	 */
	virtual double get_absolute_accuracy() = 0;

	/**
	 * Get the relative accuracy of the solver.  The contract for relative
	 * accuracy is the same as {@link #get_absolute_accuracy()}, but using
	 * relative, rather than absolute error.  If &rho; is the relative accuracy
	 * configured for a solver and {@code v} is a value returned, then a root
	 * of the function should exist somewhere in the interval
	 * ({@code v} - &rho; {@code v}, {@code v} + &rho; {@code v}).
	 *
	 * @return the relative accuracy.
	 */
	virtual double get_relative_accuracy() = 0;

	/**
	 * Get the function value accuracy of the solver.  If {@code v} is
	 * a value returned by the solver for a function {@code f}, * then by contract, {@code |f(v)|} should be less than or equal to
	 * the function value accuracy configured for the solver.
	 *
	 * @return the function value accuracy.
	 */
	virtual double get_function_value_accuracy() = 0;

	/**
	 * Solve for a zero root in the given interval.
	 * A solver may require that the interval brackets a single zero root.
	 * Solvers that do require bracketing should be able to handle the case
	 * where one of the endpoints is itself a root.
	 *
	 * @param max_eval Maximum number of evaluations.
	 * @param f Function to solve.
	 * @param min Lower bound for the interval.
	 * @param max Upper bound for the interval.
	 * @return a value where the function is zero.
	 * @
	 * if the arguments do not satisfy the requirements specified by the solver.
	 * @Math_Illegal_State_Exception if
	 * the allowed number of evaluations is exceeded.
	 */
	virtual double solve(const int& max_eval, F f, const double& min, const double& max) = 0;

	/**
	 * Solve for a zero in the given interval, start at {@code start_value}.
	 * A solver may require that the interval brackets a single zero root.
	 * Solvers that do require bracketing should be able to handle the case
	 * where one of the endpoints is itself a root.
	 *
	 * @param max_eval Maximum number of evaluations.
	 * @param f Function to solve.
	 * @param min Lower bound for the interval.
	 * @param max Upper bound for the interval.
	 * @param start_value Start value to use.
	 * @return a value where the function is zero.
	 * @
	 * if the arguments do not satisfy the requirements specified by the solver.
	 * @Math_Illegal_State_Exception if
	 * the allowed number of evaluations is exceeded.
	 */
	virtual double solve(const int& max_eval, const F& f, const double& min, const double& max, const double& start_value) = 0;

	/**
	 * Solve for a zero in the vicinity of {@code start_value}.
	 *
	 * @param f Function to solve.
	 * @param start_value Start value to use.
	 * @return a value where the function is zero.
	 * @param max_eval Maximum number of evaluations.
	 * @org.hipparchus.exception.
	 * if the arguments do not satisfy the requirements specified by the solver.
	 * @org.hipparchus.exception.Math_Illegal_State_Exception if
	 * the allowed number of evaluations is exceeded.
	 */
	virtual double solve(const int& max_eval, const F& f, const double& start_value) = 0;
};