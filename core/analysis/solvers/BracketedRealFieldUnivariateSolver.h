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

  //package org.hipparchus.analysis.solvers;

  //import org.hipparchus.Calculus_Field_Element;
  //import org.hipparchus.analysis.Calculus_Field_Univariate_Function;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.exception.Math_Runtime_Exception;
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/** Interface for {@link Univariate_Solver (univariate real) root-finding
 * algorithms} that maintain a bracketed solution. There are several advantages
 * to having such root-finding algorithms:
 * <ul>
 *  <li>The bracketed solution guarantees that the root is kept within the
 *      interval. As such, these algorithms generally also guarantee
 *      convergence.</li>
 *  <li>The bracketed solution means that we have the opportunity to only
 *      return roots that are greater than or equal to the actual root, or
 *      are less than or equal to the actual root. That is, we can control
 *      whether under-approximations and over-approximations are
 *      {@link Allowed_Solution allowed solutions}. Other root-finding
 *      algorithms can usually only guarantee that the solution (the root that
 *      was found) is around the actual root.</li>
 * </ul>
 *
 * <p>For backwards compatibility, all root-finding algorithms must have
 * {@link Allowed_Solution#ANY_SIDE ANY_SIDE} as default for the allowed
 * solutions.</p>
 *
 * @see Allowed_Solution
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Bracketed_Real_Field_Univariate_Solver
{
	/**
	 * Get the maximum number of function evaluations.
	 *
	 * @return the maximum number of function evaluations.
	 */
	virtual int get_max_evaluations();

	/**
	 * Get the number of evaluations of the objective function.
	 * The number of evaluations corresponds to the last call to the
	 * {@code optimize} method. It is 0 if the method has not been
	 * called yet.
	 *
	 * @return the number of evaluations of the objective function.
	 */
	virtual int get_evaluations();

	/**
	 * Get the absolute accuracy of the solver.  Solutions returned by the
	 * solver should be accurate to this tolerance, i.e., if &epsilon; is the
	 * absolute accuracy of the solver and {@code v} is a value returned by
	 * one of the {@code solve} methods, then a root of the function should
	 * exist somewhere in the interval ({@code v} - &epsilon;, {@code v} + &epsilon;).
	 *
	 * @return the absolute accuracy.
	 */
	virtual T get_absolute_accuracy();

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
	virtual T get_relative_accuracy();

	/**
	 * Get the function value accuracy of the solver.  If {@code v} is
	 * a value returned by the solver for a function {@code f}, * then by contract, {@code |f(v)|} should be less than or equal to
	 * the function value accuracy configured for the solver.
	 *
	 * @return the function value accuracy.
	 */
	virtual T get_function_value_accuracy();

	/**
	 * Solve for a zero in the given interval.
	 * A solver may require that the interval brackets a single zero root.
	 * Solvers that do require bracketing should be able to handle the case
	 * where one of the endpoints is itself a root.
	 *
	 * @param max_eval Maximum number of evaluations.
	 * @param f Function to solve.
	 * @param min Lower bound for the interval.
	 * @param max Upper bound for the interval.
	 * @param allowed_solution The kind of solutions that the root-finding algorithm may
	 * accept as solutions.
	 * @return A value where the function is zero.
	 * @org.hipparchus.exception.
	 * if the arguments do not satisfy the requirements specified by the solver.
	 * @org.hipparchus.exception.Math_Illegal_State_Exception if
	 * the allowed number of evaluations is exceeded.
	 */
	virtual T solve(const int& max_eval, const Calculus_Field_Univariate_Function<T>& f, const T& min, const T& max, const Allowed_Solution& allowed_solution);

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
	 * @param allowed_solution The kind of solutions that the root-finding algorithm may
	 * accept as solutions.
	 * @return A value where the function is zero.
	 * @org.hipparchus.exception.
	 * if the arguments do not satisfy the requirements specified by the solver.
	 * @org.hipparchus.exception.Math_Illegal_State_Exception if
	 * the allowed number of evaluations is exceeded.
	 */
	virtual T solve(const int& max_eval, const Calculus_Field_Univariate_Function<T>& f, const T& min, const T& max, const T& start_value, const Allowed_Solution& allowed_solution);

	/**
	 * Solve for a zero in the given interval and return a tolerance interval surrounding
	 * the root.
	 *
	 * <p> It is required that the starting interval brackets a root.
	 *
	 * @param max_eval Maximum number of evaluations.
	 * @param f       Function to solve.
	 * @param min     Lower bound for the interval. f(min) != 0.0.
	 * @param max     Upper bound for the interval. f(max) != 0.0.
	 * @return an interval [ta, tb] such that for some t in [ta, tb] f(t) == 0.0 or has a
	 * step wise discontinuity that crosses zero. Both end points also satisfy the
	 * convergence criteria so either one could be used as the root. That is the interval
	 * satisfies the condition (| tb - ta | &lt;= {@link #get_absolute_accuracy() absolute}
	 * accuracy + max(ta, tb) * {@link #get_relative_accuracy() relative} accuracy) or (
	 * max(|f(ta)|, |f(tb)|) &lt;= {@link #get_function_value_accuracy()}) or there are no
	 * numbers in the field between ta and tb. The width of the interval (tb - ta) may be
	 * zero.
	 * @ if the arguments do not satisfy the
	 *                                      requirements specified by the solver.
	 * @Math_Illegal_State_Exception    if the allowed number of evaluations is
	 *                                      exceeded.
	 */
	default Interval<T> solve_interval(const int& max_eval, const Calculus_Field_Univariate_Function<T>& f, const T& min, const T& max)
	{
		return solve_interval(max_eval, f, min, max, min.add(max.subtract(min).multiply(0.5)));
	}

	/**
	 * Solve for a zero in the given interval and return a tolerance interval surrounding
	 * the root.
	 *
	 * <p> It is required that the starting interval brackets a root.
	 *
	 * @param max_eval    Maximum number of evaluations.
	 * @param start_value start value to use.
	 * @param f          Function to solve.
	 * @param min        Lower bound for the interval. f(min) != 0.0.
	 * @param max        Upper bound for the interval. f(max) != 0.0.
	 * @return an interval [ta, tb] such that for some t in [ta, tb] f(t) == 0.0 or has a
	 * step wise discontinuity that crosses zero. Both end points also satisfy the
	 * convergence criteria so either one could be used as the root. That is the interval
	 * satisfies the condition (| tb - ta | &lt;= {@link #get_absolute_accuracy() absolute}
	 * accuracy + max(ta, tb) * {@link #get_relative_accuracy() relative} accuracy) or (
	 * max(|f(ta)|, |f(tb)|) &lt;= {@link #get_function_value_accuracy()}) or numbers in the
	 * field between ta and tb. The width of the interval (tb - ta) may be zero.
	 * @ if the arguments do not satisfy the
	 *                                      requirements specified by the solver.
	 * @Math_Illegal_State_Exception    if the allowed number of evaluations is
	 *                                      exceeded.
	 */
	Interval<T> solve_interval(const int& max_eval, const Calculus_Field_Univariate_Function<T>& f, const T& min, const T& max, const T& start_value);

	/**
	 * An interval of a function that brackets a root.
	 * <p>
	 * Contains two end points and the value of the function at the two end points.
	 *
	 * @see #solve_interval(int, Calculus_Field_Univariate_Function, Calculus_Field_Element, * Calculus_Field_Element)
	 * @param <T> the element type
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	class Interval
	{
	private:
		/** Abscissa on the left end of the interval. */
		const T my_left_abscissa;
		/** Function value at {@link #left_abscissa}. */
		const T my_left_value;
		/** Abscissa on the right end of the interval, >= {@link #left_abscissa}. */
		const T my_right_abscissa;
		/** Function value at {@link #right_abscissa}. */
		const T my_right_value;

	public:
		/**
		 * Construct a interval with the given end points.
		 *
		 * @param left_abscissa  is the abscissa value at the left side of the interval.
		 * @param left_value     is the function value at {@code left_abscissa}.
		 * @param right_abscissa is the abscissa value on the right side of the interval.
		 *                      Must be greater than or equal to {@code left_abscissa}.
		 * @param right_value    is the function value at {@code right_abscissa}.
		 */
		Interval(const T& left_abscissa, const T& left_value, const T& right_abscissa, const T& right_value)
			:
			my_left_abscissa{ left_abscissa },
			my_left_value{ left_value },
			my_right_abscissa{ right_abscissa },
			my_right_value{ right_value }
		{
		}

		/**
		 * Get the left abscissa.
		 *
		 * @return abscissa of the start of the interval.
		 */
		T get_left_abscissa() const
		{
			return my_left_abscissa;
		}

		/**
		 * Get the right abscissa.
		 *
		 * @return abscissa of the end of the interval.
		 */
		T get_right_abscissa() const
		{
			return my_right_abscissa;
		}

		/**
		 * Get the function value at {@link #get_left_abscissa()}.
		 *
		 * @return value of the function at the start of the interval.
		 */
		T get_left_value() const
		{
			return my_left_value;
		}

		/**
		 * Get the function value at {@link #get_right_abscissa()}.
		 *
		 * @return value of the function at the end of the interval.
		 */
		T get_right_value() const
		{
			return my_right_value;
		}

		/**
		 * Get the abscissa corresponding to the allowed side.
		 *
		 * @param allowed side of the root.
		 * @return the abscissa on the selected side of the root.
		 */
		T get_side(const Allowed_Solution& allowed)
		{
			const T x_a = get_left_abscissa();
			const T y_a = get_left_value();
			const T x_b = get_right_abscissa();
			switch (allowed)
			{
			case ANY_SIDE:
				const T abs_ya = get_left_value().abs();
				const T abs_y_b = get_right_value().abs();
				return abs_ya.subtract(abs_y_b).get_real() < 0
					? x_a
					: x_b;
			case LEFT_SIDE:
				return x_a;
			case RIGHT_SIDE:
				return x_b;
			case BELOW_SIDE:
				return (y_a.get_real() <= 0)
					? x_a
					: x_b;
			case ABOVE_SIDE:
				return (y_a.get_real() < 0)
					? x_b
					: x_a;
			default:
				throw std::exception("not implemented");
				// this should never happen
				//throw Math_Runtime_Exception.create_internal_error();
			}
		}
	}
};