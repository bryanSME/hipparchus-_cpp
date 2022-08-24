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
  //import org.hipparchus.Field;
  //import org.hipparchus.analysis.Calculus_Field_Univariate_Function;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Incrementor;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include <vector>
#include "../../CalculusFieldElement.hpp"

/**
 * This class : a modification of the <a
 * href="http://mathworld.wolfram.com/Brents_method.html"> Brent algorithm</a>.
 * <p>
 * The changes with respect to the original Brent algorithm are:
 * <ul>
 *   <li>the returned value is chosen in the current interval according
 *   to user specified {@link Allowed_Solution}</li>
 *   <li>the maximal order for the invert polynomial root search is
 *   user-specified instead of being invert quadratic only</li>
 * </ul><p>
 * The given interval must bracket the root.</p>
 *
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldBracketing_Nth_Order_Brent_Solver : public Bracketed_Real_Field_Univariate_Solver<T>
{
private:
	/** Maximal aging triggering an attempt to balance the bracketing interval. */
	static constexpr int MAXIMAL_AGING{ 2 };

	/** Field to which the elements belong. */
	const Field<T> my_field;

	/** Maximal order. */
	const int my_maximal_order;

	/** Function value accuracy. */
	const T my_function_value_accuracy;

	/** Absolute accuracy. */
	const T my_absolute_accuracy;

	/** Relative accuracy. */
	const T my_relative_accuracy;

	/** Evaluations counter. */
	Incrementor my_evaluations;

	/** Guess an x value by n<sup>th</sup> order inverse polynomial interpolation.
	 * <p>
	 * The x value is guessed by evaluating polynomial Q(y) at y = target_y, where Q
	 * is built such that for all considered points (x<sub>i</sub>, y<sub>i</sub>), * Q(y<sub>i</sub>) = x<sub>i</sub>.
	 * </p>
	 * @param target_y target value for y
	 * @param x reference points abscissas for interpolation, * note that this array <em>is</em> modified during computation
	 * @param y reference points ordinates for interpolation
	 * @param start start index of the points to consider (inclusive)
	 * @param end end index of the points to consider (exclusive)
	 * @return guessed root (will be a NaN if two points share the same y)
	 */
	T guess_x(const T& target_y, const std::vector<T>& x, const std::vector<T>& y, const int& start, const int& end)
	{
		// compute Q Newton coefficients by divided differences
		for (int i{ start }; i < end - 1; ++i)
		{
			const int delta = i + 1 - start;
			for (int j = end - 1; j > i; --j)
			{
				x[j] = x[j].subtract(x[j - 1]).divide(y[j].subtract(y[j - delta]));
			}
		}

		// evaluate Q(target_y)
		T x0 = field.get_zero();
		for (int j{ end - 1 }; j >= start; --j)
		{
			x0 = x[j].add(x0.multiply(target_y.subtract(y[j])));
		}

		return x0;
	}

public:
	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 * @param function_value_accuracy Function value accuracy.
	 * @param maximal_order maximal order.
	 * @exception  if maximal order is lower than 2
	 */
	FieldBracketing_Nth_Order_Brent_Solver(const T& relative_accuracy, const T& absolute_accuracy, const T& function_value_accuracy, const int& maximal_order)
		:
		my_field{ relative_accuracy.get_field() },
		my_maximal_order{ maximal_order },
		my_absolute_accuracy{ absolute_accuracy },
		my_relative_accuracy{ relative_accuracy },
		my_function_value_accuracy{ function_value_accuracy },
		my_evaluations{ Incrementor() }
	{
		if (maximal_order < 2)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, maximal_order, 2);
		}
	}

	/** Get the maximal order.
	 * @return maximal order
	 */
	int get_maximal_order() const
	{
		return my_maximal_order;
	}

	/**
	 * Get the maximal number of function evaluations.
	 *
	 * @return the maximal number of function evaluations.
	 */
	 //override
	int get_max_evaluations()
	{
		return my_evaluations.get_maximal_count();
	}

	/**
	 * Get the number of evaluations of the objective function.
	 * The number of evaluations corresponds to the last call to the
	 * {@code optimize} method. It is 0 if the method has not been
	 * called yet.
	 *
	 * @return the number of evaluations of the objective function.
	 */
	 //override
	int get_evaluations()
	{
		return my_evaluations.get_count();
	}

	/**
	 * Get the absolute accuracy.
	 * @return absolute accuracy
	 */
	 //override
	T get_absolute_accuracy() const
	{
		return my_absolute_accuracy;
	}

	/**
	 * Get the relative accuracy.
	 * @return relative accuracy
	 */
	 //override
	T get_relative_accuracy()
	{
		return my_relative_accuracy;
	}

	/**
	 * Get the function accuracy.
	 * @return function accuracy
	 */
	 //override
	T get_function_value_accuracy() const
	{
		return my_function_value_accuracy;
	}

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
	 * @return a value where the function is zero.
	 * @exception  if f is NULL.
	 * @exception  if root cannot be bracketed
	 */
	 //override
	T solve(const int& max_eval, const Calculus_Field_Univariate_Function<T>& f, const T& min, const T& max, const Allowed_Solution& allowed_solution)
	{
		return solve(max_eval, f, min, max, min.add(max).divide(2), allowed_solution);
	}

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
	 * @return a value where the function is zero.
	 * @exception  if f is NULL.
	 * @exception  if root cannot be bracketed
	 */
	 //override
	T solve(const int& max_eval, const Calculus_Field_Univariate_Function<T>& f, const T& min, const T& max, const T& start_value, const Allowed_Solution& allowed_solution)
	{
		// find interval containing root
		return solve_interval(max_eval, f, min, max, start_value).get_side(allowed_solution);
	}

	/** {@inherit_doc} */
	//override
	Interval<T> solve_interval(const int& max_eval, Calculus_Field_Univariate_Function<T>& f, const T& min, const T& max, const T& start_value)
	{
		// Checks.
		//Math_Utils::check_not_null(f);

		// Reset.
		my_evaluations = my_evaluations.with_maximal_count(max_eval);
		T zero = field.get_zero();
		T nan = zero.add(Double.NaN);

		// prepare arrays with the first points
		const auto x = Math_Arrays::build_array(field, maximal_order + 1);
		const auto y = Math_Arrays::build_array(field, maximal_order + 1);
		x[0] = min;
		x[1] = start_value;
		x[2] = max;

		// evaluate initial guess
		my_evaluations.increment();
		y[1] = f.value(x[1]);
		if (y[1].get_real() == 0.0)
		{
			// return the initial guess if it is a perfect root.
			return Interval<>(x[1], y[1], x[1], y[1]);
		}

		// evaluate first endpoint
		my_evaluations.increment();
		y[0] = f.value(x[0]);
		if (y[0].get_real() == 0.0)
		{
			// return the first endpoint if it is a perfect root.
			return Interval<>(x[0], y[0], x[0], y[0]);
		}

		int nb_points;
		int sign_change_index;
		if (y[0].multiply(y[1]).get_real() < 0)
		{
			// reduce interval if it brackets the root
			nb_points = 2;
			sign_change_index = 1;
		}
		else
		{
			// evaluate second endpoint
			evaluations.increment();
			y[2] = f.value(x[2]);
			if (y[2].get_real() == 0.0)
			{
				// return the second endpoint if it is a perfect root.
				return Interval<>(x[2], y[2], x[2], y[2]);
			}

			if (y[1].multiply(y[2]).get_real() < 0)
			{
				// use all computed point as a start sampling array for solving
				nb_points = 3;
				sign_change_index = 2;
			}
			else
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_BRACKETING_INTERVAL, x[0].get_real(), x[2].get_real(), y[0].get_real(), y[2].get_real());
			}
		}

		// prepare a work array for inverse polynomial interpolation
		const std::vector<T> tmp_x = Math_Arrays::build_array(field, x.size());

		// current tightest bracketing of the root
		T x_a = x[sign_change_index - 1];
		T y_a = y[sign_change_index - 1];
		T abs_x_a = x_a.abs();
		T abs_ya = y_a.abs();
		int aging_a{};
		T x_b = x[sign_change_index];
		T yB = y[sign_change_index];
		T abs_x_b = x_b.abs();
		T abs_y_b = yB.abs();
		int aging_b{};

		// search loop
		while (true)
		{
			// check convergence of bracketing interval
			T max_x = abs_x_a.subtract(abs_x_b).get_real() < 0 ? abs_x_b : abs_x_a;
			T max_y = abs_ya.subtract(abs_y_b).get_real() < 0 ? abs_y_b : abs_ya;
			const T x_tol = absolute_accuracy.add(relative_accuracy.multiply(max_x));
			const T midpoint = x_a.add(x_b.subtract(x_a).divide(2));
			if (x_b.subtract(x_a).subtract(x_tol).get_real() <= 0 ||
				max_y.subtract(function_value_accuracy).get_real() < 0 ||
				x_a.equals(midpoint) || x_b.equals(midpoint))
			{
				return Interval<>(x_a, y_a, x_b, yB);
			}

			// target for the next evaluation point
			T target_y;
			if (aging_a >= MAXIMAL_AGING)
			{
				// we keep updating the high bracket, try to compensate this
				target_y = yB.divide(16).negate();
			}
			else if (aging_b >= MAXIMAL_AGING)
			{
				// we keep updating the low bracket, try to compensate this
				target_y = y_a.divide(16).negate();
			}
			else
			{
				// bracketing is balanced, try to find the root itself
				target_y = zero;
			}

			// make a few attempts to guess a root, T next_x;
			int start{};
			int end = nb_points;
			do
			{
				// guess a value for current target, using inverse polynomial interpolation
				System.arraycopy(x, start, tmp_x, start, end - start);
				next_x = guess_x(target_y, tmp_x, y, start, end);

				if (!((next_x.subtract(x_a).get_real() > 0) && (next_x.subtract(x_b).get_real() < 0)))
				{
					// the guessed root is not strictly inside of the tightest bracketing interval

					// the guessed root is either not strictly inside the interval or it
					// is a NaN (which occurs when some sampling points share the same y)
					// we try again with a lower interpolation order
					if (sign_change_index - start >= end - sign_change_index)
					{
						// we have more points before the sign change, drop the lowest point
						++start;
					}
					else
					{
						// we have more points after sign change, drop the highest point
						--end;
					}

					// we need to do one more attempt
					next_x = nan;
				}
			} while (std::isnan(next_x.get_real()) && (end - start > 1));

			if (std::isnan(next_x.get_real()))
			{
				// fall back to bisection
				next_x = x_a.add(x_b.subtract(x_a).divide(2));
				start = sign_change_index - 1;
				end = sign_change_index;
			}

			// evaluate the function at the guessed root
			evaluations.increment();
			const T next_y = f.value(next_x);
			if (next_y.get_real() == 0.0)
			{
				// we have found an exact root, since it is not an approximation
				// we don't need to bother about the allowed solutions setting
				return Interval<>(next_x, next_y, next_x, next_y);
			}

			if ((nb_points > 2) && (end - start != nb_points))
			{
				// we have been forced to ignore some points to keep bracketing, // they are probably too far from the root, drop them from now on
				nb_points = end - start;
				System.arraycopy(x, start, x, 0, nb_points);
				System.arraycopy(y, start, y, 0, nb_points);
				sign_change_index -= start;
			}
			else  if (nb_points == x.size())
			{
				// we have to drop one point in order to insert the one
				nb_points--;

				// keep the tightest bracketing interval as centered as possible
				if (sign_change_index >= (x.size() + 1) / 2)
				{
					// we drop the lowest point, we have to shift the arrays and the index
					System.arraycopy(x, 1, x, 0, nb_points);
					System.arraycopy(y, 1, y, 0, nb_points);
					--sign_change_index;
				}
			}

			// insert the last computed point
			//(by construction, we know it lies inside the tightest bracketing interval)
			System.arraycopy(x, sign_change_index, x, sign_change_index + 1, nb_points - sign_change_index);
			x[sign_change_index] = next_x;
			System.arraycopy(y, sign_change_index, y, sign_change_index + 1, nb_points - sign_change_index);
			y[sign_change_index] = next_y;
			++nb_points;

			// update the bracketing interval
			if (next_y.multiply(y_a).get_real() <= 0)
			{
				// the sign change occurs before the inserted point
				x_b = next_x;
				yB = next_y;
				abs_y_b = yB.abs();
				++aging_a;
				aging_b = 0;
			}
			else
			{
				// the sign change occurs after the inserted point
				x_a = next_x;
				y_a = next_y;
				abs_ya = y_a.abs();
				aging_a = 0;
				++aging_b;

				// update the sign change index
				sign_change_index++;
			}
		}
	}
};