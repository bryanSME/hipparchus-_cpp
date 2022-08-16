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

  //import org.hipparchus.analysis.Univariate_Function;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Precision;
#include <vector>
#include "../../exception/LocalizedCoreFormats.h"
#include "AllowedSolution.h"
#include "AbstractUnivariateSolver.h"
#include "BracketedUnivariateSolver.hpp"

/**
 * This class : a modification of the <a
 * href="http://mathworld.wolfram.com/Brents_method.html"> Brent algorithm</a>.
 * <p>
 * The changes with respect to the original Brent algorithm are:
 * <ul>
 *   <li>the returned value is chosen in the current interval according
 *   to user specified {@link Allowed_Solution},</li>
 *   <li>the maximal order for the invert polynomial root search is
 *   user-specified instead of being invert quadratic only</li>
 * </ul><p>
 * The given interval must bracket the root.</p>
 *
 */
class Bracketing_Nth_Order_Brent_Solver : public Abstract_Univariate_Solver, public Bracketed_Univariate_Solver<Univariate_Function>
{
private:
	/** Default absolute accuracy. */
	static constexpr double DEFAULT_ABSOLUTE_ACCURACY{ 1e-6 };

	/** Default maximal order. */
	static constexpr int DEFAULT_MAXIMAL_ORDER{ 5 };

	/** Maximal aging triggering an attempt to balance the bracketing interval. */
	static constexpr int MAXIMAL_AGING{ 2 };

	/** Reduction factor for attempts to balance the bracketing interval. */
	static constexpr double REDUCTION_FACTOR{ 1.0 / 16.0 };

	/** Maximal order. */
	const int my_maximal_order;

	/** The kinds of solutions that the algorithm may accept. */
	Allowed_Solution my_allowed;

public:
	/**
	 * Construct a solver with default accuracy and maximal order (1e-6 and 5 respectively)
	 */
	Bracketing_Nth_Order_Brent_Solver()
	{
		Bracketing_Nth_Order_Brent_Solver(DEFAULT_ABSOLUTE_ACCURACY, DEFAULT_MAXIMAL_ORDER);
	}

	/**
	 * Construct a solver.
	 *
	 * @param absolute_accuracy Absolute accuracy.
	 * @param maximal_order maximal order.
	 * @exception  if maximal order is lower than 2
	 */
	Bracketing_Nth_Order_Brent_Solver(const double& absolute_accuracy, const int& maximal_order)
		: my_maximal_order{ maximal_order }, my_allowed{ Allowed_Solution::ANY_SIDE }
	{
		super(absolute_accuracy);
		if (maximal_order < 2)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, maximal_order, 2);
		}
	}

	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 * @param maximal_order maximal order.
	 * @exception  if maximal order is lower than 2
	 */
	Bracketing_Nth_Order_Brent_Solver(const double& relative_accuracy, const double& absolute_accuracy, const int maximal_order)
	{
		super(relative_accuracy, absolute_accuracy);
		if (maximal_order < 2)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, maximal_order, 2);
		}
		this.maximal_order = maximal_order;
		this.allowed = Allowed_Solution.ANY_SIDE;
	}

	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 * @param function_value_accuracy Function value accuracy.
	 * @param maximal_order maximal order.
	 * @exception  if maximal order is lower than 2
	 */
	Bracketing_Nth_Order_Brent_Solver(const double& relative_accuracy, const double& absolute_accuracy, const double function_value_accuracy, const int maximal_order)
	{
		super(relative_accuracy, absolute_accuracy, function_value_accuracy);
		if (maximal_order < 2)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, maximal_order, 2);
		}
		my_maximal_order = maximal_order;
		my_allowed = Allowed_Solution.ANY_SIDE;
	}

	/** Get the maximal order.
	 * @return maximal order
	 */
	int get_maximal_order() const
	{
		return my_maximal_order;
	}

	/** {@inherit_doc} */
	//override
	double solve(const int& max_eval, const Univariate_Function& f, const double& min, const double& max, Allowed_Solution& allowed_solution)
	{
		my_allowed = allowed_solution;
		return super.solve(max_eval, f, min, max);
	}

	/** {@inherit_doc} */
	//override
	double solve(const int& max_eval, const Univariate_Function& f, const double& min, const double& max, const double& start_value, Allowed_Solution& allowed_solution)
	{
		my_allowed = allowed_solution;
		return super.solve(max_eval, f, min, max, start_value);
	}

	/** {@inherit_doc} */
	//override
	Interval solve_interval(const int& max_eval, const Univariate_Function& f, const double& min, const double& max, const double& start_value)
	{
		setup(max_eval, f, min, max, start_value);
		my_allowed = NULL;
		return do_solve_interval();
	}

protected:
	/**
	 * {@inherit_doc}
	 */
	 //override
	double do_solve()
	{
		return do_solve_interval().get_side(allowed);
	}

	/**
	 * Find a root and return the containing interval.
	 *
	 * @return an interval containing the root such that both end points meet the
	 * convergence criteria.
	 */
	Interval do_solve_interval()
	{
		// prepare arrays with the first points
		auto x = std::vector<double>(my_maximal_order + 1];
		auto y = std::vector<double>(my_maximal_order + 1];
		x[0] = get_min();
		x[1] = get_start_value();
		x[2] = get_max();
		verify_interval(x[0], x[2]);
		if (x[1] < x[0] || x[2] < x[1])
		{
			throw (
				hipparchus::exception::Localized_Core_Formats_Type::START_POINT_NOT_IN_INTERVAL, x[1], x[0], x[2]);
		}

		// evaluate initial guess
		y[1] = compute_objective_value(x[1]);
		if (y[1] == 0.0)
		{
			// return the initial guess if it is a perfect root.
			return Interval(x[1], y[1], x[1], y[1]);
		}

		// evaluate first  endpoint
		y[0] = compute_objective_value(x[0]);
		if (y[0] == 0.0)
		{
			// return the first endpoint if it is a perfect root.
			return Interval(x[0], y[0], x[0], y[0]);
		}

		int nb_points;
		int sign_change_index;
		if (y[0] * y[1] < 0)
		{
			// reduce interval if it brackets the root
			nb_points = 2;
			sign_change_index = 1;
		}
		else
		{
			// evaluate second endpoint
			y[2] = compute_objective_value(x[2]);
			if (y[2] == 0.0)
			{
				// return the second endpoint if it is a perfect root.
				return Interval(x[2], y[2], x[2], y[2]);
			}

			if (y[1] * y[2] < 0)
			{
				// use all computed point as a start sampling array for solving
				nb_points = 3;
				sign_change_index = 2;
			}
			else
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_BRACKETING_INTERVAL, x[0], x[2], y[0], y[2]);
			}
		}

		// prepare a work array for inverse polynomial interpolation
		const std::vector<double> tmp_x = std::vector<double>(x.size()];

		// current tightest bracketing of the root
		double x_a = x[sign_change_index - 1];
		double y_a = y[sign_change_index - 1];
		double abs_ya = std::abs(y_a);
		int aging_a{};
		double x_b = x[sign_change_index];
		double yB = y[sign_change_index];
		double abs_y_b = std::abs(yB);
		int aging_b{};

		// search loop
		while (true)
		{
			// check convergence of bracketing interval
			const double x_tol = get_absolute_accuracy() +
				get_relative_accuracy() * std::max(std::abs(x_a), std::abs(x_b));
			if (x_b - x_a <= x_tol ||
				std::max(abs_ya, abs_y_b) < get_function_value_accuracy() ||
				Precision::equals(x_a, x_b, 1))
			{
				return Interval(x_a, y_a, x_b, yB);
			}

			// target for the next evaluation point
			double target_y;
			if (aging_a >= MAXIMAL_AGING)
			{
				// we keep updating the high bracket, try to compensate this
				const int p{ aging_a - MAXIMAL_AGING };
				const double weight_a = (1 << p) - 1;
				const double weight_b = p + 1;
				target_y = (weight_a * y_a - weight_b * REDUCTION_FACTOR * yB) / (weight_a + weight_b);
			}
			else if (aging_b >= MAXIMAL_AGING)
			{
				// we keep updating the low bracket, try to compensate this
				const int p = aging_b - MAXIMAL_AGING;
				const double weight_a{ p + 1 };
				const double weight_b{ (1 << p) - 1 };
				target_y = (weight_b * yB - weight_a * REDUCTION_FACTOR * y_a) / (weight_a + weight_b);
			}
			else
			{
				// bracketing is balanced, try to find the root itself
				target_y = 0;
			}

			// make a few attempts to guess a root, double next_x;
			int start = 0;
			int end = nb_points;
			do
			{
				// guess a value for current target, using inverse polynomial interpolation
				System.arraycopy(x, start, tmp_x, start, end - start);
				next_x = guess_x(target_y, tmp_x, y, start, end);

				if (!((next_x > x_a) && (next_x < x_b)))
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
					next_x = std::numeric_limits<double>::quiet_NaN();
				}
			} while (std::isnan(next_x) && (end - start > 1));

			if (std::isnan(next_x))
			{
				// fall back to bisection
				next_x = x_a + 0.5 * (x_b - x_a);
				start = sign_change_index - 1;
				end = sign_change_index;
			}

			// evaluate the function at the guessed root
			const double next_y = compute_objective_value(next_x);
			if (next_y == 0.0 || std::abs(next_y) < get_function_value_accuracy() && allowed == Allowed_Solution.ANY_SIDE)
			{
				// we have either:
				// - an exact root, so we don't we don't need to bother about the allowed solutions setting
				// - or an approximate root and we know allowed solutions setting if to retrieve the value closest to zero
				return Interval(next_x, next_y, next_x, next_y);
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
			if (next_y * y_a <= 0)
			{
				// the sign change occurs before the inserted point
				x_b = next_x;
				yB = next_y;
				abs_y_b = std::abs(yB);
				++aging_a;
				aging_b = 0;
			}
			else
			{
				// the sign change occurs after the inserted point
				x_a = next_x;
				y_a = next_y;
				abs_ya = std::abs(y_a);
				aging_a = 0;
				++aging_b;

				// update the sign change index
				sign_change_index++;
			}
		}
	}

private:
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
	double guess_x(const double target_y, const std::vector<double> x, const std::vector<double> y, const int start, const int end)
	{
		// compute Q Newton coefficients by divided differences
		for (int i = start; i < end - 1; ++i)
		{
			const int delta = i + 1 - start;
			for (int j = end - 1; j > i; --j)
			{
				x[j] = (x[j] - x[j - 1]) / (y[j] - y[j - delta]);
			}
		}

		// evaluate Q(target_y)
		double x0{};
		for (int j = end - 1; j >= start; --j)
		{
			x0 = x[j] + x0 * (target_y - y[j]);
		}

		return x0;
	}
};