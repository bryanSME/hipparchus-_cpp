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
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.util.FastMath;
#include "AbstractUnivariateSolver.h"
#include "BracketedUnivariateSolver.hpp"
#include "AllowedSolution.h"


/**
 * Base class for all bracketing <em>Secant</em>-based methods for root-finding
 * (approximating a zero of a univariate real function).
 *
 * <p>Implementation of the {@link Regula_Falsi_Solver <em>Regula Falsi</em>} and
 * {@link Illinois_Solver <em>Illinois</em>} methods is based on the
 * following article: M. Dowell and P. Jarratt, * <em>A modified regula falsi method for computing the root of an
 * equation</em>, BIT Numerical Mathematics, volume 11, number 2, * pages 168-174, Springer, 1971.</p>
 *
 * <p>Implementation of the {@link Pegasus_Solver <em>Pegasus</em>} method is
 * based on the following article: M. Dowell and P. Jarratt, * <em>The "Pegasus" method for computing the root of an equation</em>, * BIT Numerical Mathematics, volume 12, number 4, pages 503-508, Springer, * 1972.</p>
 *
 * <p>The {@link Secant_Solver <em>Secant</em>} method is <em>not</em> a
 * bracketing method, so it is not implemented here. It has a separate
 * implementation.</p>
 *
 */
class Base_Secant_Solver : public Abstract_Univariate_Solver, public Bracketed_Univariate_Solver<Univariate_Function>
{
private:
	/** The kinds of solutions that the algorithm may accept. */
	Allowed_Solution my_allowed;

	/** The <em>Secant</em>-based root-finding method to use. */
	const Method my_method;

protected:
	/** Default absolute accuracy. */
	static constexpr double DEFAULT_ABSOLUTE_ACCURACY{ 1e-6 };

	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 * @param method <em>Secant</em>-based root-finding method to use.
	 */
	Base_Secant_Solver(const double& relative_accuracy, const double& absolute_accuracy, const Method& method)
		:
		my_allowed{ Allowed_Solution::ANY_SIDE }, my_method{ method }
	{
		Abstract_Univariate_Solver(relative_accuracy, absolute_accuracy);
	}

	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Maximum relative error.
	 * @param absolute_accuracy Maximum absolute error.
	 * @param function_value_accuracy Maximum function value error.
	 * @param method <em>Secant</em>-based root-finding method to use
	 */
	Base_Secant_Solver(const double& relative_accuracy, const double& absolute_accuracy, const double& function_value_accuracy, const Method& method)
		:
		my_allowed{ Allowed_Solution::ANY_SIDE },
		my_method{ method }
	{
		Abstract_Univariate_Solver(relative_accuracy, absolute_accuracy, function_value_accuracy);
	}

	/**
	 * {@inherit_doc}
	 *
	 * @Math_Illegal_State_Exception if the algorithm failed due to finite
	 * precision.
	 */
	 //override
	const double do_solve()
	{
		return do_solve_interval().get_side(my_allowed);
	}

	/**
	 * Find a root and return the containing interval.
	 *
	 * @return an interval containing the root such that the selected end point meets the
	 * convergence criteria.
	 * @Math_Illegal_State_Exception if convergence fails.
	 */
	const Interval do_solve_interval()
	{
		// Get initial solution
		double x0 = get_min();
		double x1 = get_max();
		double f0 = compute_objective_value(x0);
		double f1 = compute_objective_value(x1);

		// If one of the bounds is the exact root, return it. sin_ce these are
		// not under-approximations or over-approximations, we can return them
		// regardless of the allowed solutions.
		if (f0 == 0.0)
		{
			return Interval(x0, f0, x0, f0);
		}
		if (f1 == 0.0)
		{
			return Interval(x1, f1, x1, f1);
		}

		// Verify bracketing of initial solution.
		verify_bracketing(x0, x1);

		// Get accuracies.
		const double ftol = get_function_value_accuracy();
		const double atol = get_absolute_accuracy();
		const double rtol = get_relative_accuracy();

		// Keep track of inverted intervals, meaning that the left bound is
		// larger than the right bound.
		bool inverted = false;

		// Keep finding better approximations.
		while (true)
		{
			// Calculate the next approximation.
			const double x = x1 - ((f1 * (x1 - x0)) / (f1 - f0));
			const double fx = compute_objective_value(x);

			// If the approximation is the exact root, return it. sin_ce
			// this is not an under-approximation or an over-approximation, // we can return it regardless of the allowed solutions.
			if (fx == 0.0)
			{
				return Interval(x, fx, x, fx);
			}

			// Update the bounds with the approximation.
			if (f1 * fx < 0)
			{
				// The value of x1 has switched to the other bound, thus inverting
				// the interval.
				x0 = x1;
				f0 = f1;
				inverted = !inverted;
			}
			else
			{
				switch (my_method)
				{
				case ILLINOIS:
					f0 *= 0.5;
					break;
				case PEGASUS:
					f0 *= f1 / (f1 + fx);
					break;
				case REGULA_FALSI:
					// Detect early that algorithm is stuck, instead of waiting
					// for the maximum number of iterations to be exceeded.
					if (x == x1)
					{
						throw std::exception("not implemented");
						//throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CONVERGENCE_FAILED);
					}
					break;
				default:
					// Should never happen.
					throw std::exception("not implemented");
					//throw Math_Runtime_Exception.create_internal_error();
				}
			}
			// Update from [x0, x1] to [x0, x].
			x1 = x;
			f1 = fx;

			// If the current interval is within the given accuracies, we
			// are satisfied with the current approximation.
			if (std::abs(x1 - x0) < std::max(rtol * std::abs(x1), atol) ||
				(std::abs(f1) < ftol && (my_allowed == Allowed_Solution::ANY_SIDE ||
					(inverted && my_allowed == Allowed_Solution::LEFT_SIDE) ||
					(!inverted && my_allowed == Allowed_Solution::RIGHT_SIDE) ||
					(f1 <= 0.0 && my_allowed == Allowed_Solution::BELOW_SIDE) ||
					(f1 >= 0.0 && my_allowed == Allowed_Solution::ABOVE_SIDE))))
			{
				if (inverted)
				{
					return Interval(x1, f1, x0, f0);
				}
				return Interval(x0, f0, x1, f1);
			}
		}
	}

	/** <em>Secant</em>-based root-finding methods. */
	enum Method
	{
		/**
		 * The {@link Regula_Falsi_Solver <em>Regula Falsi</em>} or
		 * <em>False Position</em> method.
		 */
		REGULA_FALSI,
		/** The {@link Illinois_Solver <em>Illinois</em>} method. */
		ILLINOIS,
		/** The {@link Pegasus_Solver <em>Pegasus</em>} method. */
		PEGASUS
	};

public:
	/**
	 * Construct a solver.
	 *
	 * @param absolute_accuracy Absolute accuracy.
	 * @param method <em>Secant</em>-based root-finding method to use.
	 */
	Base_Secant_Solver(const double& absolute_accuracy, const Method& method) 
		:
		my_allowed{ Allowed_Solution::ANY_SIDE }, 
		my_method{ method }
	{
		Abstract_Univariate_Solver(absolute_accuracy);
	}

	/** {@inherit_doc} */
	//override
	double solve(const int& max_eval, const Univariate_Function& f, const double& min, const double& max, const Allowed_Solution& allowed_solution)
	{
		return solve(max_eval, f, min, max, min + 0.5 * (max - min), allowed_solution);
	}

	/** {@inherit_doc} */
	//override
	double solve(const int& max_eval, const Univariate_Function& f, const double& min, const double& max, const double& start_value, const Allowed_Solution& allowed_solution)
	{
		my_allowed = allowed_solution;
		return Bracketed_Univariate_Solver<Univariate_Function>::solve(max_eval, f, min, max, start_value);
	}

	/** {@inherit_doc} */
	//override
	double solve(const int& max_eval, const Univariate_Function& f, const double& min, const double& max, const double& start_value)
	{
		return solve(max_eval, f, min, max, start_value, Allowed_Solution::ANY_SIDE);
	}

	/** {@inherit_doc} */
	//override
	Interval solve_interval(const int& max_eval, const Univariate_Function& f, const double& min, const double& max, const double& start_value)
	{
		setup(max_eval, f, min, max, start_value);
		my_allowed = NULL;
		return do_solve_interval();
	}
};