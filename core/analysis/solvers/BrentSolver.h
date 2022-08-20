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

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Precision;
#include "AbstractUnivariateSolver.h"
#include "../../util/Precision.h"

/**
 * This class : the <a href="http://mathworld.wolfram.com/Brents_method.html">
 * Brent algorithm</a> for finding zeros of real univariate functions.
 * The function should be continuous but not necessarily smooth.
 * The {@code solve} method returns a zero {@code x} of the function {@code f}
 * in the given interval {@code [a, b]} to within a tolerance
 * {@code 2 eps abs(x) + t} where {@code eps} is the relative accuracy and
 * {@code t} is the absolute accuracy.
 * <p>The given interval must bracket the root.</p>
 * <p>
 *  The reference implementation is given in chapter 4 of
 *  <blockquote>
 *   <b>Algorithms for Minimization Without Derivatives</b>, *   <em>Richard P. Brent</em>, *   Dover, 2002
 *  </blockquote>
 *
 * @see Base_Abstract_Univariate_Solver
 */
class Brent_Solver : public Abstract_Univariate_Solver
{
private:
	/** Default absolute accuracy. */
	static constexpr double DEFAULT_ABSOLUTE_ACCURACY{ 1e-6 };

	/**
	 * Search for a zero inside the provided interval.
	 * This implementation is based on the algorithm described at page 58 of
	 * the book
	 * <blockquote>
	 *  <b>Algorithms for Minimization Without Derivatives</b>, *  <it>Richard P. Brent</it>, *  Dover 0-486-41998-3
	 * </blockquote>
	 *
	 * @param lo Lower bound of the search interval.
	 * @param hi Higher bound of the search interval.
	 * @param f_lo Function value at the lower bound of the search interval.
	 * @param f_hi Function value at the higher bound of the search interval.
	 * @return the value where the function is zero.
	 */
	double brent(const double& lo, const double& hi, const double& f_lo, const double& f_hi)
	{
		double a = lo;
		double fa = f_lo;
		double b = hi;
		double fb = f_hi;
		double c = a;
		double fc = fa;
		double d = b - a;
		double e = d;

		const double t = get_absolute_accuracy();
		const double eps = get_relative_accuracy();

		while (true)
		{
			if (std::abs(fc) < std::abs(fb))
			{
				a = b;
				b = c;
				c = a;
				fa = fb;
				fb = fc;
				fc = fa;
			}

			const double tol = 2 * eps * std::abs(b) + t;
			const double m = 0.5 * (c - b);

			if (std::abs(m) <= tol || Precision::equals(fb, 0))
			{
				return b;
			}
			if (std::abs(e) < tol || std::abs(fa) <= std::abs(fb))
			{
				// Force bisection.
				d = m;
				e = d;
			}
			else
			{
				double s = fb / fa;
				double p;
				double q;
				// The equality test (a == c) is intentional, // it is part of the original Brent's method and
				// it should NOT be replaced by proximity test.
				if (a == c)
				{
					// Linear interpolation.
					p = 2 * m * s;
					q = 1 - s;
				}
				else
				{
					// Inverse quadratic interpolation.
					q = fa / fc;
					const double r = fb / fc;
					p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
					q = (q - 1) * (r - 1) * (s - 1);
				}
				if (p > 0)
				{
					q = -q;
				}
				else
				{
					p = -p;
				}
				s = e;
				e = d;
				if (p >= 1.5 * m * q - std::abs(tol * q) ||
					p >= std::abs(0.5 * s * q))
				{
					// Inverse quadratic interpolation gives a value
					// in the wrong direction, or progress is slow.
					// Fall back to bisection.
					d = m;
					e = d;
				}
				else
				{
					d = p / q;
				}
			}
			a = b;
			fa = fb;

			if (std::abs(d) > tol)
			{
				b += d;
			}
			else if (m > 0)
			{
				b += tol;
			}
			else
			{
				b -= tol;
			}
			fb = compute_objective_value(b);
			if ((fb > 0 && fc > 0) ||
				(fb <= 0 && fc <= 0))
			{
				c = a;
				fc = fa;
				d = b - a;
				e = d;
			}
		}
	}

public:
	/**
	 * Construct a solver with default absolute accuracy (1e-6).
	 */
	Brent_Solver()
	{
		Brent_Solver(DEFAULT_ABSOLUTE_ACCURACY);
	}
	/**
	 * Construct a solver.
	 *
	 * @param absolute_accuracy Absolute accuracy.
	 */
	Brent_Solver(double absolute_accuracy)
	{
		Abstract_Univariate_Solver(absolute_accuracy);
	}
	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 */
	Brent_Solver(double relative_accuracy, double absolute_accuracy)
	{
		Abstract_Univariate_Solver(relative_accuracy, absolute_accuracy);
	}
	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 * @param function_value_accuracy Function value accuracy.
	 *
	 * @see Base_Abstract_Univariate_Solver#Base_Abstract_Univariate_Solver(double,double,double)
	 */
	Brent_Solver(double relative_accuracy, double absolute_accuracy, double function_value_accuracy)
	{
		Abstract_Univariate_Solver(relative_accuracy, absolute_accuracy, function_value_accuracy);
	}

protected:
	/**
	 * {@inherit_doc}
	 */
	 //override
	double do_solve()
	{
		double min = get_min();
		double max = get_max();
		const double initial = get_start_value();
		const double function_value_accuracy = get_function_value_accuracy();

		verify_sequence(min, initial, max);

		// Return the initial guess if it is good enough.
		double y_initial = compute_objective_value(initial);
		if (std::abs(y_initial) <= function_value_accuracy)
		{
			return initial;
		}

		// Return the first endpoint if it is good enough.
		double y_min = compute_objective_value(min);
		if (std::abs(y_min) <= function_value_accuracy)
		{
			return min;
		}

		// Reduce interval if min and initial bracket the root.
		if (y_initial * y_min < 0)
		{
			return brent(min, initial, y_min, y_initial);
		}

		// Return the second endpoint if it is good enough.
		double y_max = compute_objective_value(max);
		if (std::abs(y_max) <= function_value_accuracy)
		{
			return max;
		}

		// Reduce interval if initial and max bracket the root.
		if (y_initial * y_max < 0)
		{
			return brent(initial, max, y_initial, y_max);
		}

		throw std::exception("not implemented");
		//throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_BRACKETING_INTERVAL, min, max, y_min, y_max);
	}
};