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
  //import org.hipparchus.analysis.Univariate_Function;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include <vector>

#include "../UnivariateFunction.h"
#include "BrentSolver.h"
#include "UnivariateSolver.h"
#include "AllowedSolution.h"
#include "BracketedUnivariateSolver.h"
#include "../UnivariateFunction.h"
#include "../../util/MathUtils.h"
#include "../../CalculusFieldElement.hpp"
#include "../../exception/LocalizedCoreFormats.h"
#include "../CalculusFieldUnivariateFunction.h"

/**
 * Utility routines for {@link Univariate_Solver} objects.
 *
 */
class Univariate_Solver_Utils
{
private:
	/**
	 * Class contains only static methods.
	 */
	Univariate_Solver_Utils() {}

	/** Compute the maximum of two values
	 * @param a first value
	 * @param b second value
	 * @param <T> type of the field elements
	 * @return b if a is lesser or equal to b, a otherwise
	 * @since 1.2
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T max(const T& a, const T& b)
	{
		return (a.subtract(b).get_real() <= 0) ? b : a;
	}

	/** Compute the minimum of two values
	 * @param a first value
	 * @param b second value
	 * @param <T> type of the field elements
	 * @return a if a is lesser or equal to b, b otherwise
	 * @since 1.2
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T min(const T& a, const T& b)
	{
		return (a.subtract(b).get_real() <= 0) ? a : b;
	}

public:
	/**
	 * Convenience method to find a zero of a univariate real function.  A default
	 * solver is used.
	 *
	 * @param function Function.
	 * @param x0 Lower bound for the interval.
	 * @param x1 Upper bound for the interval.
	 * @return a value where the function is zero.
	 * @ if the function has the same sign at the
	 * endpoints.
	 * @ if {@code function} is {@code NULL}.
	 */
	static double solve(const Univariate_Function& function, const double& x0, const double& x1)
	{
		//Math_Utils::check_not_null(function, hipparchus::exception::Localized_Core_Formats_Type::FUNCTION);
		const Univariate_Solver solver = Brent_Solver();
		return solver.solve(std::numeric_limits<int>::max(), function, x0, x1);
	}

		/**
		 * Convenience method to find a zero of a univariate real function.  A default
		 * solver is used.
		 *
		 * @param function Function.
		 * @param x0 Lower bound for the interval.
		 * @param x1 Upper bound for the interval.
		 * @param absolute_accuracy Accuracy to be used by the solver.
		 * @return a value where the function is zero.
		 * @ if the function has the same sign at the
		 * endpoints.
		 * @ if {@code function} is {@code NULL}.
		 */
		static double solve(const Univariate_Function& function, const double& x0, const double& x1, const double& absolute_accuracy)
		{
		//Math_Utils::check_not_null(function, hipparchus::exception::Localized_Core_Formats_Type::FUNCTION);
		const Univariate_Solver solver = Brent_Solver(absolute_accuracy);
		return solver.solve(std::numeric_limits<int>::max(), function, x0, x1);
	}

	/**
	 * Force a root found by a non-bracketing solver to lie on a specified side, * as if the solver were a bracketing one.
	 *
	 * @param max_eval maximal number of evaluations of the function
	 * (evaluations already done for finding the root should have already been subtracted
	 * from this number)
	 * @param f function to solve
	 * @param bracketing bracketing solver to use for shifting the root
	 * @param base_root original root found by a previous non-bracketing solver
	 * @param min minimal bound of the search interval
	 * @param max maximal bound of the search interval
	 * @param allowed_solution the kind of solutions that the root-finding algorithm may
	 * accept as solutions.
	 * @return a root approximation, on the specified side of the exact root
	 * @ if the function has the same sign at the
	 * endpoints.
	 */
	static double force_side(const int& max_eval, const Univariate_Function& f, const Bracketed_Univariate_Solver<Univariate_Function>& bracketing, const double& base_root, const double& min, const double& max, const Allowed_Solution allowed_solution)
	{
		if (allowed_solution == Allowed_Solution::ANY_SIDE)
		{
			// no further bracketing required
			return base_root;
		}

		// find a very small interval bracketing the root
		const double step = std::max(bracketing.get_absolute_accuracy(), std::abs(base_root * bracketing.get_relative_accuracy()));
		double x_lo = std::max(min, base_root - step);
		double f_lo = f.value(x_lo);
		double x_hi = std::min(max, base_root + step);
		double f_hi = f.value(x_hi);
		int remaining_eval = max_eval - 2;
		while (remaining_eval > 0)
		{
			if ((f_lo >= 0 && f_hi <= 0) || (f_lo <= 0 && f_hi >= 0))
			{
				// compute the root on the selected side
				return bracketing.solve(remaining_eval, f, x_lo, x_hi, base_root, allowed_solution);
			}

			// try increasing the interval
			bool change_lo = false;
			bool change_hi = false;
			if (f_lo < f_hi)
			{
				// increasing function
				if (f_lo >= 0)
				{
					change_lo = true;
				}
				else
				{
					change_hi = true;
				}
			}
			else if (f_lo > f_hi)
			{
				// decreasing function
				if (f_lo <= 0)
				{
					change_lo = true;
				}
				else
				{
					change_hi = true;
				}
			}
			else
			{
				// unknown variation
				change_lo = true;
				change_hi = true;
			}

			// update the lower bound
			if (change_lo)
			{
				x_lo = std::max(min, x_lo - step);
				f_lo = f.value(x_lo);
				remaining_eval--;
			}

			// update the higher bound
			if (change_hi)
			{
				x_hi = std::min(max, x_hi + step);
				f_hi = f.value(x_hi);
				remaining_eval--;
			}
		}

		throw std::exception("not implemented");
		//throw (Localized_Core_Formats::FAILED_BRACKETING, x_lo, x_hi, f_lo, f_hi, max_eval - remaining_eval, max_eval, base_root, min, max);
	}

	/**
	 * This method simply calls {@link #bracket(Univariate_Function, double, double, double, * double, double, int) bracket(function, initial, lower_bound, upper_bound, q, r, maximum_iterations)}
	 * with {@code q} and {@code r} set to 1.0 and {@code maximum_iterations} set to {@code std::numeric_limits<int>::max()}.
	 * <p>
	 * <strong>Note: </strong> this method can take {@code std::numeric_limits<int>::max()}
	 * iterations to throw a {@code Math_Illegal_State_Exception.}  Unless you are
	 * confident that there is a root between {@code lower_bound} and
	 * {@code upper_bound} near {@code initial}, it is better to use
	 * {@link #bracket(Univariate_Function, double, double, double, double,double, int)
	 * bracket(function, initial, lower_bound, upper_bound, q, r, maximum_iterations)}, * explicitly specifying the maximum number of iterations.</p>
	 *
	 * @param function Function.
	 * @param initial Initial midpoint of interval being expanded to
	 * bracket a root.
	 * @param lower_bound Lower bound (a is never lower than this value)
	 * @param upper_bound Upper bound (b never is greater than this
	 * value).
	 * @return a two-element array holding a and b.
	 * @ if a root cannot be bracketed.
	 * @ if {@code maximum_iterations <= 0}.
	 * @ if {@code function} is {@code NULL}.
	 */
	static std::vector<double> bracket(const Univariate_Function& function, const double& initial, const double& lower_bound, const double& upper_bound)
	{
		return bracket(function, initial, lower_bound, upper_bound, 1.0, 1.0, std::numeric_limits<int>::max());
	}

	/**
	* This method simply calls {@link #bracket(Univariate_Function, double, double, double, * double, double, int) bracket(function, initial, lower_bound, upper_bound, q, r, maximum_iterations)}
	* with {@code q} and {@code r} set to 1.0.
	* @param function Function.
	* @param initial Initial midpoint of interval being expanded to
	* bracket a root.
	* @param lower_bound Lower bound (a is never lower than this value).
	* @param upper_bound Upper bound (b never is greater than this
	* value).
	* @param maximum_iterations Maximum number of iterations to perform
	* @return a two element array holding a and b.
	* @ if the algorithm fails to find a and b
	* satisfying the desired conditions.
	* @ if {@code maximum_iterations <= 0}.
	* @ if {@code function} is {@code NULL}.
	*/
	static std::vector<double> bracket(const Univariate_Function& function, const double& initial, const double& lower_bound, const double& upper_bound, const int& maximum_iterations)
	{
		return bracket(function, initial, lower_bound, upper_bound, 1.0, 1.0, maximum_iterations);
	}

	/**
	 * This method attempts to find two values a and b satisfying <ul>
	 * <li> {@code lower_bound <= a < initial < b <= upper_bound} </li>
	 * <li> {@code f(a) * f(b) <= 0} </li>
	 * </ul>
	 * If {@code f} is continuous on {@code [a,b]}, this means that {@code a}
	 * and {@code b} bracket a root of {@code f}.
	 * <p>
	 * The algorithm checks the sign of \\( f(l_k) \\) and \\( f(u_k) \\) for increasing
	 * values of k, where \\( l_k = max(lower, initial - \\delta_k) \\), * \\( u_k = min(upper, initial + \\delta_k) \\), using recurrence
	 * \\( \\delta_{k+1} = r \\delta_k + q, \\delta_0 = 0\\) and starting search with \\( k=1 \\).
	 * The algorithm stops when one of the following happens: <ul>
	 * <li> at least one positive and one negative value have been found --  success!</li>
	 * <li> both endpoints have reached their respective limits --  </li>
	 * <li> {@code maximum_iterations} iterations elapse --  </li></ul>
	 * <p>
	 * If different signs are found at first iteration ({@code k=1}), then the returned
	 * interval will be \\( [a, b] = [l_1, u_1] \\). If different signs are found at a later
	 * iteration {@code k>1}, then the returned interval will be either
	 * \\( [a, b] = [l_{k+1}, l_{k}] \\) or \\( [a, b] = [u_{k}, u_{k+1}] \\). A root solver called
	 * with these parameters will therefore start with the smallest bracketing interval known
	 * at this step.
	 * </p>
	 * <p>
	 * Interval expansion rate is tuned by changing the recurrence parameters {@code r} and
	 * {@code q}. When the multiplicative factor {@code r} is set to 1, the sequence is a
	 * simple arithmetic sequence with linear increase. When the multiplicative factor {@code r}
	 * is larger than 1, the sequence has an asymptotically exponential rate. Note than the
	 * additive parameter {@code q} should never be set to zero, otherwise the interval would
	 * degenerate to the single initial point for all values of {@code k}.
	 * </p>
	 * <p>
	 * As a rule of thumb, when the location of the root is expected to be approximately known
	 * within some error margin, {@code r} should be set to 1 and {@code q} should be set to the
	 * order of magnitude of the error margin. When the location of the root is really a wild guess, * then {@code r} should be set to a value larger than 1 (typically 2 to double the interval
	 * length at each iteration) and {@code q} should be set according to half the initial
	 * search interval length.
	 * </p>
	 * <p>
	 * As an example, if we consider the trivial function {@code f(x) = 1 - x} and use
	 * {@code initial = 4}, {@code r = 1}, {@code q = 2}, the algorithm will compute
	 * {@code f(4-2) = f(2) = -1} and {@code f(4+2) = f(6) = -5} for {@code k = 1}, then
	 * {@code f(4-4) = f(0) = +1} and {@code f(4+4) = f(8) = -7} for {@code k = 2}. Then it will
	 * return the interval {@code [0, 2]} as the smallest one known to be bracketing the root.
	 * As shown by this example, the initial value (here {@code 4}) may lie outside of the returned
	 * bracketing interval.
	 * </p>
	 * @param function function to check
	 * @param initial Initial midpoint of interval being expanded to
	 * bracket a root.
	 * @param lower_bound Lower bound (a is never lower than this value).
	 * @param upper_bound Upper bound (b never is greater than this
	 * value).
	 * @param q additive offset used to compute bounds sequence (must be strictly positive)
	 * @param r multiplicative factor used to compute bounds sequence
	 * @param maximum_iterations Maximum number of iterations to perform
	 * @return a two element array holding the bracketing values.
	 * @exception  if function cannot be bracketed in the search interval
	 */
	static std::vector<double> bracket(const Univariate_Function& function, const double& initial, const double& lower_bound, const double& upper_bound, const double& q, const double& r, const int& maximum_iterations)
	{
		//Math_Utils::check_not_null(function, hipparchus::exception::Localized_Core_Formats_Type::FUNCTION
		if (q <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, q, 0);
		}
		if (maximum_iterations <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INVALID_MAX_ITERATIONS, maximum_iterations);
		}
		verify_sequence(lower_bound, initial, upper_bound);

		// initialize the recurrence
		double a = initial;
		double b = initial;
		double fa = std::numeric_limits<double>::quiet_NaN();
		double fb = std::numeric_limits<double>::quiet_NaN();
		double delta = 0;

		for (int num_iterations = 0; (num_iterations < maximum_iterations) && (a > lower_bound || b < upper_bound); ++num_iterations)
		{
			const double previous_a = a;
			const double previous_fa = fa;
			const double previous_b = b;
			const double previous_fb = fb;

			delta = r * delta + q;
			a = std::max(initial - delta, lower_bound);
			b = std::min(initial + delta, upper_bound);
			fa = function.value(a);
			fb = function.value(b);

			if (num_iterations == 0)
			{
				// at first iteration, we don't have a previous interval
				// we simply compare both sides of the initial interval
				if (fa * fb <= 0)
				{
					// the first interval already brackets a root
					return std::vector<double> { a, b };
				}
			}
			else
			{
				// we have a previous interval with constant sign and expand it, // we expect sign changes to occur at boundaries
				if (fa * previous_fa <= 0)
				{
					// sign change detected at near lower bound
					return std::vector<double> { a, previous_a };
				}
				else if (fb * previous_fb <= 0)
				{
					// sign change detected at near upper bound
					return std::vector<double> { previous_b, b };
				}
			}
		}

		// no bracketing found
		throw std::exception("not implemented");
		//throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_BRACKETING_INTERVAL, a, b, fa, fb);
	}

	/**
	 * This method simply calls {@link #bracket(Calculus_Field_Univariate_Function, * Calculus_Field_Element, Calculus_Field_Element, Calculus_Field_Element, Calculus_Field_Element, * Calculus_Field_Element, int) bracket(function, initial, lower_bound, upper_bound, q, r, maximum_iterations)}
	 * with {@code q} and {@code r} set to 1.0 and {@code maximum_iterations} set to {@code std::numeric_limits<int>::max()}.
	 * <p>
	 * <strong>Note: </strong> this method can take {@code std::numeric_limits<int>::max()}
	 * iterations to throw a {@code Math_Illegal_State_Exception.}  Unless you are
	 * confident that there is a root between {@code lower_bound} and
	 * {@code upper_bound} near {@code initial}, it is better to use
	 * {@link #bracket(Univariate_Function, double, double, double, double,double, int)
	 * bracket(function, initial, lower_bound, upper_bound, q, r, maximum_iterations)}, * explicitly specifying the maximum number of iterations.</p>
	 *
	 * @param function Function.
	 * @param initial Initial midpoint of interval being expanded to
	 * bracket a root.
	 * @param lower_bound Lower bound (a is never lower than this value)
	 * @param upper_bound Upper bound (b never is greater than this
	 * value).
	 * @param <T> type of the field elements
	 * @return a two-element array holding a and b.
	 * @ if a root cannot be bracketed.
	 * @ if {@code maximum_iterations <= 0}.
	 * @ if {@code function} is {@code NULL}.
	 * @since 1.2
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static std::vector<T> bracket(const Calculus_Field_Univariate_Function<T>& function, const T& initial, const T& lower_bound, const T& upper_bound)
	{
		return bracket(function, initial, lower_bound, upper_bound, initial.get_field().get_one(), initial.get_field().get_one(), std::numeric_limits<int>::max());
	}

	/**
	 * This method simply calls {@link #bracket(Calculus_Field_Univariate_Function, * Calculus_Field_Element, Calculus_Field_Element, Calculus_Field_Element, Calculus_Field_Element, * Calculus_Field_Element, int) bracket(function, initial, lower_bound, upper_bound, q, r, maximum_iterations)}
	 * with {@code q} and {@code r} set to 1.0.
	 * @param function Function.
	 * @param initial Initial midpoint of interval being expanded to
	 * bracket a root.
	 * @param lower_bound Lower bound (a is never lower than this value).
	 * @param upper_bound Upper bound (b never is greater than this
	 * value).
	 * @param maximum_iterations Maximum number of iterations to perform
	 * @param <T> type of the field elements
	 * @return a two element array holding a and b.
	 * @ if the algorithm fails to find a and b
	 * satisfying the desired conditions.
	 * @ if {@code maximum_iterations <= 0}.
	 * @ if {@code function} is {@code NULL}.
	 * @since 1.2
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static std::vector<T> bracket(const Calculus_Field_Univariate_Function<T>& function, const T& initial, const T& lower_bound, const T& upper_bound, const int& maximum_iterations)
	{
		return bracket(function, initial, lower_bound, upper_bound, initial.get_field().get_one(), initial.get_field().get_one(), maximum_iterations);
	}

	/**
	* This method attempts to find two values a and b satisfying <ul>
	* <li> {@code lower_bound <= a < initial < b <= upper_bound} </li>
	* <li> {@code f(a) * f(b) <= 0} </li>
	* </ul>
	* If {@code f} is continuous on {@code [a,b]}, this means that {@code a}
	* and {@code b} bracket a root of {@code f}.
	* <p>
	* The algorithm checks the sign of \\( f(l_k) \\) and \\( f(u_k) \\) for increasing
	* values of k, where \\( l_k = max(lower, initial - \\delta_k) \\), * \\( u_k = min(upper, initial + \\delta_k) \\), using recurrence
	* \\( \\delta_{k+1} = r \\delta_k + q, \\delta_0 = 0\\) and starting search with \\( k=1 \\).
	* The algorithm stops when one of the following happens: <ul>
	* <li> at least one positive and one negative value have been found --  success!</li>
	* <li> both endpoints have reached their respective limits --  </li>
	* <li> {@code maximum_iterations} iterations elapse --  </li></ul>
	* <p>
	* If different signs are found at first iteration ({@code k=1}), then the returned
	* interval will be \\( [a, b] = [l_1, u_1] \\). If different signs are found at a later
	* iteration {@code k>1}, then the returned interval will be either
	* \\( [a, b] = [l_{k+1}, l_{k}] \\) or \\( [a, b] = [u_{k}, u_{k+1}] \\). A root solver called
	* with these parameters will therefore start with the smallest bracketing interval known
	* at this step.
	* </p>
	* <p>
	* Interval expansion rate is tuned by changing the recurrence parameters {@code r} and
	* {@code q}. When the multiplicative factor {@code r} is set to 1, the sequence is a
	* simple arithmetic sequence with linear increase. When the multiplicative factor {@code r}
	* is larger than 1, the sequence has an asymptotically exponential rate. Note than the
	* additive parameter {@code q} should never be set to zero, otherwise the interval would
	* degenerate to the single initial point for all values of {@code k}.
	* </p>
	* <p>
	* As a rule of thumb, when the location of the root is expected to be approximately known
	* within some error margin, {@code r} should be set to 1 and {@code q} should be set to the
	* order of magnitude of the error margin. When the location of the root is really a wild guess, * then {@code r} should be set to a value larger than 1 (typically 2 to double the interval
	* length at each iteration) and {@code q} should be set according to half the initial
	* search interval length.
	* </p>
	* <p>
	* As an example, if we consider the trivial function {@code f(x) = 1 - x} and use
	* {@code initial = 4}, {@code r = 1}, {@code q = 2}, the algorithm will compute
	* {@code f(4-2) = f(2) = -1} and {@code f(4+2) = f(6) = -5} for {@code k = 1}, then
	* {@code f(4-4) = f(0) = +1} and {@code f(4+4) = f(8) = -7} for {@code k = 2}. Then it will
	* return the interval {@code [0, 2]} as the smallest one known to be bracketing the root.
	* As shown by this example, the initial value (here {@code 4}) may lie outside of the returned
	* bracketing interval.
	* </p>
	* @param function function to check
	* @param initial Initial midpoint of interval being expanded to
	* bracket a root.
	* @param lower_bound Lower bound (a is never lower than this value).
	* @param upper_bound Upper bound (b never is greater than this
	* value).
	* @param q additive offset used to compute bounds sequence (must be strictly positive)
	* @param r multiplicative factor used to compute bounds sequence
	* @param maximum_iterations Maximum number of iterations to perform
	* @param <T> type of the field elements
	* @return a two element array holding the bracketing values.
	* @exception  if function cannot be bracketed in the search interval
	* @since 1.2
	*/
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static std::vector<T> bracket(const Calculus_Field_Univariate_Function<T>& function, const T& initial, const T& lower_bound, const T& upper_bound, const T& q, const T& r, const int& maximum_iterations)
	{
		//Math_Utils::check_not_null(function, hipparchus::exception::Localized_Core_Formats_Type::FUNCTION);

		if (q.get_real() <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, q, 0);
		}
		if (maximum_iterations <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INVALID_MAX_ITERATIONS, maximum_iterations);
		}
		verify_sequence(lower_bound.get_real(), initial.get_real(), upper_bound.get_real());

		// initialize the recurrence
		T a = initial;
		T b = initial;
		T fa = NULL;
		T fb = NULL;
		T delta = initial.get_field().get_zero();

		for (const int& num_iterations = 0;
			(num_iterations < maximum_iterations) &&
			(a.get_real() > lower_bound.get_real() || b.get_real() < upper_bound.get_real());
			++num_iterations)
		{
			const T previous_a = a;
			const T previous_fa = fa;
			const T previous_b = b;
			const T previous_fb = fb;

			delta = r.multiply(delta).add(q);
			a = max(initial.subtract(delta), lower_bound);
			b = min(initial.add(delta), upper_bound);
			fa = function.value(a);
			fb = function.value(b);

			if (num_iterations == 0)
			{
				// at first iteration, we don't have a previous interval
				// we simply compare both sides of the initial interval
				if (fa.multiply(fb).get_real() <= 0)
				{
					// the first interval already brackets a root
					const std::vector<T> interval = Math_Arrays::build_array(initial.get_field(), 2);
					interval[0] = a;
					interval[1] = b;
					return interval;
				}
			}
			else
			{
				// we have a previous interval with constant sign and expand it, // we expect sign changes to occur at boundaries
				if (fa.multiply(previous_fa).get_real() <= 0)
				{
					// sign change detected at near lower bound
					const std::vector<T> interval = Math_Arrays::build_array(initial.get_field(), 2);
					interval[0] = a;
					interval[1] = previous_a;
					return interval;
				}
				if (fb.multiply(previous_fb).get_real() <= 0)
				{
					// sign change detected at near upper bound
					const std::vector<T> interval = Math_Arrays::build_array(initial.get_field(), 2);
					interval[0] = previous_b;
					interval[1] = b;
					return interval;
				}
			}
		}

		// no bracketing found
		throw std::exception("not implemented");
		//throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_BRACKETING_INTERVAL, a.get_real(), b.get_real(), fa.get_real(), fb.get_real());
	}

	/**
	 * Compute the midpoint of two values.
	 *
	 * @param a first value.
	 * @param b second value.
	 * @return the midpoint.
	 */
	static double midpoint(const double& a, const double& b)
	{
		return (a + b) * 0.5;
	}

	/**
	 * Check whether the interval bounds bracket a root. That is, if the
	 * values at the endpoints are not equal to zero, then the function takes
	 * opposite signs at the endpoints.
	 *
	 * @param function Function.
	 * @param lower Lower endpoint.
	 * @param upper Upper endpoint.
	 * @return {@code true} if the function values have opposite signs at the
	 * given points.
	 * @ if {@code function} is {@code NULL}.
	 */
	static bool is_bracketing(const Univariate_Function& function, const double& lower, const double& upper)
	{
		//Math_Utils::check_not_null(function, hipparchus::exception::Localized_Core_Formats_Type::FUNCTION);
		const double f_lo = function.value(lower);
		const double f_hi = function.value(upper);
		return (f_lo >= 0 && f_hi <= 0) || (f_lo <= 0 && f_hi >= 0);
	}

	/**
	 * Check whether the arguments form a (strictly) increasing sequence.
	 *
	 * @param start First number.
	 * @param mid Second number.
	 * @param end Third number.
	 * @return {@code true} if the arguments form an increasing sequence.
	 */
	static bool is_sequence(const double& start, const double& mid, const double& end)
	{
		return (start < mid) && (mid < end);
	}

	/**
	 * Check that the endpoints specify an interval.
	 *
	 * @param lower Lower endpoint.
	 * @param upper Upper endpoint.
	 * @ if {@code lower >= upper}.
	 */
	static void verify_interval(const double& lower, const double& upper)

	{
		if (lower >= upper)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::ENDPOINTS_NOT_AN_INTERVAL, lower, upper, false);
		}
	}

	/**
	 * Check that {@code lower < initial < upper}.
	 *
	 * @param lower Lower endpoint.
	 * @param initial Initial value.
	 * @param upper Upper endpoint.
	 * @ if {@code lower >= initial} or
	 * {@code initial >= upper}.
	 */
	static void verify_sequence(const double& lower, const double& initial, const double& upper)
	{
		verify_interval(lower, initial);
		verify_interval(initial, upper);
	}

	/**
	 * Check that the endpoints specify an interval and the end points
	 * bracket a root.
	 *
	 * @param function Function.
	 * @param lower Lower endpoint.
	 * @param upper Upper endpoint.
	 * @ if the function has the same sign at the
	 * endpoints.
	 * @ if {@code function} is {@code NULL}.
	 */
	static void verify_bracketing(const Univariate_Function& function, const double& lower, const double& upper)
	{
		//Math_Utils::check_not_null(function, hipparchus::exception::Localized_Core_Formats_Type::FUNCTION);
		verify_interval(lower, upper);
		if (!is_bracketing(function, lower, upper))
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_BRACKETING_INTERVAL, lower, upper, function.value(lower), function.value(upper));
		}
	}
};