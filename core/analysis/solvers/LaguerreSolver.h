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

  //import org.hipparchus.analysis.polynomials.Polynomial_Function;
  //import org.hipparchus.complex.std::complex<double>;
  //import org.hipparchus.complex.Complex_Utils;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
#include <limits>
#include <complex>
#include <vector>
#include "AbstractPolynomialSolver.h"
#include "../../exception/LocalizedCoreFormats.h"
#include "../../complex/ComplexUtils.hpp"

/**
 * Implements the <a href="http://mathworld.wolfram.com/Laguerres_method.html">
 * Laguerre's Method</a> for root finding of real coefficient polynomials.
 * For reference, see
 * <blockquote>
 *  <b>A First Course in Numerical Analysis</b>, *  ISBN 048641454X, chapter 8.
 * </blockquote>
 * Laguerre's method is global in the sense that it can start with any initial
 * approximation and be able to solve all roots from that point.
 * The algorithm requires a bracketing condition.
 *
 */
class Laguerre_Solver : public Abstract_Polynomial_Solver
{
private:
	/** Default absolute accuracy. */
	static constexpr double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;
	/** std::complex<double> solver. */
	const auto complex_solver = Complex_Solver();

	/**
	 * Find a real root in the given interval.
	 *
	 * Despite the bracketing condition, the root returned by
	 * {@link Laguerre_Solver.Complex_Solver#solve(std::complex<double>[],std::complex<double>)} may
	 * not be a real zero inside {@code [min, max]}.
	 * For example, <code> p(x) = x<sup>3</sup> + 1, </code>
	 * with {@code min = -2}, {@code max = 2}, {@code initial = 0}.
	 * When it occurs, this code calls
	 * {@link Laguerre_Solver.Complex_Solver#solve_all(std::complex<double>[],std::complex<double>)}
	 * in order to obtain all roots and picks up one real root.
	 *
	 * @param lo Lower bound of the search interval.
	 * @param hi Higher bound of the search interval.
	 * @return the point at which the function value is zero.
	 */
	double laguerre(const double& lo, const double& hi)
	{
		auto c = Complex_Utils<double>::convert_to_complex(get_coefficients());

		const auto initial = std::complex<double>(0.5 * (lo + hi), 0);
		const auto z = complex_solver.solve(c, initial);
		if (complex_solver.is_root(lo, hi, z))
		{
			return z.get_real();
		}
		double r = std::numeric_limits<double>::quiet_NaN();
		// Solve all roots and select the one we are seeking.
		auto root = complex_solver.solve_all(c, initial);
		for (int i{}; i < root.size(); i++)
		{
			if (complex_solver.is_root(lo, hi, root[i]))
			{
				r = root[i].get_real();
				break;
			}
		}
		return r;
	}

public:
	/**
	 * Construct a solver with default accuracy (1e-6).
	 */
	Laguerre_Solver()
	{
		Laguerre_Solver(DEFAULT_ABSOLUTE_ACCURACY);
	}
	/**
	 * Construct a solver.
	 *
	 * @param absolute_accuracy Absolute accuracy.
	 */
	Laguerre_Solver(double absolute_accuracy)
	{
		super(absolute_accuracy);
	}
	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 */
	Laguerre_Solver(double relative_accuracy, double absolute_accuracy)
	{
		super(relative_accuracy, absolute_accuracy);
	}
	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 * @param function_value_accuracy Function value accuracy.
	 */
	Laguerre_Solver(double relative_accuracy, double absolute_accuracy, double function_value_accuracy)
	{
		super(relative_accuracy, absolute_accuracy, function_value_accuracy);
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	double do_solve()
	{
		const double min = get_min();
		const double max = get_max();
		const double initial = get_start_value();
		const double function_value_accuracy = get_function_value_accuracy();

		verify_sequence(min, initial, max);

		// Return the initial guess if it is good enough.
		const double y_initial = compute_objective_value(initial);
		if (std::abs(y_initial) <= function_value_accuracy)
		{
			return initial;
		}

		// Return the first endpoint if it is good enough.
		const double y_min = compute_objective_value(min);
		if (std::abs(y_min) <= function_value_accuracy)
		{
			return min;
		}

		// Reduce interval if min and initial bracket the root.
		if (y_initial * y_min < 0)
		{
			return laguerre(min, initial);
		}

		// Return the second endpoint if it is good enough.
		const double y_max = compute_objective_value(max);
		if (std::abs(y_max) <= function_value_accuracy)
		{
			return max;
		}

		// Reduce interval if initial and max bracket the root.
		if (y_initial * y_max < 0)
		{
			return laguerre(initial, max);
		}

		throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_BRACKETING_INTERVAL, min, max, y_min, y_max);
	}

	/**
	 * Find all complex roots for the polynomial with the given
	 * coefficients, starting from the given initial value.
	 * <p>
	 * Note: This method is not part of the API of {@link Base_Univariate_Solver}.</p>
	 *
	 * @param coefficients Polynomial coefficients.
	 * @param initial Start value.
	 * @return the point at which the function value is zero.
	 * @org.hipparchus.exception.Math_Illegal_State_Exception
	 * if the maximum number of evaluations is exceeded.
	 * @ if the {@code coefficients} is
	 * {@code NULL}.
	 * @ if the {@code coefficients} array is empty.
	 */
	std::vector<std::complex<double>>solve_all_complex(const std::vector<double>& coefficients, const double& initial)
	{
		setup(std::numeric_limits<int>::max(), Polynomial_Function(coefficients), -INFINITY, INFINITY, initial);
		return complex_solver.solve_all(Complex_Utils<double>::convert_to_complex(coefficients), std::complex<double>(initial, 0.0));
	}

	/**
	 * Find a complex root for the polynomial with the given coefficients, * starting from the given initial value.
	 * <p>
	 * Note: This method is not part of the API of {@link Base_Univariate_Solver}.</p>
	 *
	 * @param coefficients Polynomial coefficients.
	 * @param initial Start value.
	 * @return the point at which the function value is zero.
	 * @org.hipparchus.exception.Math_Illegal_State_Exception
	 * if the maximum number of evaluations is exceeded.
	 * @ if the {@code coefficients} is
	 * {@code NULL}.
	 * @ if the {@code coefficients} array is empty.
	 */
	std::complex<double> solve_complex(const std::vector<double>& coefficients, const double& initial)
	{
		setup(std::numeric_limits<int>::max(), Polynomial_Function(coefficients), -INFINITY, INFINITY, initial);
		return complex_solver.solve(Complex_Utils<double>::convert_to_complex(coefficients), std::complex<double>(initial, 0.0));
	}

	/**
	 * Class for searching all (complex) roots.
	 */
	class Complex_Solver
	{
	public:
		/**
		 * Check whether the given complex root is actually a real zero
		 * in the given interval, within the solver tolerance level.
		 *
		 * @param min Lower bound for the interval.
		 * @param max Upper bound for the interval.
		 * @param z std::complex<double> root.
		 * @return {@code true} if z is a real zero.
		 */
		bool is_root(double min, const double& max, std::complex<double> z)
		{
			if (is_sequence(min, z.get_real(), max))
			{
				const double z_abs = z.norm();
				double tolerance = std::max(get_relative_accuracy() * z_abs, get_absolute_accuracy());
				return (std::abs(z.get_imaginary()) <= tolerance) ||
					(z_abs <= get_function_value_accuracy());
			}
			return false;
		}

		/**
		 * Find all complex roots for the polynomial with the given
		 * coefficients, starting from the given initial value.
		 *
		 * @param coefficients Polynomial coefficients.
		 * @param initial Start value.
		 * @return the point at which the function value is zero.
		 * @org.hipparchus.exception.Math_Illegal_State_Exception
		 * if the maximum number of evaluations is exceeded.
		 * @ if the {@code coefficients} is
		 * {@code NULL}.
		 * @ if the {@code coefficients} array is empty.
		 */
		std::vector<std::complex<double>>solve_all(std::complex<double> coefficients[], std::complex<double> initial)
		{
			if (coefficients == NULL)
			{
				throw std::exception("not implemented");
				//throw ();
			}
			const int n = coefficients.size() - 1;
			if (n == 0)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::POLYNOMIAL);
			}
			// Coefficients for deflated polynomial.
			const std::complex<double> c[] = std::complex<double>[n + 1];
			for (int i{}; i <= n; i++)
			{
				c[i] = coefficients[i];
			}

			// Solve individual roots successively.
			const std::complex<double> root[] = std::complex<double>[n];
			for (int i{}; i < n; i++)
			{
				const std::complex<double> subarray[] = std::complex<double>[n - i + 1];
				System.arraycopy(c, 0, subarray, 0, subarray.size());
				root[i] = solve(subarray, initial);
				// Polynomial deflation using synthetic division.
				std::complex<double> newc = c[n - i];
				std::complex<double> oldc;
				for (int j = n - i - 1; j >= 0; j--)
				{
					oldc = c[j];
					c[j] = newc;
					newc = oldc.add(newc.multiply(root[i]));
				}
			}

			return root;
		};

		/**
		 * Find a complex root for the polynomial with the given coefficients, * starting from the given initial value.
		 *
		 * @param coefficients Polynomial coefficients.
		 * @param initial Start value.
		 * @return the point at which the function value is zero.
		 * @org.hipparchus.exception.Math_Illegal_State_Exception
		 * if the maximum number of evaluations is exceeded.
		 * @ if the {@code coefficients} is
		 * {@code NULL}.
		 * @ if the {@code coefficients} array is empty.
		 */
		std::complex<double> solve(std::complex<double> coefficients[], std::complex<double> initial)
		{
			if (coefficients == NULL)
			{
				throw ();
			}

			const int n = coefficients.size() - 1;
			if (n == 0)
			{
				throw (hipparchus::exception::Localized_Core_Formats_Type::POLYNOMIAL);
			}

			const double& absolute_accuracy = get_absolute_accuracy();
			const double relative_accuracy = get_relative_accuracy();
			const double function_value_accuracy = get_function_value_accuracy();

			const std::complex<double> n_c = std::complex<double>(n, 0);
			const std::complex<double> n1C = std::complex<double>(n - 1, 0);

			std::complex<double> z = initial;
			std::complex<double> oldz = std::complex<double>(INFINITY, INFINITY);
			while (true)
			{
				// Compute pv (polynomial value), dv (derivative value), and
				// d2v (second derivative value) simultaneously.
				std::complex<double> pv = coefficients[n];
				std::complex<double> dv = std::complex<double>.ZERO;
				std::complex<double> d2v = std::complex<double>.ZERO;
				for (int j = n - 1; j >= 0; j--)
				{
					d2v = dv.add(z.multiply(d2v));
					dv = pv.add(z.multiply(dv));
					pv = coefficients[j].add(z.multiply(pv));
				}
				d2v = d2v.multiply(new std::complex<double>(2.0, 0.0));

				// Check for convergence.
				const double& tolerance = std::max(relative_accuracy * z.norm(), absolute_accuracy);
				if ((z.subtract(oldz)).norm() <= tolerance)
				{
					return z;
				}
				if (pv.norm() <= function_value_accuracy)
				{
					return z;
				}

				// Now pv != 0, calculate the approximation.
				const auto G = dv.divide(pv);
				const auto G2 = G.multiply(G);
				const auto H = G2.subtract(d2v.divide(pv));
				const auto delta = n1C.multiply((n_c.multiply(H)).subtract(G2));
				// Choose a denominator larger in magnitude.
				const auto delta_sqrt = delta.sqrt();
				const auto dplus = G.add(delta_sqrt);
				const auto dminus = G.subtract(delta_sqrt);
				const auto denominator = dplus.norm() > dminus.norm() ? dplus : dminus;
				// Perturb z if denominator is zero, for instance, // p(x) = x^3 + 1, z = 0.
				if (denominator.is_zero())
				{
					z = z.add(new std::complex<double>(absolute_accuracy, absolute_accuracy));
					oldz = std::complex<double>(INFINITY, INFINITY);
				}
				else
				{
					oldz = z;
					z = z.subtract(n_c.divide(denominator));
				}
				increment_evaluation_count();
			}
		}
	};
};