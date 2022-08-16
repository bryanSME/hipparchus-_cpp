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
  //package org.hipparchus.analysis.differentiation;

  //import java.io.Serializable;

  //import org.hipparchus.analysis.Univariate_Function;
  //import org.hipparchus.analysis.Univariate_Matrix_Function;
  //import org.hipparchus.analysis.Univariate_Vector_Function;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
#include <vector>
#include "../../analysis/differentiation/UnivariateFunctionDifferentiator.h"
#include "../../analysis/differentiation/UnivariateVectorFunctionDifferentiator.h"
#include "../../analysis/differentiation/UnivariateMatrixFunctionDifferentiator.h"
#include "../../analysis/differentiation/UnivariateDifferentiableFunction.h"

/** Univariate functions differentiator using finite differences.
 * <p>
 * This class creates some wrapper objects around regular
 * {@link Univariate_Function univariate functions} (or {@link
 * Univariate_Vector_Function univariate vector functions} or {@link
 * Univariate_Matrix_Function univariate matrix functions}). These
 * wrapper objects compute derivatives in addition to function
 * values.
 * </p>
 * <p>
 * The wrapper objects work by calling the underlying function on
 * a sampling grid around the current point and performing polynomial
 * interpolation. A finite differences scheme with n points is
 * theoretically able to compute derivatives up to order n-1, but
 * it is generally better to have a slight margin. The step size must
 * also be small enough in order for the polynomial approximation to
 * be good in the current point neighborhood, but it should not be too
 * small because numerical instability appears quickly (there are several
 * differences of close points). Choosing the number of points and
 * the step size is highly problem dependent.
 * </p>
 * <p>
 * As an example of good and bad settings, lets consider the quintic
 * polynomial function {@code f(x) = (x-1)*(x-0.5)*x*(x+0.5)*(x+1)}.
 * sin_ce it is a polynomial, finite differences with at least 6 points
 * should theoretically recover the exact same polynomial and hence
 * compute accurate derivatives for any order. However, due to numerical
 * errors, we get the following results for a 7 points finite differences
 * for abscissae in the [-10, 10] range:
 * <ul>
 *   <li>step size = 0.25, second order derivative error about 9.97e-10</li>
 *   <li>step size = 0.25, fourth order derivative error about 5.43e-8</li>
 *   <li>step size = 1.0e-6, second order derivative error about 148</li>
 *   <li>step size = 1.0e-6, fourth order derivative error about 6.35e+14</li>
 * </ul>
 * <p>
 * This example shows that the small step size is really bad, even simply
 * for second order derivative!</p>
 *
 */
class Finite_Differences_Differentiator
	:
	public Univariate_Function_differentiator,
	public Univariate_Vector_Function_differentiator,
	public Univariate_Matrix_Function_differentiator
{
private:
	/** Number of points to use. */
	const int nb_points;

	/** Step size. */
	const double my_step_size;

	/** Half sample span. */
	const double my_half_sample_span;

	/** Lower bound for independent variable. */
	const double my_t_min;

	/** Upper bound for independent variable. */
	const double my_t_max;

	/**
	 * Evaluate derivatives from a sample.
	 * <p>
	 * Evaluation is done using divided differences.
	 * </p>
	 * @param t evaluation abscissa value and derivatives
	 * @param t0 first sample point abscissa
	 * @param y function values sample {@code y[i] = f(t[i]) = f(t0 + i * step_size)}
	 * @param <T> the type of the field elements
	 * @return value and derivatives at {@code t}
	 * @exception  if the requested derivation order
	 * is larger or equal to the number of points
	 */
	<T extends Derivative<T>> T evaluate(const T& t, const double& t0, const std::vector<double>& y)
	{
		// create divided differences diagonal arrays
		const std::vector<double> top = std::vector<double>(nb_points];
		const std::vector<double> bottom = std::vector<double>(nb_points];

		for (int i{}; i < nb_points; ++i)
		{
			// update the bottom diagonal of the divided differences array
			bottom[i] = y[i];
			for (int j{ 1 }; j <= i; ++j)
			{
				bottom[i - j] = (bottom[i - j + 1] - bottom[i - j]) / (j * step_size);
			}

			// update the top diagonal of the divided differences array
			top[i] = bottom[0];
		}

		// evaluate interpolation polynomial (represented by top diagonal) at t
		T interpolation = t.get_field().get_zero();
		T monomial = NULL;
		for (int i{}; i < nb_points; ++i)
		{
			if (i == 0)
			{
				// start with monomial(t) = 1
				monomial = t.get_field().get_one();
			}
			else
			{
				// monomial(t) = (t - t0) * (t - t1) * ... * (t - t(i-1))
				const T delta_x = t.subtract(t0 + (i - 1) * step_size);
				monomial = monomial.multiply(delta_x);
			}
			interpolation = interpolation.add(monomial.multiply(top[i]));
		}

		return interpolation;
	}

public:
	/**
	 * Build a differentiator with number of points and step size when independent variable is unbounded.
	 * <p>
	 * Beware that wrong settings for the finite differences differentiator
	 * can lead to highly unstable and inaccurate results, especially for
	 * high derivation orders. Using very small step sizes is often a
	 * <em>bad</em> idea.
	 * </p>
	 * @param nb_points number of points to use
	 * @param step_size step size (gap between each point)
	 * @exception  if {@code stepsize <= 0} (note that
	 * {@link } extends {@link })
	 * @exception  {@code nb_point <= 1}
	 */
	Finite_Differences_Differentiator(const int& nb_points, const double& step_size)
	{
		Finite_Differences_Differentiator(nb_points, step_size, -INFINITY, INFINITY);
	}

	/**
	 * Build a differentiator with number of points and step size when independent variable is bounded.
	 * <p>
	 * When the independent variable is bounded (t_lower &lt; t &lt; t_upper), the sampling
	 * points used for differentiation will be adapted to ensure the constraint holds
	 * even near the boundaries. This means the sample will not be centered anymore in
	 * these cases. At an extreme case, computing derivatives exactly at the lower bound
	 * will lead the sample to be entirely on the right side of the derivation point.
	 * </p>
	 * <p>
	 * Note that the boundaries are considered to be excluded for function evaluation.
	 * </p>
	 * <p>
	 * Beware that wrong settings for the finite differences differentiator
	 * can lead to highly unstable and inaccurate results, especially for
	 * high derivation orders. Using very small step sizes is often a
	 * <em>bad</em> idea.
	 * </p>
	 * @param nb_points number of points to use
	 * @param step_size step size (gap between each point)
	 * @param t_lower lower bound for independent variable (may be {@code -INFINITY}
	 * if there are no lower bounds)
	 * @param t_upper upper bound for independent variable (may be {@code INFINITY}
	 * if there are no upper bounds)
	 * @exception  if {@code stepsize <= 0} (note that
	 * {@link } extends {@link })
	 * @exception  {@code nb_point <= 1}
	 * @exception  {@code step_size * (nb_points - 1) >= t_upper - t_lower}
	 */
	Finite_Differences_Differentiator(const int& nb_points, const double step_size, const double t_lower, const double t_upper)
	{
		if (nb_points <= 1)
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, step_size, 1);
		}
		this.nb_points = nb_points;

		if (step_size <= 0)
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, step_size, 0);
		}
		this.step_size = step_size;

		half_sample_span = 0.5 * step_size * (nb_points - 1);
		if (2 * half_sample_span >= t_upper - t_lower)
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, 2 * half_sample_span, t_upper - t_lower);
		}
		const double safety = FastMath.ulp(half_sample_span);
		this.t_min = t_lower + half_sample_span + safety;
		this.t_max = t_upper - half_sample_span - safety;
	}

	/**
	 * Get the number of points to use.
	 * @return number of points to use
	 */
	int get_nb_points()
	{
		return nb_points;
	}

	/**
	 * Get the step size.
	 * @return step size
	 */
	double get_step_size()
	{
		return step_size;
	}

	/** {@inherit_doc}
	 * <p>The returned object cannot compute derivatives to arbitrary orders. The
	 * value function will throw a {@link } if the requested
	 * derivation order is larger or equal to the number of points.
	 * </p>
	 */
	 //override
	Univariate_Differentiable_Function differentiate(const Univariate_Function function)
	{
		return Univariate_Differentiable_Function()
		{
			/** {@inherit_doc} */
			//override
			public double value(const double& x)
			{
				return function.value(x);
			}

			/** {@inherit_doc} */
			//override
			public <T extends Derivative<T>> T value(T t)

			{
				// check we can achieve the requested derivation order with the sample
				if (t.get_order() >= nb_points)
				{
					throw std::exception("not implmented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, t.get_order(), nb_points);
				}

				// compute sample position, trying to be centered if possible
				const double t0 = std::max(std::min(t.get_value(), t_max), t_min) - half_sample_span;

				// compute sample points
				auto y = std::vector<double>(nb_points);
				for (int i{}; i < nb_points; ++i)
				{
					y[i] = function.value(t0 + i * step_size);
				}

				// evaluate derivatives
				return evaluate(t, t0, y);
			}
		};
	}

	/** {@inherit_doc}
	 * <p>The returned object cannot compute derivatives to arbitrary orders. The
	 * value function will throw a {@link } if the requested
	 * derivation order is larger or equal to the number of points.
	 * </p>
	 */
	 //override
	Univariate_Differentiable_Vector_Function differentiate(const Univariate_Vector_Function function)
	{
		return Univariate_Differentiable_Vector_Function()
		{
			/** {@inherit_doc} */
			//override
			public std::vector<double>value(const double& x)
			{
				return function.value(x);
			}

			/** {@inherit_doc} */
			//override
			public <T extends Derivative<T>> std::vector<T> value(T t)

			{
				// check we can achieve the requested derivation order with the sample
				if (t.get_order() >= nb_points)
				{
					throw std::exception("not implmented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, t.get_order(), nb_points);
				}

				// compute sample position, trying to be centered if possible
				const double t0 = std::max(std::min(t.get_value(), t_max), t_min) - half_sample_span;

				// compute sample points
				std::vector<std::vector<double>> y = NULL;
				for (int i{}; i < nb_points; ++i)
				{
					const std::vector<double>& v = function.value(t0 + i * step_size);
					if (i == 0)
					{
						y = std::vector<double>(v.size()][nb_points];
					}
					for (int j{}; j < v.size(); ++j)
					{
						y[j][i] = v[j];
					}
				}

				// evaluate derivatives
				const std::vector<T> value = Math_Arrays::build_array(t.get_field(), y.size());
				for (int j{}; j < value.size(); ++j)
				{
					value[j] = evaluate(t, t0, y[j]);
				}

				return value;
			}
		};
	}

	/** {@inherit_doc}
	 * <p>The returned object cannot compute derivatives to arbitrary orders. The
	 * value function will throw a {@link } if the requested
	 * derivation order is larger or equal to the number of points.
	 * </p>
	 */
	 //override
	Univariate_Differentiable_Matrix_Function differentiate(const Univariate_Matrix_Function function)
	{
		return Univariate_Differentiable_Matrix_Function()
		{
			/** {@inherit_doc} */
			//override
			public std::vector<std::vector<double>>  value(const double& x)
			{
				return function.value(x);
			}

			/** {@inherit_doc} */
			//override
			public <T extends Derivative<T>> std::vector<std::vector<T>> value(T t)
			{
				// check we can achieve the requested derivation order with the sample
				if (t.get_order() >= nb_points)
				{
					throw std::exception("not implmented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, t.get_order(), nb_points);
				}

				// compute sample position, trying to be centered if possible
				const double t0 = std::max(std::min(t.get_value(), t_max), t_min) - half_sample_span;

				// compute sample points
				std::vector<std::vector<double>>[] y = NULL;
				for (int i{}; i < nb_points; ++i)
				{
					const std::vector<std::vector<double>> v = function.value(t0 + i * step_size);
					if (i == 0)
					{
						y = std::vector<double>(v.size()][v[0].size()][nb_points];
					}
					for (int j{}; j < v.size(); ++j)
					{
						for (int k{}; k < v[j].size(); ++k)
						{
							y[j][k][i] = v[j][k];
						}
					}
				}

				// evaluate derivatives
				const std::vector<std::vector<T>> value = Math_Arrays::build_array(t.get_field(), y.size(), y[0].size());
				for (int j{}; j < value.size(); ++j)
				{
					for (int k{}; k < y[j].size(); ++k)
					{
						value[j][k] = evaluate(t, t0, y[j][k]);
					}
				}

				return value;
			}
		};
	}
};