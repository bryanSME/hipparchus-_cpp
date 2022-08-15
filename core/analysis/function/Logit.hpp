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

  //package org.hipparchus.analysis.function;

#include <type_traits>
#include "../differentiation/UnivariateDifferentiableFunction.h"
#include "../ParametricUnivariateFunction.h"
#include "../../util/MathUtils.h"
//import org.hipparchus.analysis.Parametric_Univariate_Function ;
//import org.hipparchus.analysis.differentiation.Derivative;
//import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Function;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/**
 * <a href="http://en.wikipedia.org/wiki/Logit">
 *  Logit</a> function.
 * It is the inverse of the {@link Sigmoid sigmoid} function.
 *
 */
class Logit : public Univariate_Differentiable_Function
{
private:
	/** Lower bound. */
	const double my_lo;
	/** Higher bound. */
	const double my_hi;

	/**
	 * @param x Value at which to compute the logit.
	 * @param lo Lower bound.
	 * @param hi Higher bound.
	 * @return the value of the logit function at {@code x}.
	 * @ if {@code x < lo} or {@code x > hi}.
	 */
	static double value(const double& x, const double& lo, const double& hi)

	{
		Math_Utils::check_range_inclusive(x, lo, hi);
		return std::log((x - lo) / (hi - x));
	}

public:
	/**
	 * Usual logit function, where the lower bound is 0 and the higher
	 * bound is 1.
	 */
	Logit()
	{
		this(0, 1);
	}

	/**
	 * Logit function.
	 *
	 * @param lo Lower bound of the function domain.
	 * @param hi Higher bound of the function domain.
	 */
	Logit(const double& lo, const double& hi) : my_lo{ lo }, my_hi{ hi } {};

	/** {@inherit_doc} */
	//override
	double value(const double& x) const
	{
		return value(x, my_lo, my_hi);
	}

	/** {@inherit_doc}
	 * @exception  if parameter is outside of function domain
	 */
	 //override
	template<typename T, typename std::enable_if<std::is_base_of<Derivative<T>, T>::value>::type* = nullptr>
	T value(T t)

	{
		const double x = t.get_value();
		Math_Utils::check_range_inclusive(x, lo, hi);
		std::vector<double> f = std::vector<double>(t.get_order() + 1];

		// function value
		f[0] = std::log((x - lo) / (hi - x));

		if (std::isinf(f[0]))
		{
			if (f.size() > 1)
			{
				f[1] = INFINITY;
			}
			// fill the array with infinities
			// (for x close to lo the signs will flip between -inf and +inf, //  for x close to hi the signs will always be +inf)
			// this is probably overkill, since the call to compose at the end
			// of the method will transform most infinities into NaN ...
			for (int i{ 2 }; i < f.size(); ++i)
			{
				f[i] = f[i - 2];
			}
		}
		else
		{
			// function derivatives
			const double inv_l = 1.0 / (x - lo);
			double x_l = inv_l;
			const double inv_h = 1.0 / (hi - x);
			double xH = inv_h;
			for (int i{ 1 }; i < f.size(); ++i)
			{
				f[i] = x_l + xH;
				x_l *= -i * inv_l;
				xH *= i * inv_h;
			}
		}

		return t.compose(f);
	}

	/**
	 * Parametric function where the input array contains the parameters of
	 * the logit function, ordered as follows:
	 * <ul>
	 *  <li>Lower bound</li>
	 *  <li>Higher bound</li>
	 * </ul>
	 */
	static class Parametric : Parametric_Univariate_Function
	{
	public:
		/**
		 * Computes the value of the logit at {@code x}.
		 *
		 * @param x Value for which the function must be computed.
		 * @param param Values of lower bound and higher bounds.
		 * @return the value of the function.
		 * @Null_Argument_Exception if {@code param} is {@code NULL}.
		 * @ if the size of {@code param} is
		 * not 2.
		 */
		 //override
		double value(const double& x, double ... param)
		{
			validate_parameters(param);
			return Logit::value(x, param[0], param[1]);
		}

		/**
		 * Computes the value of the gradient at {@code x}.
		 * The components of the gradient vector are the partial
		 * derivatives of the function with respect to each of the
		 * <em>parameters</em> (lower bound and higher bound).
		 *
		 * @param x Value at which the gradient must be computed.
		 * @param param Values for lower and higher bounds.
		 * @return the gradient vector at {@code x}.
		 * @Null_Argument_Exception if {@code param} is {@code NULL}.
		 * @ if the size of {@code param} is
		 * not 2.
		 */
		 //override
		std::vector<double> gradient(const double& x, double ... param)
		{
			validate_parameters(param);

			const double lo = param[0];
			const double hi = param[1];

			return std::vector<double> { 1 / (lo - x), 1 / (hi - x) };
		}

	private:
		/**
		 * Validates parameters to ensure they are appropriate for the evaluation of
		 * the {@link #value(double,std::vector<double>)} and {@link #gradient(double,std::vector<double>)}
		 * methods.
		 *
		 * @param param Values for lower and higher bounds.
		 * @Null_Argument_Exception if {@code param} is {@code NULL}.
		 * @ if the size of {@code param} is
		 * not 2.
		 */
		void validate_parameters(const std::vector<double>& param)
		{
			//Math_Utils::check_not_null(param);
			Math_Utils::check_dimension(param.size(), 2);
		}
	};
};