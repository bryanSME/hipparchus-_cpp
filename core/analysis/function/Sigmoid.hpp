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

  //import java.util.Arrays;

  //import org.hipparchus.analysis.Parametric_Univariate_Function ;
  //import org.hipparchus.analysis.differentiation.Derivative;
  //import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Function;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
#include <vector>
#include "../../util/MathUtils.h"
#include "../ParametricUnivariateFunction.h"
#include "../differentiation/UnivariateDifferentiableFunction.h"

/**
 * <a href="http://en.wikipedia.org/wiki/Sigmoid_function">
 *  Sigmoid</a> function.
 * It is the inverse of the {@link Logit logit} function.
 * A more flexible version, the generalised logistic, is implemented
 * by the {@link Logistic} class.
 *
 */
class Sigmoid : public Univariate_Differentiable_Function
{
private:
	/** Lower asymptote. */
	const double my_lo;
	/** Higher asymptote. */
	const double my_hi;

	/**
	 * @param x Value at which to compute the sigmoid.
	 * @param lo Lower asymptote.
	 * @param hi Higher asymptote.
	 * @return the value of the sigmoid function at {@code x}.
	 */
	static double value(const double& x, const double& lo, const double& hi)
	{
		return lo + (hi - lo) / (1 + std::exp(-x));
	}

public:
	/**
	 * Usual sigmoid function, where the lower asymptote is 0 and the higher
	 * asymptote is 1.
	 */
	Sigmoid()
	{
		Sigmoid(0, 1);
	}

	/**
	 * Sigmoid function.
	 *
	 * @param lo Lower asymptote.
	 * @param hi Higher asymptote.
	 */
	Sigmoid(const double& lo, const double& hi) : my_lo{ lo }, my_hi{ hi } {};

	/** {@inherit_doc} */
	//override
	double value(const double& x)
	{
		return value(x, my_lo, my_hi);
	}

	/** {@inherit_doc}
	 */
	 //override
	template<typename T, typename std::enable_if<std::is_base_of<Derivative<T>, T>::value>::type* = nullptr>
	T value(const T& t)
	{
		auto f = std::vector<double>(t.get_order() + 1];
		const double exp = std::exp(-t.get_value());
		if (std::isinf(exp))
		{
			// special handling near lower boundary, to avoid NaN
			f[0] = lo;
			Arrays.fill(f, 1, f.size(), 0.0);
		}
		else
		{
			// the nth order derivative of sigmoid has the form:
			// dn(sigmoid(x)/dxn = P_n(exp(-x)) / (1+exp(-x))^(n+1)
			// where P_n(t) is a degree n polynomial with normalized higher term
			// P_0(t) = 1, P_1(t) = t, P_2(t) = t^2 - t, P_3(t) = t^3 - 4 t^2 + t...
			// the general recurrence relation for P_n is:
			// P_n(x) = n t P_(n-1)(t) - t (1 + t) P_(n-1)'(t)
			auto p = std::vector<double>(f.size());

			const double inv = 1 / (1 + exp);
			double coeff = my_hi - my_lo;
			for (int n{}; n < f.size(); ++n)
			{
				// update and evaluate polynomial P_n(t)
				double v{};
				p[n] = 1;
				for (int k = n; k >= 0; --k)
				{
					v = v * exp + p[k];
					if (k > 1)
					{
						p[k - 1] = (n - k + 2) * p[k - 2] - (k - 1) * p[k - 1];
					}
					else
					{
						p[0] = 0;
					}
				}

				coeff *= inv;
				f[n] = coeff * v;
			}

			// fix function value
			f[0] += lo;
		}

		return t.compose(f);
	}

	/**
	 * Parametric function where the input array contains the parameters of
	 * the {@link Sigmoid#Sigmoid(double,double) sigmoid function}, ordered
	 * as follows:
	 * <ul>
	 *  <li>Lower asymptote</li>
	 *  <li>Higher asymptote</li>
	 * </ul>
	 */
	class Parametric : Parametric_Univariate_Function
	{
	public:
		/**
		 * Computes the value of the sigmoid at {@code x}.
		 *
		 * @param x Value for which the function must be computed.
		 * @param param Values of lower asymptote and higher asymptote.
		 * @return the value of the function.
		 * @ if {@code param} is {@code NULL}.
		 * @ if the size of {@code param} is
		 * not 2.
		 */
		 //override
		double value(const double& x, double ... param)
		{
			validate_parameters(param);
			return Sigmoid::value(x, param[0], param[1]);
		}

		/**
		 * Computes the value of the gradient at {@code x}.
		 * The components of the gradient vector are the partial
		 * derivatives of the function with respect to each of the
		 * <em>parameters</em> (lower asymptote and higher asymptote).
		 *
		 * @param x Value at which the gradient must be computed.
		 * @param param Values for lower asymptote and higher asymptote.
		 * @return the gradient vector at {@code x}.
		 * @ if {@code param} is {@code NULL}.
		 * @ if the size of {@code param} is
		 * not 2.
		 */
		 //override
		std::vector<double> gradient(const double& x, double ... param)
		{
			validate_parameters(param);
			const double inv_exp1 = 1 / (1 + std::exp(-x));
			return std::vector<double> { 1 - inv_exp1, inv_exp1 };
		}

	private:
		/**
		 * Validates parameters to ensure they are appropriate for the evaluation of
		 * the {@link #value(double,std::vector<double>)} and {@link #gradient(double,std::vector<double>)}
		 * methods.
		 *
		 * @param param Values for lower and higher asymptotes.
		 * @ if {@code param} is {@code NULL}.
		 * @ if the size of {@code param} is
		 * not 2.
		 */
		void validate_parameters(const std::vector<double>& param)
		{
			//Math_Utils::check_not_null(param);
			Math_Utils::check_dimension(param.size(), 2);
		}
	}
};