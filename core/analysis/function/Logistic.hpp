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

  //import org.hipparchus.analysis.Parametric_Univariate_Function ;
  //import org.hipparchus.analysis.differentiation.Derivative;
  //import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Function;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include <cmath>
#include <vector>
#include "../differentiation/UnivariateDifferentiableFunction.h"
#include "../differentiation/Derivative.h"
#include "../ParametricUnivariateFunction.h"
#include "../../util/MathUtils.h"

/**
 * <a href="http://en.wikipedia.org/wiki/Generalised_logistic_function">
 *  Generalised logistic</a> function.
 *
 */
class Logistic : public Univariate_Differentiable_Function
{
private:
	/** Lower asymptote. */
	const double my_a;
	/** Upper asymptote. */
	const double my_k;
	/** Growth rate. */
	const double my_b;
	/** Parameter that affects near which asymptote maximum growth occurs. */
	const double my_one_over_n;
	/** Parameter that affects the position of the curve along the ordinate axis. */
	const double my_q;
	/** Abscissa of maximum growth. */
	const double my_m;

	/**
	 * @param m_minus_x {@code m - x}.
	 * @param k {@code k}.
	 * @param b {@code b}.
	 * @param q {@code q}.
	 * @param a {@code a}.
	 * @param one_over_n {@code 1 / n}.
	 * @return the value of the function.
	 */
	static double value(double m_minus_x, double k, double b, double q, const double& a, double one_over_n)
	{
		return a + (k - a) / std::pow(1 + q * std::exp(b * m_minus_x), one_over_n);
	}

public:
	/**
	 * @param k If {@code b > 0}, value of the function for x going towards +&infin;.
	 * If {@code b < 0}, value of the function for x going towards -&infin;.
	 * @param m Abscissa of maximum growth.
	 * @param b Growth rate.
	 * @param q Parameter that affects the position of the curve along the
	 * ordinate axis.
	 * @param a If {@code b > 0}, value of the function for x going towards -&infin;.
	 * If {@code b < 0}, value of the function for x going towards +&infin;.
	 * @param n Parameter that affects near which asymptote the maximum
	 * growth occurs.
	 * @ if {@code n <= 0}.
	 */
	Logistic(const double& k, const double& m, const double& b, const double& q, const double& a, const double& n)
		: my_k{ k }, my_m{ m }, my_b{ b }, my_q{ q }, my_a{ a }, my_one_over_n{ 1 / n }
	{
		if (n <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, n, 0);
		}
	}

	/** {@inherit_doc} */
	//override
	double value(double x)
	{
		return value(my_m - x, my_k, my_b, my_q, my_a, my_one_over_n);
	}

	/** {@inherit_doc}
	 */
	 //override
	template<typename T, typename std::enable_if<std::is_base_of<Derivative<T>, T>::value>::type* = nullptr>
	T value(T t)
	{
		return t.negate().add(m).multiply(b).exp().multiply(q).add(1).pow(one_over_n).reciprocal().multiply(k - a).add(a);
	}

	/**
	 * Parametric function where the input array contains the parameters of
	 * the {@link Logistic#Logistic(double,double,double,double,double,double)
	 * logistic function}, ordered as follows:
	 * <ul>
	 *  <li>k</li>
	 *  <li>m</li>
	 *  <li>b</li>
	 *  <li>q</li>
	 *  <li>a</li>
	 *  <li>n</li>
	 * </ul>
	 */
	class Parametric : public Parametric_Univariate_Function
	{
	public:
		/**
		 * Computes the value of the sigmoid at {@code x}.
		 *
		 * @param x Value for which the function must be computed.
		 * @param param Values for {@code k}, {@code m}, {@code b}, {@code q}, * {@code a} and  {@code n}.
		 * @return the value of the function.
		 * @ if {@code param} is {@code NULL}.
		 * @ if the size of {@code param} is
		 * not 6.
		 * @ if {@code param[5] <= 0}.
		 */
		 //override
		double value(const double& x, double ... param)
		{
			validate_parameters(param);
			return Logistic::value(param[1] - x, param[0], param[2], param[3], param[4], 1 / param[5]);
		}

		/**
		 * Computes the value of the gradient at {@code x}.
		 * The components of the gradient vector are the partial
		 * derivatives of the function with respect to each of the
		 * <em>parameters</em>.
		 *
		 * @param x Value at which the gradient must be computed.
		 * @param param Values for {@code k}, {@code m}, {@code b}, {@code q}, * {@code a} and  {@code n}.
		 * @return the gradient vector at {@code x}.
		 * @ if {@code param} is {@code NULL}.
		 * @ if the size of {@code param} is
		 * not 6.
		 * @ if {@code param[5] <= 0}.
		 */
		 //override
		std::vector<double> gradient(const double& x, double ... param)
		{
			validate_parameters(param);

			const double b = param[2];
			const double q = param[3];

			const double m_minus_x = param[1] - x;
			const double one_over_n = 1 / param[5];
			const double exp = std::exp(b * m_minus_x);
			const double q_exp = q * exp;
			const double q_exp1 = q_exp + 1;
			const double factor1 = (param[0] - param[4]) * one_over_n / std::pow(q_exp1, one_over_n);
			const double factor2 = -factor1 / q_exp1;

			// Components of the gradient.
			const double gk = Logistic::value(m_minus_x, 1, b, q, 0, one_over_n);
			const double gm = factor2 * b * q_exp;
			const double gb = factor2 * m_minus_x * q_exp;
			const double gq = factor2 * exp;
			const double ga = Logistic::value(m_minus_x, 0, b, q, 1, one_over_n);
			const double gn = factor1 * std::log(q_exp1) * one_over_n;

			return std::vector<double> { gk, gm, gb, gq, ga, gn };
		}
	private:
		/**
		 * Validates parameters to ensure they are appropriate for the evaluation of
		 * the {@link #value(double,std::vector<double>)} and {@link #gradient(double,std::vector<double>)}
		 * methods.
		 *
		 * @param param Values for {@code k}, {@code m}, {@code b}, {@code q}, * {@code a} and {@code n}.
		 * @ if {@code param} is {@code NULL}.
		 * @ if the size of {@code param} is
		 * not 6.
		 * @ if {@code param[5] <= 0}.
		 */
		void validate_parameters(const std::vector<double>& param)
		{
			//Math_Utils::check_not_null(param);
			Math_Utils::check_dimension(param.size(), 6);
			if (param[5] <= 0)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, param[5], 0);
			}
		}
	}
};