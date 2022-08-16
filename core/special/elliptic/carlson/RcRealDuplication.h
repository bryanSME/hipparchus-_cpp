#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
 //package org.hipparchus.special.elliptic.carlson;

 //import org.hipparchus.util.FastMath;
 //import org.hipparchus.util.Math_Arrays;
#include <vector>
#include "RealDuplication.h"
#include "../../../../core/util/MathArrays.h"

/** Duplication algorithm for Carlson R<sub>C</sub> elliptic integral.
 * @since 2.0
 */
class Rc_Real_Duplication : public Real_Duplication
{
public:
	/** Constant term in R<sub>C</sub> polynomial. */
	static constexpr double S0{ 80080 };

	/** Coefficient of s² in R<sub>C</sub> polynomial. */
	static constexpr double S2{ 24024 };

	/** Coefficient of s³ in R<sub>C</sub> polynomial. */
	static constexpr double S3{ 11440 };

	/** Coefficient of s⁴ in R<sub>C</sub> polynomial. */
	static constexpr double S4{ 30030 };

	/** Coefficient of s⁵ in R<sub>C</sub> polynomial. */
	static constexpr double S5{ 32760 };

	/** Coefficient of s⁶ in R<sub>C</sub> polynomial. */
	static constexpr double S6{ 61215 };

	/** Coefficient of s⁷ in R<sub>C</sub> polynomial. */
	static constexpr double S7{ 90090 };

	/** Denominator in R<sub>C</sub> polynomial. */
	static constexpr double DENOMINATOR{ 80080 };

	/** Simple constructor.
	 * @param x first symmetric variable of the integral
	 * @param y second symmetric variable of the integral
	 */
	Rc_Real_Duplication(const double& x, const double& y)
	{
		Real_Duplication({ x, y });
		//super(x, y);
	}

protected:
	/** {@inherit_doc} */
	//override
	void initial_mean_point(std::vector<double>& va)
	{
		va[2] = (va[0] + va[1] * 2) / 3.0;
	}

	/** {@inherit_doc} */
	//override
	double convergence_criterion(const double r, const double max)
	{
		return max / std::sqrt(std::sqrt(std::sqrt(r * 3.0)));
	}

	/** {@inherit_doc} */
	//override
	void update(const int& m, std::vector<double>& va_m, const std::vector<double>& sqrt_m, const  double& four_m)
	{
		const auto lambda_a = sqrt_m[0] * sqrt_m[1] * 2;
		auto lambda_b = va_m[1];
		va_m[0] = Math_Arrays::linear_combination(0.25, va_m[0], 0.25, lambda_a, 0.25, lambda_b); // xₘ
		va_m[1] = Math_Arrays::linear_combination(0.25, va_m[1], 0.25, lambda_a, 0.25, lambda_b); // yₘ
		va_m[2] = Math_Arrays::linear_combination(0.25, va_m[2], 0.25, lambda_a, 0.25, lambda_b); // aₘ
	}

	/** {@inherit_doc} */
	//override
	double evaluate(const std::vector<double>& va0, const double& a_m, const double& four_m)
	{
		// compute the single polynomial independent variable
		const double s = (va0[1] - va0[2]) / (a_m * four_m);

		// evaluate integral using equation 2.13 in Carlson[1995]
		const double poly = ((((((S7 * s + S6) * s + S5) * s + S4) * s + S3) * s + S2) * s * s + S0) / DENOMINATOR;
		return poly / std::sqrt(a_m);
	}
};