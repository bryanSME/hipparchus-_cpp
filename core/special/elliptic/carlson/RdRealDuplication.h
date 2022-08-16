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
#include "RealDuplication.h"
#include "../../../../core/util/MathArrays.h"

/** Duplication algorithm for Carlson R<sub>D</sub> elliptic integral.
 * @since 2.0
 */
class Rd_Real_Duplication : public Real_Duplication
{
public:
	/** Constant term in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double CONSTANT{ 4084080 };

	/** Coefficient of E₂ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E2{ -875160 };

	/** Coefficient of E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E3{ 680680 };

	/** Coefficient of E₂² in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E2_E2{ 417690 };

	/** Coefficient of E₄ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E4{ -556920 };

	/** Coefficient of E₂E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E2_E3{ -706860 };

	/** Coefficient of E₅ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E5{ 471240 };

	/** Coefficient of E₂³ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E2_E2_E2{ -255255 };

	/** Coefficient of E₃² in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E3_E3{ 306306 };

	/** Coefficient of E₂E₄ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E2_E4{ 612612 };

	/** Coefficient of E₂²E₃ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E2_E2_E3{ 675675 };

	/** Coefficient of E₃E₄+E₂E₅ in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double E3_E4_P_E2_E5{ -540540 };

	/** Denominator in R<sub>J</sub> and R<sub>D</sub> polynomials. */
	static constexpr double DENOMINATOR{ 4084080 };

private:
	/** Partial sum. */
	double my_sum;

public:
	/** Simple constructor.
	 * @param x first symmetric variable of the integral
	 * @param y second symmetric variable of the integral
	 * @param z third symmetric variable of the integral
	 */
	Rd_Real_Duplication(const double& x, const double& y, const double& z)
	{
		Real_Duplication({ x, y, z });
		//super(x, y, z);
		my_sum = 0;
	}

protected:

	/** {@inherit_doc} */
	//override
	void initial_mean_point(std::vector<double>& va)
	{
		va[3] = (va[0] + va[1] + va[2] * 3.0) / 5.0;
	}

	/** {@inherit_doc} */
	//override
	double convergence_criterion(const double& r, const double& max) const
	{
		return max / (std::sqrt(std::sqrt(std::sqrt(r * 0.25))));
	}

	/** {@inherit_doc} */
	//override
	void update(const int m, std::vector<double>& va_m, const std::vector<double> sqrt_m, const  double four_m)
	{
		// equation 2.29 in Carlson[1995]
		const double lambda_a = sqrt_m[0] * sqrt_m[1];
		const double lambda_b = sqrt_m[0] * sqrt_m[2];
		const double lambda_c = sqrt_m[1] * sqrt_m[2];

		// running sum in equation 2.34 in Carlson[1995]
		const double lambda = lambda_a + lambda_b + lambda_c;
		my_sum += 1.0 / ((va_m[2] + lambda) * sqrt_m[2] * four_m);

		// equations 2.29 and 2.30 in Carlson[1995]
		va_m[0] = Math_Arrays::linear_combination(0.25, va_m[0], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // xₘ
		va_m[1] = Math_Arrays::linear_combination(0.25, va_m[1], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // yₘ
		va_m[2] = Math_Arrays::linear_combination(0.25, va_m[2], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // zₘ
		va_m[3] = Math_Arrays::linear_combination(0.25, va_m[3], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // aₘ
	}

	/** {@inherit_doc} */
	//override
	double evaluate(const std::vector<double>& va0, const double& a_m, const double& four_m)
	{
		// compute symmetric differences
		const double inv = 1.0 / (a_m * four_m);
		const double big_x = (va0[3] - va0[0]) * inv;
		const double big_y = (va0[3] - va0[1]) * inv;
		const double big_z = (big_x + big_y) / -3;
		const double big_x_y = big_x * big_y;
		const double big_z2 = big_z * big_z;

		// compute elementary symmetric functions (we already know e1 = 0 by construction)
		const double e2 = big_x_y - big_z2 * 6;
		const double e3 = (big_x_y * 3 - big_z2 * 8) * big_z;
		const double e4 = (big_x_y - big_z2) * 3 * big_z2;
		const double e5 = big_x_y * big_z2 * big_z;

		const double e2e2 = e2 * e2;
		const double e2e3 = e2 * e3;
		const double e2e4 = e2 * e4;
		const double e2e5 = e2 * e5;
		const double e3e3 = e3 * e3;
		const double e3e4 = e3 * e4;
		const double e2e2e2 = e2e2 * e2;
		const double e2e2e3 = e2e2 * e3;

		// evaluate integral using equation 19.36.1 in DLMF
		// (which add more terms than equation 2.7 in Carlson[1995])
		const double poly = ((e3e4 + e2e5) * E3_E4_P_E2_E5 +
			e2e2e3 * E2_E2_E3 +
			e2e4 * E2_E4 +
			e3e3 * E3_E3 +
			e2e2e2 * E2_E2_E2 +
			e5 * E5 +
			e2e3 * E2_E3 +
			e4 * E4 +
			e2e2 * E2_E2 +
			e3 * E3 +
			e2 * E2 +
			CONSTANT) /
			DENOMINATOR;
		const double poly_term = poly / (a_m * std::sqrt(a_m) * four_m);

		return poly_term + my_sum * 3;
	}
};