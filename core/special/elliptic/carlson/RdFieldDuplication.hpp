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

 //import org.hipparchus.Calculus_Field_Element;
 //import org.hipparchus.complex.std::complex<double>;
 //import org.hipparchus.complex.Field_Complex<double>;
 //import org.hipparchus.util.FastMath;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/** Duplication algorithm for Carlson R<sub>D</sub> elliptic integral.
 * @param <T> type of the field elements (really {@link std::complex<double>} or {@link Field_Complex<double>})
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Rd_Field_Duplication : Field_Duplication<T>
{
private:
	/** Partial sum. */
	T my_sum;

	/** Simple constructor.
	 * @param x first symmetric variable of the integral
	 * @param y second symmetric variable of the integral
	 * @param z third symmetric variable of the integral
	 */
	Rd_Field_Duplication(const T& x, const T y, const T z)
	{
		super(x, y, z);
		my_sum = x.get_field().get_zero();
	}

protected:
	/** {@inherit_doc} */
	//override
	void initial_mean_point(const std::vector<T> va)
	{
		va[3] = va[0].add(va[1]).add(va[2].multiply(3)).divide(5.0);
	}

	/** {@inherit_doc} */
	//override
	T convergence_criterion(const T r, const T max)
	{
		return max.divide(std::sqrt(std::sqrt(std::sqrt(r.multiply(0.25)))));
	}

	/** {@inherit_doc} */
	//override
	void update(const int m, const std::vector<T> va_m, const std::vector<T> sqrt_m, const  double four_m)
	{
		// equation 2.29 in Carlson[1995]
		const T lambda_a = sqrt_m[0].multiply(sqrt_m[1]);
		const T lambda_b = sqrt_m[0].multiply(sqrt_m[2]);
		const T lambda_c = sqrt_m[1].multiply(sqrt_m[2]);

		// running sum in equation 2.34 in Carlson[1995]
		const T lambda = lambda_a.add(lambda_b).add(lambda_c);
		my_sum = my_sum.add(va_m[2].add(lambda).multiply(sqrt_m[2]).multiply(four_m).reciprocal());

		// equations 2.29 and 2.30 in Carlson[1995]
		va_m[0] = va_m[0].linear_combination(0.25, va_m[0], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // xₘ
		va_m[1] = va_m[1].linear_combination(0.25, va_m[1], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // yₘ
		va_m[2] = va_m[2].linear_combination(0.25, va_m[2], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // zₘ
		va_m[3] = va_m[3].linear_combination(0.25, va_m[3], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // aₘ
	}

	/** {@inherit_doc} */
	//override
	T evaluate(const std::vector<T> va0, const T a_m, const  double four_m)
	{
		// compute symmetric differences
		const T inv = a_m.multiply(four_m).reciprocal();
		const T& big_x = va0[3].subtract(va0[0]).multiply(inv);
		const T& big_y = va0[3].subtract(va0[1]).multiply(inv);
		const T& big_z = big_x.add(big_y).divide(-3);
		const T& big_x_y = big_x.multiply(big_y);
		const T& big_z2 = big_z.multiply(big_z);

		// compute elementary symmetric functions (we already know e1 = 0 by construction)
		const T e2 = big_x_y.subtract(big_z2.multiply(6));
		const T e3 = big_x_y.multiply(3).subtract(big_z2.multiply(8)).multiply(big_z);
		const T e4 = big_x_y.subtract(big_z2).multiply(3).multiply(big_z2);
		const T e5 = big_x_y.multiply(big_z2).multiply(big_z);

		const T e2e2 = e2.multiply(e2);
		const T e2e3 = e2.multiply(e3);
		const T e2e4 = e2.multiply(e4);
		const T e2e5 = e2.multiply(e5);
		const T e3e3 = e3.multiply(e3);
		const T e3e4 = e3.multiply(e4);
		const T e2e2e2 = e2e2.multiply(e2);
		const T e2e2e3 = e2e2.multiply(e3);

		// evaluate integral using equation 19.36.1 in DLMF
		// (which add more terms than equation 2.7 in Carlson[1995])
		const T poly = e3e4.add(e2e5).multiply(Rd_Real_Duplication.E3_E4_P_E2_E5).
			add(e2e2e3.multiply(Rd_Real_Duplication.E2_E2_E3)).
			add(e2e4.multiply(Rd_Real_Duplication.E2_E4)).
			add(e3e3.multiply(Rd_Real_Duplication.E3_E3)).
			add(e2e2e2.multiply(Rd_Real_Duplication.E2_E2_E2)).
			add(e5.multiply(Rd_Real_Duplication.E5)).
			add(e2e3.multiply(Rd_Real_Duplication.E2_E3)).
			add(e4.multiply(Rd_Real_Duplication.E4)).
			add(e2e2.multiply(Rd_Real_Duplication.E2_E2)).
			add(e3.multiply(Rd_Real_Duplication.E3)).
			add(e2.multiply(Rd_Real_Duplication.E2)).
			add(Rd_Real_Duplication.CONSTANT).
			divide(Rd_Real_Duplication.DENOMINATOR);
		const T poly_term = poly.divide(a_m.multiply(std::sqrt(a_m)).multiply(four_m));

		return poly_term.add(my_sum.multiply(3));
	}
};