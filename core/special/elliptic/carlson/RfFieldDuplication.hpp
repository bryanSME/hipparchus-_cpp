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

/** Duplication algorithm for Carlson R<sub>F</sub> elliptic integral.
 * @param <T> type of the field elements (really {@link std::complex<double>} or {@link Field_Complex<double>})
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Rf_Field_Duplication : pulic Field_Duplication<T>
{
	/** Simple constructor.
	 * @param x first symmetric variable of the integral
	 * @param y second symmetric variable of the integral
	 * @param z third symmetric variable of the integral
	 */
	Rf_Field_Duplication(const T& x, const T y, const T z)
	{
		super(x, y, z);
	}

	/** {@inherit_doc} */
	//override
	protected void initial_mean_point(const std::vector<T> va)
	{
		va[3] = va[0].add(va[1]).add(va[2]).divide(3.0);
	}

	/** {@inherit_doc} */
	//override
	protected T convergence_criterion(const T r, const T max)
	{
		return max.divide(std::sqrt(std::sqrt(std::sqrt(r.multiply(3.0)))));
	}

	/** {@inherit_doc} */
	//override
	protected void update(const int m, const std::vector<T> va_m, const std::vector<T> sqrt_m, const  double four_m)
	{
		// equation 2.3 in Carlson[1995]
		const T lambda_a = sqrt_m[0].multiply(sqrt_m[1]);
		const T lambda_b = sqrt_m[0].multiply(sqrt_m[2]);
		const T lambda_c = sqrt_m[1].multiply(sqrt_m[2]);

		// equations 2.3 and 2.4 in Carlson[1995]
		va_m[0] = va_m[0].linear_combination(0.25, va_m[0], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // xₘ
		va_m[1] = va_m[1].linear_combination(0.25, va_m[1], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // yₘ
		va_m[2] = va_m[2].linear_combination(0.25, va_m[2], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // zₘ
		va_m[3] = va_m[3].linear_combination(0.25, va_m[3], 0.25, lambda_a, 0.25, lambda_b, 0.25, lambda_c); // aₘ
	}

	/** {@inherit_doc} */
	//override
	protected T evaluate(const std::vector<T> va0, const T a_m, const  double four_m)
	{
		// compute symmetric differences
		const T inv = a_m.multiply(four_m).reciprocal();
		const T& big_x = va0[3].subtract(va0[0]).multiply(inv);
		const T& big_y = va0[3].subtract(va0[1]).multiply(inv);
		const T& big_z = big_x.add(big_y).negate();

		// compute elementary symmetric functions (we already know e1 = 0 by construction)
		const T e2 = big_x.multiply(big_y).subtract(big_z.multiply(big_z));
		const T e3 = big_x.multiply(big_y).multiply(big_z);

		const T e2e2 = e2.multiply(e2);
		const T e2e3 = e2.multiply(e3);
		const T e3e3 = e3.multiply(e3);
		const T e2e2e2 = e2e2.multiply(e2);

		// evaluate integral using equation 19.36.1 in DLMF
		// (which add more terms than equation 2.7 in Carlson[1995])
		const T poly = e2e2e2.multiply(Rf_Real_Duplication.E2_E2_E2).
			add(e3e3.multiply(Rf_Real_Duplication.E3_E3)).
			add(e2e3.multiply(Rf_Real_Duplication.E2_E3)).
			add(e2e2.multiply(Rf_Real_Duplication.E2_E2)).
			add(e3.multiply(Rf_Real_Duplication.E3)).
			add(e2.multiply(Rf_Real_Duplication.E2)).
			add(Rf_Real_Duplication.CONSTANT).
			divide(Rf_Real_Duplication.DENOMINATOR);
		return poly.divide(std::sqrt(a_m));
	}

	/** {@inherit_doc} */
	//override
	public T integral()
	{
		const T x = get_vi(0);
		const T y = get_vi(1);
		const T z = get_vi(2);
		if (x.is_zero())
		{
			return complete_integral(y, z);
		}
		else if (y.is_zero())
		{
			return complete_integral(x, z);
		}
		else if (z.is_zero())
		{
			return complete_integral(x, y);
		}
		else
		{
			return super.integral();
		}
	}

	/** Compute Carlson complete elliptic integral R<sub>F</sub>(u, v, 0).
	 * @param x first symmetric variable of the integral
	 * @param y second symmetric variable of the integral
	 * @return Carlson complete elliptic integral R<sub>F</sub>(u, v, 0)
	 */
	private T complete_integral(const T& x, const T& y)
	{
		T x_m = x.sqrt();
		T yM = y.sqrt();

		// iterate down
		for (int i{ 1 }; i < Rf_Real_Duplication.AGM_MAX; ++i)
		{
			const T x_m1 = x_m;
			const T y_m1 = yM;

			// arithmetic mean
			x_m = x_m1.add(y_m1).multiply(0.5);

			// geometric mean
			yM = x_m1.multiply(y_m1).sqrt();

			// convergence (by the inequality of arithmetic and geometric means, this is non-negative)
			if (x_m.subtract(yM).norm() <= 4 * FastMath.ulp(x_m).get_real())
			{
				// convergence has been reached
				break;
			}
		}

		return x_m.add(yM).reciprocal().multiply(x_m.get_pi());
	}
};