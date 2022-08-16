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
#include <vector>
#include <cmath>

/** Duplication algorithm for Carlson symmetric forms.
 * <p>
 * The algorithms are described in B. C. Carlson 1995 paper
 * "Numerical computation of real or complex elliptic integrals", with
 * improvements described in the appendix of B. C. Carlson and James FitzSimons
 * 2000 paper "Reduction theorems for elliptic integrands with the square root
 * of two quadratic factors". They are also described in
 * <a href="https://dlmf.nist.gov/19.36#i">section 19.36(i)</a>
 * of Digital Library of Mathematical Functions.
 * </p>
 * @since 2.0
 */
class Real_Duplication
{
private:
	/** Max number of iterations. */
	static constexpr int M_MAX{ 16 };

	/** Symmetric variables of the integral, plus mean point. */
	std::vector<double> my_initial_va;

	/** Convergence criterion. */
	double my_q;

public:
	Real_Duplication() = default;

	/** Constructor.
	 * @param v symmetric variables of the integral
	 */
	Real_Duplication(const std::vector<double>& v)
	{
		my_initial_va = v;
		initial_mean_point(my_initial_va);

		double max{};
		const auto a0 = my_initial_va[v.size()];
		for (const auto& vi : v)
		{
			max = std::max(max, std::abs(a0 - vi));
		}
		throw std::exception("not fully implemented - Real_Duplication");
		//my_q = convergence_criterion(std::ulp(1.0), max);
	}

	/** Compute Carlson elliptic integral.
	 * @return Carlson elliptic integral
	 */
	double integral()
	{
		// duplication iterations
		const auto n = my_initial_va.size() - 1;
		const auto va_m = my_initial_va;
		auto sqrt_m = std::vector<double>(n);
		double four_m{ 1.0 };
		for (int m{}; m < M_MAX; ++m)
		{
			if (m > 0 && my_q < four_m * std::abs(va_m[n]))
			{
				// convergence reached
				break;
			}

			// apply duplication once more
			// (we know that {Field}std::complex<double>.sqrt() returns the root with nonnegative real part)
			for (int i{}; i < n; ++i)
			{
				sqrt_m[i] = std::sqrt(va_m[i]);
			}
			update(m, va_m, sqrt_m, four_m);

			four_m *= 4;
		}

		return evaluate(my_initial_va, va_m[n], four_m);
	}

protected:
	/** Get the i<sup>th</sup> symmetric variable.
	 * @param i index of the variable
	 * @return i<sup>th</sup> symmetric variable
	 */
	double get_vi(const int& i)
	{
		return my_initial_va[i];
	}

	/** Compute initial mean point.
	 * <p>
	 * The initial mean point is put as the last array element
	 * </>
	 * @param va symmetric variables of the integral (plus placeholder for initial mean point)
	 */
	virtual void initial_mean_point(std::vector<double> va);

	/** Compute convergence criterion.
	 * @param r relative tolerance
	 * @param max max(|a0-v[i]|)
	 * @return convergence criterion
	 */
	virtual double convergence_criterion(double r, double max);

	/** Update reduced variables in place.
	 * <ul>
	 *  <li>vₘ₊₁|i] ← (vₘ[i] + λₘ) / 4</li>
	 *  <li>aₘ₊₁ ← (aₘ + λₘ) / 4</li>
	 * </ul>
	 * @param m iteration index
	 * @param va_m reduced variables and mean point (updated in place)
	 * @param sqrt_m square roots of reduced variables
	 * @param four_m 4<sup>m</sup>
	 */
	virtual void update(const int& m, std::vector<double> va_m, std::vector<double> sqrt_m, double four_m);

	/** Evaluate integral.
	 * @param va0 initial symmetric variables and mean point of the integral
	 * @param a_m reduced mean point
	 * @param four_m 4<sup>m</sup>
	 * @return integral value
	 */
	virtual double evaluate(std::vector<double> va0, double a_m, double four_m);
};