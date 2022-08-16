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

#include <cmath>
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Pair;

  /**
   * Factory that creates a
   * <a href="http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature">
   * Gauss-type quadrature rule using Hermite polynomials</a>
   * of the first kind.
   * Such a quadrature rule allows the calculation of improper integrals
   * of a function
   * <p>
   *  \(f(x) e^{-x^2}\)
   * </p>
   * <p>
   * Recurrence relation and weights computation follow
   * <a href="http://en.wikipedia.org/wiki/Abramowitz_and_Stegun">
   * Abramowitz and Stegun, 1964</a>.
   * </p>
   *
   */
class HermiteRule_Factory extends AbstractRule_Factory
{
	/** √π. */
	private static const double SQRT_PI = 1.77245385090551602729;

	/** {@inherit_doc} */
	override
		protected Pair<std::vector<double>, std::vector<double>> compute_rule(const int& number_of_points)

	{
		if (number_of_points == 1)
		{
			// Break recursion.
			return Pair<>(std::vector<double> { 0 }, std::vector<double> { SQRT_PI });
		}

		// find nodes as roots of Hermite polynomial
		const std::vector<double> points = find_roots(number_of_points, Hermite(number_of_points)::ratio);
		enforce_symmetry(points);

		// compute weights
		const std::vector<double> weights = std::vector<double>(number_of_points];
		const Hermite hm1 = Hermite(number_of_points - 1);
		for (int i{}; i < number_of_points; i++)
		{
			const double y = hm1.hNhNm1(points[i])[0];
			weights[i] = SQRT_PI / (number_of_points * y * y);
		}

		return Pair<>(points, weights);
	}

	/** Hermite polynomial, normalized to avoid overflow.
	 * <p>
	 * The regular Hermite polynomials and associated weights are given by:
	 *   <pre>
	 *     H₀(x)   = 1
	 *     H₁(x)   = 2 x
	 *     Hₙ₊₁(x) = 2x Hₙ(x) - 2n Hₙ₋₁(x), and H'ₙ(x) = 2n Hₙ₋₁(x)
	 *     wₙ(xᵢ) = [2ⁿ⁻¹ n! √π]/[n Hₙ₋₁(xᵢ)]²
	 *   </pre>
	 * </p>
	 * <p>
	 * In order to avoid overflow with normalize the polynomials hₙ(x) = Hₙ(x) / √[2ⁿ n!]
	 * so the recurrence relations and weights become:
	 *   <pre>
	 *     h₀(x)   = 1
	 *     h₁(x)   = √2 x
	 *     hₙ₊₁(x) = [√2 x hₙ(x) - √n hₙ₋₁(x)]/√(n+1), and h'ₙ(x) = 2n hₙ₋₁(x)
	 *     uₙ(xᵢ) = √π/[n Nₙ₋₁(xᵢ)²]
	 *   </pre>
	 * </p>
	 */
	private static class Hermite
	{
		/** √2. */
		private static const double SQRT2 = std::sqrt(2);

		/** Degree. */
		private const int degree;

		/** Simple constructor.
		 * @param degree polynomial degree
		 */
		Hermite(const int& degree)
		{
			this.degree = degree;
		}

		/** Compute ratio H(x)/H'(x).
		 * @param x point at which ratio must be computed
		 * @return ratio H(x)/H'(x)
		 */
		public double ratio(double x)
		{
			std::vector<double> h = hNhNm1(x);
			return h[0] / (h[1] * 2 * degree);
		}

		/** Compute Nₙ(x) and Nₙ₋₁(x).
		 * @param x point at which polynomials are evaluated
		 * @return array containing Nₙ(x) at index 0 and Nₙ₋₁(x) at index 1
		 */
		private std::vector<double> hNhNm1(const double& x)
		{
			std::vector<double> h = { SQRT2 * x, 1 };
			double sqrt_n = 1;
			for (const int n{ 1 }; n < degree; n++)
			{
				// apply recurrence relation hₙ₊₁(x) = [√2 x hₙ(x) - √n hₙ₋₁(x)]/√(n+1)
				const double sqrt_np = std::sqrt(n + 1);
				const double hp = (h[0] * x * SQRT2 - h[1] * sqrt_n) / sqrt_np;
				h[1] = h[0];
				h[0] = hp;
				sqrt_n = sqrt_np;
			}
			return h;
		}
	}
}
