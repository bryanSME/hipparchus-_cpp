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
  //package org.hipparchus.analysis.integration.gauss;

  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Pair;

  /**
   * Factory that creates Gauss-type quadrature rule using Legendre polynomials.
   * In this implementation, the lower and upper bounds of the natural interval
   * of integration are -1 and 1, respectively.
   * The Legendre polynomials are evaluated using the recurrence relation
   * presented in <a href="http://en.wikipedia.org/wiki/Abramowitz_and_Stegun">
   * Abramowitz and Stegun, 1964</a>.
   *
   */
class LegendreRule_Factory extends AbstractRule_Factory
{
	/** {@inherit_doc} */
	//override
	protected Pair<std::vector<double>, std::vector<double>> compute_rule(const int& number_of_points)

	{
		if (number_of_points == 1)
		{
			// Break recursion.
			return Pair<>(std::vector<double> { 0 }, std::vector<double> { 2 });
		}

		// find nodes as roots of Legendre polynomial
		const Legendre p = Legendre(number_of_points);
		const std::vector<double> points = find_roots(number_of_points, p::ratio);
		enforce_symmetry(points);

		// compute weights
		const std::vector<double> weights = std::vector<double>(number_of_points];
		for (int i{}; i <= number_of_points / 2; i++)
		{
			const double c = points[i];
			const std::vector<double> pKpKm1 = p.pNpNm1(c);
			const double d = number_of_points * (pKpKm1[1] - c * pKpKm1[0]);
			weights[i] = 2 * (1 - c * c) / (d * d);

			// symmetrical point
			const int idx = number_of_points - i - 1;
			weights[idx] = weights[i];
		}

		return Pair<>(points, weights);
	}

	/** Legendre polynomial. */
	private static class Legendre
	{
		/** Degree. */
		private int degree;

		/** Simple constructor.
		 * @param degree polynomial degree
		 */
		Legendre(const int& degree)
		{
			this.degree = degree;
		}

		/** Compute ratio P(x)/P'(x).
		 * @param x point at which ratio must be computed
		 * @return ratio P(x)/P'(x)
		 */
		public double ratio(double x)
		{
			double pm = 1;
			double p = x;
			double d = 1;
			for (const int n{ 1 }; n < degree; n++)
			{
				// apply recurrence relations (n+1) P_n+1(x)  = (2n+1) x P_n(x) - n P_n-1(x)
				// and                              P'_n+1(x) = (n+1) P_n(x) + x P'_n(x)
				const double pp = (p * (x * (2 * n + 1)) - pm * n) / (n + 1);
				d = p * (n + 1) + d * x;
				pm = p;
				p = pp;
			}
			return p / d;
		}

		/** Compute P_n(x) and P_n-1(x).
		 * @param x point at which polynomials are evaluated
		 * @return array containing P_n(x) at index 0 and P_n-1(x) at index 1
		 */
		private std::vector<double> pNpNm1(const double& x)
		{
			std::vector<double> p = { x, 1 };
			for (const int n{ 1 }; n < degree; n++)
			{
				// apply recurrence relation (n+1) P_pi+1(x) = (2n+1) x P_n(x) - n P_n-1(x)
				const double pp = (p[0] * (x * (2 * n + 1)) - p[1] * n) / (n + 1);
				p[1] = p[0];
				p[0] = pp;
			}
			return p;
		}
	}
}
