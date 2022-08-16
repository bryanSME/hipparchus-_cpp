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

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
#include <vector>
#include "Gamma.h"
#include "../util/MathUtils.h"
#include "../util/ContinuedFraction.h"

  /**
   * <p>
   * This is a utility class that provides computation methods related to the
   * Beta family of functions.
   * </p>
   * <p>
   * Implementation of {@link #log_beta(double, double)} is based on the
   * algorithms described in
   * <ul>
   * <li><a href="http://dx.doi.org/10.1145/22721.23109">Didonato and Morris
   *     (1986)</a>, <em>Computation of the Incomplete Gamma Function Ratios
   *     and their Inverse</em>, TOMS 12(4), 377-393,</li>
   * <li><a href="http://dx.doi.org/10.1145/131766.131776">Didonato and Morris
   *     (1992)</a>, <em>Algorithm 708: Significant Digit Computation of the
   *     Incomplete Beta Function Ratios</em>, TOMS 18(3), 360-373,</li>
   * </ul>
   * and implemented in the
   * <a href="http://www.dtic.mil/docs/citations/ADA476840">NSWC Library of Mathematical Functions</a>, * available
   * <a href="http://www.ualberta.ca/CNS/RESEARCH/Software/NumericalNSWC/site.html">here</a>.
   * This library is "approved for public release", and the
   * <a href="http://www.dtic.mil/dtic/pdf/announcements/CopyrightGuidance.pdf">Copyright guidance</a>
   * indicates that unless otherwise stated in the code, all FORTRAN functions in
   * this library are license free. sin_ce no such notice appears in the code these
   * functions can safely be ported to Hipparchus.
   * </p>
   */
class Beta
{
private:
	/** Maximum allowed numerical error. */
	static constexpr double DEFAULT_EPSILON{ 1E-14 };

	/** The constant value of ½log 2π. */
	static constexpr double HALF_LOG_TWO_PI{ .9189385332046727 };

	/**
	 * <p>
	 * The coefficients of the series expansion of the Δ function. This function
	 * is defined as follows
	 * </p>
	 * <center>Δ(x) = log Γ(x) - (x - 0.5) log a + a - 0.5 log 2π,</center>
	 * <p>
	 * see equation (23) in Didonato and Morris (1992). The series expansion, * which applies for x ≥ 10, reads
	 * </p>
	 * <pre>
	 *                 14
	 *                ====
	 *             1  \                2 n
	 *     Δ(x) = ---  >    d  (10 / x)
	 *             x  /      n
	 *                ====
	 *                n = 0
	 * <pre>
	 */
	static const std::vector<double> DELTA
	{
		.833333333333333333333333333333E-01,
		-.277777777777777777777777752282E-04,
		.793650793650793650791732130419E-07,
		-.595238095238095232389839236182E-09,
		.841750841750832853294451671990E-11,
		-.191752691751854612334149171243E-12,
		.641025640510325475730918472625E-14,
		-.295506514125338232839867823991E-15,
		.179643716359402238723287696452E-16,
		-.139228964661627791231203060395E-17,
		.133802855014020915603275339093E-18,
		-.154246009867966094273710216533E-19,
		.197701992980957427278370133333E-20,
		-.234065664793997056856992426667E-21,
		.171348014966398575409015466667E-22
	};

	/**
	 * Default constructor.  Prohibit instantiation.
	 */
	Beta() = default;

	/**
	 * Returns the value of log Γ(a + b) for 1 ≤ a, b ≤ 2. Based on the
	 * <em>NSWC Library of Mathematics Subroutines</em> double precision
	 * implementation, {@code DGSMLN}. In {@code Beta_Test.test_log_gamma_sum()}, * this private method is accessed through reflection.
	 *
	 * @param a First argument.
	 * @param b Second argument.
	 * @return the value of {@code log(Gamma(a + b))}.
	 * @ if {@code a} or {@code b} is lower than
	 * {@code 1.0} or greater than {@code 2.0}.
	 */
	static double log_gamma_sum(const double& a, const double& b)
	{
		Math_Utils::check_range_inclusive(a, 1, 2);
		Math_Utils::check_range_inclusive(b, 1, 2);

		const double x = (a - 1.0) + (b - 1.0);
		if (x <= 0.5)
		{
			return Gamma::log_gamma1p(1.0 + x);
		}
		if (x <= 1.5)
		{
			return Gamma::log_gamma1p(x) + std::log1p(x);
		}
		return Gamma::log_gamma1p(x - 1.0) + std::log(x * (1.0 + x));
	}

	/**
	 * Returns the value of log[Γ(b) / Γ(a + b)] for a ≥ 0 and b ≥ 10. Based on
	 * the <em>NSWC Library of Mathematics Subroutines</em> double precision
	 * implementation, {@code DLGDIV}. In
	 * {@code Beta_Test.test_log_gamma_minus_log_gamma_sum()}, this private method is
	 * accessed through reflection.
	 *
	 * @param a First argument.
	 * @param b Second argument.
	 * @return the value of {@code log(Gamma(b) / Gamma(a + b))}.
	 * @ if {@code a < 0.0} or {@code b < 10.0}.
	 */
	static double log_gamma_minus_log_gamma_sum(const double& a, const double& b)
	{
		if (a < 0.0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, a, 0.0);
		}
		if (b < 10.0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, b, 10.0);
		}

		/*
		 * d = a + b - 0.5
		 */
		double d;
		double w;
		if (a <= b)
		{
			d = b + (a - 0.5);
			w = delta_minus_delta_sum(a, b);
		}
		else
		{
			d = a + (b - 0.5);
			w = delta_minus_delta_sum(b, a);
		}

		const double u = d * std::log1p(a / b);
		const double v = a * (std::log(b) - 1.0);

		return u <= v ? (w - u) - v : (w - v) - u;
	}

	/**
	 * Returns the value of Δ(b) - Δ(a + b), with 0 ≤ a ≤ b and b ≥ 10. Based
	 * on equations (26), (27) and (28) in Didonato and Morris (1992).
	 *
	 * @param a First argument.
	 * @param b Second argument.
	 * @return the value of {@code Delta(b) - Delta(a + b)}
	 * @ if {@code a < 0} or {@code a > b}
	 * @ if {@code b < 10}
	 */
	static double delta_minus_delta_sum(const double& a, const double& b)
	{
		Math_Utils::check_range_inclusive(a, 0, b);
		if (b < 10)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, b, 10);
		}

		const double h = a / b;
		const double p = h / (1.0 + h);
		const double q = 1.0 / (1.0 + h);
		const double q2 = q * q;
		/*
		 * s[i] = 1 + q + ... - q**(2 * i)
		 */
		auto s = std::vector<double>(DELTA.size());
		s[0] = 1.0;
		for (int i{ 1 }; i < s.size(); i++)
		{
			s[i] = 1.0 + (q + q2 * s[i - 1]);
		}
		/*
		 * w = Delta(b) - Delta(a + b)
		 */
		const double sqrt_t = 10.0 / b;
		const double t = sqrt_t * sqrt_t;
		double w = DELTA[DELTA.size() - 1] * s[s.size() - 1];
		for (int i = DELTA.size() - 2; i >= 0; i--)
		{
			w = t * w + DELTA[i] * s[i];
		}
		return w * p / b;
	}

	/**
	 * Returns the value of Δ(p) + Δ(q) - Δ(p + q), with p, q ≥ 10. Based on
	 * the <em>NSWC Library of Mathematics Subroutines</em> double precision
	 * implementation, {@code DBCORR}. In
	 * {@code Beta_Test.test_sum_delta_minus_delta_sum()}, this private method is
	 * accessed through reflection.
	 *
	 * @param p First argument.
	 * @param q Second argument.
	 * @return the value of {@code Delta(p) + Delta(q) - Delta(p + q)}.
	 * @ if {@code p < 10.0} or {@code q < 10.0}.
	 */
	static double sum_delta_minus_delta_sum(const double p, const double q)
	{
		if (p < 10.0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, p, 10.0);
		}
		if (q < 10.0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, q, 10.0);
		}

		const double& a = std::min(p, q);
		const double b = std::max(p, q);
		const double sqrt_t = 10.0 / a;
		const double t = sqrt_t * sqrt_t;
		double z = DELTA[DELTA.size() - 1];
		for (int i = DELTA.size() - 2; i >= 0; i--)
		{
			z = t * z + DELTA[i];
		}
		return z / a + delta_minus_delta_sum(a, b);
	}

public:
	/**
	 * Returns the
	 * <a href="http://mathworld.wolfram.com/RegularizedBetaFunction.html">
	 * regularized beta function</a> I(x, a, b).
	 *
	 * @param x Value.
	 * @param a Parameter {@code a}.
	 * @param b Parameter {@code b}.
	 * @return the regularized beta function I(x, a, b).
	 * @org.hipparchus.exception.Math_Illegal_State_Exception
	 * if the algorithm fails to converge.
	 */
	static double regularized_beta(const double& x, const double& a, double b)
	{
		return regularized_beta(x, a, b, DEFAULT_EPSILON, std::numeric_limits<int>::max());
	}

	/**
	 * Returns the
	 * <a href="http://mathworld.wolfram.com/RegularizedBetaFunction.html">
	 * regularized beta function</a> I(x, a, b).
	 *
	 * @param x Value.
	 * @param a Parameter {@code a}.
	 * @param b Parameter {@code b}.
	 * @param epsilon When the absolute value of the nth item in the
	 * series is less than epsilon the approximation ceases to calculate
	 * further elements in the series.
	 * @return the regularized beta function I(x, a, b)
	 * @org.hipparchus.exception.Math_Illegal_State_Exception
	 * if the algorithm fails to converge.
	 */
	static double regularized_beta(const double& x, const double& a, double b, double epsilon)
	{
		return regularized_beta(x, a, b, epsilon, std::numeric_limits<int>::max());
	}

	/**
	 * Returns the regularized beta function I(x, a, b).
	 *
	 * @param x the value.
	 * @param a Parameter {@code a}.
	 * @param b Parameter {@code b}.
	 * @param max_iterations Maximum number of "iterations" to complete.
	 * @return the regularized beta function I(x, a, b)
	 * @org.hipparchus.exception.Math_Illegal_State_Exception
	 * if the algorithm fails to converge.
	 */
	static double regularized_beta(const double& x, const double& a, double b, int max_iterations)
	{
		return regularized_beta(x, a, b, DEFAULT_EPSILON, max_iterations);
	}

	/**
	 * Returns the regularized beta function I(x, a, b).
	 *
	 * The implementation of this method is based on:
	 * <ul>
	 * <li>
	 * <a href="http://mathworld.wolfram.com/RegularizedBetaFunction.html">
	 * Regularized Beta Function</a>.</li>
	 * <li>
	 * <a href="http://functions.wolfram.com/06.21.10.0001.01">
	 * Regularized Beta Function</a>.</li>
	 * </ul>
	 *
	 * @param x the value.
	 * @param a Parameter {@code a}.
	 * @param b Parameter {@code b}.
	 * @param epsilon When the absolute value of the nth item in the
	 * series is less than epsilon the approximation ceases to calculate
	 * further elements in the series.
	 * @param max_iterations Maximum number of "iterations" to complete.
	 * @return the regularized beta function I(x, a, b)
	 * @org.hipparchus.exception.Math_Illegal_State_Exception
	 * if the algorithm fails to converge.
	 */
	static double regularized_beta(const double& x, const double& a, const double& b, const double& epsilon, const int& max_iterations)
	{
		if (std::isnan(x) ||
			std::isnan(a) ||
			std::isnan(b) ||
			x < 0 ||
			x > 1 ||
			a <= 0 ||
			b <= 0)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (x > (a + 1) / (2 + b + a) &&
			1 - x <= (b + 1) / (2 + b + a))
		{
			return 1 - regularized_beta(1 - x, b, a, epsilon, max_iterations);
		}
		auto fraction = Continued_Fraction()
		{
		protected:
			/** {@inherit_doc} */
			//override
			double get_b(const int& n, const double& x)
			{
				if (n % 2 == 0)
				{ // even
					const double m = n / 2.0;
					return (m * (b - m) * x) /
						((a + (2 * m) - 1) * (a + (2 * m)));
				}
				const double m = (n - 1.0) / 2.0;
				return -((a + m) * (a + b + m) * x) /
					((a + (2 * m)) * (a + (2 * m) + 1.0));
			}

			/** {@inherit_doc} */
			//override
			protected double get_a(const int& n, double x)
			{
				return 1.0;
			}
		};
		return std::exp((a * std::log(x)) + (b * std::log1p(-x)) -
			std::log(a) - log_beta(a, b)) *
			1.0 / fraction.evaluate(x, epsilon, max_iterations);
	}

	/**
	 * Returns the value of log B(p, q) for 0 ≤ x ≤ 1 and p, q &gt; 0. Based on the
	 * <em>NSWC Library of Mathematics Subroutines</em> implementation, * {@code DBETLN}.
	 *
	 * @param p First argument.
	 * @param q Second argument.
	 * @return the value of {@code log(Beta(p, q))}, {@code NaN} if
	 * {@code p <= 0} or {@code q <= 0}.
	 */
	static double log_beta(const double& p, const double& q)
	{
		if (std::isnan(p) || std::isnan(q) || (p <= 0.0) || (q <= 0.0))
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		const double& a = std::min(p, q);
		const double b = std::max(p, q);
		if (a >= 10.0)
		{
			const double w = sum_delta_minus_delta_sum(a, b);
			const double h = a / b;
			const double c = h / (1.0 + h);
			const double u = -(a - 0.5) * std::log(c);
			const double v = b * std::log1p(h);
			if (u <= v)
			{
				return (((-0.5 * std::log(b) + HALF_LOG_TWO_PI) + w) - u) - v;
			}
			else
			{
				return (((-0.5 * std::log(b) + HALF_LOG_TWO_PI) + w) - v) - u;
			}
		}
		else if (a > 2.0)
		{
			if (b > 1000.0)
			{
				const int n = static_cast<int>(std::floor(a - 1.0));
				double prod = 1.0;
				double ared = a;
				for (int i{}; i < n; i++)
				{
					ared -= 1.0;
					prod *= ared / (1.0 + ared / b);
				}
				return (std::log(prod) - n * std::log(b)) +
					(Gamma::log_gamma(ared) +
						log_gamma_minus_log_gamma_sum(ared, b));
			}
			else
			{
				double prod1 = 1.0;
				double ared = a;
				while (ared > 2.0)
				{
					ared -= 1.0;
					const double h = ared / b;
					prod1 *= h / (1.0 + h);
				}
				if (b < 10.0)
				{
					double prod2 = 1.0;
					double bred = b;
					while (bred > 2.0)
					{
						bred -= 1.0;
						prod2 *= bred / (ared + bred);
					}
					return std::log(prod1) +
						std::log(prod2) +
						(Gamma::log_gamma(ared) +
							(Gamma::log_gamma(bred) -
								log_gamma_sum(ared, bred)));
				}
				else
				{
					return std::log(prod1) +
						Gamma::log_gamma(ared) +
						log_gamma_minus_log_gamma_sum(ared, b);
				}
			}
		}
		else if (a >= 1.0)
		{
			if (b > 2.0)
			{
				if (b < 10.0)
				{
					double prod = 1.0;
					double bred = b;
					while (bred > 2.0)
					{
						bred -= 1.0;
						prod *= bred / (a + bred);
					}
					return std::log(prod) +
						(Gamma::log_gamma(a) +
							(Gamma::log_gamma(bred) -
								log_gamma_sum(a, bred)));
				}
				else
				{
					return Gamma::log_gamma(a) +
						log_gamma_minus_log_gamma_sum(a, b);
				}
			}
			else
			{
				return Gamma::log_gamma(a) +
					Gamma::log_gamma(b) -
					log_gamma_sum(a, b);
			}
		}
		else
		{
			if (b >= 10.0)
			{
				return Gamma::log_gamma(a) +
					log_gamma_minus_log_gamma_sum(a, b);
			}
			else
			{
				// The following command is the original NSWC implementation.
				// return Gamma::log_gamma(a) +
				// (Gamma::log_gamma(b) - Gamma::log_gamma(a + b));
				// The following command turns out to be more accurate.
				return std::log(Gamma::gamma(a) * Gamma::gamma(b) /
					Gamma::gamma(a + b));
			}
		}
	}
};