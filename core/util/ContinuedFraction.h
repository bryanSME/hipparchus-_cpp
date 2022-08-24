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
  //package org.hipparchus.util;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;

  /**
   * Provides a generic means to evaluate continued fractions.  Subclasses simply
   * provided the a and b coefficients to evaluate the continued fraction.
   * <p>
   * References:
   * <ul>
   * <li><a href="http://mathworld.wolfram.com/Continued_Fraction.html">
   * Continued Fraction</a></li>
   * </ul>
   */
class Continued_Fraction
{
private:
	/** Maximum allowed numerical error. */
	static constexpr double DEFAULT_EPSILON{ 10e-9 };

protected:
	/**
	 * Default constructor.
	 */
	Continued_Fraction()
	{
		//super();
	}

	/**
	 * Access the n-th a coefficient of the continued fraction.  sin_ce a can be
	 * a function of the evaluation point, x, that is passed in as well.
	 * @param n the coefficient index to retrieve.
	 * @param x the evaluation point.
	 * @return the n-th a coefficient.
	 */
	virtual double get_a(const int& n, const double& x);

	/**
	 * Access the n-th b coefficient of the continued fraction.  sin_ce b can be
	 * a function of the evaluation point, x, that is passed in as well.
	 * @param n the coefficient index to retrieve.
	 * @param x the evaluation point.
	 * @return the n-th b coefficient.
	 */
	virtual double get_b(const int& n, double x);

public:
	/**
	 * Evaluates the continued fraction at the value x.
	 * @param x the evaluation point.
	 * @return the value of the continued fraction evaluated at x.
	 * @Math_Illegal_State_Exception if the algorithm fails to converge.
	 */
	double evaluate(double x)
	{
		return evaluate(x, DEFAULT_EPSILON, std::numeric_limits<int>::max());
	}

	/**
	 * Evaluates the continued fraction at the value x.
	 * @param x the evaluation point.
	 * @param epsilon maximum error allowed.
	 * @return the value of the continued fraction evaluated at x.
	 * @Math_Illegal_State_Exception if the algorithm fails to converge.
	 */
	double evaluate(const double& x, double epsilon) Math_Illegal_State_Exception
	{
		return evaluate(x, epsilon, std::numeric_limits<int>::max());
	}

	/**
	 * Evaluates the continued fraction at the value x.
	 * @param x the evaluation point.
	 * @param max_iterations maximum number of convergents
	 * @return the value of the continued fraction evaluated at x.
	 * @Math_Illegal_State_Exception if the algorithm fails to converge.
	 * @Math_Illegal_State_Exception if maximal number of iterations is reached
	 */
	double evaluate(const double& x, int max_iterations)
	{
		return evaluate(x, DEFAULT_EPSILON, max_iterations);
	}

	/**
	 * Evaluates the continued fraction at the value x.
	 * <p>
	 * The implementation of this method is based on the modified Lentz algorithm as described
	 * on page 18 ff. in:
	 * <ul>
	 *   <li>
	 *   I. J. Thompson,  A. R. Barnett. "Coulomb and Bessel Functions of std::complex<double> Arguments and Order."
	 *   <a target="_blank" href="http://www.fresco.org.uk/papers/Thompson-JCP64p490.pdf">
	 *   http://www.fresco.org.uk/papers/Thompson-JCP64p490.pdf</a>
	 *   </li>
	 * </ul>
	 * <b>Note:</b> the implementation uses the terms a<sub>i</sub> and b<sub>i</sub> as defined in
	 * <a href="http://mathworld.wolfram.com/Continued_Fraction.html">Continued Fraction @ MathWorld</a>.
	 * </p>
	 *
	 * @param x the evaluation point.
	 * @param epsilon maximum error allowed.
	 * @param max_iterations maximum number of convergents
	 * @return the value of the continued fraction evaluated at x.
	 * @Math_Illegal_State_Exception if the algorithm fails to converge.
	 * @Math_Illegal_State_Exception if maximal number of iterations is reached
	 */
	double evaluate(const double& x, double epsilon, int max_iterations)
	{
		const double small{ 1e-50 };
		double h_prev = get_a(0, x);

		// use the value of small as epsilon criteria for zero checks
		if (Precision::equals(h_prev, 0.0, small))
		{
			h_prev = small;
		}

		int n{ 1 };
		double d_prev{};
		double c_prev = h_prev;
		double hN = h_prev;

		while (n < max_iterations)
		{
			const double a = get_a(n, x);
			const double b = get_b(n, x);

			double dN = a + b * d_prev;
			if (Precision::equals(dN, 0.0, small))
			{
				dN = small;
			}
			double cN = a + b / c_prev;
			if (Precision::equals(cN, 0.0, small))
			{
				cN = small;
			}

			dN = 1 / dN;
			const double delta_n = cN * dN;
			hN = h_prev * delta_n;

			if (std::isinf(hN))
			{
				throw std::exception("not implemented");
				//throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CONTINUED_FRACTION_INFINITY_DIVERGENCE, x);
			}
			if (std::isnan(hN))
			{
				throw std::exception("not implemented");
				//throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CONTINUED_FRACTION_NAN_DIVERGENCE, x);
			}

			if (std::abs(delta_n - 1.0) < epsilon)
			{
				break;
			}

			d_prev = dN;
			c_prev = cN;
			h_prev = hN;
			n++;
		}

		if (n >= max_iterations)
		{
			throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::NON_CONVERGENT_CONTINUED_FRACTION, max_iterations, x);
		}

		return hN;
	}
};