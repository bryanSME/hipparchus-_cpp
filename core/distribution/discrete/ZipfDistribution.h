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

  //package org.hipparchus.distribution.discrete;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;

  /**
   * Implementation of the Zipf distribution.
   * <p>
   * <strong>Parameters:</strong>
   * For a random variable {@code X} whose values are distributed according to this
   * distribution, the probability mass function is given by
   * <pre>
   *   P(X = k) = H(N,s) * 1 / k^s    for {@code k = 1,2,...,N}.
   * </pre>
   * {@code H(N,s)} is the normalizing constant
   * which corresponds to the generalized harmonic number of order N of s.
   * <p>
   * <ul>
   * <li>{@code N} is the number of elements</li>
   * <li>{@code s} is the exponent</li>
   * </ul>
   *
   * @see <a href="https://en.wikipedia.org/wiki/Zipf's_law">Zipf's law (Wikipedia)</a>
   * @see <a href="https://en.wikipedia.org/wiki/Harmonic_number#Generalized_harmonic_numbers">Generalized harmonic numbers</a>
   */
class Zipf_Distribution : Abstract_Integer_Distribution
{
protected:
	/**
		* Used by {@link #get_numerical_mean()}.
		*
		* @return the mean of this distribution
		*/
	double calculate_numerical_mean() const
	{
		const int N = get_number_of_elements();
		const double s = get_exponent();

		const double Hs1 = generalized_harmonic(N, s - 1);
		const double Hs = my_nth_harmonic;

		return Hs1 / Hs;
	}
	/**
	 * Used by {@link #get_numerical_variance()}.
	 *
	 * @return the variance of this distribution
	 */
	double calculate_numerical_variance()
	{
		const int N = get_number_of_elements();
		const double s = get_exponent();

		const double Hs2 = generalized_harmonic(N, s - 2);
		const double Hs1 = generalized_harmonic(N, s - 1);
		const double Hs = nth_harmonic;

		return (Hs2 / Hs) - ((Hs1 * Hs1) / (Hs * Hs));
	}

private:
	/** Number of elements. */
	const int number_of_elements;
	/** Exponent parameter of the distribution. */
	const double exponent;
	/** Cached values of the nth generalized harmonic. */
	const double nth_harmonic;
	/** Cached numerical mean */
	double numerical_mean = std::numeric_limits<double>::quiet_NaN();
	/** Whether or not the numerical mean has been calculated */
	bool numerical_mean_is_calculated;
	/** Cached numerical variance */
	double numerical_variance = std::numeric_limits<double>::quiet_NaN();
	/** Whether or not the numerical variance has been calculated */
	bool numerical_variance_is_calculated;

	/**
	 * Calculates the Nth generalized harmonic number. See
	 * <a href="http://mathworld.wolfram.com/HarmonicSeries.html">Harmonic
	 * Series</a>.
	 *
	 * @param n Term in the series to calculate (must be larger than 1)
	 * @param m Exponent (special case {@code m = 1} is the harmonic series).
	 * @return the n<sup>th</sup> generalized harmonic number.
	 */
	double generalized_harmonic(const int& n, const double& m) const
	{
		double value{};
		for (int k{ n }; k > 0; --k)
		{
			value += 1.0 / std::pow(k, m);
		}
		return value;
	}

public:
	/**
	 * Create a Zipf distribution with the given number of elements and
	 * exponent.
	 *
	 * @param number_of_elements Number of elements.
	 * @param exponent Exponent.
	 * @exception  if {@code number_of_elements <= 0}
	 * or {@code exponent <= 0}.
	 */
	Zipf_Distribution(const int& number_of_elements, const double& exponent)
		:
		my_number_of_elements{ number_of_elements },
		my_exponent{ exponent }
		my_nth_harmonic{ generalized_harmonic(number_of_elements, exponent) }
	{
		if (number_of_elements <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSION, number_of_elements);
		}
		if (exponent <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::EXPONENT, exponent);
		}
	}

	/**
	 * Get the number of elements (e.g. corpus size) for the distribution.
	 *
	 * @return the number of elements
	 */
	int get_number_of_elements() const
	{
		return my_number_of_elements;
	}

	/**
	 * Get the exponent characterizing the distribution.
	 *
	 * @return the exponent
	 */
	double get_exponent() const
	{
		return my_exponent;
	}

	/** {@inherit_doc} */
	//override
	double probability(const int& x) const
	{
		if (x <= 0 || x > my_number_of_elements)
		{
			return 0.0;
		}

		return (1.0 / std::pow(x, my_exponent)) / my_nth_harmonic;
	}

	/** {@inherit_doc} */
	//override
	double log_probability(const int& x) const
	{
		if (x <= 0 || x > my_number_of_elements)
		{
			return -INFINITY;
		}

		return -std::log(x) * my_exponent - std::log(my_nth_harmonic);
	}

	/** {@inherit_doc} */
	//override
	double cumulative_probability(const int& x) const
	{
		if (x <= 0)
		{
			return 0.0;
		}
		if (x >= my_number_of_elements)
		{
			return 1.0;
		}

		return generalized_harmonic(x, my_exponent) / my_nth_harmonic;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For number of elements {@code N} and exponent {@code s}, the mean is
	 * {@code Hs1 / Hs}, where
	 * <ul>
	 *  <li>{@code Hs1 = generalized_harmonic(N, s - 1)},</li>
	 *  <li>{@code Hs = generalized_harmonic(N, s)}.</li>
	 * </ul>
	 */
	 //override
	double get_numerical_mean()
	{
		if (!my_numerical_mean_is_calculated)
		{
			my_numerical_mean = calculate_numerical_mean();
			my_numerical_mean_is_calculated = true;
		}
		return my_numerical_mean;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For number of elements {@code N} and exponent {@code s}, the mean is
	 * {@code (Hs2 / Hs) - (Hs1^2 / Hs^2)}, where
	 * <ul>
	 *  <li>{@code Hs2 = generalized_harmonic(N, s - 2)},</li>
	 *  <li>{@code Hs1 = generalized_harmonic(N, s - 1)},</li>
	 *  <li>{@code Hs = generalized_harmonic(N, s)}.</li>
	 * </ul>
	 */
	 //override
	double get_numerical_variance() const
	{
		if (!numerical_variance_is_calculated)
		{
			numerical_variance = calculate_numerical_variance();
			numerical_variance_is_calculated = true;
		}
		return numerical_variance;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always 1 no matter the parameters.
	 *
	 * @return lower bound of the support (always 1)
	 */
	 //override
	int get_support_lower_bound() const
	{
		return 1;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is the number of elements.
	 *
	 * @return upper bound of the support
	 */
	 //override
	int get_support_upper_bound() const
	{
		return get_number_of_elements();
	}

	/**
	 * {@inherit_doc}
	 *
	 * The support of this distribution is connected.
	 *
	 * @return {@code true}
	 */
	 //override
	bool is_support_connected() const
	{
		return true;
	}
};