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

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Implementation of the triangular real distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Triangular_distribution">
   * Triangular distribution (Wikipedia)</a>
   */
class Triangular_Distribution : Abstract_Real_Distribution
{
private:
	/** Lower limit of this distribution (inclusive). */
	const double my_a;
	/** Upper limit of this distribution (inclusive). */
	const double my_b;
	/** Mode of this distribution. */
	const double my_c;

public:

	/**
	 * Creates a triangular real distribution using the given lower limit, * upper limit, and mode.
	 *
	 * @param a Lower limit of this distribution (inclusive).
	 * @param b Upper limit of this distribution (inclusive).
	 * @param c Mode of this distribution.
	 * @ if {@code a >= b} or if {@code c > b}.
	 * @ if {@code c < a}.
	 */
	Triangular_Distribution(const double& a, const double& c, double b) : my_a{ a }, my_b{ b }, my_c{ c }
	{
		super();

		if (a >= b)
		{
			throw (
				hipparchus::exception::Localized_Core_Formats_Type::LOWER_BOUND_NOT_BELOW_UPPER_BOUND, a, b, false);
		}
		if (c < a)
		{
			throw (
				hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, c, a, true);
		}
		if (c > b)
		{
			throw (
				hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, c, b, true);
		}
	}

	/**
	 * Returns the mode {@code c} of this distribution.
	 *
	 * @return the mode {@code c} of this distribution
	 */
	double get_mode() const
	{
		return my_c;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For lower limit {@code a}, upper limit {@code b} and mode {@code c}, the
	 * PDF is given by
	 * <ul>
	 * <li>{@code 2 * (x - a) / [(b - a) * (c - a)]} if {@code a <= x < c},</li>
	 * <li>{@code 2 / (b - a)} if {@code x = c},</li>
	 * <li>{@code 2 * (b - x) / [(b - a) * (b - c)]} if {@code c < x <= b},</li>
	 * <li>{@code 0} otherwise.
	 * </ul>
	 */
	 //override
	double density(double x)
	{
		if (x < my_a)
		{
			return 0;
		}
		if (my_a <= x && x < my_c)
		{
			double divident = 2 * (x - a);
			double divisor = (b - a) * (c - a);
			return divident / divisor;
		}
		if (x == my_c)
		{
			return 2 / (my_b - my_a);
		}
		if (my_c < x && x <= my_b)
		{
			const double divident = 2 * (b - x);
			const double divisor = (b - a) * (b - c);
			return divident / divisor;
		}
		return 0;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For lower limit {@code a}, upper limit {@code b} and mode {@code c}, the
	 * CDF is given by
	 * <ul>
	 * <li>{@code 0} if {@code x < a},</li>
	 * <li>{@code (x - a)^2 / [(b - a) * (c - a)]} if {@code a <= x < c},</li>
	 * <li>{@code (c - a) / (b - a)} if {@code x = c},</li>
	 * <li>{@code 1 - (b - x)^2 / [(b - a) * (b - c)]} if {@code c < x <= b},</li>
	 * <li>{@code 1} if {@code x > b}.</li>
	 * </ul>
	 */
	 //override
	double cumulative_probability(const double& x)
	{
		if (x < a)
		{
			return 0;
		}
		if (a <= x && x < c)
		{
			const double divident = (x - a) * (x - a);
			const double divisor = (b - a) * (c - a);
			return divident / divisor;
		}
		if (x == c)
		{
			return (c - a) / (b - a);
		}
		if (c < x && x <= b)
		{
			const double divident = (b - x) * (b - x);
			const double divisor = (b - a) * (b - c);
			return 1 - (divident / divisor);
		}
		return 1;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For lower limit {@code a}, upper limit {@code b}, and mode {@code c}, * the mean is {@code (a + b + c) / 3}.
	 */
	 //override
	double get_numerical_mean() const
	{
		return (a + b + c) / 3;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For lower limit {@code a}, upper limit {@code b}, and mode {@code c}, * the variance is {@code (a^2 + b^2 + c^2 - a * b - a * c - b * c) / 18}.
	 */
	 //override
	double get_numerical_variance() const
	{
		return (a * a + b * b + c * c - a * b - a * c - b * c) / 18;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is equal to the lower limit parameter
	 * {@code a} of the distribution.
	 *
	 * @return lower bound of the support
	 */
	 //override
	double get_support_lower_bound() const
	{
		return my_a;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is equal to the upper limit parameter
	 * {@code b} of the distribution.
	 *
	 * @return upper bound of the support
	 */
	 //override
	double get_support_upper_bound() const
	{
		return my_b;
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

	/** {@inherit_doc} */
	//override
	double inverse_cumulative_probability(const double& p) const
	{
		Math_Utils::check_range_inclusive(p, 0, 1);

		if (p == 0)
		{
			return my_a;
		}
		if (p == 1)
		{
			return my_b;
		}
		if (p < (c - a) / (b - a))
		{
			return my_a + std::sqrt(p * (my_b - my_a) * (my_c - my_a));
		}
		return my_b - std::sqrt((1 - p) * (my_b - my_a) * (my_b - my_c));
	}
}
