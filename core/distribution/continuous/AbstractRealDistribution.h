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

#include <exception>
#include <cmath>
#include "../../util/MathUtils.h"

  //import java.io.Serializable;

  //import org.hipparchus.analysis.Univariate_Function;
  //import org.hipparchus.analysis.solvers.Univariate_Solver_Utils;
  //import org.hipparchus.distribution.Real_Distribution;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
#include "../RealDistribution.h"

/**
 * Base class for probability distributions on the reals.
 * <p>
 * Default implementations are provided for some of the methods
 * that do not vary from distribution to distribution.
 */
class Abstract_Real_Distribution : public Real_Distribution
{
private:

	/** Inverse cumulative probability accuracy. */
	const double my_solver_absolute_accuracy;

protected:
	/** Default absolute accuracy for inverse cumulative computation. */
	static constexpr double DEFAULT_SOLVER_ABSOLUTE_ACCURACY = 1e-9;

	/**
	 * @param solver_absolute_accuracy the absolute accuracy to use when
	 * computing the inverse cumulative probability.
	 */
	Abstract_Real_Distribution(const double& solver_absolute_accuracy) : my_solver_absolute_accuracy{ solver_absolute_accuracy } {};

	/**
	 * Create a real distribution with default solver absolute accuracy.
	 */
	Abstract_Real_Distribution() : my_solver_absolute_accuracy{ DEFAULT_SOLVER_ABSOLUTE_ACCURACY } {};

	/**
	 * Returns the solver absolute accuracy for inverse cumulative computation.
	 * You can //override this method in order to use a Brent solver with an
	 * absolute accuracy different from the default.
	 *
	 * @return the maximum absolute error in inverse cumulative probability estimates
	 */
	double get_solver_absolute_accuracy() const
	{
		return my_solver_absolute_accuracy;
	}

public:

	/**
	 * For a random variable {@code X} whose values are distributed according
	 * to this distribution, this method returns {@code P(x0 < X <= x1)}.
	 *
	 * @param x0 Lower bound (excluded).
	 * @param x1 Upper bound (included).
	 * @return the probability that a random variable with this distribution
	 * takes a value between {@code x0} and {@code x1}, excluding the lower
	 * and including the upper endpoint.
	 * @ if {@code x0 > x1}.
	 *
	 * The default implementation uses the identity
	 * {@code P(x0 < X <= x1) = P(X <= x1) - P(X <= x0)}
	 */
	 //override
	double probability(double x0, double x1)
	{
		if (x0 > x1)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT, x0, x1, true);
		}
		return cumulative_probability(x1) - cumulative_probability(x0);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The default implementation returns
	 * <ul>
	 * <li>{@link #get_support_lower_bound()} for {@code p = 0},</li>
	 * <li>{@link #get_support_upper_bound()} for {@code p = 1}.</li>
	 * </ul>
	 */
	 //override
	double inverse_cumulative_probability(const double& p)
	{
		/*
		 * IMPLEMENTATION NOTES
		 * --------------------
		 * Where applicable, use is made of the one-sided Chebyshev inequality
		 * to bracket the root. This inequality states that
		 * P(X - mu >= k * sig) <= 1 / (1 + k^2), * mu: mean, sig: standard deviation. Equivalently
		 * 1 - P(X < mu + k * sig) <= 1 / (1 + k^2), * F(mu + k * sig) >= k^2 / (1 + k^2).
		 *
		 * For k = sqrt(p / (1 - p)), we find
		 * F(mu + k * sig) >= p, * and (mu + k * sig) is an upper-bound for the root.
		 *
		 * Then, introducing Y = -X, mean(Y) = -mu, sd(Y) = sig, and
		 * P(Y >= -mu + k * sig) <= 1 / (1 + k^2), * P(-X >= -mu + k * sig) <= 1 / (1 + k^2), * P(X <= mu - k * sig) <= 1 / (1 + k^2), * F(mu - k * sig) <= 1 / (1 + k^2).
		 *
		 * For k = sqrt((1 - p) / p), we find
		 * F(mu - k * sig) <= p, * and (mu - k * sig) is a lower-bound for the root.
		 *
		 * In cases where the Chebyshev inequality does not apply, geometric
		 * progressions 1, 2, 4, ... and -1, -2, -4, ... are used to bracket
		 * the root.
		 */

		Math_Utils::check_range_inclusive(p, 0, 1);

		double lower_bound = get_support_lower_bound();
		if (p == 0.0)
		{
			return lower_bound;
		}

		double upper_bound = get_support_upper_bound();
		if (p == 1.0)
		{
			return upper_bound;
		}

		const double mu = get_numerical_mean();
		const double sig = std::sqrt(get_numerical_variance());
		const bool chebyshev_applies;
		chebyshev_applies = !(std::isinf(mu) || std::isnan(mu) ||
			std::isinf(sig) || std::isnan(sig));

		if (lower_bound == -INFINITY)
		{
			if (chebyshev_applies)
			{
				lower_bound = mu - sig * std::sqrt((1. - p) / p);
			}
			else
			{
				lower_bound = -1.0;
				while (cumulative_probability(lower_bound) >= p)
				{
					lower_bound *= 2.0;
				}
			}
		}

		if (upper_bound == INFINITY)
		{
			if (chebyshev_applies)
			{
				upper_bound = mu + sig * std::sqrt(p / (1. - p));
			}
			else
			{
				upper_bound = 1.0;
				while (cumulative_probability(upper_bound) < p)
				{
					upper_bound *= 2.0;
				}
			}
		}

		const auto to_solve = Univariate_Function()
		{
			/** {@inherit_doc} */
			//override
			public double value(const double& x)
			{
				return cumulative_probability(x) - p;
			}
		};

		double x = Univariate_Solver_Utils.solve(to_solve, lower_bound, upper_bound, get_solver_absolute_accuracy());

		if (!is_support_connected())
		{
			/* Test for plateau. */
			const double dx = get_solver_absolute_accuracy();
			if (x - dx >= get_support_lower_bound())
			{
				double px = cumulative_probability(x);
				if (cumulative_probability(x - dx) == px)
				{
					upper_bound = x;
					while (upper_bound - lower_bound > dx)
					{
						const double mid_point = 0.5 * (lower_bound + upper_bound);
						if (cumulative_probability(mid_point) < px)
						{
							lower_bound = mid_point;
						}
						else
						{
							upper_bound = mid_point;
						}
					}
					return upper_bound;
				}
			}
		}
		return x;
	}

	/**
	 * {@inherit_doc}
	 * <p>
	 * The default implementation simply computes the logarithm of {@code density(x)}.
	 */
	 //override
	double log_density(const double& x) const
	{
		return std::log(density(x));
	}
};