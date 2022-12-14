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
  //package org.hipparchus.random;

#include <type_traits>
#include <exception>
#include <cmath>
#include <string>
#include <vector>
  //import java.io.Serializable;
  //import java.util.Array_list;
  //import java.util.Arrays;
  //import java.util.Collection;
  //import java.util.List;
  //import java.util.Map;
  //import java.util.concurrent.ConcurrentHash_Map;

  //import org.hipparchus.distribution.Enumerated_Distribution;
  //import org.hipparchus.distribution.Integer_Distribution;
  //import org.hipparchus.distribution.Real_Distribution;
  //import org.hipparchus.distribution.continuous.Beta_Distribution;
  //import org.hipparchus.distribution.continuous.Enumerated_Real_Distribution;
  //import org.hipparchus.distribution.continuous.Exponential_Distribution;
  //import org.hipparchus.distribution.continuous.Gamma_Distribution;
  //import org.hipparchus.distribution.continuous.Log_Normal_Distribution;
  //import org.hipparchus.distribution.continuous.Normal_Distribution;
  //import org.hipparchus.distribution.continuous.UniformReal_Distribution;
  //import org.hipparchus.distribution.discrete.Enumerated_Integer_Distribution;
  //import org.hipparchus.distribution.discrete.Poisson_Distribution;
  //import org.hipparchus.distribution.discrete.UniformInteger_Distribution;
  //import org.hipparchus.distribution.discrete.Zipf_Distribution;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Combinatorics_Utils;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Pair;
  //import org.hipparchus.util.Precision;
  //import org.hipparchus.util.Resizable_Double_Array;
#include "RandomGenerator.h"
#include "ForwardingRandomGenerator.h"

  /**
   * A class for generating random data.
   */
class Random_Data_Generator : public Forwarding_Random_Generator, public Random_Generator
{
private:
	/**
	 * Used when generating Exponential samples.
	 * Table containing the constants
	 * q_i = sum_{j=1}^i (ln 2)^j/j! = ln 2 + (ln 2)^2/2 + ... + (ln 2)^i/i!
	 * until the largest representable fraction below 1 is exceeded.
	 *
	 * Note that
	 * 1 = 2 - 1 = exp(ln 2) - 1 = sum_{n=1}^infty (ln 2)^n / n!
	 * thus q_i -> 1 as i -> +inf, * so the higher i, the closer to one we get (the series is not alternating).
	 *
	 * By trying, n = 16 in Java is enough to reach 1.0.
	 */
	static const std::vector<double> EXPONENTIAL_SA_QI;

	/** Map of <classname, switch constant> for continuous distributions */
	static const Map<Class< ? extends Real_Distribution>, Real_Distribution_Sampler> CONTINUOUS_SAMPLERS = ConcurrentHash_Map<>();

	/** Map of <classname, switch constant> for discrete distributions */
	static const Map<Class< ? extends Integer_Distribution>, Integer_DistributionSampler> DISCRETE_SAMPLERS = ConcurrentHash_Map<>();

	/** The default sampler for continuous distributions using the inversion technique. */
	static const Real_Distribution_Sampler DEFAULT_REAL_SAMPLER = (generator, dist)->dist.inverse_cumulative_probability(generator.next_double());

	/** The default sampler for discrete distributions using the inversion technique. */
	static const Integer_DistributionSampler DEFAULT_INTEGER_SAMPLER = (generator, dist)->dist.inverse_cumulative_probability(generator.next_double());

	/** Source of random data */
	const Random_Generator random_generator;

	/** The sampler to be used for the next_zip_f method */
	transient Zipf_Rejection_Inversion_Sampler zipf_sampler;

	/**
	 * Interface for samplers of continuous distributions.
	 */
	 //@Functional_Interface
	interface Real_Distribution_Sampler
	{
		/**
		 * Return the next sample following the given distribution.
		 *
		 * @param generator the random data generator to use
		 * @param distribution the distribution to use
		 * @return the next sample
		 */
		double next_sample(Random_Data_Generator generator, Real_Distribution distribution);
	}

	/**
	 * Interface for samplers of discrete distributions.
	 */
	 //@Functional_Interface
	private interface Integer_DistributionSampler
	{
		/**
		 * Return the next sample following the given distribution.
		 *
		 * @param generator the random data generator to use
		 * @param distribution the distribution to use
		 * @return the next sample
		 */
		int next_sample(Random_Data_Generator generator, Integer_Distribution distribution);
	}

	/**
	 * Initialize tables.
	 */
	static
	{
		/**
		 * Filling EXPONENTIAL_SA_QI table.
		 * Note that we don't want qi = 0 in the table.
		 */
		const double LN2 = std::log(2);
		double qi = 0;
		int i = 1;

		/**
		 * Arithmetic_Utils provides factorials up to 20, so let's use that
		 * limit together with Precision.EPSILON to generate the following
		 * code (a priori, we know that there will be 16 elements, but it is
		 * better to not hardcode it).
		 */
		const Resizable_Double_Array ra = Resizable_Double_Array(20);

		while (qi < 1)
		{
			qi += std::pow(LN2, i) / Combinatorics_Utils.factorial(i);
			ra.add_element(qi);
			++i;
		}

		EXPONENTIAL_SA_QI = ra.get_elements();

		// Continuous samplers

		CONTINUOUS_SAMPLERS.put(Beta_Distribution.class, (generator, dist) ->
		{
			Beta_Distribution beta = (Beta_Distribution)dist;
			return generator.next_beta(beta.get_alpha(), beta.get_beta());
		});

		CONTINUOUS_SAMPLERS.put(Exponential_Distribution.class, (generator, dist)->generator.next_exponential(dist.get_numerical_mean()));

		CONTINUOUS_SAMPLERS.put(Gamma_Distribution.class, (generator, dist) ->
		{
			Gamma_Distribution gamma = (Gamma_Distribution)dist;
			return generator.next_gamma(gamma.get_shape(), gamma.get_scale());
		});

		CONTINUOUS_SAMPLERS.put(Normal_Distribution.class, (generator, dist) ->
		{
			Normal_Distribution normal = (Normal_Distribution)dist;
			return generator.next_normal(normal.get_mean(), normal.get_standard_deviation());
		});

		CONTINUOUS_SAMPLERS.put(Log_Normal_Distribution.class, (generator, dist) ->
		{
			Log_Normal_Distribution log_normal = (Log_Normal_Distribution)dist;
			return generator.next_log_normal(log_normal.get_shape(), log_normal.get_location());
		});

		CONTINUOUS_SAMPLERS.put(UniformReal_Distribution.class, (generator, dist)->generator.next_uniform(dist.get_support_lower_bound(), dist.get_support_upper_bound()));

		CONTINUOUS_SAMPLERS.put(Enumerated_Real_Distribution.class, (generator, dist) ->
		{
			const Enumerated_Real_Distribution edist =
				(Enumerated_Real_Distribution)dist;
			Enumerated_DistributionSampler<Double> sampler =
				generator.new Enumerated_DistributionSampler<Double>(edist.get_pmf());
			return sampler.sample();
		});

		// Discrete samplers

		DISCRETE_SAMPLERS.put(Poisson_Distribution.class, (generator, dist)->generator.next_poisson(dist.get_numerical_mean()));

		DISCRETE_SAMPLERS.put(UniformInteger_Distribution.class, (generator, dist)->generator.next_int(dist.get_support_lower_bound(), dist.get_support_upper_bound()));
		DISCRETE_SAMPLERS.put(Zipf_Distribution.class, (generator, dist) ->
		{
			Zipf_Distribution zipf_dist = (Zipf_Distribution)dist;
			return generator.next_zipf(zipf_dist.get_number_of_elements(), zipf_dist.get_exponent());
		});

		DISCRETE_SAMPLERS.put(Enumerated_Integer_Distribution.class, (generator, dist) ->
		{
			const Enumerated_Integer_Distribution edist =
				(Enumerated_Integer_Distribution)dist;
			Enumerated_DistributionSampler<Integer> sampler =
				generator.new Enumerated_DistributionSampler<Integer>(edist.get_pmf());
			return sampler.sample();
		});
	}

	/**
	 * Construct a Random_Data_Generator with a default Random_Generator as its source of random data.
	 */
	public Random_Data_Generator()
	{
		this(new Well19937c());
	}

	/**
	 * Construct a Random_Data_Generator with a default Random_Generator as its source of random data, initialized
	 * with the given seed value.
	 *
	 * @param seed seed value
	 */
	public Random_Data_Generator(long seed)
	{
		this(new Well19937c(seed));
	}

	/**
	 * Construct a Random_Data_Generator using the given Random_Generator as its source of random data.
	 *
	 * @param random_generator the underlying PRNG
	 * @ if random_generator is NULL
	 */
	private Random_Data_Generator(Random_Generator random_generator)
	{
		//Math_Utils::check_not_null(random_generator);
		this.random_generator = random_generator;
	}

	/**
	 * Factory method to create a {@code Random_data} instance using the supplied
	 * {@code Random_Generator}.
	 *
	 * @param random_generator source of random bits
	 * @return a Random_data using the given Random_Generator to source bits
	 * @ if random_generator is NULL
	 */
	public static Random_Data_Generator of(Random_Generator random_generator)
	{
		return Random_Data_Generator(random_generator);
	}

	/** {@inherit_doc} */
	//override
	protected Random_Generator delegate()
	{
		return random_generator;
	}

	/**
	 * Returns the next pseudo-random beta-distributed value with the given
	 * shape and scale parameters.
	 *
	 * @param alpha First shape parameter (must be positive).
	 * @param beta Second shape parameter (must be positive).
	 * @return beta-distributed random deviate
	 */
	public double next_beta(double alpha, double beta)
	{
		return Cheng_Beta_Sampler.sample(random_generator, alpha, beta);
	}

	/**
	 * Returns the next pseudo-random, exponentially distributed deviate.
	 *
	 * @param mean mean of the exponential distribution
	 * @return exponentially distributed deviate about the given mean
	 */
	public double next_exponential(double mean)
	{
		if (mean <= 0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::MEAN, mean);
		}
		// Step 1:
		double a = 0;
		double u = random_generator.next_double();

		// Step 2 and 3:
		while (u < 0.5)
		{
			a += EXPONENTIAL_SA_QI[0];
			u *= 2;
		}

		// Step 4 (now u >= 0.5):
		u += u - 1;

		// Step 5:
		if (u <= EXPONENTIAL_SA_QI[0])
		{
			return mean * (a + u);
		}

		// Step 6:
		int i = 0; // Should be 1, be we iterate before it in while using 0
		double u2 = random_generator.next_double();
		double umin = u2;

		// Step 7 and 8:
		do
		{
			++i;
			u2 = random_generator.next_double();

			if (u2 < umin)
			{
				umin = u2;
			}

			// Step 8:
		} while (u > EXPONENTIAL_SA_QI[i]); // Ensured to exit since EXPONENTIAL_SA_QI[MAX] = 1

		return mean * (a + umin * EXPONENTIAL_SA_QI[0]);
	}

	/**
	 * Returns the next pseudo-random gamma-distributed value with the given shape and scale parameters.
	 *
	 * @param shape shape parameter of the distribution
	 * @param scale scale parameter of the distribution
	 * @return gamma-distributed random deviate
	 */
	public double next_gamma(double shape, double scale)
	{
		if (shape < 1)
		{
			// [1]: p. 228, Algorithm GS

			while (true)
			{
				// Step 1:
				const double u = random_generator.next_double();
				const double b_g_s = 1 + shape / FastMath.E;
				const double p = b_g_s * u;

				if (p <= 1)
				{
					// Step 2:

					const double x = std::pow(p, 1 / shape);
					const double u2 = random_generator.next_double();

					if (u2 > std::exp(-x))
					{
						// Reject
						continue;
					}
					else
					{
						return scale * x;
					}
				}
				else
				{
					// Step 3:

					const double x = -1 * std::log((b_g_s - p) / shape);
					const double u2 = random_generator.next_double();

					if (u2 > std::pow(x, shape - 1))
					{
						// Reject
						continue;
					}
					else
					{
						return scale * x;
					}
				}
			}
		}

		// Now shape >= 1

		const double d = shape - 0.333333333333333333;
		const double c = 1 / (3 * std::sqrt(d));

		while (true)
		{
			const double x = random_generator.next_gaussian();
			const double v = (1 + c * x) * (1 + c * x) * (1 + c * x);

			if (v <= 0)
			{
				continue;
			}

			const double x2 = x * x;
			const double u = random_generator.next_double();

			// Squeeze
			if (u < 1 - 0.0331 * x2 * x2)
			{
				return scale * d * v;
			}

			if (std::log(u) < 0.5 * x2 + d * (1 - v + std::log(v)))
			{
				return scale * d * v;
			}
		}
	}

	/**
	 * Returns the next normally-distributed pseudo-random deviate.
	 *
	 * @param mean mean of the normal distribution
	 * @param standard_deviation standard deviation of the normal distribution
	 * @return a random value, normally distributed with the given mean and standard deviation
	 */
	public double next_normal(double mean, double standard_deviation)
	{
		if (standard_deviation <= 0)
		{
			throw  (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, standard_deviation, 0);
		}
		return standard_deviation * next_gaussian() + mean;
	}

	/**
	 * Returns the next log-normally-distributed pseudo-random deviate.
	 *
	 * @param shape shape parameter of the log-normal distribution
	 * @param scale scale parameter of the log-normal distribution
	 * @return a random value, normally distributed with the given mean and standard deviation
	 */
	public double next_log_normal(double shape, double scale)
	{
		if (shape <= 0)
		{
			throw  (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, shape, 0);
		}
		return std::exp(scale + shape * next_gaussian());
	}

	/**
	 * Returns a poisson-distributed deviate with the given mean.
	 *
	 * @param mean expected value
	 * @return poisson deviate
	 * @ if mean is not strictly positive
	 */
	public int next_poisson(double mean)
	{
		if (mean <= 0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, mean, 0);
		}
		const double pivot = 40.0;
		if (mean < pivot)
		{
			double p = std::exp(-mean);
			long n = 0;
			double r = 1.0;
			double rnd;

			while (n < 1000 * mean)
			{
				rnd = random_generator.next_double();
				r *= rnd;
				if (r >= p)
				{
					n++;
				}
				else
				{
					return static_cast<int>(std::min(n, std::numeric_limits<int>::max());
				}
			}
			return static_cast<int>(std::min(n, std::numeric_limits<int>::max());
		}
		else
		{
			const double lambda = std::floor(mean);
			const double lambda_fractional = mean - lambda;
			const double log_lambda = std::log(lambda);
			const double log_lambda_factorial = Combinatorics_Utils.factorial_log(static_cast<int>(lambda);
			const long y2 = lambda_fractional < Double.MIN_VALUE ? 0 : next_poisson(lambda_fractional);
			const double delta = std::sqrt(lambda * std::log(32 * lambda / std::numbers::pi + 1));
			const double half_delta = delta / 2;
			const double twolpd = 2 * lambda + delta;
			const double& a1 = std::sqrt(std::numbers::pi * twolpd) * std::exp(1 / (8 * lambda));
			const double& a2 = (twolpd / delta) * std::exp(-delta * (1 + delta) / twolpd);
			const double& a_sum = a1 + a2 + 1;
			const double p1 = a1 / a_sum;
			const double p2 = a2 / a_sum;
			const double c1 = 1 / (8 * lambda);

			double x;
			double y = 0;
			double v;
			int a;
			double t;
			double qr;
			double qa;
			for (;;)
			{
				const double u = random_generator.next_double();
				if (u <= p1)
				{
					const double n = random_generator.next_gaussian();
					x = n * std::sqrt(lambda + half_delta) - 0.5d;
					if (x > delta || x < -lambda)
					{
						continue;
					}
					y = x < 0 ? std::floor(x) : std::ceil(x);
					const double e = next_exponential(1);
					v = -e - (n * n / 2) + c1;
				}
				else
				{
					if (u > p1 + p2)
					{
						y = lambda;
						break;
					}
					else
					{
						x = delta + (twolpd / delta) * next_exponential(1);
						y = std::ceil(x);
						v = -next_exponential(1) - delta * (x + 1) / twolpd;
					}
				}
				a = x < 0 ? 1 : 0;
				t = y * (y + 1) / (2 * lambda);
				if (v < -t && a == 0)
				{
					y = lambda + y;
					break;
				}
				qr = t * ((2 * y + 1) / (6 * lambda) - 1);
				qa = qr - (t * t) / (3 * (lambda + a * (y + 1)));
				if (v < qa)
				{
					y = lambda + y;
					break;
				}
				if (v > qr)
				{
					continue;
				}
				if (v < y * log_lambda - Combinatorics_Utils.factorial_log(static_cast<int>((y + lambda)) + log_lambda_factorial)
				{
					y = lambda + y;
						break;
				}
			}
			return static_cast<int>(std::min(y2 + static_cast<long>(y, std::numeric_limits<int>::max());
		}
	}

	/**
	 * Returns a random deviate from the given distribution.
	 *
	 * @param dist the distribution to sample from
	 * @return a random value following the given distribution
	 */
	public double next_deviate(Real_Distribution dist)
	{
		return get_sampler(dist).next_sample(this, dist);
	}

	/**
	 * Returns an array of random deviates from the given distribution.
	 *
	 * @param dist the distribution to sample from
	 * @param size the number of values to return
	 *
	 * @return an array of {@code size} values following the given distribution
	 */
	public std::vector<double> next_deviates(Real_Distribution dist, int size)
	{
		//TODO: check parameters

		Real_Distribution_Sampler sampler = get_sampler(dist);
		std::vector<double> out = std::vector<double>(size];
		for (int i{}; i < size; i++)
		{
			out[i] = sampler.next_sample(this, dist);
		}
		return out;
	}

	/**
	 * Returns a random deviate from the given distribution.
	 *
	 * @param dist the distribution to sample from
	 * @return a random value following the given distribution
	 */
	public int next_deviate(Integer_Distribution dist)
	{
		return get_sampler(dist).next_sample(this, dist);
	}

	/**
	 * Returns an array of random deviates from the given distribution.
	 *
	 * @param dist the distribution to sample from
	 * @param size the number of values to return
	 *
	 * @return an array of {@code size }values following the given distribution
	 */
	public std::vector<int> next_deviates(Integer_Distribution dist, int size)
	{
		//TODO: check parameters

		Integer_DistributionSampler sampler = get_sampler(dist);
		std::vector<int> out = int[size];
		for (int i{}; i < size; i++)
		{
			out[i] = sampler.next_sample(this, dist);
		}
		return out;
	}

	/**
	 * Returns a sampler for the given continuous distribution.
	 * @param dist the distribution
	 * @return a sampler for the distribution
	 */
	private Real_Distribution_Sampler get_sampler(Real_Distribution dist)
	{
		Real_Distribution_Sampler sampler = CONTINUOUS_SAMPLERS.get(dist.get_class());
		if (sampler != NULL)
		{
			return sampler;
		}
		return DEFAULT_REAL_SAMPLER;
	}

	/**
	 * Returns a sampler for the given discrete distribution.
	 * @param dist the distribution
	 * @return a sampler for the distribution
	 */
	private Integer_DistributionSampler get_sampler(Integer_Distribution dist)
	{
		Integer_DistributionSampler sampler = DISCRETE_SAMPLERS.get(dist.get_class());
		if (sampler != NULL)
		{
			return sampler;
		}
		return DEFAULT_INTEGER_SAMPLER;
	}

	/**
	 * Returns a uniformly distributed random integer between lower and upper (inclusive).
	 *
	 * @param lower lower bound for the generated value
	 * @param upper upper bound for the generated value
	 * @return a random integer value within the given bounds
	 * @ if lower is not strictly less than or equal to upper
	 */
	public int next_int(const int& lower, int upper)
	{
		if (lower >= upper)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::LOWER_BOUND_NOT_BELOW_UPPER_BOUND, lower, upper);
		}
		const int max = (upper - lower) + 1;
		if (max <= 0)
		{
			// The range is too wide to fit in a positive int (larger
			// than 2^31); as it covers more than half the integer range, // we use a simple rejection method.
			while (true)
			{
				const int r = next_int();
				if (r >= lower &&
					r <= upper)
				{
					return r;
				}
			}
		}
		else
		{
			// We can shift the range and directly generate a positive int.
			return lower + next_int(max);
		}
	}

	/**
	 * Returns a uniformly distributed random long integer between lower and upper (inclusive).
	 *
	 * @param lower lower bound for the generated value
	 * @param upper upper bound for the generated value
	 * @return a random long integer value within the given bounds
	 * @ if lower is not strictly less than or equal to upper
	 */
	public long next_long(const long lower, const long upper)
	{
		if (lower >= upper)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::LOWER_BOUND_NOT_BELOW_UPPER_BOUND, lower, upper);
		}
		const long max = (upper - lower) + 1;
		if (max <= 0)
		{
			// the range is too wide to fit in a positive long (larger than 2^63); as it covers
			// more than half the long range, we use directly a simple rejection method
			while (true)
			{
				const long r = random_generator.next_long();
				if (r >= lower && r <= upper)
				{
					return r;
				}
			}
		}
		else if (max < std::numeric_limits<int>::max())
		{
			// we can shift the range and generate directly a positive int
			return lower + random_generator.next_int(static_cast<int>(max);
		}
		else
		{
			// we can shift the range and generate directly a positive long
			return lower + next_long(max);
		}
	}

	/**
	 * Returns a double value uniformly distributed over [lower, upper]
	 * @param lower lower bound
	 * @param upper upper bound
	 * @return uniform deviate
	 * @ if upper is less than or equal to upper
	 */
	public double next_uniform(double lower, double upper)
	{
		if (upper <= lower)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::LOWER_BOUND_NOT_BELOW_UPPER_BOUND, lower, upper);
		}
		if (std::isinf(lower) || std::isinf(upper))
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::INFINITE_BOUND);
		}
		if (std::isnan(lower) || std::isnan(upper))
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NAN_NOT_ALLOWED);
		}
		const double u = random_generator.next_double();
		return u * upper + (1 - u) * lower;
	}

	/**
	 * Returns an integer value following a Zipf distribution with the given parameter.
	 *
	 * @param number_of_elements number of elements of the distribution
	 * @param exponent exponent of the distribution
	 * @return random Zipf value
	 */
	public int next_zipf(const int& number_of_elements, double exponent)
	{
		if (zipf_sampler == NULL || zipf_sampler.get_exponent() != exponent || zipf_sampler.get_number_of_elements() != number_of_elements)
		{
			zipf_sampler = Zipf_Rejection_Inversion_Sampler(number_of_elements, exponent);
		}
		return zipf_sampler.sample(random_generator);
	}

	/**
	 * Generates a random string of hex characters of length {@code len}.
	 * <p>
	 * The generated string will be random, but not cryptographically secure.
	 * <p>
	 * <strong>Algorithm Description:</strong> hex strings are generated using a
	 * 2-step process.
	 * <ol>
	 * <li>{@code len / 2 + 1} binary bytes are generated using the underlying
	 * Random</li>
	 * <li>Each binary std::byte is translated into 2 hex digits</li>
	 * </ol>
	 * </p>
	 *
	 * @param len the desired string length.
	 * @return the random string.
	 * @ if {@code len <= 0}.
	 */
	public std::string next_hex_string(const int& len)
	{
		if (len <= 0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::size(), len);
		}

		// Initialize output buffer
		std::stringBuilder out_buffer = std::stringstream();

		// Get int(len/2)+1 random bytes
		std::vector<std::byte>random_bytes = byte[(len / 2) + 1];
		random_generator.next_bytes(random_bytes);

		// Convert each std::byte to 2 hex digits
		for (int i{}; i < random_bytes.size(); i++)
		{
			Integer c = static_cast<int>(random_bytes[i]);

			/*
			 * Add 128 to std::byte value to make interval 0-255 before doing hex
			 * conversion. This guarantees <= 2 hex digits from to_hex_string()
			 * to_hex_string would otherwise add 2^32 to negative arguments.
			 */
			std::string hex = Integer.to_hex_string(c.int_value() + 128);

			// Make sure we add 2 hex digits for each byte
			if (hex.size()() == 1)
			{
				out_buffer.append('0');
			}
			out_buffer.append(hex);
		}
		return out_buffer.to_string().substring(0, len);
	}

	/**
	 * Generates an integer array of length {@code k} whose entries are selected
	 * randomly, without repetition, from the integers {@code 0, ..., n - 1}
	 * (inclusive).
	 * <p>
	 * Generated arrays represent permutations of {@code n} taken {@code k} at a
	 * time.</p>
	 * This method calls {@link Math_Arrays#shuffle(std::vector<int>,Random_Generator)
	 * Math_Arrays::shuffle} in order to create a random shuffle of the set
	 * of natural numbers {@code { 0, 1, ..., n - 1 }}.
	 *
	 *
	 * @param n the domain of the permutation
	 * @param k the size of the permutation
	 * @return a random {@code k}-permutation of {@code n}, as an array of
	 * integers
	 * @ if {@code k > n}.
	 * @ if {@code k <= 0}.
	 */
	public std::vector<int> next_permutation(const int& n, const int& k)

	{
		if (k > n)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::PERMUTATION_EXCEEDS_N, k, n, true);
		}
		if (k <= 0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::PERMUTATION_SIZE, k);
		}

		const std::vector<int> index = Math_Arrays::natural(n);
		Math_Arrays::shuffle(index, random_generator);

		// Return a array containing the first "k" entries of "index".
		return Arrays.copy_of(index, k);
	}

	/**
	 * Returns an array of {@code k} objects selected randomly from the
	 * Collection {@code c}.
	 * <p>
	 * Sampling from {@code c} is without replacement; but if {@code c} contains
	 * identical objects, the sample may include repeats.  If all elements of
	 * {@code c} are distinct, the resulting object array represents a
	 * <a href="http://rkb.home.cern.ch/rkb/AN16pp/node250.html#SECTION0002500000000000000000">
	 * Simple Random Sample</a> of size {@code k} from the elements of
	 * {@code c}.</p>
	 * <p>This method calls {@link #next_permutation(int,int) next_permutation(c.size(), k)}
	 * in order to sample the collection.
	 * </p>
	 *
	 * @param c the collection to be sampled
	 * @param k the size of the sample
	 * @return a random sample of {@code k} elements from {@code c}
	 * @ if {@code k > c.size()}.
	 * @ if {@code k <= 0}.
	 */
	public std::vector<Object> next_sample(Collection< ? > c, const int& k)
	{
		int len = c.size();
		if (k > len)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::SAMPLE_SIZE_EXCEEDS_COLLECTION_SIZE, k, len, true);
		}
		if (k <= 0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_SAMPLES, k);
		}

		std::vector<Object> objects = c.to_array();
		std::vector<int> index = next_permutation(len, k);
		std::vector<Object> result = Object[k];
		for (int i{}; i < k; i++)
		{
			result[i] = objects[index[i]];
		}
		return result;
	}

	/**
	 * Returns an array of {@code k} double values selected randomly from the
	 * double array {@code a}.
	 * <p>
	 * Sampling from {@code a} is without replacement; but if {@code a} contains
	 * identical objects, the sample may include repeats.  If all elements of
	 * {@code a} are distinct, the resulting object array represents a
	 * <a href="http://rkb.home.cern.ch/rkb/AN16pp/node250.html#SECTION0002500000000000000000">
	 * Simple Random Sample</a> of size {@code k} from the elements of
	 * {@code a}.</p>
	 *
	 * @param a the array to be sampled
	 * @param k the size of the sample
	 * @return a random sample of {@code k} elements from {@code a}
	 * @ if {@code k > c.size()}.
	 * @ if {@code k <= 0}.
	 */
	public std::vector<double> next_sample(std::vector<double> a, const int& k)
	{
		int len = a.size();
		if (k > len)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::SAMPLE_SIZE_EXCEEDS_COLLECTION_SIZE, k, len, true);
		}
		if (k <= 0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_OF_SAMPLES, k);
		}
		std::vector<int> index = next_permutation(len, k);
		std::vector<double> result = std::vector<double>(k];
		for (int i{}; i < k; i++)
		{
			result[i] = a[index[i]];
		}
		return result;
	}

	/**
	 * Generates a random sample of size sample_size from {0, 1, ... , weights.size() - 1}, * using weights as probabilities.
	 * <p>
	 * For 0 &lt; i &lt; weights.size(), the probability that i is selected (on any draw) is weights[i].
	 * If necessary, the weights array is normalized to sum to 1 so that weights[i] is a probability
	 * and the array sums to 1.
	 * <p>
	 * Weights can be 0, but must not be negative, infinite or NaN.
	 * At least one weight must be positive.
	 *
	 * @param sample_size size of sample to generate
	 * @param weights probability sampling weights
	 * @return an array of integers between 0 and weights.size() - 1
	 * @ if weights contains negative, NaN or infinite values or only 0s or sample_size is less than 0
	 */
	public std::vector<int> next_sample_with_replacement(const int& sample_size, std::vector<double> weights)
	{
		// Check sample size
		if (sample_size < 0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_NUMBER_OF_SAMPLES);
		}

		// Check and normalize weights
		std::vector<double> norm_wt = Enumerated_Distribution.check_and_normalize(weights);

		// Generate sample values by dividing [0,1] into subintervals corresponding to weights.
		const std::vector<int> out = int[sample_size];
		const int len = norm_wt.size();
		for (int i{}; i < sample_size; i++)
		{
			const double u = random_generator.next_double();
			double cum = norm_wt[0];
			int j = 1;
			while (cum < u && j < len)
			{
				cum += norm_wt[j++];
			}
			out[i] = --j;
		}
		return out;
	}

	/**
	 * Utility class implementing Cheng's algorithms for beta distribution sampling.
	 *
	 * <blockquote>
	 * <pre>
	 * R. C. H. Cheng, * "Generating beta variates with nonintegral shape parameters", * Communications of the ACM, 21, 317-322, 1978.
	 * </pre>
	 * </blockquote>
	 */
	private static class Cheng_Beta_Sampler
	{
		/** Private constructor for utility class. */
		private Cheng_Beta_Sampler() { // NOPMD - PMD fails to detect this is a utility class
			// not called
		}

		/**
		 * Returns the next sample following a beta distribution
		 * with given alpha and beta parameters.
		 *
		 * @param generator the random generator to use
		 * @param alpha the alpha parameter
		 * @param beta the beta parameter
		 * @return the next sample
		 */
		public static double sample(Random_Generator generator, double alpha, double beta)
		{
			// TODO: validate parameters
			const double& a = std::min(alpha, beta);
			const double b = std::max(alpha, beta);

			if (a > 1)
			{
				return algorithm_b_b(generator, alpha, a, b);
			}
			else
			{
				return algorithm_b_c(generator, alpha, b, a);
			}
		}

		/**
		 * Returns one Beta sample using Cheng's BB algorithm, * when both &alpha; and &beta; are greater than 1.
		 *
		 * @param generator the random generator to use
		 * @param a0 distribution first shape parameter (&alpha;)
		 * @param a min(&alpha;, &beta;) where &alpha;, &beta; are the two distribution shape parameters
		 * @param b max(&alpha;, &beta;) where &alpha;, &beta; are the two distribution shape parameters
		 * @return sampled value
		 */
		private static double algorithm_b_b(const Random_Generator generator, const double& a0, const double& a, const double& b)
		{
			const double& alpha = a + b;
			const double beta = std::sqrt((alpha - 2.) / (2. * a * b - alpha));
			const double gamma = a + 1. / beta;

			double r;
			double w;
			double t;
			do
			{
				const double u1 = generator.next_double();
				const double u2 = generator.next_double();
				const double v = beta * (std::log(u1) - std::log1p(-u1));
				w = a * std::exp(v);
				const double z = u1 * u1 * u2;
				r = gamma * v - 1.3862944;
				const double s = a + r - w;
				if (s + 2.609438 >= 5 * z)
				{
					break;
				}

				t = std::log(z);
				if (s >= t)
				{
					break;
				}
			} while (r + alpha * (std::log(alpha) - std::log(b + w)) < t);

			w = std::min(w, Double.MAX_VALUE);
			return Precision::equals(a, a0) ? w / (b + w) : b / (b + w);
		}

		/**
		 * Returns a Beta-distribute value using Cheng's BC algorithm, * when at least one of &alpha; and &beta; is smaller than 1.
		 *
		 * @param generator the random generator to use
		 * @param a0 distribution first shape parameter (&alpha;)
		 * @param a max(&alpha;, &beta;) where &alpha;, &beta; are the two distribution shape parameters
		 * @param b min(&alpha;, &beta;) where &alpha;, &beta; are the two distribution shape parameters
		 * @return sampled value
		 */
		private static double algorithm_b_c(const Random_Generator generator, const double& a0, const double& a, const double& b)
		{
			const double& alpha = a + b;
			const double beta = 1. / b;
			const double delta = 1. + a - b;
			const double k1 = delta * (0.0138889 + 0.0416667 * b) / (a * beta - 0.777778);
			const double k2 = 0.25 + (0.5 + 0.25 / delta) * b;

			double w;
			for (;;)
			{
				const double u1 = generator.next_double();
				const double u2 = generator.next_double();
				const double y = u1 * u2;
				const double z = u1 * y;
				if (u1 < 0.5)
				{
					if (0.25 * u2 + z - y >= k1)
					{
						continue;
					}
				}
				else
				{
					if (z <= 0.25)
					{
						const double v = beta * (std::log(u1) - std::log1p(-u1));
						w = a * std::exp(v);
						break;
					}

					if (z >= k2)
					{
						continue;
					}
				}

				const double v = beta * (std::log(u1) - std::log1p(-u1));
				w = a * std::exp(v);
				if (alpha * (std::log(alpha) - std::log(b + w) + v) - 1.3862944 >= std::log(z))
				{
					break;
				}
			}

			w = std::min(w, Double.MAX_VALUE);
			return Precision::equals(a, a0) ? w / (b + w) : b / (b + w);
		}
	}

	/**
	 * Utility class implementing a rejection inversion sampling method for a discrete, * bounded Zipf distribution that is based on the method described in
	 * <p>
	 * Wolfgang H??rmann and Gerhard Derflinger
	 * "Rejection-inversion to generate variates from monotone discrete distributions."
	 * ACM Transactions on Modeling and Computer Simulation (TOMACS) 6.3 (1996): 169-184.
	 * <p>
	 * The paper describes an algorithm for exponents larger than 1 (Algorithm ZRI).
	 * The original method uses {@code H(x) := (v + x)^(1 - q) / (1 - q)}
	 * as the integral of the hat function. This function is undefined for
	 * q = 1, which is the reason for the limitation of the exponent.
	 * If instead the integral function
	 * {@code H(x) := ((v + x)^(1 - q) - 1) / (1 - q)} is used, * for which a meaningful limit exists for q = 1, * the method works for all positive exponents.
	 * <p>
	 * The following implementation uses v := 0 and generates integral numbers
	 * in the range [1, number_of_elements]. This is different to the original method
	 * where v is defined to be positive and numbers are taken from [0, i_max].
	 * This explains why the implementation looks slightly different.
	 *
	 */
	static const class Zipf_Rejection_Inversion_Sampler
	{
		/** Exponent parameter of the distribution. */
		private const double exponent;
		/** Number of elements. */
		private const int& number_of_elements;
		/** Constant equal to {@code h_integral(1.5) - 1}. */
		private const double h_integral_x1;
		/** Constant equal to {@code h_integral(number_of_elements + 0.5)}. */
		private const double h_integral_number_of_elements;
		/** Constant equal to {@code 2 - h_integral_inverse(h_integral(2.5) - h(2)}. */
		private const double s;

		/** Simple constructor.
		 * @param number_of_elements number of elements
		 * @param exponent exponent parameter of the distribution
		 */
		Zipf_Rejection_Inversion_Sampler(const int& number_of_elements, const double exponent)
		{
			this.exponent = exponent;
			this.number_of_elements = number_of_elements;
			this.h_integral_x1 = h_integral(1.5) - 1d;
			this.h_integral_number_of_elements = h_integral(number_of_elements + 0.5);
			this.s = 2d - h_integral_inverse(h_integral(2.5) - h(2));
		}

		/** Generate one integral number in the range [1, number_of_elements].
		 * @param random random generator to use
		 * @return generated integral number in the range [1, number_of_elements]
		 */
		int sample(const Random_Generator random)
		{
			while (true)
			{
				const double u = h_integral_number_of_elements + random.next_double() * (h_integral_x1 - h_integral_number_of_elements);
				// u is uniformly distributed in (h_integral_x1, h_integral_number_of_elements]

				double x = h_integral_inverse(u);

				int k = static_cast<int>((x + 0.5);

				// Limit k to the range [1, number_of_elements]
				// (k could be outside due to numerical inaccuracies)
				if (k < 1)
				{
					k = 1;
				}
				else if (k > number_of_elements)
				{
					k = number_of_elements;
				}

				// Here, the distribution of k is given by:
				//
				//   P(k = 1) = C * (h_integral(1.5) - h_integral_x1) = C
				//   P(k = m) = C * (h_integral(m + 1/2) - h_integral(m - 1/2)) for m >= 2
				//
				//   where C := 1 / (h_integral_number_of_elements - h_integral_x1)

				if (k - x <= s || u >= h_integral(k + 0.5) - h(k))
				{
					// Case k = 1:
					//
					//   The right inequality is always true, because replacing k by 1 gives
					//   u >= h_integral(1.5) - h(1) = h_integral_x1 and u is taken from
					//   (h_integral_x1, h_integral_number_of_elements].
					//
					//   Therefore, the acceptance rate for k = 1 is P(accepted | k = 1) = 1
					//   and the probability that 1 is returned as random value is
					//   P(k = 1 and accepted) = P(accepted | k = 1) * P(k = 1) = C = C / 1^exponent
					//
					// Case k >= 2:
					//
					//   The left inequality (k - x <= s) is just a short cut
					//   to avoid the more expensive evaluation of the right inequality
					//   (u >= h_integral(k + 0.5) - h(k)) in many cases.
					//
					//   If the left inequality is true, the right inequality is also true:
					//     Theorem 2 in the paper is valid for all positive exponents, because
					//     the requirements h'(x) = -exponent/x^(exponent + 1) < 0 and
					//     (-1/h_inverse'(x))'' = (1+1/exponent) * x^(1/exponent-1) >= 0
					//     are both fulfilled.
					//     Therefore, f(x) := x - h_integral_inverse(h_integral(x + 0.5) - h(x))
					//     is a non-decreasing function. If k - x <= s holds, //     k - x <= s + f(k) - f(2) is obviously also true which is equivalent to
					//     -x <= -h_integral_inverse(h_integral(k + 0.5) - h(k)), //     -h_integral_inverse(u) <= -h_integral_inverse(h_integral(k + 0.5) - h(k)), //     and constly u >= h_integral(k + 0.5) - h(k).
					//
					//   Hence, the right inequality determines the acceptance rate:
					//   P(accepted | k = m) = h(m) / (h_integrated(m+1/2) - h_integrated(m-1/2))
					//   The probability that m is returned is given by
					//   P(k = m and accepted) = P(accepted | k = m) * P(k = m) = C * h(m) = C / m^exponent.
					//
					// In both cases the probabilities are proportional to the probability mass function
					// of the Zipf distribution.

					return k;
				}
			}
		}

		/**
		 * {@code H(x) :=}
		 * <ul>
		 * <li>{@code (x^(1-exponent) - 1)/(1 - exponent)}, if {@code exponent != 1}</li>
		 * <li>{@code log(x)}, if {@code exponent == 1}</li>
		 * </ul>
		 * H(x) is an integral function of h(x), * the derivative of H(x) is h(x).
		 *
		 * @param x free parameter
		 * @return {@code H(x)}
		 */
		private double h_integral(const double& x)
		{
			const double log_x = std::log(x);
			return helper2((1d - exponent) * log_x) * log_x;
		}

		/**
		 * {@code h(x) := 1/x^exponent}
		 *
		 * @param x free parameter
		 * @return h(x)
		 */
		private double h(const double& x)
		{
			return std::exp(-exponent * std::log(x));
		}

		/**
		 * The inverse function of H(x).
		 *
		 * @param x free parameter
		 * @return y for which {@code H(y) = x}
		 */
		private double h_integral_inverse(const double& x)
		{
			double t = x * (1d - exponent);
			if (t < -1d)
			{
				// Limit value to the range [-1, +inf).
				// t could be smaller than -1 in some rare cases due to numerical errors.
				t = -1;
			}
			return std::exp(helper1(t) * x);
		}

		/**
		 * @return the exponent of the distribution being sampled
		 */
		public double get_exponent()
		{
			return exponent;
		}

		/**
		 * @return the number of elements of the distribution being sampled
		 */
		public int get_number_of_elements()
		{
			return number_of_elements;
		}

		/**
		 * Helper function that calculates {@code log(1+x)/x}.
		 * <p>
		 * A Taylor series expansion is used, if x is close to 0.
		 *
		 * @param x a value larger than or equal to -1
		 * @return {@code log(1+x)/x}
		 */
		static double helper1(const double& x)
		{
			if (std::abs(x) > 1e-8)
			{
				return std::log1p(x) / x;
			}
			else
			{
				return 1. - x * ((1. / 2.) - x * ((1. / 3.) - x * (1. / 4.)));
			}
		}

		/**
		 * Helper function to calculate {@code (exp(x)-1)/x}.
		 * <p>
		 * A Taylor series expansion is used, if x is close to 0.
		 *
		 * @param x free parameter
		 * @return {@code (exp(x)-1)/x} if x is non-zero, or 1 if x=0
		 */
		static double helper2(const double& x)
		{
			if (std::abs(x) > 1e-8)
			{
				return std::expm1(x) / x;
			}
			else
			{
				return 1. + x * (1. / 2.) * (1. + x * (1. / 3.) * (1. + x * (1. / 4.)));
			}
		}
	}

	/**
	 * Sampler for enumerated distributions.
	 *
	 * @param <T> type of sample space objects
	 */
	private const class Enumerated_DistributionSampler<T>
	{
		/** Probabilities */
		private const std::vector<double> weights;
		/** Values */
		private const List<T> values;
		/**
		 * Create an Enumerated_DistributionSampler from the provided pmf.
		 *
		 * @param pmf probability mass function describing the distribution
		 */
		Enumerated_DistributionSampler(List<Pair<T, Double>> pmf)
		{
			const int& num_masses = pmf.size();
			weights = std::vector<double>(num_masses];
			values = Array_list<>();
			for (int i{}; i < num_masses; i++)
			{
				weights[i] = pmf.get(i).get_second();
				values.add(pmf.get(i).get_first());
			}
		}
		/**
		 * @return a random value from the distribution
		 */
		public T sample()
		{
			std::vector<int> chosen = next_sample_with_replacement(1, weights);
			return values.get(chosen[0]);
		}
	}
};