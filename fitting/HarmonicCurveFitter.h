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
  //package org.hipparchus.fitting;

  //import java.util.Array_list;
  //import java.util.Collection;
  //import java.util.List;

  //import org.hipparchus.analysis.function.Harmonic_Oscillator;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.linear.Diagonal_Matrix;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Builder;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Sin_Cos;
#include <vector>
#include <numbers>
  /**
   * Fits points to a {@link
   * org.hipparchus.analysis.function.Harmonic_Oscillator.Parametric harmonic oscillator}
   * function.
   * <br/>
   * The {@link #with_start_point(std::vector<double>) initial guess values} must be passed
   * in the following order:
   * <ul>
   *  <li>Amplitude</li>
   *  <li>Angular frequency</li>
   *  <li>phase</li>
   * </ul>
   * The optimal values will be returned in the same order.
   *
   */
class Harmonic_Curve_Fitter : public Abstract_Curve_Fitter
{
	/** Parametric function to be fitted. */
	private static const Harmonic_Oscillator.Parametric FUNCTION = Harmonic_Oscillator.Parametric();
	/** Initial guess. */
	private const std::vector<double> initial_guess;
	/** Maximum number of iterations of the optimization algorithm. */
	private const int max_iter;

	/**
	 * Constructor used by the factory methods.
	 *
	 * @param initial_guess Initial guess. If set to {@code NULL}, the initial guess
	 * will be estimated using the {@link Parameter_Guesser}.
	 * @param max_iter Maximum number of iterations of the optimization algorithm.
	 */
	private Harmonic_Curve_Fitter(std::vector<double> initial_guess, int max_iter)
	{
		this.initial_guess = initial_guess == NULL ? NULL : initial_guess.clone();
		this.max_iter = max_iter;
	}

	/**
	 * Creates a default curve fitter.
	 * The initial guess for the parameters will be {@link Parameter_Guesser}
	 * computed automatically, and the maximum number of iterations of the
	 * optimization algorithm is set to {@link Integer#MAX_VALUE}.
	 *
	 * @return a curve fitter.
	 *
	 * @see #with_start_point(std::vector<double>)
	 * @see #with_max_iterationsstatic_cast<int>(
	 */
	public static Harmonic_Curve_Fitter create()
	{
		return Harmonic_Curve_Fitter(null, std::numeric_limits<int>::max());
	}

	/**
	 * Configure the start point (initial guess).
	 * @param new_start start point (initial guess)
	 * @return a instance.
	 */
	public Harmonic_Curve_Fitter with_start_point(std::vector<double> new_start)
	{
		return Harmonic_Curve_Fitter(new_start.clone(), max_iter);
	}

	/**
	 * Configure the maximum number of iterations.
	 * @param new_max_iter maximum number of iterations
	 * @return a instance.
	 */
	public Harmonic_Curve_Fitter with_max_iterations(const int& new_max_iter)
	{
		return Harmonic_Curve_Fitter(initial_guess, new_max_iter);
	}

	/** {@inherit_doc} */
 //override
	protected Least_Squares_Problem get_problem(Collection<Weighted_Observed_Point> observations)
	{
		// Prepare least-squares problem.
		const int len = observations.size();
		const std::vector<double> target = std::vector<double>(len];
		const std::vector<double> weights = std::vector<double>(len];

		int i = 0;
		for (Weighted_Observed_Point obs : observations)
		{
			target[i] = obs.get_y();
			weights[i] = obs.get_weight();
			++i;
		}

		const Abstract_Curve_Fitter.Theoretical_Values_Function model
			= Abstract_Curve_Fitter.Theoretical_Values_Function(FUNCTION, observations);

		const std::vector<double> start_point = initial_guess != NULL ?
			initial_guess :
			// Compute estimation.
			Parameter_Guesser(observations).guess();

		// Return a optimizer set up to fit a Gaussian curve to the
		// observed points.
		return Least_Squares_Builder().
			max_evaluations(std::numeric_limits<int>::max()).
			max_iterations(max_iter).
			start(start_point).
			target(target).
			weight(new Diagonal_Matrix(weights)).
			model(model.get_model_function(), model.get_model_function_jacobian()).
			build();
	}

	/**
	 * This class guesses harmonic coefficients from a sample.
	 * <p>The algorithm used to guess the coefficients is as follows:</p>
	 *
	 * <p>We know \( f(t) \) at some sampling points \( t_i \) and want
	 * to find \( a \), \( \omega \) and \( \phi \) such that
	 * \( f(t) = a \cos (\omega t + \phi) \).
	 * </p>
	 *
	 * <p>From the analytical expression, we can compute two primitives :
	 * \[
	 *     If2(t) = \int f^2 dt  = a^2 (t + S(t)) / 2
	 * \]
	 * \[
	 *     If'2(t) = \int f'^2 dt = a^2 \omega^2 (t - S(t)) / 2
	 * \]
	 * where \(S(t) = \frac{\sin(2 (\omega t + \phi))}{2\omega}\)
	 * </p>
	 *
	 * <p>We can remove \(S\) between these expressions :
	 * \[
	 *     If'2(t) = a^2 \omega^2 t - \omega^2 If2(t)
	 * \]
	 * </p>
	 *
	 * <p>The preceding expression shows that \(If'2 (t)\) is a linear
	 * combination of both \(t\) and \(If2(t)\):
	 * \[
	 *   If'2(t) = A t + B If2(t)
	 * \]
	 * </p>
	 *
	 * <p>From the primitive, we can deduce the same form for definite
	 * integrals between \(t_1\) and \(t_i\) for each \(t_i\) :
	 * \[
	 *   If2(t_i) - If2(t_1) = A (t_i - t_1) + B (If2 (t_i) - If2(t_1))
	 * \]
	 * </p>
	 *
	 * <p>We can find the coefficients \(A\) and \(B\) that best fit the sample
	 * to this linear expression by computing the definite integrals for
	 * each sample points.
	 * </p>
	 *
	 * <p>For a bilinear expression \(z(x_i, y_i) = A x_i + B y_i\), the
	 * coefficients \(A\) and \(B\) that minimize a least-squares criterion
	 * \(\sum (z_i - z(x_i, y_i))^2\) are given by these expressions:</p>
	 * \[
	 *   A = \frac{\sum y_i y_i \sum x_i z_i - \sum x_i y_i \sum y_i z_i}
	 *            {\sum x_i x_i \sum y_i y_i - \sum x_i y_i \sum x_i y_i}
	 * \]
	 * \[
	 *   B = \frac{\sum x_i x_i \sum y_i z_i - \sum x_i y_i \sum x_i z_i}
	 *            {\sum x_i x_i \sum y_i y_i - \sum x_i y_i \sum x_i y_i}
	 *
	 * \]
	 *
	 * <p>In fact, we can assume that both \(a\) and \(\omega\) are positive and
	 * compute them directly, knowing that \(A = a^2 \omega^2\) and that
	 * \(B = -\omega^2\). The complete algorithm is therefore:</p>
	 *
	 * For each \(t_i\) from \(t_1\) to \(t_{n-1}\), compute:
	 * \[ f(t_i) \]
	 * \[ f'(t_i) = \frac{f (t_{i+1}) - f(t_{i-1})}{t_{i+1} - t_{i-1}} \]
	 * \[ x_i = t_i  - t_1 \]
	 * \[ y_i = \int_{t_1}^{t_i} f^2(t) dt \]
	 * \[ z_i = \int_{t_1}^{t_i} f'^2(t) dt \]
	 * and update the sums:
	 * \[ \sum x_i x_i, \sum y_i y_i, \sum x_i y_i, \sum x_i z_i, \sum y_i z_i \]
	 *
	 * Then:
	 * \[
	 *  a = \sqrt{\frac{\sum y_i y_i  \sum x_i z_i - \sum x_i y_i \sum y_i z_i }
	 *                 {\sum x_i y_i  \sum x_i z_i - \sum x_i x_i \sum y_i z_i }}
	 * \]
	 * \[
	 *  \omega = \sqrt{\frac{\sum x_i y_i \sum x_i z_i - \sum x_i x_i \sum y_i z_i}
	 *                      {\sum x_i x_i \sum y_i y_i - \sum x_i y_i \sum x_i y_i}}
	 * \]
	 *
	 * <p>Once we know \(\omega\) we can compute:
	 * \[
	 *    fc = \omega f(t) \cos(\omega t) - f'(t) \sin(\omega t)
	 * \]
	 * \[
	 *    fs = \omega f(t) \sin(\omega t) + f'(t) \cos(\omega t)
	 * \]
	 * </p>
	 *
	 * <p>It appears that \(fc = a \omega \cos(\phi)\) and
	 * \(fs = -a \omega \sin(\phi)\), so we can use these
	 * expressions to compute \(\phi\). The best estimate over the sample is
	 * given by averaging these expressions.
	 * </p>
	 *
	 * <p>sin_ce integrals and means are involved in the preceding
	 * estimations, these operations run in \(O(n)\) time, where \(n\) is the
	 * number of measurements.</p>
	 */
	public static class Parameter_Guesser
	{
	private:
		/** Amplitude. */
		double my_a;
		/** Angular frequency. */
		double my_omega;
		/** Phase. */
		double my_phi;

		/**
		 * Simple constructor.
		 *
		 * @param observations Sampled observations.
		 * @ if the sample is too short.
		 * @ if the abscissa range is zero.
		 * @Math_Illegal_State_Exception when the guessing procedure cannot
		 * produce sensible results.
		 */
		public Parameter_Guesser(Collection<Weighted_Observed_Point>& observations)
		{
			if (observations.size() < 4)
			{
				throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_OBSERVED_POINTS_IN_SAMPLE, observations.size(), 4, true);
			}

			const Weighted_Observed_Point sorted = sort_observations(observations).to_array(Weighted_Observed_Point[0]);

			auto a_omega = guess_a_omega(sorted);
			my_a = a_omega[0];
			my_omega = a_omega[1];

			my_phi = guess_phi(sorted);
		}

		/**
		 * Gets an estimation of the parameters.
		 *
		 * @return the guessed parameters, in the following order:
		 * <ul>
		 *  <li>Amplitude</li>
		 *  <li>Angular frequency</li>
		 *  <li>Phase</li>
		 * </ul>
		 */
		public std::vector<double> guess()
		{
			return std::vector<double> { my_a, my_omega, my_phi };
		}

		/**
		 * Sort the observations with respect to the abscissa.
		 *
		 * @param unsorted Input observations.
		 * @return the input observations, sorted.
		 */
		private List<Weighted_Observed_Point> sort_observations(Collection<Weighted_Observed_Point> unsorted)
		{
			const List<Weighted_Observed_Point> observations = Array_list<>(unsorted);

			// sin_ce the samples are almost always already sorted, this
			// method is implemented as an insertion sort that reorders the
			// elements in place. Insertion sort is very efficient in this case.
			Weighted_Observed_Point curr = observations.get(0);
			const int len = observations.size();
			for (int j{ 1 }; j < len; j++)
			{
				Weighted_Observed_Point prec = curr;
				curr = observations.get(j);
				if (curr.get_x() < prec.get_x())
				{
					// the current element should be inserted closer to the beginning
					int i = j - 1;
					Weighted_Observed_Point mI = observations.get(i);
					while ((i >= 0) && (curr.get_x() < mI.get_x()))
					{
						observations.set(i + 1, mI);
						if (i != 0)
						{
							mI = observations.get(i - 1);
						}
						--i;
					}
					observations.set(i + 1, curr);
					curr = observations.get(j);
				}
			}

			return observations;
		}

		/**
		 * Estimate a first guess of the amplitude and angular frequency.
		 *
		 * @param observations Observations, sorted w.r.t. abscissa.
		 * @ if the abscissa range is zero.
		 * @Math_Illegal_State_Exception when the guessing procedure cannot
		 * produce sensible results.
		 * @return the guessed amplitude (at index 0) and circular frequency
		 * (at index 1).
		 */
		private std::vector<double> guess_a_omega(Weighted_Observed_Point observations)
		{
			const std::vector<double> a_omega = std::vector<double>(2);

			// initialize the sums for the linear model between the two integrals
			double sx2 = 0;
			double sy2 = 0;
			double sxy = 0;
			double sxz = 0;
			double syz = 0;

			double current_x = observations[0].get_x();
			double current_y = observations[0].get_y();
			double f2_integral = 0;
			double f_prime2_integral = 0;
			const double start_x = current_x;
			for (int i{ 1 }; i < observations.size(); ++i)
			{
				// one step forward
				const double previous_x = current_x;
				const double previous_y = current_y;
				current_x = observations[i].get_x();
				current_y = observations[i].get_y();

				// update the integrals of f<sup>2</sup> and f'<sup>2</sup>
				// considering a linear model for f (and therefore constant f')
				const double dx = current_x - previous_x;
				const double dy = current_y - previous_y;
				const double f2_step_integral =
					dx * (previous_y * previous_y + previous_y * current_y + current_y * current_y) / 3;
				const double f_prime2_step_integral = dy * dy / dx;

				const double x = current_x - start_x;
				f2_integral += f2_step_integral;
				f_prime2_integral += f_prime2_step_integral;

				sx2 += x * x;
				sy2 += f2_integral * f2_integral;
				sxy += x * f2_integral;
				sxz += x * f_prime2_integral;
				syz += f2_integral * f_prime2_integral;
			}

			// compute the amplitude and pulsation coefficients
			double c1 = sy2 * sxz - sxy * syz;
			double c2 = sxy * sxz - sx2 * syz;
			double c3 = sx2 * sy2 - sxy * sxy;
			if ((c1 / c2 < 0) || (c2 / c3 < 0))
			{
				const int last = observations.size() - 1;
				// Range of the observations, assuming that the
				// observations are sorted.
				const double x_range = observations[last].get_x() - observations[0].get_x();
				if (x_range == 0)
				{
					throw (hipparchus::exception::Localized_Core_Formats_Type::ZERO_NOT_ALLOWED);
				}
				a_omega[1] = 2 * std::numbers::pi / x_range;

				double y_min = INFINITY;
				double y_max = -INFINITY;
				for (int i{ 1 }; i < observations.size(); ++i)
				{
					const double y = observations[i].get_y();
					if (y < y_min)
					{
						y_min = y;
					}
					if (y > y_max)
					{
						y_max = y;
					}
				}
				a_omega[0] = 0.5 * (y_max - y_min);
			}
			else
			{
				if (c2 == 0)
				{
					// In some ill-conditioned cases (cf. MATH-844), the guesser
					// procedure cannot produce sensible results.
					throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR);
				}

				a_omega[0] = std::sqrt(c1 / c2);
				a_omega[1] = std::sqrt(c2 / c3);
			}

			return a_omega;
		}

		/**
		 * Estimate a first guess of the phase.
		 *
		 * @param observations Observations, sorted w.r.t. abscissa.
		 * @return the guessed phase.
		 */
		private double guess_phi(Weighted_Observed_Point observations)
		{
			// initialize the means
			double fc_mean = 0;
			double fs_mean = 0;

			double current_x = observations[0].get_x();
			double current_y = observations[0].get_y();
			for (int i{ 1 }; i < observations.size(); ++i)
			{
				// one step forward
				const double previous_x = current_x;
				const double previous_y = current_y;
				current_x = observations[i].get_x();
				current_y = observations[i].get_y();
				const double current_y_prime = (current_y - previous_y) / (current_x - previous_x);

				double omega_x = omega * current_x;
				Sin_Cos sc = Sin_Cos(omega_x);
				fc_mean += omega * current_y * sc.cos() - current_y_prime * sc.sin();
				fs_mean += omega * current_y * sc.sin() + current_y_prime * sc.cos();
			}

			return std::atan2(-fs_mean, fc_mean);
		}
	}
}