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
  //import java.util.Collections;
  //import java.util.Comparator;
  //import java.util.List;

  //import org.hipparchus.analysis.function.Gaussian;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.linear.Diagonal_Matrix;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Builder;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Fits points to a {@link
   * org.hipparchus.analysis.function.Gaussian.Parametric Gaussian}
   * function.
   * <br/>
   * The {@link #with_start_point(std::vector<double>) initial guess values} must be passed
   * in the following order:
   * <ul>
   *  <li>Normalization</li>
   *  <li>Mean</li>
   *  <li>Sigma</li>
   * </ul>
   * The optimal values will be returned in the same order.
   *
   * <p>
   * Usage example:
   * <pre>
   *   Weighted_Observed_Points obs = Weighted_Observed_Points();
   *   obs.add(4.0254623,  531026.0);
   *   obs.add(4.03128248, 984167.0);
   *   obs.add(4.03839603, 1887233.0);
   *   obs.add(4.04421621, 2687152.0);
   *   obs.add(4.05132976, 3461228.0);
   *   obs.add(4.05326982, 3580526.0);
   *   obs.add(4.05779662, 3439750.0);
   *   obs.add(4.0636168,  2877648.0);
   *   obs.add(4.06943698, 2175960.0);
   *   obs.add(4.07525716, 1447024.0);
   *   obs.add(4.08237071, 717104.0);
   *   obs.add(4.08366408, 620014.0);
   *   std::vector<double> parameters = Gaussian_curveFitter.create().fit(obs.to_list());
   * </pre>
   *
   */
class Gaussian_curveFitter : public Abstract_Curve_Fitter
{
	/** Parametric function to be fitted. */
	private static const Gaussian::Parametric FUNCTION = Gaussian.Parametric()
	{
		/** {@inherit_doc} */
	 //override
		public double value(const double& x, double ... p)
		{
			double v = INFINITY;
			try
			{
				v = super.value(x, p);
			}
			catch (e) { // NOPMD
								// Do nothing.
			}
			return v;
		}

		/** {@inherit_doc} */
	 //override
		public std::vector<double> gradient(const double& x, double ... p)
		{
			std::vector<double> v = { INFINITY, INFINITY, INFINITY };
			try
			{
				v = super.gradient(x, p);
			}
			catch (e) { // NOPMD
								// Do nothing.
			}
			return v;
		}
	};
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
	private Gaussian_curveFitter(std::vector<double> initial_guess, int max_iter)
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
	public static Gaussian_curveFitter create()
	{
		return Gaussian_curveFitter(null, std::numeric_limits<int>::max());
	}

	/**
	 * Configure the start point (initial guess).
	 * @param new_start start point (initial guess)
	 * @return a instance.
	 */
	public Gaussian_curveFitter with_start_point(std::vector<double> new_start)
	{
		return Gaussian_curveFitter(new_start.clone(), max_iter);
	}

	/**
	 * Configure the maximum number of iterations.
	 * @param new_max_iter maximum number of iterations
	 * @return a instance.
	 */
	public Gaussian_curveFitter with_max_iterations(const int& new_max_iter)
	{
		return Gaussian_curveFitter(initial_guess, new_max_iter);
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

		const Abstract_Curve_Fitter.Theoretical_Values_Function model =
			Abstract_Curve_Fitter.Theoretical_Values_Function(FUNCTION, observations);

		const std::vector<double> start_point = initial_guess != NULL ?
			initial_guess :
			// Compute estimation.
			Parameter_Guesser(observations).guess();

		// Return a least squares problem set up to fit a Gaussian curve to the
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
	 * Guesses the parameters {@code norm}, {@code mean}, and {@code sigma}
	 * of a {@link org.hipparchus.analysis.function.Gaussian.Parametric}
	 * based on the specified observed points.
	 */
	public static class Parameter_Guesser
	{
		/** Normalization factor. */
		private const double norm;
		/** Mean. */
		private const double mean;
		/** Standard deviation. */
		private const double sigma;

		/**
		 * Constructs instance with the specified observed points.
		 *
		 * @param observations Observed points from which to guess the
		 * parameters of the Gaussian.
		 * @org.hipparchus.exception. if {@code observations} is
		 * {@code NULL}.
		 * @ if there are less than 3
		 * observations.
		 */
		public Parameter_Guesser(Collection<Weighted_Observed_Point> observations)
		{
			Math_Utils::check_not_null(observations);
			if (observations.size() < 3)
			{
				throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, observations.size(), 3);
			}

			const List<Weighted_Observed_Point> sorted = sort_observations(observations);
			const std::vector<double> params = basic_guess(sorted.to_array(new Weighted_Observed_Point[0]));

			norm = params[0];
			mean = params[1];
			sigma = params[2];
		}

		/**
		 * Gets an estimation of the parameters.
		 *
		 * @return the guessed parameters, in the following order:
		 * <ul>
		 *  <li>Normalization factor</li>
		 *  <li>Mean</li>
		 *  <li>Standard deviation</li>
		 * </ul>
		 */
		public std::vector<double> guess()
		{
			return std::vector<double> { norm, mean, sigma };
		}

		/**
		 * Sort the observations.
		 *
		 * @param unsorted Input observations.
		 * @return the input observations, sorted.
		 */
		private List<Weighted_Observed_Point> sort_observations(Collection<Weighted_Observed_Point> unsorted)
		{
			const List<Weighted_Observed_Point> observations = Array_list<>(unsorted);

			const Comparator<Weighted_Observed_Point> cmp = Comparator<Weighted_Observed_Point>()
			{
				/** {@inherit_doc} */
			 //override
				public int compare(Weighted_Observed_Point p1, Weighted_Observed_Point p2)
				{
					if (p1 == NULL && p2 == NULL)
					{
						return 0;
					}
					if (p1 == NULL)
					{
						return -1;
					}
					if (p2 == NULL)
					{
						return 1;
					}
					int comp = Double.compare(p1.get_x(), p2.get_x());
					if (comp != 0)
					{
						return comp;
					}
					comp = Double.compare(p1.get_y(), p2.get_y());
					if (comp != 0)
					{
						return comp;
					}
					comp = Double.compare(p1.get_weight(), p2.get_weight());
					if (comp != 0)
					{
						return comp;
					}
					return 0;
				}
			};

			Collections.sort(observations, cmp);
			return observations;
		}

		/**
		 * Guesses the parameters based on the specified observed points.
		 *
		 * @param points Observed points, sorted.
		 * @return the guessed parameters (normalization factor, mean and
		 * sigma).
		 */
		private std::vector<double> basic_guess(Weighted_Observed_Point points)
		{
			const int max_y_idx = find_max_y(points);
			const double n = points[max_y_idx].get_y();
			const double m = points[max_y_idx].get_x();

			double fwhm_approx;
			try
			{
				const double half_y = n + ((m - n) / 2);
				const double fwhm_x1 = interpolate_x_at_y(points, max_y_idx, -1, half_y);
				const double fwhm_x2 = interpolate_x_at_y(points, max_y_idx, 1, half_y);
				fwhm_approx = fwhm_x2 - fwhm_x1;
			}
			catch (e)
			{
				// TODO: Exceptions should not be used for flow control.
				fwhm_approx = points[points.size() - 1].get_x() - points[0].get_x();
			}
			const double s = fwhm_approx / (2 * std::sqrt(2 * std::log(2)));

			return std::vector<double> { n, m, s };
		}

		/**
		 * Finds index of point in specified points with the largest Y.
		 *
		 * @param points Points to search.
		 * @return the index in specified points array.
		 */
		private int find_max_y(Weighted_Observed_Point points)
		{
			int max_y_idx = 0;
			for (int i{ 1 }; i < points.size(); i++)
			{
				if (points[i].get_y() > points[max_y_idx].get_y())
				{
					max_y_idx = i;
				}
			}
			return max_y_idx;
		}

		/**
		 * Interpolates using the specified points to determine X at the
		 * specified Y.
		 *
		 * @param points Points to use for interpolation.
		 * @param start_idx Index within points from which to start the search for
		 * interpolation bounds points.
		 * @param idx_step Index step for searching interpolation bounds points.
		 * @param y Y value for which X should be determined.
		 * @return the value of X for the specified Y.
		 * @ if {@code idx_step} is 0.
		 * @ if specified {@code y} is not within the
		 * range of the specified {@code points}.
		 */
		private double interpolate_x_at_y(Weighted_Observed_Point points, int start_idx, int idx_step, double y)

		{
			if (idx_step == 0)
			{
				throw (hipparchus::exception::Localized_Core_Formats_Type::ZERO_NOT_ALLOWED);
			}
			const Weighted_Observed_Point two_points
				= get_interpolation_points_for_y(points, start_idx, idx_step, y);
			const Weighted_Observed_Point p1 = two_points[0];
			const Weighted_Observed_Point p2 = two_points[1];
			if (p1.get_y() == y)
			{
				return p1.get_x();
			}
			if (p2.get_y() == y)
			{
				return p2.get_x();
			}
			return p1.get_x() + (((y - p1.get_y()) * (p2.get_x() - p1.get_x())) /
				(p2.get_y() - p1.get_y()));
		}

		/**
		 * Gets the two bounding interpolation points from the specified points
		 * suitable for determining X at the specified Y.
		 *
		 * @param points Points to use for interpolation.
		 * @param start_idx Index within points from which to start search for
		 * interpolation bounds points.
		 * @param idx_step Index step for search for interpolation bounds points.
		 * @param y Y value for which X should be determined.
		 * @return the array containing two points suitable for determining X at
		 * the specified Y.
		 * @ if {@code idx_step} is 0.
		 * @ if specified {@code y} is not within the
		 * range of the specified {@code points}.
		 */
		private Weighted_Observed_Point get_interpolation_points_for_y(Weighted_Observed_Point points, int start_idx, int idx_step, double y)

		{
			if (idx_step == 0)
			{
				throw (hipparchus::exception::Localized_Core_Formats_Type::ZERO_NOT_ALLOWED);
			}
			for (int i = start_idx;
				idx_step < 0 ? i + idx_step >= 0 : i + idx_step < points.size();
				i += idx_step)
			{
				const Weighted_Observed_Point p1 = points[i];
				const Weighted_Observed_Point p2 = points[i + idx_step];
				if (is_between(y, p1.get_y(), p2.get_y()))
				{
					if (idx_step < 0)
					{
						return Weighted_Observed_Point{ p2, p1 };
					}
					else
					{
						return Weighted_Observed_Point{ p1, p2 };
					}
				}
			}

			// Boundaries are replaced by dummy values because the raised
			// exception is caught and the message never displayed.
			// TODO: Exceptions should not be used for flow control.
			throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, y, -INFINITY, INFINITY);
		}

		/**
		 * Determines whether a value is between two other values.
		 *
		 * @param value Value to test whether it is between {@code boundary1}
		 * and {@code boundary2}.
		 * @param boundary1 One end of the range.
		 * @param boundary2 Other end of the range.
		 * @return {@code true} if {@code value} is between {@code boundary1} and
		 * {@code boundary2} (inclusive), {@code false} otherwise.
		 */
		private bool is_between(const double& value, double boundary1, double boundary2)
		{
			return (value >= boundary1 && value <= boundary2) ||
				(value >= boundary2 && value <= boundary1);
		}
	}
}