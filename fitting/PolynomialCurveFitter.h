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

  //import java.util.Collection;

  //import org.hipparchus.analysis.polynomials.Polynomial_Function;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.linear.Diagonal_Matrix;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Builder;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem;

  /**
   * Fits points to a {@link
   * org.hipparchus.analysis.polynomials.Polynomial_Function.Parametric polynomial}
   * function.
   * <br/>
   * The size of the {@link #with_start_point(std::vector<double>) initial guess} array defines the
   * degree of the polynomial to be fitted.
   * They must be sorted in increasing order of the polynomial's degree.
   * The optimal values of the coefficients will be returned in the same order.
   *
   */
class Polynomial_Curve_Fitter extends Abstract_Curve_Fitter
{
	/** Parametric function to be fitted. */
	private static const Polynomial_Function.Parametric FUNCTION = Polynomial_Function.Parametric();
	/** Initial guess. */
	private const std::vector<double> initial_guess;
	/** Maximum number of iterations of the optimization algorithm. */
	private const int max_iter;

	/**
	 * Constructor used by the factory methods.
	 *
	 * @param initial_guess Initial guess.
	 * @param max_iter Maximum number of iterations of the optimization algorithm.
	 * @Math_Runtime_Exception if {@code initial_guess} is {@code NULL}.
	 */
	private Polynomial_Curve_Fitter(std::vector<double> initial_guess, int max_iter)
	{
		this.initial_guess = initial_guess.clone();
		this.max_iter = max_iter;
	}

	/**
	 * Creates a default curve fitter.
	 * Zero will be used as initial guess for the coefficients, and the maximum
	 * number of iterations of the optimization algorithm is set to
	 * {@link Integer#MAX_VALUE}.
	 *
	 * @param degree Degree of the polynomial to be fitted.
	 * @return a curve fitter.
	 *
	 * @see #with_start_point(std::vector<double>)
	 * @see #with_max_iterationsstatic_cast<int>(
	 */
	public static Polynomial_Curve_Fitter create(const int& degree)
	{
		return Polynomial_Curve_Fitter(std::vector<double>(degree + 1], std::numeric_limits<int>::max());
	}

	/**
	 * Configure the start point (initial guess).
	 * @param new_start start point (initial guess)
	 * @return a instance.
	 */
	public Polynomial_Curve_Fitter with_start_point(std::vector<double> new_start)
	{
		return Polynomial_Curve_Fitter(new_start.clone(), max_iter);
	}

	/**
	 * Configure the maximum number of iterations.
	 * @param new_max_iter maximum number of iterations
	 * @return a instance.
	 */
	public Polynomial_Curve_Fitter with_max_iterations(const int& new_max_iter)
	{
		return Polynomial_Curve_Fitter(initial_guess, new_max_iter);
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

		if (initial_guess == NULL)
		{
			throw Math_Runtime_Exception.create_internal_error();
		}

		// Return a least squares problem set up to fit a polynomial curve to the
		// observed points.
		return Least_Squares_Builder().
			max_evaluations(std::numeric_limits<int>::max()).
			max_iterations(max_iter).
			start(initial_guess).
			target(target).
			weight(new Diagonal_Matrix(weights)).
			model(model.get_model_function(), model.get_model_function_jacobian()).
			build();
	}
}