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
#include <vector>
#include <limits>
#include "AbstractCurveFitter.h"

  //import org.hipparchus.analysis.Parametric_Univariate_Function ;
  //import org.hipparchus.linear.Diagonal_Matrix;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Builder;
  //import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem;

  /**
   * Fits points to a user-defined {@link Parametric_Univariate_Function  function}.
   *
   */
class Simple_Curve_Fitter : public Abstract_Curve_Fitter
{
private:
	/** Function to fit. */
	const Parametric_Univariate_Function  function;
	/** Initial guess for the parameters. */
	const std::vector<double> initial_guess;
	/** Maximum number of iterations of the optimization algorithm. */
	const int max_iter;

	/**
	 * Constructor used by the factory methods.
	 *
	 * @param function Function to fit.
	 * @param initial_guess Initial guess. Cannot be {@code NULL}. Its length must
	 * be consistent with the number of parameters of the {@code function} to fit.
	 * @param max_iter Maximum number of iterations of the optimization algorithm.
	 */
	Simple_Curve_Fitter(Parametric_Univariate_Function  function, std::vector<double> initial_guess, int max_iter)
	{
		this.function = function;
		this.initial_guess = initial_guess.clone();
		this.max_iter = max_iter;
	}

public:
	/**
	 * Creates a curve fitter.
	 * The maximum number of iterations of the optimization algorithm is set
	 * to {@link Integer#MAX_VALUE}.
	 *
	 * @param f Function to fit.
	 * @param start Initial guess for the parameters.  Cannot be {@code NULL}.
	 * Its length must be consistent with the number of parameters of the
	 * function to fit.
	 * @return a curve fitter.
	 *
	 * @see #with_start_point(std::vector<double>)
	 * @see #with_max_iterationsstatic_cast<int>(
	 */
	static Simple_Curve_Fitter create(Parametricconst Univariate_Function& f, std::vector<double> start)
	{
		return Simple_Curve_Fitter(f, start, std::numeric_limits<int>::max());
	}

	/**
	 * Configure the start point (initial guess).
	 * @param new_start start point (initial guess)
	 * @return a instance.
	 */
	Simple_Curve_Fitter with_start_point(std::vector<double> new_start)
	{
		return Simple_Curve_Fitter(function, new_start.clone(), max_iter);
	}

	/**
	 * Configure the maximum number of iterations.
	 * @param new_max_iter maximum number of iterations
	 * @return a instance.
	 */
	Simple_Curve_Fitter with_max_iterations(const int& new_max_iter)
	{
		return Simple_Curve_Fitter(function, initial_guess, new_max_iter);
	}

protected:
	/** {@inherit_doc} */
	//override
	Least_Squares_Problem get_problem(Collection<Weighted_Observed_Point> observations)
	{
		// Prepare least-squares problem.
		const int len = observations.size();
		auto target = std::vector<double>(len);
		auto weights = std::vector<double>(len);

		int count{};
		for (Weighted_Observed_Point obs : observations)
		{
			target[count] = obs.get_y();
			weights[count] = obs.get_weight();
			++count;
		}

		const Abstract_Curve_Fitter.Theoretical_Values_Function model
			= Abstract_Curve_Fitter.Theoretical_Values_Function(function, observations);

		// Create an optimizer for fitting the curve to the observed points.
		return Least_Squares_Builder().
			max_evaluations(std::numeric_limits<int>::max()).
			max_iterations(max_iter).
			start(initial_guess).
			target(target).
			weight(new Diagonal_Matrix(weights)).
			model(model.get_model_function(), model.get_model_function_jacobian()).
			build();
	}
};