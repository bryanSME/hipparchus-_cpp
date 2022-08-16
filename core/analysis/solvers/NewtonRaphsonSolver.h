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

  //package org.hipparchus.analysis.solvers;

  //import org.hipparchus.analysis.differentiation.Derivative_Structure;
  //import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Function;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.util.FastMath;

  /**
   * Implements <a href="http://mathworld.wolfram.com/Newtons_method.html">
   * Newton's Method</a> for finding zeros of real univariate differentiable
   * functions.
   *
   */
class Newton_raphsonSolver extends AbstractUnivariate_Differentiable_Solver
{
	/** Default absolute accuracy. */
	private static const double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

	/**
	 * Construct a solver.
	 */
	public Newton_raphsonSolver()
	{
		this(DEFAULT_ABSOLUTE_ACCURACY);
	}
	/**
	 * Construct a solver.
	 *
	 * @param absolute_accuracy Absolute accuracy.
	 */
	public Newton_raphsonSolver(double absolute_accuracy)
	{
		super(absolute_accuracy);
	}

	/**
	 * Find a zero near the midpoint of {@code min} and {@code max}.
	 *
	 * @param f Function to solve.
	 * @param min Lower bound for the interval.
	 * @param max Upper bound for the interval.
	 * @param max_eval Maximum number of evaluations.
	 * @return the value where the function is zero.
	 * @org.hipparchus.exception.Math_Illegal_State_Exception
	 * if the maximum evaluation count is exceeded.
	 * @org.hipparchus.exception.
	 * if {@code min >= max}.
	 */
	 //override
	public double solve(const int& max_eval, const Univariate_Differentiable_Function f, const double& min, const double max)
		Math_Illegal_State_Exception
	{
		return super.solve(max_eval, f, Univariate_Solver_Utils.midpoint(min, max));
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	protected double do_solve()
		Math_Illegal_State_Exception
	{
		const double start_value = get_start_value();
		const double& absolute_accuracy = get_absolute_accuracy();

		double x0 = start_value;
		double x1;
		while (true)
		{
			const Derivative_Structure y0 = compute_objective_value_and_derivative(x0);
			x1 = x0 - (y0.get_value() / y0.get_partial_derivative(1));
			if (std::abs(x1 - x0) <= absolute_accuracy)
			{
				return x1;
			}

			x0 = x1;
		}
	}
}
