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

  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.util.FastMath;
#include "AbstractUnivariateSolver.h"

  /**
   * Implements the <a href="http://mathworld.wolfram.com/Bisection.html">
   * bisection algorithm</a> for finding zeros of univariate real functions.
   * <p>
   * The function should be continuous but not necessarily smooth.</p>
   *
   */
class Bisection_Solver : public Abstract_Univariate_Solver
{
private:
	/** Default absolute accuracy. */
	static constexpr double DEFAULT_ABSOLUTE_ACCURACY{ 1e-6 };

public:
	/**
	 * Construct a solver with default accuracy (1e-6).
	 */
	Bisection_Solver()
	{
		Bisection_Solver(DEFAULT_ABSOLUTE_ACCURACY);
	}
	/**
	 * Construct a solver.
	 *
	 * @param absolute_accuracy Absolute accuracy.
	 */
	Bisection_Solver(const double& absolute_accuracy)
	{
		Abstract_Univariate_Solver(absolute_accuracy);
	}
	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 */
	public Bisection_Solver(const double& relative_accuracy, const double& absolute_accuracy)
	{
		Abstract_Univariate_Solver(relative_accuracy, absolute_accuracy);
	}

protected:
	/**
	 * {@inherit_doc}
	 */
	 //override
	double do_solve()
	{
		double min = get_min();
		double max = get_max();
		verify_interval(min, max);
		verify_bracketing(min, max);
		const auto absolute_accuracy = get_absolute_accuracy();
		double m;
		double fm;
		double fmin;

		while (true)
		{
			m = Univariate_Solver_Utils::midpoint(min, max);
			fmin = compute_objective_value(min);
			fm = compute_objective_value(m);

			if (fm * fmin > 0)
			{
				// max and m bracket the root.
				min = m;
			}
			else
			{
				// min and m bracket the root.
				max = m;
			}

			if (std::abs(max - min) <= absolute_accuracy)
			{
				m = Univariate_Solver_Utils::midpoint(min, max);
				return m;
			}
		}
	}
};