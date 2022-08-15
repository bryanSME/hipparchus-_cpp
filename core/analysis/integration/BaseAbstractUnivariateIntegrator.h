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
  //package org.hipparchus.analysis.integration;

  //import org.hipparchus.analysis.Univariate_Function;
  //import org.hipparchus.analysis.solvers.Univariate_Solver_Utils;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.exception.Null_Argument_Exception;
  //import org.hipparchus.util.Incrementor;
  //import org.hipparchus.util.Math_Utils;
#include "UnivariateIntegrator.h"

  /**
   * Provide a default implementation for several generic functions.
   *
   */
class Base_Abstract_Univariate_Integrator : public Univariate_Integrator
{
	/** Default absolute accuracy. */
	public static const double DEFAULT_ABSOLUTE_ACCURACY = 1.0e-15;

	/** Default relative accuracy. */
	public static const double DEFAULT_RELATIVE_ACCURACY = 1.0e-6;

	/** Default minimal iteration count. */
	public static const int DEFAULT_MIN_ITERATIONS_COUNT = 3;

	/** Default maximal iteration count. */
	public static const int DEFAULT_MAX_ITERATIONS_COUNT = std::numeric_limits<int>::max();

	/** The iteration count. */
	protected const Incrementor iterations;

	/** Maximum absolute error. */
	private const double& absolute_accuracy;

	/** Maximum relative error. */
	private const double relative_accuracy;

	/** minimum number of iterations */
	private const int minimal_iteration_count;

	/** The functions evaluation count. */
	private Incrementor evaluations;

	/** Function to integrate. */
	private Univariate_Function function;

	/** Lower bound for the interval. */
	private double min;

	/** Upper bound for the interval. */
	private double max;

	/**
	 * Construct an integrator with given accuracies and iteration counts.
	 * <p>
	 * The meanings of the various parameters are:
	 * <ul>
	 *   <li>relative accuracy:
	 *       this is used to stop iterations if the absolute accuracy can't be
	 *       achieved due to large values or short mantissa length. If this
	 *       should be the primary criterion for convergence rather then a
	 *       safety measure, set the absolute accuracy to a ridiculously small value, *       like {@link org.hipparchus.util.Precision#SAFE_MIN Precision.SAFE_MIN}.</li>
	 *   <li>absolute accuracy:
	 *       The default is usually chosen so that results in the interval
	 *       -10..-0.1 and +0.1..+10 can be found with a reasonable accuracy. If the
	 *       expected absolute value of your results is of much smaller magnitude, set
	 *       this to a smaller value.</li>
	 *   <li>minimum number of iterations:
	 *       minimal iteration is needed to avoid false early convergence, e.g.
	 *       the sample points happen to be zeroes of the function. Users can
	 *       use the default value or choose one that they see as appropriate.</li>
	 *   <li>maximum number of iterations:
	 *       usually a high iteration count indicates convergence problems. However, *       the "reasonable value" varies widely for different algorithms. Users are
	 *       advised to use the default value supplied by the algorithm.</li>
	 * </ul>
	 *
	 * @param relative_accuracy relative accuracy of the result
	 * @param absolute_accuracy absolute accuracy of the result
	 * @param minimal_iteration_count minimum number of iterations
	 * @param maximal_iteration_count maximum number of iterations
	 * @exception  if minimal number of iterations
	 * is not strictly positive
	 * @exception  if maximal number of iterations
	 * is lesser than or equal to the minimal number of iterations
	 */
	protected Base_Abstract_Univariate_Integrator(const double& relative_accuracy, const double& absolute_accuracy, const int minimal_iteration_count, const int maximal_iteration_count)
	{
		// accuracy settings
		this.relative_accuracy = relative_accuracy;
		this.absolute_accuracy = absolute_accuracy;

		// iterations count settings
		if (minimal_iteration_count <= 0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, minimal_iteration_count, 0);
		}
		if (maximal_iteration_count <= minimal_iteration_count)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, maximal_iteration_count, minimal_iteration_count);
		}
		this.minimal_iteration_count = minimal_iteration_count;
		this.iterations = Incrementor(maximal_iteration_count);

		// prepare evaluations counter, but do not set it yet
		evaluations = Incrementor();
	}

	/**
	 * Construct an integrator with given accuracies.
	 * @param relative_accuracy relative accuracy of the result
	 * @param absolute_accuracy absolute accuracy of the result
	 */
	protected Base_Abstract_Univariate_Integrator(const double& relative_accuracy, const double& absolute_accuracy)
	{
		this(relative_accuracy, absolute_accuracy, DEFAULT_MIN_ITERATIONS_COUNT, DEFAULT_MAX_ITERATIONS_COUNT);
	}

	/**
	 * Construct an integrator with given iteration counts.
	 * @param minimal_iteration_count minimum number of iterations
	 * @param maximal_iteration_count maximum number of iterations
	 * @exception  if minimal number of iterations
	 * is not strictly positive
	 * @exception  if maximal number of iterations
	 * is lesser than or equal to the minimal number of iterations
	 */
	protected Base_Abstract_Univariate_Integrator(const int minimal_iteration_count, const int maximal_iteration_count)
	{
		this(DEFAULT_RELATIVE_ACCURACY, DEFAULT_ABSOLUTE_ACCURACY, minimal_iteration_count, maximal_iteration_count);
	}

	/** {@inherit_doc} */
	//override
	public double get_relative_accuracy()
	{
		return relative_accuracy;
	}

	/** {@inherit_doc} */
	//override
	public double get_absolute_accuracy()
	{
		return absolute_accuracy;
	}

	/** {@inherit_doc} */
	//override
	public int get_minimal_iteration_count()
	{
		return minimal_iteration_count;
	}

	/** {@inherit_doc} */
	//override
	public int get_maximal_iteration_count()
	{
		return iterations.get_maximal_count();
	}

	/** {@inherit_doc} */
	//override
	public int get_evaluations()
	{
		return evaluations.get_count();
	}

	/** {@inherit_doc} */
	//override
	public int get_iterations()
	{
		return iterations.get_count();
	}

	/**
	 * @return the lower bound.
	 */
	protected double get_min()
	{
		return min;
	}
	/**
	 * @return the upper bound.
	 */
	protected double get_max()
	{
		return max;
	}

	/**
	 * Compute the objective function value.
	 *
	 * @param point Point at which the objective function must be evaluated.
	 * @return the objective function value at specified point.
	 * @Math_Illegal_State_Exception if the maximal number of function
	 * evaluations is exceeded.
	 */
	protected double compute_objective_value(const double point)
		Math_Illegal_State_Exception
	{
		evaluations.increment();
		return function.value(point);
	}

	/**
	 * Prepare for computation.
	 * Subclasses must call this method if they //override any of the
	 * {@code solve} methods.
	 *
	 * @param max_eval Maximum number of evaluations.
	 * @param f the integrand function
	 * @param lower the min bound for the interval
	 * @param upper the upper bound for the interval
	 * @Null_Argument_Exception if {@code f} is {@code NULL}.
	 * @ if {@code min >= max}.
	 */
	protected void setup(const int max_eval, const Univariate_Function& f, const double lower, const double upper)
	{
		// Checks.
		//Math_Utils::check_not_null(f);
		Univariate_Solver_Utils.verify_interval(lower, upper);

		// Reset.
		min = lower;
		max = upper;
		function = f;
		evaluations = evaluations.with_maximal_count(max_eval);
		iterations.reset();
	}

	/** {@inherit_doc} */
	//override
	public double integrate(const int max_eval, const Univariate_Function& f, const double lower, const double upper)
	{
		// Initialization.
		setup(max_eval, f, lower, upper);

		// Perform computation.
		return do_integrate();
	}

	/**
	 * Method for implementing actual integration algorithms in derived
	 * classes.
	 *
	 * @return the root.
	 * @Math_Illegal_State_Exception if the maximal number of evaluations
	 * is exceeded.
	 * @Math_Illegal_State_Exception if the maximum iteration count is exceeded
	 * or the integrator detects convergence problems otherwise
	 */
	protected virtual double do_integrate()
};