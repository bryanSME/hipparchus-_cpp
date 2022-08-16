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

  //import org.hipparchus.Calculus_Field_Element;
  //import org.hipparchus.Field;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.util.FastMath;
#include <type_traits>
#include "../../CalculusFieldElement.hpp"

/**
 * Implements the <a href="http://en.wikipedia.org/wiki/Midpoint_method">
 * Midpoint Rule</a> for integration of real univariate functions. For
 * reference, see <b>Numerical Mathematics</b>, ISBN 0387989595, * chapter 9.2.
 * <p>
 * The function should be integrable.</p>
 * @param <T> Type of the field elements.
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldMid_pointIntegrator : public BaseAbstractField_Univariate_Integrator<T>
{
	/** Maximum number of iterations for midpoint. */
	public static const int MIDPOINT_MAX_ITERATIONS_COUNT = 64;

	/**
	 * Build a midpoint integrator with given accuracies and iterations counts.
	 * @param field field to which function argument and value belong
	 * @param relative_accuracy relative accuracy of the result
	 * @param absolute_accuracy absolute accuracy of the result
	 * @param minimal_iteration_count minimum number of iterations
	 * @param maximal_iteration_count maximum number of iterations
	 * (must be less than or equal to {@link #MIDPOINT_MAX_ITERATIONS_COUNT}
	 * @exception  if minimal number of iterations
	 * is not strictly positive
	 * @exception  if maximal number of iterations
	 * is lesser than or equal to the minimal number of iterations
	 * @exception  if maximal number of iterations
	 * is greater than {@link #MIDPOINT_MAX_ITERATIONS_COUNT}
	 */
	public FieldMid_pointIntegrator(const Field<T> field, const double relative_accuracy, const double& absolute_accuracy, const int minimal_iteration_count, const int maximal_iteration_count)

	{
		super(field, relative_accuracy, absolute_accuracy, minimal_iteration_count, maximal_iteration_count);
		if (maximal_iteration_count > MIDPOINT_MAX_ITERATIONS_COUNT)
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, maximal_iteration_count, MIDPOINT_MAX_ITERATIONS_COUNT);
		}
	}

	/**
	 * Build a midpoint integrator with given iteration counts.
	 * @param field field to which function argument and value belong
	 * @param minimal_iteration_count minimum number of iterations
	 * @param maximal_iteration_count maximum number of iterations
	 * (must be less than or equal to {@link #MIDPOINT_MAX_ITERATIONS_COUNT}
	 * @exception  if minimal number of iterations
	 * is not strictly positive
	 * @exception  if maximal number of iterations
	 * is lesser than or equal to the minimal number of iterations
	 * @exception  if maximal number of iterations
	 * is greater than {@link #MIDPOINT_MAX_ITERATIONS_COUNT}
	 */
	public FieldMid_pointIntegrator(const Field<T> field, const int minimal_iteration_count, const int maximal_iteration_count)

	{
		super(field, minimal_iteration_count, maximal_iteration_count);
		if (maximal_iteration_count > MIDPOINT_MAX_ITERATIONS_COUNT)
		{
			throw std::exception("not implmented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE_BOUND_EXCLUDED, maximal_iteration_count, MIDPOINT_MAX_ITERATIONS_COUNT);
		}
	}

	/**
	 * Construct a midpoint integrator with default settings.
	 * @param field field to which function argument and value belong
	 * (max iteration count set to {@link #MIDPOINT_MAX_ITERATIONS_COUNT})
	 */
	public FieldMid_pointIntegrator(const Field<T> field)
	{
		super(field, DEFAULT_MIN_ITERATIONS_COUNT, MIDPOINT_MAX_ITERATIONS_COUNT);
	}

	/**
	 * Compute the n-th stage integral of midpoint rule.
	 * This function should only be called by API <code>integrate()</code> in the //package.
	 * To save time it does not verify arguments - caller does.
	 * <p>
	 * The interval is divided equally into 2^n sections rather than an
	 * arbitrary m sections because this configuration can best utilize the
	 * already computed values.</p>
	 *
	 * @param n the stage of 1/2 refinement. Must be larger than 0.
	 * @param previous_stage_result Result from the previous call to the
	 * {@code stage} method.
	 * @param min Lower bound of the integration interval.
	 * @param diff_max_min Difference between the lower bound and upper bound
	 * of the integration interval.
	 * @return the value of n-th stage integral
	 * @Math_Illegal_State_Exception if the maximal number of evaluations
	 * is exceeded.
	 */
	private T stage(const int& n, T previous_stage_result, T min, T diff_max_min)
		Math_Illegal_State_Exception
	{
		// number of points in this stage
		const long np = 1L << (n - 1);
		T sum = get_field().get_zero();

		// spacing between adjacent points
		const T spacing = diff_max_min.divide(np);

		// the first point
		T x = min.add(spacing.multiply(0.5));
		for (long i = 0; i < np; i++)
		{
			sum = sum.add(compute_objective_value(x));
			x = x.add(spacing);
		}
		// add the sum to previously calculated result
		return previous_stage_result.add(sum.multiply(spacing)).multiply(0.5);
	}

	/** {@inherit_doc} */
	//override
	protected T do_integrate()

	{
		const T min = get_min();
		const T diff = get_max().subtract(min);
		const T mid_point = min.add(diff.multiply(0.5));

		T oldt = diff.multiply(compute_objective_value(mid_point));

		while (true)
		{
			iterations.increment();
			const int i = iterations.get_count();
			const T t = stage(i, oldt, min, diff);
			if (i >= get_minimal_iteration_count())
			{
				const double delta = std::abs(t.subtract(oldt)).get_real();
				const double r_limit = std::abs(oldt).add(std::abs(t)).multiply(0.5 * get_relative_accuracy()).get_real();
				if ((delta <= r_limit) || (delta <= get_absolute_accuracy()))
				{
					return t;
				}
			}
			oldt = t;
		}
	}
};