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

#include <type_traits>
#include "../../CalculusFieldElement.hpp"
#include "../CalculusFieldUnivariateFunction.hpp"
#include "FieldUnivariateIntegrator.h"

  /**
   * Provide a default implementation for several generic functions.
   *
   * @param <T> Type of the field elements.
   * @since 2.0
   */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class BaseAbstractField_Univariate_Integrator : Field_Univariate_Integrator<T>
{
private:
	/** Maximum absolute error. */
	const double my_absolute_accuracy;

	/** Maximum relative error. */
	const double my_relative_accuracy;

	/** minimum number of iterations */
	const int my_minimal_iteration_count;

	/** The functions evaluation count. */
	Incrementor my_evaluations;

	/** Field to which function argument and value belong. */
	const Field<T> my_field;

	/** Function to integrate. */
	Calculus_Field_Univariate_Function<T> my_function;

	/** Lower bound for the interval. */
	T my_min;

	/** Upper bound for the interval. */
	T my_max;

protected:
	/** The iteration count. */
	const Incrementor my_iterations;

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
	  * @param field field to which function argument and value belong
	  * @param relative_accuracy relative accuracy of the result
	  * @param absolute_accuracy absolute accuracy of the result
	  * @param minimal_iteration_count minimum number of iterations
	  * @param maximal_iteration_count maximum number of iterations
	  * @exception  if minimal number of iterations
	  * is not strictly positive
	  * @exception  if maximal number of iterations
	  * is lesser than or equal to the minimal number of iterations
	  */
	BaseAbstractField_Univariate_Integrator(const Field<T>& field, const double& relative_accuracy, const double& absolute_accuracy, const int& minimal_iteration_count, const int& maximal_iteration_count)
		:
		my_relative_accuracy{ relative_accuracy },
		my_absolute_accuracy{ absolute_accuracy },
		my_minimal_iteration_count{ minimal_iteration_count },
		my_iterations{ Incrementor(maximal_iteration_count) },
		my_evaluations{ Incrementor() },
		my_field{ field }
	{
		// iterations count settings
		if (minimal_iteration_count <= 0)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, minimal_iteration_count, 0);
		}
		if (maximal_iteration_count <= minimal_iteration_count)
		{
			throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, maximal_iteration_count, minimal_iteration_count);
		}
	}

	/**
	 * Construct an integrator with given accuracies.
	 * @param field field to which function argument and value belong
	 * @param relative_accuracy relative accuracy of the result
	 * @param absolute_accuracy absolute accuracy of the result
	 */
	BaseAbstractField_Univariate_Integrator(const Field<T>& field, const double& relative_accuracy, const double& absolute_accuracy)
	{
		BaseAbstractField_Univariate_Integrator(field, relative_accuracy, absolute_accuracy, DEFAULT_MIN_ITERATIONS_COUNT, DEFAULT_MAX_ITERATIONS_COUNT);
	}

	/**
	 * Construct an integrator with given iteration counts.
	 * @param field field to which function argument and value belong
	 * @param minimal_iteration_count minimum number of iterations
	 * @param maximal_iteration_count maximum number of iterations
	 * @exception  if minimal number of iterations
	 * is not strictly positive
	 * @exception  if maximal number of iterations
	 * is lesser than or equal to the minimal number of iterations
	 */
	BaseAbstractField_Univariate_Integrator(const Field<T>& field, const int& minimal_iteration_count, const int& maximal_iteration_count)
	{
		BaseAbstractField_Univariate_Integrator(field, DEFAULT_RELATIVE_ACCURACY, DEFAULT_ABSOLUTE_ACCURACY, minimal_iteration_count, maximal_iteration_count);
	}

	/**
	 * @return the lower bound.
	 */
	T get_min() const
	{
		return my_min;
	}
	/**
	 * @return the upper bound.
	 */
	T get_max() const
	{
		return my_max;
	}

	/**
	 * Compute the objective function value.
	 *
	 * @param point Point at which the objective function must be evaluated.
	 * @return the objective function value at specified point.
	 * @Math_Illegal_State_Exception if the maximal number of function
	 * evaluations is exceeded.
	 */
	T compute_objective_value(const T& point)
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
	 * @ if {@code f} is {@code NULL}.
	 * @ if {@code min >= max}.
	 */
	void setup(const int& max_eval, const Calculus_Field_Univariate_Function<T>& f, const T& lower, const T& upper)
	{
		// Checks.
		//Math_Utils::check_not_null(f);
		Univariate_Solver_Utils.verify_interval(lower.get_real(), upper.get_real());

		// Reset.
		my_min = lower;
		my_max = upper;
		my_function = f;
		my_evaluations = evaluations.with_maximal_count(max_eval);
		my_iterations.reset();
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
	virtual T do_integrate()
	{
	}

public:

	/** Default absolute accuracy. */
	static constexpr double DEFAULT_ABSOLUTE_ACCURACY{ 1.0e-15 };

	/** Default relative accuracy. */
	static constexpr double DEFAULT_RELATIVE_ACCURACY{ 1.0e-6 };

	/** Default minimal iteration count. */
	static constexpr int DEFAULT_MIN_ITERATIONS_COUNT{ 3 };

	/** Default maximal iteration count. */
	static constexpr int DEFAULT_MAX_ITERATIONS_COUNT{ std::numeric_limits<int>::max() };

	/** Get the field to which function argument and value belong.
	 * @return field to which function argument and value belong
	 */
	Field<T> get_field() const
	{
		return field;
	}

	/** {@inherit_doc} */
	//override
	double get_relative_accuracy() const
	{
		return my_relative_accuracy;
	}

	/** {@inherit_doc} */
	//override
	double get_absolute_accuracy() const
	{
		return my_absolute_accuracy;
	}

	/** {@inherit_doc} */
	//override
	int get_minimal_iteration_count() const
	{
		return my_minimal_iteration_count;
	}

	/** {@inherit_doc} */
	//override
	int get_maximal_iteration_count()
	{
		return my_iterations.get_maximal_count();
	}

	/** {@inherit_doc} */
	//override
	int get_evaluations() const
	{
		return my_evaluations.get_count();
	}

	/** {@inherit_doc} */
	//override
	int get_iterations()
	{
		return my_iterations.get_count();
	}

	/** {@inherit_doc} */
	//override
	T integrate(const int& max_eval, const Calculus_Field_Univariate_Function<T>& f, const T& lower, const T& upper)
	{
		// Initialization.
		setup(max_eval, f, lower, upper);

		// Perform computation.
		return do_integrate();
	}
};