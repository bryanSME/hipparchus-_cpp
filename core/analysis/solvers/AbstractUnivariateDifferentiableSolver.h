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

  //import org.hipparchus.analysis.differentiation.DS_Factory;
  //import org.hipparchus.analysis.differentiation.Derivative_Structure;
  //import org.hipparchus.analysis.differentiation.Univariate_Differentiable_Function;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
#include "BaseAbstractUnivariateSolver.hpp"
#include "UnivariateDifferentiableSolver.h"
#include "../differentiation/UnivariateDifferentiableFunction.h"
#include "../differentiation/DSFactory.h"
#include "../differentiation/DerivativeStructure.h"

/**
 * Provide a default implementation for several functions useful to generic
 * solvers.
 *
 */
class AbstractUnivariate_Differentiable_Solver
	:
	public Base_Abstract_Univariate_Solver<Univariate_Differentiable_Function>,
	public Univariate_Differentiable_Solver
{
private:

	/** Function to solve. */
	Univariate_Differentiable_Function my_function;

	/** Factory for Derivative_Structure instances. */
	const DS_Factory my_factory;

protected:
	/**
	 * Construct a solver with given absolute accuracy.
	 *
	 * @param absolute_accuracy Maximum absolute error.
	 */
	AbstractUnivariate_Differentiable_Solver(const double& absolute_accuracy) : my_factory{ DS_Factory(1, 1) }
	{
		throw std::exception("not implemented");
		//super(absolute_accuracy);
	}

	/**
	 * Construct a solver with given accuracies.
	 *
	 * @param relative_accuracy Maximum relative error.
	 * @param absolute_accuracy Maximum absolute error.
	 * @param function_value_accuracy Maximum function value error.
	 */
	AbstractUnivariate_Differentiable_Solver(const double& relative_accuracy, const double& absolute_accuracy, const double& function_value_accuracy)
		:
		my_factory{ DS_Factory(1, 1) }
	{
		throw std::exception("not implemented");
		//super(relative_accuracy, absolute_accuracy, function_value_accuracy);
	}

	/**
	 * Compute the objective function value.
	 *
	 * @param point Point at which the objective function must be evaluated.
	 * @return the objective function value and derivative at specified point.
	 * @Math_Illegal_State_Exception
	 * if the maximal number of evaluations is exceeded.
	 */
	Derivative_Structure compute_objective_value_and_derivative(const double& point)
	{
		throw std::exception("not implemented");
		//increment_evaluation_count();
		//return my_function.value(my_factory.variable(0, point));
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	void setup(const int& max_eval, const Univariate_Differentiable_Function& f, const double& min, const double& max, const double& start_value)
	{
		throw std::exception("not implemented");
		//super.setup(max_eval, f, min, max, start_value);
		//my_function = f;
	}
};