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

#include <vector>
#include "PolynomialSolver.h"
#include "BaseAbstractUnivariateSolver.hpp"
#include "PolynomialSolver.h"
#include "../polynomials/PolynomialFunction.hpp"
#include "BaseAbstractUnivariateSolver.hpp"
  //import org.hipparchus.analysis.polynomials.Polynomial_Function;

  /**
   * Base class for solvers.
   *
   */
class Abstract_Polynomial_Solver : public Base_Abstract_Univariate_Solver<Polynomial_Function>, public Polynomial_Solver
{
private:
	/** Function. */
	Polynomial_Function polynomial_function;

protected:
	/**
	 * Construct a solver with given absolute accuracy.
	 *
	 * @param absolute_accuracy Maximum absolute error.
	 */
	Abstract_Polynomial_Solver(const double& absolute_accuracy)
	{
		super(absolute_accuracy);
	}
	/**
	 * Construct a solver with given accuracies.
	 *
	 * @param relative_accuracy Maximum relative error.
	 * @param absolute_accuracy Maximum absolute error.
	 */
	Abstract_Polynomial_Solver(const double& relative_accuracy, const double& absolute_accuracy)
	{
		super(relative_accuracy, absolute_accuracy);
	}
	/**
	 * Construct a solver with given accuracies.
	 *
	 * @param relative_accuracy Maximum relative error.
	 * @param absolute_accuracy Maximum absolute error.
	 * @param function_value_accuracy Maximum function value error.
	 */
	Abstract_Polynomial_Solver(const double& relative_accuracy, const double& absolute_accuracy, const double& function_value_accuracy)
	{
		super(relative_accuracy, absolute_accuracy, function_value_accuracy);
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	void setup(const int& max_eval, Polynomial_Function& f, const double& min, const double& max, const double& start_value)
	{
		super.setup(max_eval, f, min, max, start_value);
		polynomial_function = f;
	}

	/**
	 * @return the coefficients of the polynomial function.
	 */
	std::vector<double> get_coefficients()
	{
		return polynomial_function.get_coefficients();
	}
};