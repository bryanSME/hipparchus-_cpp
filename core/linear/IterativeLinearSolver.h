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
  //package org.hipparchus.linear;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Iteration_Manager;
  //import org.hipparchus.util.Math_Utils;
#include "../util/IterationManager.h"
#include "RealVector.h"

/**
 * This virtual class defines an iterative solver for the linear system A
 * &middot; x = b. In what follows, the <em>residual</em> r is defined as r = b
 * - A &middot; x, where A is the linear operator of the linear system, b is the
 * right-hand side vector, and x the current estimate of the solution.
 *
 */
class Iterative_Linear_Solver
{
private:
	/** The object in charge of managing the iterations. */
	const Iteration_Manager my_manager;

protected:
	/**
 * Performs all dimension checks on the parameters of
 * {@link #solve(Real_Linear_Operator, Real_Vector, Real_Vector) solve} and
 * {@link #solve_in_place(Real_Linear_Operator, Real_Vector, Real_Vector) solve_in_place}, * and an exception if one of the checks fails.
 *
 * @param a the linear operator A of the system
 * @param b the right-hand side vector
 * @param x0 the initial guess of the solution
 * @ if one of the parameters is {@code NULL}
 * @ if {@code a} is not square
 * @ if {@code b} or {@code x0} have
 * dimensions inconsistent with {@code a}
 */
	protected static void check_parameters(const Real_Linear_Operator& a, const Real_Vector& b, const Real_Vector& x0)
	{
		throw std::exception("Not Implemented");
		//Math_Utils::check_not_null(a);
		//Math_Utils::check_not_null(b);
		//Math_Utils::check_not_null(x0);
		/*if (a.get_row_dimension() != a.get_column_dimension())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_OPERATOR, a.get_row_dimension(), a.get_column_dimension());
		}
		if (b.get_dimension() != a.get_row_dimension())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_dimension(), a.get_row_dimension());
		}
		if (x0.get_dimension() != a.get_column_dimension())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, x0.get_dimension(), a.get_column_dimension());
		}*/
	}

public:
	/**
	 * Creates a instance of this class, with default iteration manager.
	 *
	 * @param max_iterations the maximum number of iterations
	 */
	Iterative_Linear_Solver(const int& max_iterations)
	{
		my_manager = Iteration_Manager(max_iterations);
	}

	/**
	 * Creates a instance of this class, with custom iteration manager.
	 *
	 * @param manager the custom iteration manager
	 * @ if {@code manager} is {@code NULL}
	 */
	Iterative_Linear_Solver(const Iteration_Manager& manager)
	{
		//Math_Utils::check_not_null(manager);
		my_manager = manager;
	}

	/**
	 * Returns the iteration manager attached to this solver.
	 *
	 * @return the manager
	 */
	Iteration_Manager get_iteration_manager() const
	{
		return my_manager;
	}

	/**
	 * Returns an estimate of the solution to the linear system A &middot; x =
	 * b.
	 *
	 * @param a the linear operator A of the system
	 * @param b the right-hand side vector
	 * @return a vector containing the solution
	 * @ if one of the parameters is {@code NULL}
	 * @ if {@code a} is not square
	 * @ if {@code b} has dimensions
	 * inconsistent with {@code a}
	 * @Math_Illegal_State_Exception at exhaustion of the iteration count, * unless a custom
	 * {@link org.hipparchus.util.Incrementor.Max_Count_Exceeded_Callback callback}
	 * has been set at construction of the {@link Iteration_Manager}
	 */
	Real_Vector solve(const Real_Linear_Operator& a, const Real_Vector& b)
	{
		//Math_Utils::check_not_null(a);
		const Real_Vector x = Array_Real_Vector(a.get_column_dimension());
		x.set(0.);
		return solve_in_place(a, b, x);
	}

	/**
	 * Returns an estimate of the solution to the linear system A &middot; x =
	 * b.
	 *
	 * @param a the linear operator A of the system
	 * @param b the right-hand side vector
	 * @param x0 the initial guess of the solution
	 * @return a vector containing the solution
	 * @ if one of the parameters is {@code NULL}
	 * @ if {@code a} is not square
	 * @ if {@code b} or {@code x0} have
	 * dimensions inconsistent with {@code a}
	 * @Math_Illegal_State_Exception at exhaustion of the iteration count, * unless a custom
	 * {@link org.hipparchus.util.Incrementor.Max_Count_Exceeded_Callback callback}
	 * has been set at construction of the {@link Iteration_Manager}
	 */
	Real_Vector solve(const Real_Linear_Operator& a, const Real_Vector& b, const Real_Vector& x0)
	{
		//Math_Utils::check_not_null(x0);
		return solve_in_place(a, b, x0.copy());
	}

	/**
	 * Returns an estimate of the solution to the linear system A &middot; x =
	 * b. The solution is computed in-place (initial guess is modified).
	 *
	 * @param a the linear operator A of the system
	 * @param b the right-hand side vector
	 * @param x0 initial guess of the solution
	 * @return a reference to {@code x0} (shallow copy) updated with the
	 * solution
	 * @ if one of the parameters is {@code NULL}
	 * @ if {@code a} is not square
	 * @ if {@code b} or {@code x0} have
	 * dimensions inconsistent with {@code a}
	 * @Math_Illegal_State_Exception at exhaustion of the iteration count, * unless a custom
	 * {@link org.hipparchus.util.Incrementor.Max_Count_Exceeded_Callback callback}
	 * has been set at construction of the {@link Iteration_Manager}
	 */
	virtual Real_Vector solve_in_place(const Real_Linear_Operator& a, const Real_Vector& b, const Real_Vector& x0);
};