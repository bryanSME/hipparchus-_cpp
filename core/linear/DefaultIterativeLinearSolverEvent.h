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
  //import org.hipparchus.exception.Math_Runtime_Exception;

  /**
   * A default concrete implementation of the virtual class
   * {@link Iterative_Linear_SolverEvent}.
   *
   */
class Default_iterativeLinearSolverEvent : public Iterative_Linear_SolverEvent
{
private:

	/** The right-hand side vector. */
	const Real_Vector my_b;

	/** The current estimate of the residual. */
	const Real_Vector my_r;

	/** The current estimate of the norm of the residual. */
	const double my_rnorm;

	/** The current estimate of the solution. */
	const Real_Vector my_x;

public:
	/**
	 * Creates a instance of this class. This implementation does
	 * <em>not</em> deep copy the specified vectors {@code x}, {@code b}, * {@code r}. Therefore the user must make sure that these vectors are
	 * either unmodifiable views or deep copies of the same vectors actually
	 * used by the {@code source}. Failure to do so may compromise subsequent
	 * iterations of the {@code source}. If the residual vector {@code r} is
	 * {@code NULL}, then {@link #get_residual()} a
	 * {@link Math_Runtime_Exception}, and
	 * {@link #provides_residual()} returns {@code false}.
	 *
	 * @param source the iterative solver which fired this event
	 * @param iterations the number of iterations performed at the time
	 * {@code this} event is created
	 * @param x the current estimate of the solution
	 * @param b the right-hand side vector
	 * @param r the current estimate of the residual (can be {@code NULL})
	 * @param rnorm the norm of the current estimate of the residual
	 */
	Default_iterativeLinearSolverEvent(const Object& source, const int& iterations, const Real_Vector& x, const Real_Vector& b, const Real_Vector& r, const double& rnorm)
		: 
		my_x{ x },
		my_b{ b },
		my_r{ r },
		my_rnorm{ rnorm }
	{
		Iterative_Linear_SolverEvent(source, iterations);
	}

	/**
	 * Creates a instance of this class. This implementation does
	 * <em>not</em> deep copy the specified vectors {@code x}, {@code b}.
	 * Therefore the user must make sure that these vectors are either
	 * unmodifiable views or deep copies of the same vectors actually used by
	 * the {@code source}. Failure to do so may compromise subsequent iterations
	 * of the {@code source}. Callling {@link #get_residual()} on instances
	 * returned by this constructor a
	 * {@link Math_Runtime_Exception}, while
	 * {@link #provides_residual()} returns {@code false}.
	 *
	 * @param source the iterative solver which fired this event
	 * @param iterations the number of iterations performed at the time
	 * {@code this} event is created
	 * @param x the current estimate of the solution
	 * @param b the right-hand side vector
	 * @param rnorm the norm of the current estimate of the residual
	 */
	Default_iterativeLinearSolverEvent(const Object& source, const int& iterations, const Real_Vector& x, const Real_Vector& b, const double& rnorm)
		:
		my_x{ x },
		my_b{ b },
		my_r{ NULL },
		my_rnorm{ rnorm }
	{
		Iterative_Linear_SolverEvent(source, iterations);
	}

	/** {@inherit_doc} */
	//override
	double get_norm_of_residual()
	{
		return rnorm;
	}

	/**
	 * {@inherit_doc}
	 *
	 * This implementation a {@link Math_Runtime_Exception}
	 * if no residual vector {@code r} was provided at construction time.
	 */
	 //override
	Real_Vector get_residual() const
	{
		return my_r;
	}

	/** {@inherit_doc} */
	//override
	Real_Vector get_right_hand_side_vector() const
	{
		return my_b;
	}

	/** {@inherit_doc} */
	//override
	Real_Vector get_solution() const
	{
		return my_x;
	}

	/**
	 * {@inherit_doc}
	 *
	 * This implementation returns {@code true} if a non-{@code NULL} value was
	 * specified for the residual vector {@code r} at construction time.
	 *
	 * @return {@code true} if {@code r != NULL}
	 */
	 //override
	bool provides_residual()
	{
		return r != NULL;
	}
};