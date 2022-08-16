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

#include "BaseSecantSolver.h"

  /**
   * Implements the <em>Pegasus</em> method for root-finding (approximating
   * a zero of a univariate real function). It is a modified
   * {@link Regula_Falsi_Solver <em>Regula Falsi</em>} method.
   *
   * <p>Like the <em>Regula Falsi</em> method, convergence is guaranteed by
   * maintaining a bracketed solution. The <em>Pegasus</em> method however, * should converge much faster than the original <em>Regula Falsi</em>
   * method. Furthermore, this implementation of the <em>Pegasus</em> method
   * should not suffer from the same implementation issues as the <em>Regula
   * Falsi</em> method, which may fail to convergence in certain cases. Also, * the <em>Pegasus</em> method should converge faster than the
   * {@link Illinois_Solver <em>Illinois</em>} method, another <em>Regula
   * Falsi</em>-based method.</p>
   *
   * <p>The <em>Pegasus</em> method assumes that the function is continuous, * but not necessarily smooth.</p>
   *
   * <p>Implementation based on the following article: M. Dowell and P. Jarratt, * <em>The "Pegasus" method for computing the root of an equation</em>, * BIT Numerical Mathematics, volume 12, number 4, pages 503-508, Springer, * 1972.</p>
   *
   */
class Pegasus_Solver : public Base_Secant_Solver
{
public:
	/** Construct a solver with default accuracy (1e-6). */
	Pegasus_Solver()
	{
		super(DEFAULT_ABSOLUTE_ACCURACY, Method::PEGASUS);
	}

	/**
	 * Construct a solver.
	 *
	 * @param absolute_accuracy Absolute accuracy.
	 */
	Pegasus_Solver(const double& absolute_accuracy)
	{
		super(absolute_accuracy, Method::PEGASUS);
	}

	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 */
	Pegasus_Solver(const double& relative_accuracy, const double& absolute_accuracy)
	{
		super(relative_accuracy, absolute_accuracy, Method::PEGASUS);
	}

	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 * @param function_value_accuracy Maximum function value error.
	 */
	Pegasus_Solver(const double& relative_accuracy, const double& absolute_accuracy, const double function_value_accuracy)
	{
		super(relative_accuracy, absolute_accuracy, function_value_accuracy, Method::PEGASUS);
	}
};