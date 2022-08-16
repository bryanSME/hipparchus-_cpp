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
   * Implements the <em>Illinois</em> method for root-finding (approximating
   * a zero of a univariate real function). It is a modified
   * {@link Regula_Falsi_Solver <em>Regula Falsi</em>} method.
   *
   * <p>Like the <em>Regula Falsi</em> method, convergence is guaranteed by
   * maintaining a bracketed solution. The <em>Illinois</em> method however, * should converge much faster than the original <em>Regula Falsi</em>
   * method. Furthermore, this implementation of the <em>Illinois</em> method
   * should not suffer from the same implementation issues as the <em>Regula
   * Falsi</em> method, which may fail to convergence in certain cases.</p>
   *
   * <p>The <em>Illinois</em> method assumes that the function is continuous, * but not necessarily smooth.</p>
   *
   * <p>Implementation based on the following article: M. Dowell and P. Jarratt, * <em>A modified regula falsi method for computing the root of an
   * equation</em>, BIT Numerical Mathematics, volume 11, number 2, * pages 168-174, Springer, 1971.</p>
   *
   */
class Illinois_Solver : public Base_Secant_Solver
{
public:
	/** Construct a solver with default accuracy (1e-6). */
	Illinois_Solver()
	{
		super(DEFAULT_ABSOLUTE_ACCURACY, Method::ILLINOIS);
	}

	/**
	 * Construct a solver.
	 *
	 * @param absolute_accuracy Absolute accuracy.
	 */
	Illinois_Solver(const double& absolute_accuracy)
	{
		super(absolute_accuracy, Method::ILLINOIS);
	}

	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 */
	Illinois_Solver(const double& relative_accuracy, const double& absolute_accuracy)
	{
		super(relative_accuracy, absolute_accuracy, Method::ILLINOIS);
	}

	/**
	 * Construct a solver.
	 *
	 * @param relative_accuracy Relative accuracy.
	 * @param absolute_accuracy Absolute accuracy.
	 * @param function_value_accuracy Maximum function value error.
	 */
	Illinois_Solver(const double& relative_accuracy, const double& absolute_accuracy, const double& function_value_accuracy)
	{
		super(relative_accuracy, absolute_accuracy, function_value_accuracy, Method::ILLINOIS);
	}
};