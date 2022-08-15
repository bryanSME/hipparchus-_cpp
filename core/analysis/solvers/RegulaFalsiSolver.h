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
 * Implements the <em>Regula Falsi</em> or <em>False position</em> method for
 * root-finding (approximating a zero of a univariate real function). It is a
 * modified {@link Secant_Solver <em>Secant</em>} method.
 *
 * <p>The <em>Regula Falsi</em> method is included for completeness, for
 * testing purposes, for educational purposes, for comparison to other
 * algorithms, etc. It is however <strong>not</strong> intended to be used
 * for actual problems, as one of the bounds often remains fixed, resulting
 * in very slow convergence. Instead, one of the well-known modified
 * <em>Regula Falsi</em> algorithms can be used ({@link Illinois_Solver
 * <em>Illinois</em>} or {@link Pegasus_Solver <em>Pegasus</em>}). These two
 * algorithms solve the fundamental issues of the original <em>Regula
 * Falsi</em> algorithm, and greatly out-performs it for most, if not all, * (practical) functions.
 *
 * <p>Unlike the <em>Secant</em> method, the <em>Regula Falsi</em> guarantees
 * convergence, by maintaining a bracketed solution. Note however, that due to
 * the finite/limited precision of Java's {@link Double double} type, which is
 * used in this implementation, the algorithm may get stuck in a situation
 * where it no longer makes any progress. Such cases are detected and result
 * in a {@code Math_Illegal_State_Exception} exception being thrown. In other words, * the algorithm theoretically guarantees convergence, but the implementation
 * does not.</p>
 *
 * <p>The <em>Regula Falsi</em> method assumes that the function is continuous, * but not necessarily smooth.</p>
 *
 * <p>Implementation based on the following article: M. Dowell and P. Jarratt, * <em>A modified regula falsi method for computing the root of an
 * equation</em>, BIT Numerical Mathematics, volume 11, number 2, * pages 168-174, Springer, 1971.</p>
 *
 */
class Regula_Falsi_Solver : public Base_Secant_Solver
{
public:
    /** Construct a solver with default accuracy (1e-6). */
    Regula_Falsi_Solver()
    {
        super(DEFAULT_ABSOLUTE_ACCURACY, Method::REGULA_FALSI);
    }

    /**
     * Construct a solver.
     *
     * @param absolute_accuracy Absolute accuracy.
     */
    Regula_Falsi_Solver(const double& absolute_accuracy)
    {
        super(absolute_accuracy, Method::REGULA_FALSI);
    }

    /**
     * Construct a solver.
     *
     * @param relative_accuracy Relative accuracy.
     * @param absolute_accuracy Absolute accuracy.
     */
    Regula_Falsi_Solver(const double& relative_accuracy, const double& absolute_accuracy)
    {
        super(relative_accuracy, absolute_accuracy, Method::REGULA_FALSI);
    }

    /**
     * Construct a solver.
     *
     * @param relative_accuracy Relative accuracy.
     * @param absolute_accuracy Absolute accuracy.
     * @param function_value_accuracy Maximum function value error.
     */
    Regula_Falsi_Solver(const double& relative_accuracy, const double& absolute_accuracy, const double function_value_accuracy)
    {
        super(relative_accuracy, absolute_accuracy, function_value_accuracy, Method::REGULA_FALSI);
    }
};