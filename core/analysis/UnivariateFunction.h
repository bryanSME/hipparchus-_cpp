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
  //package org.hipparchus.analysis;

  /**
   * An interface representing a univariate real function.
   * <p>
   * When a <em>user-defined</em> function encounters an error during
   * evaluation, the {@link #valuestatic_cast<double>( value} method should throw a
   * <em>user-defined</em> unchecked exception.</p>
   * <p>
   * The following code excerpt shows the recommended way to do that using
   * a root solver as an example, but the same construct is applicable to
   * ODE integrators or optimizers.</p>
   *
   * <pre>
   * private static class Local_Exception extends Runtime_Exception
   {
   *     // The x value that caused the problem.
   *     private const double x;
   *
   *     public Local_Exception(double x)
   {
   *         this.x = x;
   *     }
   *
   *     public double get_x()
   {
   *         return x;
   *     }
   * }
   *
   * private static class My_Function : Univariate_Function
   {
   *     public double value(double x)
   {
   *         double y = huge_formula(x);
   *         if (something_bad_happens)
   {
   *           throw Local_Exception(x);
   *         }
   *         return y;
   *     }
   * }
   *
   * public void compute()
   {
   *     try
   {
   *         solver.solve(max_eval, My_Function(a, b, c), min, max);
   *     }
  catch (Local_Exception le)
   {
   *         // Retrieve the x value.
   *     }
   * }
   * </pre>
   *
   * As shown, the exception is local to the user's code and it is guaranteed
   * that Hipparchus will not catch it.
   *
   */
class Univariate_Function
{
	/**
	 * Compute the value of the function.
	 *
	 * @param x Point at which the function value should be computed.
	 * @return the value of the function.
	 * @Illegal_Argument_Exception when the activated method itself can
	 * ascertain that a precondition, specified in the API expressed at the
	 * level of the activated method, has been violated.
	 * When Hipparchus an {@code Illegal_Argument_Exception}, it is
	 * usually the consequence of checking the actual parameters passed to
	 * the method.
	 */
	virtual double value(const double& x) = 0;
};