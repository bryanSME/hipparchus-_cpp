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
  //package org.hipparchus.migration.ode;

  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.ode.Ordinary_Differential_Equation;

  /** Interface expanding {@link First_Order_Differential_Equations first order
   *  differential equations} in order to compute exactly the main state jacobian
   *  matrix for {@link Jacobian_Matrices partial derivatives equations}.
   * @deprecated as of 1.0, replaced with {@link org.hipparchus.ode.ODE_Jacobians_Provider}
   */
@Deprecated
class Main_State_Jacobian_Provider extends Ordinary_Differential_Equation
{
	/** Compute the jacobian matrix of ODE with respect to main state.
	 * @param t current value of the independent <I>time</I> variable
	 * @param y array containing the current value of the main state vector
	 * @param y_dot array containing the current value of the time derivative of the main state vector
	 * @return Jacobian matrix of the ODE w.r.t. the main state vector
	 * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
	 * @exception  if arrays dimensions do not match equations settings
	 */
	std::vector<std::vector<double>> compute_main_state_jacobian(double t, std::vector<double> y, std::vector<double> y_dot)
		;
}
