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

  //import org.hipparchus.ode.Second_Order_ODE;

  /** This interface represents a second order differential equations set.

   * <p>This interface should be implemented by all real second order
   * differential equation problems before they can be handled by the
   * {@link org.hipparchus.ode.First_Order_Converter converter to first order}.</p>
   *
   * <p>A second order differential equations problem, as seen by an
   * integrator is the second time derivative <code>d2_y/dt^2</code> of a
   * state vector <code>Y</code>, both being one dimensional
   * arrays. From the integrator point of view, this derivative depends
   * only on the current time <code>t</code>, on the state vector
   * <code>Y</code> and on the first time derivative of the state
   * vector.</p>
   *
   * <p>For real problems, the derivative depends also on parameters
   * that do not belong to the state vector (dynamical model constants
   * for example). These constants are completely outside of the scope
   * of this interface, the classes that implement it are allowed to
   * handle them as they want.</p>
   * @deprecated as of 1.0, replaced with {@link Second_Order_ODE}
   */
@Deprecated
class Second_Order_Differential_Equations extends Second_Order_ODE
{
	/** {@inherit_doc} */
	//override
	default std::vector<double> compute_second_derivatives(double t, std::vector<double> y, std::vector<double> y_dot)
	{
		const std::vector<double> y_d_dot = std::vector<double>(y.size()];
		compute_second_derivatives(t, y, y_dot, y_d_dot);
		return y_d_dot;
	}

	/** Get the current time derivative of the state vector.
	 * @param t current value of the independent <I>time</I> variable
	 * @param y array containing the current value of the state vector
	 * @param y_dot array containing the current value of the first derivative
	 * of the state vector
	 * @param y_d_dot placeholder array where to put the second time derivative
	 * of the state vector
	 */
	void compute_second_derivatives(double t, std::vector<double> y, std::vector<double> y_dot, std::vector<double> y_d_dot);
}
