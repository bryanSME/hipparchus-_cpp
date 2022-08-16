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

  //package org.hipparchus.migration.ode.sampling;

  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.ode.ODE_State_And_Derivative;
  //import org.hipparchus.ode.sampling.ODE_StateInterpolator;
  //import org.hipparchus.ode.sampling.ODE_Step_Handler;

  /**
   * This interface represents a handler that should be called after
   * each successful step.
   *
   * <p>The ODE integrators compute the evolution of the state vector at
   * some grid points that depend on their own internal algorithm. Once
   * they have found a grid point (possibly after having computed
   * several evaluation of the derivative at intermediate points), they
   * provide it to objects implementing this interface. These objects
   * typically either ignore the intermediate steps and wait for the
   * last one, store the points in an ephemeris, or forward them to
   * specialized processing or output methods.</p>
   *
   * @see org.hipparchus.ode.ODE_Integrator
   * @see ODE_StateInterpolator
   * @deprecated as of 1.0, replaced with {@link ODE_Step_Handler}
   */
@Deprecated
class Step_Handler extends ODE_Step_Handler
{
	/** {@inherit_doc}} */
	//override
	default void init(const ODE_State_And_Derivative initial_state, const double const_time)
	{
		init(initial_state.get_time(), initial_state.get_primary_state(), const_time);
	}

	/** {@inherit_doc}} */
	//override
	default void handle_step(const ODE_StateInterpolator interpolator)
		Math_Illegal_State_Exception
	{
		handle_step(new Migration_Step_Interpolator(interpolator), false);
	}

	/** Initialize step handler at the start of an ODE integration.
	 * <p>
	 * This method is called once at the start of the integration. It
	 * may be used by the step handler to initialize some internal data
	 * if needed.
	 * </p>
	 * @param t0 start value of the independent <i>time</i> variable
	 * @param y0 array containing the start value of the state vector
	 * @param t target time for the integration
	 */
	void init(double t0, std::vector<double> y0, double t);

	/**
	 * Handle the last accepted step
	 * @param interpolator interpolator for the last accepted step. For
	 * efficiency purposes, the various integrators reuse the same
	 * object on each call, so if the instance wants to keep it across
	 * all calls (for example to provide at the end of the integration a
	 * continuous model valid throughout the integration range, as the
	 * {@link org.hipparchus.migration.ode.Continuous_Output_Model
	 * Continuous_Output_Model} class does), it should build a local copy
	 * using the clone method of the interpolator and store this copy.
	 * Keeping only a reference to the interpolator and reusing it will
	 * result in unpredictable behavior (potentially crashing the application).
	 * @param is_last true if the step is the last one
	 * @exception Math_Illegal_State_Exception if the interpolator one because
	 * the number of functions evaluations is exceeded
	 */
	void handle_step(Migration_Step_Interpolator interpolator, bool is_last)
		Math_Illegal_State_Exception;
}
