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

  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.ode.Dense_Output_Model;

  /**
   * This class stores all information provided by an ODE integrator
   * during the integration process and build a continuous model of the
   * solution from this.
   *
   * <p>This class act as a step handler from the integrator point of
   * view. It is called iteratively during the integration process and
   * stores a copy of all steps information in a sorted collection for
   * later use. Once the integration process is over, the user can use
   * the {@link #set_interpolated_time set_interpolated_time} and {@link
   * #get_interpolated_state get_interpolated_state} to retrieve this
   * information at any time. It is important to wait for the
   * integration to be over before attempting to call {@link
   * #set_interpolated_time set_interpolated_time} because some internal
   * variables are set only once the last step has been handled.</p>
   *
   * <p>This is useful for example if the main loop of the user
   * application should remain independent from the integration process
   * or if one needs to mimic the behaviour of an analytical model
   * despite a numerical model is used (i.e. one needs the ability to
   * get the model value at any time or to navigate through the
   * data).</p>
   *
   * <p>If problem modeling is done with several separate
   * integration phases for contiguous intervals, the same
   * Continuous_Output_Model can be used as step handler for all
   * integration phases as long as they are performed in order and in
   * the same direction. As an example, one can extrapolate the
   * trajectory of a satellite with one model (i.e. one set of
   * differential equations) up to the beginning of a maneuver, use
   * another more complex model including thrusters modeling and
   * accurate attitude control during the maneuver, and revert to the
   * first model after the end of the maneuver. If the same continuous
   * output model handles the steps of all integration phases, the user
   * do not need to bother when the maneuver begins or ends, he has all
   * the data available in a transparent manner.</p>
   *
   * <p>An important feature of this class is that it : the
   * <code>Serializable</code> interface. This means that the result of
   * an integration can be serialized and reused later (if stored into a
   * persistent medium like a filesystem or a database) or elsewhere (if
   * sent to another application). Only the result of the integration is
   * stored, there is no reference to the integrated problem by
   * itself.</p>
   *
   * <p>One should be aware that the amount of data stored in a
   * Continuous_Output_Model instance can be important if the state vector
   * is large, if the integration interval is long or if the steps are
   * small (which can result from small tolerance settings in {@link
   * org.hipparchus.ode.nonstiff.Adaptive_Stepsize_Integrator adaptive
   * step size integrators}).</p>
   *
   * @see org.hipparchus.migration.ode.sampling.Step_Handler
   * @see org.hipparchus.ode.sampling.ODE_StateInterpolator
   * @deprecated as of 1.0, replaced with {@link Dense_Output_Model}
   */
@Deprecated
class Continuous_Output_Model extends Dense_Output_Model
{
	/** Serializable version identifier */
	20160403L;

	/** Interpolation time. */
	private double interpolated_time;

	/** Set the time of the interpolated point.
	 * <p>This method should <strong>not</strong> be called before the
	 * integration is over because some internal variables are set only
	 * once the last step has been handled.</p>
	 * <p>Setting the time outside of the integration interval is now
	 * allowed, but should be used with care since the accuracy of the
	 * interpolator will probably be very poor far from this interval.
	 * This allowance has been added to simplify implementation of search
	 * algorithms near the interval endpoints.</p>
	 * <p>Note that each time this method is called, the internal arrays
	 * returned in {@link #get_interpolated_state()}, {@link
	 * #get_interpolated_derivatives()} and {@link #get_interpolated_secondary_statestatic_cast<int>(}
	 * <em>will</em> be overwritten. So if their content must be preserved
	 * across several calls, user must copy them.</p>
	 * @param time time of the interpolated point
	 * @see #get_interpolated_state()
	 * @see #get_interpolated_derivatives()
	 * @see #get_interpolated_secondary_statestatic_cast<int>(
	 */
	public void set_interpolated_time(const double time)
	{
		this.interpolated_time = time;
	}

	/**
	 * Get the time of the interpolated point.
	 * If {@link #set_interpolated_time} has not been called, it returns
	 * the const integration time.
	 * @return interpolation point time
	 */
	public double get_interpolated_time()
	{
		return interpolated_time;
	}

	/**
	 * Get the state vector of the interpolated point.
	 * <p>The returned vector is a reference to a reused array, so
	 * it should not be modified and it should be copied if it needs
	 * to be preserved across several calls to the associated
	 * {@link #set_interpolated_timestatic_cast<double>(} method.</p>
	 * @return state vector at time {@link #get_interpolated_time}
	 * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
	 * @see #set_interpolated_timestatic_cast<double>(
	 * @see #get_interpolated_derivatives()
	 * @see #get_interpolated_secondary_statestatic_cast<int>(
	 * @see #get_interpolated_secondary_derivativesstatic_cast<int>(
	 */
	public std::vector<double> get_interpolated_state() Math_Illegal_State_Exception
	{
		return get_interpolated_state(get_interpolated_time()).get_primary_state();
	}

	/**
	 * Get the derivatives of the state vector of the interpolated point.
	 * <p>The returned vector is a reference to a reused array, so
	 * it should not be modified and it should be copied if it needs
	 * to be preserved across several calls to the associated
	 * {@link #set_interpolated_timestatic_cast<double>(} method.</p>
	 * @return derivatives of the state vector at time {@link #get_interpolated_time}
	 * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
	 * @see #set_interpolated_timestatic_cast<double>(
	 * @see #get_interpolated_state()
	 * @see #get_interpolated_secondary_statestatic_cast<int>(
	 * @see #get_interpolated_secondary_derivativesstatic_cast<int>(
	 */
	public std::vector<double> get_interpolated_derivatives() Math_Illegal_State_Exception
	{
		return get_interpolated_state(get_interpolated_time()).get_primary_derivative();
	}

	/** Get the interpolated secondary state corresponding to the secondary equations.
	 * <p>The returned vector is a reference to a reused array, so
	 * it should not be modified and it should be copied if it needs
	 * to be preserved across several calls to the associated
	 * {@link #set_interpolated_timestatic_cast<double>(} method.</p>
	 * @param secondary_state_index index of the secondary set, as returned by {@link
	 * org.hipparchus.ode.Expandable_ODE#add_secondary_equations(org.hipparchus.ode.Secondary_ODE)
	 * Expandable_ODE.add_secondary_equations(secondary)}
	 * @return interpolated secondary state at the current interpolation date
	 * @see #set_interpolated_timestatic_cast<double>(
	 * @see #get_interpolated_state()
	 * @see #get_interpolated_derivatives()
	 * @see #get_interpolated_secondary_derivativesstatic_cast<int>(
	 * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
	 */
	public std::vector<double> get_interpolated_secondary_state(const int secondary_state_index)
		Math_Illegal_State_Exception
	{
		return get_interpolated_state(get_interpolated_time()).get_secondary_state(secondary_state_index);
	}

	/** Get the interpolated secondary derivatives corresponding to the secondary equations.
	 * <p>The returned vector is a reference to a reused array, so
	 * it should not be modified and it should be copied if it needs
	 * to be preserved across several calls to the associated
	 * {@link #set_interpolated_timestatic_cast<double>(} method.</p>
	 * @param secondary_state_index index of the secondary set, as returned by {@link
	 * org.hipparchus.ode.Expandable_ODE#add_secondary_equations(org.hipparchus.ode.Secondary_ODE)
	 * Expandable_ODE.add_secondary_equations(secondary)}
	 * @return interpolated secondary derivatives at the current interpolation date
	 * @see #set_interpolated_timestatic_cast<double>(
	 * @see #get_interpolated_state()
	 * @see #get_interpolated_derivatives()
	 * @see #get_interpolated_secondary_statestatic_cast<int>(
	 * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
	 */
	public std::vector<double> get_interpolated_secondary_derivatives(const int secondary_state_index)
		Math_Illegal_State_Exception
	{
		return get_interpolated_state(get_interpolated_time()).get_secondary_derivative(secondary_state_index);
	}
}
