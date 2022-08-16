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
  //import org.hipparchus.ode.sampling.ODE_StateInterpolator;

  /** This interface represents an interpolator over the last step
   * during an ODE integration.
   *
   * <p>The various ODE integrators provide objects implementing this
   * interface to the step handlers. These objects are often custom
   * objects tightly bound to the integrator internal algorithms. The
   * handlers can use these objects to retrieve the state vector at
   * intermediate times between the previous and the current grid points
   * (this feature is often called dense output).</p>
   * <p>One important thing to note is that the step handlers may be so
   * tightly bound to the integrators that they often share some internal
   * state arrays. This imply that one should <em>never</em> use a direct
   * reference to a step interpolator outside of the step handler, either
   * for future use or for use in another thread. If such a need arise, the
   * step interpolator <em>must</em> be copied using the dedicated
   * {@link #copy()} method.
   * </p>
   *
   * @see org.hipparchus.ode.ODE_Integrator
   * @see Step_Handler
   * @deprecated as of 1.0, replaced with {@link ODE_StateInterpolator}
   */
@Deprecated
class Step_Interpolator extends ODE_StateInterpolator
{
	/**
	 * Get the previous grid point time.
	 * @return previous grid point time
	 * @deprecated as of 1.0, replaced with {@link #get_previous_state()}/{@link org.hipparchus.ode.ODE_State#get_time()}
	 */
	@Deprecated
		double get_previous_time();

	/**
	 * Get the current grid point time.
	 * @return current grid point time
	 * @deprecated as of 1.0, replaced with {@link #get_current_state()}/{@link org.hipparchus.ode.ODE_State#get_time()}
	 */
	@Deprecated
		double get_current_time();

	/**
	 * Get the time of the interpolated point.
	 * If {@link #set_interpolated_time} has not been called, it returns
	 * the current grid point time.
	 * @return interpolation point time
	 * @deprecated as of 1.0, replaced with {@link #get_interpolated_statestatic_cast<double>(}/{@link org.hipparchus.ode.ODE_State#get_time()}
	 */
	@Deprecated
		double get_interpolated_time();

	/**
	 * Set the time of the interpolated point.
	 * <p>Setting the time outside of the current step is now allowed, but
	 * should be used with care since the accuracy of the interpolator will
	 * probably be very poor far from this step. This allowance has been
	 * added to simplify implementation of search algorithms near the
	 * step endpoints.</p>
	 * <p>Setting the time changes the instance internal state. This includes
	 * the internal arrays returned in {@link #get_interpolated_state()}, * {@link #get_interpolated_derivatives()}, {@link
	 * #get_interpolated_secondary_statestatic_cast<int>(} and {@link
	 * #get_interpolated_secondary_derivativesstatic_cast<int>(}. So if their content must be preserved
	 * across several calls, user must copy them.</p>
	 * @param time time of the interpolated point
	 * @see #get_interpolated_state()
	 * @see #get_interpolated_derivatives()
	 * @see #get_interpolated_secondary_statestatic_cast<int>(
	 * @see #get_interpolated_secondary_derivativesstatic_cast<int>(
	 * @deprecated as of 1.0, replaced with {@link #get_interpolated_statestatic_cast<double>(}
	 */
	@Deprecated
		void set_interpolated_time(double time);

	/**
	 * Get the state vector of the interpolated point.
	 * <p>The returned vector is a reference to a reused array, so
	 * it should not be modified and it should be copied if it needs
	 * to be preserved across several calls to the associated
	 * {@link #set_interpolated_timestatic_cast<double>(} method.</p>
	 * @return state vector at time {@link #get_interpolated_time}
	 * @see #get_interpolated_derivatives()
	 * @see #get_interpolated_secondary_statestatic_cast<int>(
	 * @see #get_interpolated_secondary_derivativesstatic_cast<int>(
	 * @see #set_interpolated_timestatic_cast<double>(
	 * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
	 * @deprecated as of 1.0, replaced with {@link #get_interpolated_statestatic_cast<double>(}.{@link org.hipparchus.ode.ODE_State#get_primary_state()}
	 */
	@Deprecated
		std::vector<double> get_interpolated_state() Math_Illegal_State_Exception;

	/**
	 * Get the derivatives of the state vector of the interpolated point.
	 * <p>The returned vector is a reference to a reused array, so
	 * it should not be modified and it should be copied if it needs
	 * to be preserved across several calls to the associated
	 * {@link #set_interpolated_timestatic_cast<double>(} method.</p>
	 * @return derivatives of the state vector at time {@link #get_interpolated_time}
	 * @see #get_interpolated_state()
	 * @see #get_interpolated_secondary_statestatic_cast<int>(
	 * @see #get_interpolated_secondary_derivativesstatic_cast<int>(
	 * @see #set_interpolated_timestatic_cast<double>(
	 * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
	 * @deprecated as of 1.0, replaced with {@link #get_interpolated_statestatic_cast<double>(}.{@link org.hipparchus.ode.ODE_State_And_Derivative#get_primary_derivative()}
	 */
	@Deprecated
		std::vector<double> get_interpolated_derivatives() Math_Illegal_State_Exception;

	/** Get the interpolated secondary state corresponding to the secondary equations.
	 * <p>The returned vector is a reference to a reused array, so
	 * it should not be modified and it should be copied if it needs
	 * to be preserved across several calls to the associated
	 * {@link #set_interpolated_timestatic_cast<double>(} method.</p>
	 * @param index index of the secondary set, as returned by {@link
	 * org.hipparchus.ode.Expandable_ODE#add_secondary_equations(org.hipparchus.ode.Secondary_ODE)
	 * Expandable_ODE.add_secondary_equations(secondary)}
	 * @return interpolated secondary state at the current interpolation date
	 * @see #get_interpolated_state()
	 * @see #get_interpolated_derivatives()
	 * @see #get_interpolated_secondary_derivativesstatic_cast<int>(
	 * @see #set_interpolated_timestatic_cast<double>(
	 * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
	 * @deprecated as of 1.0, replaced with {@link #get_interpolated_statestatic_cast<double>(}.{@link org.hipparchus.ode.ODE_State#get_secondary_statestatic_cast<int>(}
	 */
	@Deprecated
		std::vector<double> get_interpolated_secondary_state(const int& index) Math_Illegal_State_Exception;

	/** Get the interpolated secondary derivatives corresponding to the secondary equations.
	 * <p>The returned vector is a reference to a reused array, so
	 * it should not be modified and it should be copied if it needs
	 * to be preserved across several calls.</p>
	 * @param index index of the secondary set, as returned by {@link
	 * org.hipparchus.ode.Expandable_ODE#add_secondary_equations(org.hipparchus.ode.Secondary_ODE)
	 * Expandable_ODE.add_secondary_equations(secondary)}
	 * @return interpolated secondary derivatives at the current interpolation date
	 * @see #get_interpolated_state()
	 * @see #get_interpolated_derivatives()
	 * @see #get_interpolated_secondary_statestatic_cast<int>(
	 * @see #set_interpolated_timestatic_cast<double>(
	 * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
	 * @deprecated as of 1.0, replaced with {@link #get_interpolated_statestatic_cast<double>(}.{@link org.hipparchus.ode.ODE_State_And_Derivative#get_secondary_derivativestatic_cast<int>(}
	 */
	@Deprecated
		std::vector<double> get_interpolated_secondary_derivatives(const int& index) Math_Illegal_State_Exception;

	/** Check if the natural integration direction is forward.
	 * <p>This method provides the integration direction as specified by
	 * the integrator itself, it avoid some nasty problems in
	 * degenerated cases like NULL steps due to cancellation at step
	 * initialization, step control or discrete events
	 * triggering.</p>
	 * @return true if the integration variable (time) increases during
	 * integration
	 */
	 //override
	bool is_forward();

	/** Copy the instance.
	 * <p>The copied instance is guaranteed to be independent from the
	 * original one. Both can be used with different settings for
	 * interpolated time without any side effect.</p>
	 * @return a deep copy of the instance, which can be used independently.
	 * @see #set_interpolated_timestatic_cast<double>(
	 * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
	 * during step constization
	 */
	Step_Interpolator copy() Math_Illegal_State_Exception;
}
