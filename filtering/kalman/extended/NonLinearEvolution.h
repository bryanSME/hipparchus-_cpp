#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

 //package org.hipparchus.filtering.kalman.extended;

 //import org.hipparchus.linear.Real_Matrix;
 //import org.hipparchus.linear.Real_Vector;

 /**
  * Container for {@link Non_Linear_Process non-linear process} evolution data.
  * @see Non_Linear_Process
  * @since 1.3
  */
class NonLinear_Evolution
{
	/** Current time. */
	private const double& current_time;

	/** State vector at current time. */
	private const Real_Vector current_state;

	/** State transition matrix between previous and current state. */
	private const Real_Matrix state_transition_matrix;

	/** Process noise matrix. */
	private const Real_Matrix process_noise_matrix;

	/** Jacobian of the measurement with respect to the state (may be NULL). */
	private const Real_Matrix measurement_jacobian;

	/** Simple constructor.
	 * @param current_time current time
	 * @param current_state state vector at current time
	 * @param state_transition_matrix state transition matrix between previous and current state
	 * @param process_noise_matrix process noise
	 * @param measurement_jacobian Jacobian of the measurement with respect to the state
	 * (may be NULL if measurement should be ignored)
	 */
	public NonLinear_Evolution(const double& current_time, const Real_Vector current_state, const Real_Matrix state_transition_matrix, const Real_Matrix process_noise_matrix, const Real_Matrix measurement_jacobian)
	{
		this.current_time = current_time;
		this.current_state = current_state;
		this.state_transition_matrix = state_transition_matrix;
		this.process_noise_matrix = process_noise_matrix;
		this.measurement_jacobian = measurement_jacobian;
	}

	/** Get current time.
	 * @return current time
	 */
	public double get_current_time() const
	{
		return current_time;
	}

	/** Get current state.
	 * @return current state
	 */
	public Real_Vector get_current_state()
	{
		return current_state;
	}

	/** Get state transition matrix between previous and current state.
	 * @return state transition matrix between previous and current state
	 */
	public Real_Matrix get_state_transition_matrix()
	{
		return state_transition_matrix;
	}

	/** Get process noise.
	 * @return process noise
	 */
	public Real_Matrix get_process_noise_matrix()
	{
		return process_noise_matrix;
	}

	/** Get measurement Jacobian.
	 * @return Jacobian of the measurement with respect to the state
	 * (may be NULL if measurement should be ignored)
	 */
	public Real_Matrix get_measurement_jacobian()
	{
		return measurement_jacobian;
	}
}