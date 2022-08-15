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

 //package org.hipparchus.filtering.kalman.linear;

 //import org.hipparchus.linear.Real_Matrix;
 //import org.hipparchus.linear.Real_Vector;

 /**
  * Container for {@link Linear_Process linear process} evolution data.
  * @see Linear_Process
  * @since 1.3
  */
class Linear_Evolution
{
	/** State transition matrix A<sub>k-1</sub>. */
	private const Real_Matrix state_transition_matrix;

	/** Control matrix B<sub>k-1</sub> (can be NULL if the process is not controlled). */
	private const Real_Matrix control_matrix;

	/** Command u<sub>k-1</sub>. (can be NULL if the process is not controlled). */
	private const Real_Vector command;

	/** Process noise matrix Q<sub>k-1</sub>. */
	private const Real_Matrix process_noise_matrix;

	/** Jacobian of the measurement with respect to the state (may be NULL). */
	private const Real_Matrix measurement_jacobian;

	/** Simple constructor.
	 * @param state_transition_matrix state transition matrix A<sub>k-1</sub>
	 * @param control_matrix control matrix B<sub>k-1</sub> (can be NULL if the process is not controlled)
	 * @param command u<sub>k-1</sub>. (can be NULL if the process is not controlled)
	 * @param process_noise_matrix process noise matrix Q<sub>k-1</sub>
	 * @param measurement_jacobian Jacobian of the measurement with respect to the state
	 * (may be NULL if measurement should be ignored)
	 */
	public Linear_Evolution(const Real_Matrix state_transition_matrix, const Real_Matrix control_matrix, const Real_Vector command, const Real_Matrix process_noise_matrix, const Real_Matrix measurement_jacobian)
	{
		this.state_transition_matrix = state_transition_matrix;
		this.control_matrix = control_matrix;
		this.command = command;
		this.process_noise_matrix = process_noise_matrix;
		this.measurement_jacobian = measurement_jacobian;
	}

	/** Get the state transition matrix A<sub>k-1</sub>.
	 * @return state transition matrix A<sub>k-1</sub>
	 */
	public Real_Matrix get_state_transition_matrix()
	{
		return state_transition_matrix;
	}

	/** Get the control matrix B<sub>k-1</sub>.
	 * @return control matrix B<sub>k-1</sub> (can be NULL if there is no control)
	 */
	public Real_Matrix get_control_matrix()
	{
		return control_matrix;
	}

	/** Get the command u<sub>k-1</sub>.
	 * @return command vector u<sub>k-1</sub> (can be NULL if there is no control)
	 */
	public Real_Vector get_command()
	{
		return command;
	}

	/** Get the process noise matrix Q<sub>k-1</sub>.
	 * @return process noise matrix<sub>k-1</sub>
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