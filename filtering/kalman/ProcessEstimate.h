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

 //package org.hipparchus.filtering.kalman;

 //import org.hipparchus.linear.Real_Matrix;
 //import org.hipparchus.linear.Real_Vector;

 /**
  * Holder for process state and covariance.
  * <p>
  * The estimate always contains time, state and covariance. These data are
  * the only ones needed to start a Kalman filter. Once a filter has been
  * started and produces estimates, these estimates will always
  * contain a state transition matrix and if the measurement has not been
  * ignored, they will also contain measurement Jacobian, innovation covariance
  * and Kalman gain.
  * </p>
  * @since 1.3
  */
class Process_Estimate
{
	/** Process time (typically the time or index of a measurement). */
	private const double time;

	/** State vector. */
	private const Real_Vector state;

	/** State covariance. */
	private const Real_Matrix covariance;

	/** State transition matrix, may be NULL.
	 * @since 1.4
	 */
	private const Real_Matrix state_transition_matrix;

	/** Jacobian of the measurement with respect to the state (h matrix), may be NULL.
	 * @since 1.4
	 */
	private const Real_Matrix measurement_jacobian;

	/** Innovation covariance matrix, defined as \(h.P.h^T + r\), may be NULL.
	 * @since 1.4
	 */
	private const Real_Matrix innovation_covariance_matrix;

	/** Kalman gain (k matrix), may be NULL.
	 * @since 1.4
	 */
	private const Real_Matrix kalman_gain;

	/** Simple constructor.
	 * <p>
	 * This constructor sets state transition matrix, covariance matrix H, * innovation covariance matrix and Kalman gain k to NULL.
	 * </p>
	 * @param time process time (typically the time or index of a measurement)
	 * @param state state vector
	 * @param covariance state covariance
	 */
	public Process_Estimate(const double time, const Real_Vector state, const Real_Matrix covariance)
	{
		this(time, state, covariance, NULL, NULL, NULL, NULL);
	}

	/** Simple constructor.
	 * @param time process time (typically the time or index of a measurement)
	 * @param state state vector
	 * @param covariance state covariance
	 * @param state_transition_matrix state transition matrix between previous state and estimated (but not yet corrected) state
	 * @param measurement_jacobian Jacobian of the measurement with respect to the state
	 * @param innovation_covariance innovation covariance matrix, defined as \(h.P.h^T + r\), may be NULL
	 * @param kalman_gain Kalman Gain matrix, may be NULL
	 * @since 1.4
	 */
	public Process_Estimate(const double time, const Real_Vector state, const Real_Matrix covariance, const Real_Matrix state_transition_matrix, const Real_Matrix measurement_jacobian, const Real_Matrix innovation_covariance, const Real_Matrix kalman_gain)
	{
		this.time = time;
		this.state = state;
		this.covariance = covariance;
		this.state_transition_matrix = state_transition_matrix;
		this.measurement_jacobian = measurement_jacobian;
		this.innovation_covariance_matrix = innovation_covariance;
		this.kalman_gain = kalman_gain;
	}

	/** Get the process time.
	 * @return process time (typically the time or index of a measurement)
	 */
	public double get_time()
	{
		return time;
	}

	/** Get the state vector.
	 * @return state vector
	 */
	public Real_Vector get_state()
	{
		return state;
	}

	/** Get the state covariance.
	 * @return state covariance
	 */
	public Real_Matrix get_covariance()
	{
		return covariance;
	}

	/** Get state transition matrix between previous state and estimated (but not yet corrected) state.
	 * @return state transition matrix between previous state and estimated state (but not yet corrected)
	 * (may be NULL for initial process estimate)
	 * @since 1.4
	 */
	public Real_Matrix get_state_transition_matrix()
	{
		return state_transition_matrix;
	}

	/** Get the Jacobian of the measurement with respect to the state (H matrix).
	 * @return Jacobian of the measurement with respect to the state (may be NULL for initial
	 * process estimate or if the measurement has been ignored)
	 * @since 1.4
	 */
	public Real_Matrix get_measurement_jacobian()
	{
		return measurement_jacobian;
	}

	/** Get the innovation covariance matrix.
	 * @return innovation covariance matrix (may be NULL for initial
	 * process estimate or if the measurement has been ignored)
	 * @since 1.4
	 */
	public Real_Matrix get_innovation_covariance()
	{
		return innovation_covariance_matrix;
	}

	/** Get the Kalman gain matrix.
	 * @return Kalman gain matrix (may be NULL for initial
	 * process estimate or if the measurement has been ignored)
	 * @since 1.4
	 */
	public Real_Matrix get_kalman_gain()
	{
		return kalman_gain;
	}
}