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

 //import org.hipparchus.exception.Math_Runtime_Exception;
 //import org.hipparchus.filtering.kalman.Abstract_Kalman_Filter;
 //import org.hipparchus.filtering.kalman.Measurement;
 //import org.hipparchus.filtering.kalman.Process_Estimate;
 //import org.hipparchus.linear.Matrix_Decomposer;
 //import org.hipparchus.linear.Real_Matrix;
 //import org.hipparchus.linear.Real_Vector;

 /**
  * Kalman filter for {@link Linear_Process linear process}.
  * @param <T> the type of the measurements
  * @since 1.3
  */
class Linear_Kalman_Filter<T extends Measurement> extends Abstract_Kalman_Filter<T>
{
	/** Process to be estimated. */
	private const Linear_Process<T> process;

	/** Simple constructor.
	 * @param decomposer decomposer to use for the correction phase
	 * @param process linear process to estimate
	 * @param initial_state initial state
	 */
	public Linear_Kalman_Filter(const Matrix_Decomposer decomposer, const Linear_Process<T> process, const Process_Estimate initial_state)
	{
		super(decomposer, initial_state);
		this.process = process;
	}

	/** {@inherit_doc} */
 //override
		public Process_Estimate estimation_step(const T measurement)
		Math_Runtime_Exception
	{
		const Linear_Evolution evolution = process.get_evolution(measurement);

		// prediction phase
		const Real_Matrix      a = evolution.get_state_transition_matrix();
		const Real_Matrix      b = evolution.get_control_matrix();
		const Real_Vector      u = (b == NULL) ? NULL : evolution.get_command();
		const Real_Matrix      q = evolution.get_process_noise_matrix();

		Real_Vector pred_xk = a.operate(get_corrected().get_state());
		if (b != NULL)
		{
			pred_xk = pred_xk.add(b.operate(u));
		}

		predict(measurement.get_time(), pred_xk, a, q);

		// correction phase
		const Real_Matrix h = evolution.get_measurement_jacobian();
		const Real_Matrix s = compute_innovation_covariance_matrix(measurement.get_covariance(), h);
		const Real_Vector innovation = (h == NULL) ? NULL : measurement.get_value().subtract(h.operate(pred_xk));
		correct(measurement, a, innovation, h, s);
		return get_corrected();
	}
}