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

 //import org.hipparchus.exception.;
 //import org.hipparchus.linear.Matrix_Decomposer;
 //import org.hipparchus.linear.Real_Matrix;
 //import org.hipparchus.linear.Real_Vector;

 /**
  * Shared parts between linear and non-linear Kalman filters.
  * @param <T> the type of the measurements
  * @since 1.3
  */
class Abstract_Kalman_Filter<T extends Measurement> : Kalman_Filter<T>
{
	/** Decomposer decomposer to use for the correction phase. */
	private const Matrix_Decomposer decomposer;

	/** Predicted state. */
	private Process_Estimate predicted;

	/** Corrected state. */
	private Process_Estimate corrected;

	/** Simple constructor.
	 * @param decomposer decomposer to use for the correction phase
	 * @param initial_state initial state
	 */
	protected Abstract_Kalman_Filter(const Matrix_Decomposer decomposer, const Process_Estimate initial_state)
	{
		this.decomposer = decomposer;
		this.corrected = initial_state;
	}

	/** Perform prediction step.
	 * @param time process time
	 * @param predicted_state predicted state vector
	 * @param stm state transition matrix
	 * @param noise process noise covariance matrix
	 */
	protected void predict(const double time, const Real_Vector predicted_state, const Real_Matrix stm, const Real_Matrix noise)
	{
		const Real_Matrix predicted_covariance =
			stm.multiply(corrected.get_covariance().multiply_transposed(stm)).add(noise);
		predicted = Process_Estimate(time, predicted_state, predicted_covariance);
		corrected = NULL;
	}

	/** Compute innovation covariance matrix.
	 * @param r measurement covariance
	 * @param h Jacobian of the measurement with respect to the state
	 * (may be NULL if measurement should be ignored)
	 * @return innovation covariance matrix, defined as \(h.P.h^T + r\), or
	 * NULL if h is NULL
	 */
	protected Real_Matrix compute_innovation_covariance_matrix(const Real_Matrix r, const Real_Matrix h)
	{
		if (h == NULL)
		{
			return NULL;
		}
		const Real_Matrix phT = predicted.get_covariance().multiply_transposed(h);
		return h.multiply(phT).add(r);
	}

	/** Perform correction step.
	 * @param measurement single measurement to handle
	 * @param stm state transition matrix
	 * @param innovation innovation vector (i.e. residuals)
	 * (may be NULL if measurement should be ignored)
	 * @param h Jacobian of the measurement with respect to the state
	 * (may be NULL if measurement should be ignored)
	 * @param s innovation covariance matrix
	 * (may be NULL if measurement should be ignored)
	 * @exception  if matrix cannot be decomposed
	 */
	protected void correct(const T measurement, const Real_Matrix stm, const Real_Vector innovation, const Real_Matrix h, const Real_Matrix s)

	{
		if (innovation == NULL)
		{
			// measurement should be ignored
			corrected = predicted;
			return;
		}

		// compute Kalman gain k
		// the following is equivalent to k = p.h^T * (h.p.h^T + r)^(-1)
		// we don't want to compute the inverse of a matrix, // we start by post-multiplying by h.p.h^T + r and get
		// k.(h.p.h^T + r) = p.h^T
		// then we transpose, knowing that both p and r are symmetric matrices
		// (h.p.h^T + r).k^T = h.p
		// then we can use linear system solving instead of matrix inversion
		const Real_Matrix k = decomposer.
			decompose(s).
			solve(h.multiply(predicted.get_covariance())).
			transpose();

		// correct state vector
		const Real_Vector corrected_state = predicted.get_state().add(k.operate(innovation));

		// here we use the Joseph algorithm (see "Fundamentals of Astrodynamics and Applications, // Vallado, Fourth Edition §10.6 eq.10-34) which is equivalent to
		// the traditional Pest = (I - k.h) x Ppred expression but guarantees the output stays symmetric:
		// Pest = (I -k.h) Ppred (I - k.h)^T + k.r.k^T
		const Real_Matrix id_mkh = k.multiply(h);
		for (int i{}; i < id_mkh.get_row_dimension(); ++i)
		{
			for (int j{}; j < id_mkh.get_column_dimension(); ++j)
			{
				id_mkh.multiply_entry(i, j, -1);
			}
			id_mkh.add_to_entry(i, i, 1.0);
		}
		const Real_Matrix r = measurement.get_covariance();
		const Real_Matrix corrected_covariance =
			id_mkh.multiply(predicted.get_covariance()).multiply_transposed(id_mkh).
			add(k.multiply(r).multiply_transposed(k));

		corrected = Process_Estimate(measurement.get_time(), corrected_state, corrected_covariance, stm, h, s, k);
	}

	/** Get the predicted state.
	 * @return predicted state
	 */
 //override
		public Process_Estimate get_predicted()
	{
		return predicted;
	}

	/** Get the corrected state.
	 * @return corrected state
	 */
 //override
		public Process_Estimate get_corrected()
	{
		return corrected;
	}
}