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
  //package org.hipparchus.stat.regression;

  //import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
  //import org.hipparchus.linear.LU_Decomposition;
  //import org.hipparchus.linear.Real_Matrix;
  //import org.hipparchus.linear.Real_Vector;
#include <vector>
#include "AbstractMultipleLinearRegression.h"
#include "../../core/linear/RealMatrix.h"

/**
 * The GLS implementation of multiple linear regression.
 *
 * GLS assumes a general covariance matrix Omega of the error
 * <pre>
 * u ~ N(0, Omega)
 * </pre>
 *
 * Estimated by GLS, * <pre>
 * b=(X' Omega^-1 X)^-1X'Omega^-1 y
 * </pre>
 * whose variance is
 * <pre>
 * Var(b)=(X' Omega^-1 X)^-1
 * </pre>
 */
class GLSMultiple_Linear_Regression : public Abstract_Multiple_Linear_Regression
{
private:
	/** Covariance matrix. */
	Real_Matrix my_omega;

	/** Inverse of covariance matrix. */
	Real_Matrix my_omega_inverse;

public:
	/** Replace sample data, overriding any previous sample.
	 * @param y y values of the sample
	 * @param x x values of the sample
	 * @param covariance array representing the covariance matrix
	 */
	void new_sample_data(const std::vector<double>& y, const std::vector<std::vector<double>>& x, const std::vector<std::vector<double>>& covariance)
	{
		validate_sample_data(x, y);
		new_y_sample_data(y);
		new_x_sample_data(x);
		validate_covariance_data(x, covariance);
		new_covariance_data(covariance);
	}

protected:
	/**
	 * Add the covariance data.
	 *
	 * @param omega the [n,n] array representing the covariance
	 */
	void new_covariance_data(std::vector<std::vector<double>> omega)
	{
		my_omega = Array_2D_Row_Real_Matrix(omega);
		my_omega_inverse = NULL;
	}

	/**
	 * Get the inverse of the covariance.
	 * <p>The inverse of the covariance matrix is lazily evaluated and cached.</p>
	 * @return inverse of the covariance
	 */
	Real_Matrix get_omega_inverse()
	{
		if (my_omega_inverse == NULL)
		{
			my_omega_inverse = LU_Decomposition(my_omega).get_solver().get_inverse();
		}
		return my_omega_inverse;
	}

	/**
	 * Calculates beta by GLS.
	 * <pre>
	 *  b=(X' Omega^-1 X)^-1X'Omega^-1 y
	 * </pre>
	 * @return beta
	 */
	 //override
	Real_Vector calculate_beta()
	{
		Real_Matrix OI = get_omega_inverse();
		Real_Matrix XT = get_x().transpose();
		Real_Matrix XTOIX = XT.multiply(OI).multiply(get_x());
		Real_Matrix inverse = LU_Decomposition(XTOIX).get_solver().get_inverse();
		return inverse.multiply(XT).multiply(OI).operate(get_y());
	}

	/**
	 * Calculates the variance on the beta.
	 * <pre>
	 *  Var(b)=(X' Omega^-1 X)^-1
	 * </pre>
	 * @return The beta variance matrix
	 */
	 //override
	Real_Matrix calculate_beta_variance()
	{
		Real_Matrix OI = get_omega_inverse();
		Real_Matrix XTOIX = get_x().transpose_multiply(OI).multiply(get_x());
		return LU_Decomposition(XTOIX).get_solver().get_inverse();
	}

	/**
	 * Calculates the estimated variance of the error term using the formula
	 * <pre>
	 *  Var(u) = Tr(u' Omega^-1 u)/(n-k)
	 * </pre>
	 * where n and k are the row and column dimensions of the design
	 * matrix X.
	 *
	 * @return error variance
	 */
	 //override
	double calculate_error_variance()
	{
		Real_Vector residuals = calculate_residuals();
		double t = residuals.dot_product(get_omega_inverse().operate(residuals));
		return t / (get_x().get_row_dimension() - get_x().get_column_dimension());
	}
};