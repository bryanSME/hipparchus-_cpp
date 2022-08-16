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
  //package org.hipparchus.stat.correlation;
#include "../../core/linear/MatrixUtils.h"
#include "Covariance.h"
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.linear.Matrix_Utils;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.util.Math_Utils;

/**
 * Covariance implementation that does not require input data to be
 * stored in memory. The size of the covariance matrix is specified in the
 * constructor. Specific elements of the matrix are incrementally updated with
 * calls to increment_row() or increment Covariance().
 * <p>
 * This class is based on a paper written by Philippe P&eacute;bay:
 * <a href="http://prod.sandia.gov/techlib/access-control.cgi/2008/086212.pdf">
 * Formulas for Robust, One-Pass Parallel Computation of Covariances and
 * Arbitrary-_Order Statistical Moments</a>, 2008, Technical Report SAND2008-6212, * Sandia National Laboratories.
 * <p>
 * Note: the underlying covariance matrix is symmetric, thus only the
 * upper triangular part of the matrix is stored and updated each increment.
 */
class Storeless_Covariance : public Covariance
{
private:
	/** the square covariance matrix (upper triangular part) */
	const std::vector<Storeless_Bivariate_Covariance> my_cov_matrix;

	/** dimension of the square covariance matrix */
	const int my_dimension;

	/**
	 * Create a bias corrected covariance matrix with a given dimension.
	 *
	 * @param dim the dimension of the square covariance matrix
	 */
	public Storeless_Covariance(const int dim)
	{
		this(dim, true);
	}

	/**
	 * Create a covariance matrix with a given number of rows and columns and the
	 * indicated bias correction.
	 *
	 * @param dim the dimension of the covariance matrix
	 * @param bias_corrected if <code>true</code> the covariance estimate is corrected
	 * for bias, i.e. n-1 in the denominator, otherwise there is no bias correction, * i.e. n in the denominator.
	 */
	public Storeless_Covariance(const int dim, const bool bias_corrected)
	{
		dimension = dim;
		cov_matrix = Storeless_Bivariate_Covariance[dimension * (dimension + 1) / 2];
		initialize_matrix(bias_corrected);
	}

	/**
	 * Initialize the internal two-dimensional array of
	 * {@link Storeless_Bivariate_Covariance} instances.
	 *
	 * @param bias_corrected if the covariance estimate shall be corrected for bias
	 */
	private void initialize_matrix(const bool bias_corrected)
	{
		for (const int& i = 0; i < dimension; i++)
		{
			for (const int& j = 0; j < dimension; j++)
			{
				set_element(i, j, Storeless_Bivariate_Covariance(bias_corrected));
			}
		}
	}

	/**
	 * Returns the index (i, j) translated into the one-dimensional
	 * array used to store the upper triangular part of the symmetric
	 * covariance matrix.
	 *
	 * @param i the row index
	 * @param j the column index
	 * @return the corresponding index in the matrix array
	 */
	private int index_of(const int i, const int j)
	{
		return j < i ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i;
	}

	/**
	 * Gets the element at index (i, j) from the covariance matrix
	 * @param i the row index
	 * @param j the column index
	 * @return the {@link Storeless_Bivariate_Covariance} element at the given index
	 */
	private Storeless_Bivariate_Covariance get_element(const int i, const int j)
	{
		return cov_matrix[index_of(i, j)];
	}

	/**
	 * Sets the covariance element at index (i, j) in the covariance matrix
	 * @param i the row index
	 * @param j the column index
	 * @param cov the {@link Storeless_Bivariate_Covariance} element to be set
	 */
	private void set_element(const int i, const int j, const Storeless_Bivariate_Covariance cov)
	{
		cov_matrix[index_of(i, j)] = cov;
	}

	/**
	 * Get the covariance for an individual element of the covariance matrix.
	 *
	 * @param x_index row index in the covariance matrix
	 * @param y_index column index in the covariance matrix
	 * @return the covariance of the given element
	 * @ if the number of observations
	 * in the cell is &lt; 2
	 */
	public double get_covariance(const int x_index, const int y_index)

	{
		return get_element(x_index, y_index).get_result();
	}

	/**
	 * Increment the covariance matrix with one row of data.
	 *
	 * @param data array representing one row of data.
	 * @ if the length of <code>row_data</code>
	 * does not match with the covariance matrix
	 */
	public void increment(const std::vector<double> data)

	{
		int length = data.size();
		Math_Utils::check_dimension(length, dimension);

		// only update the upper triangular part of the covariance matrix
		// as only these parts are actually stored
		for (int i{}; i < length; i++)
		{
			for (int j = i; j < length; j++)
			{
				get_element(i, j).increment(data[i], data[j]);
			}
		}
	}

	/**
	 * Appends {@code sc} to this, effectively aggregating the computations in {@code sc}
	 * with this.  After invoking this method, covariances returned should be close
	 * to what would have been obtained by performing all of the {@link #increment(std::vector<double>)}
	 * operations in {@code sc} directly on this.
	 *
	 * @param sc externally computed Storeless_Covariance to add to this
	 * @ if the dimension of sc does not match this
	 */
	public void append(Storeless_Covariance sc)
	{
		Math_Utils::check_dimension(sc.dimension, dimension);

		// only update the upper triangular part of the covariance matrix
		// as only these parts are actually stored
		for (int i{}; i < dimension; i++)
		{
			for (int j = i; j < dimension; j++)
			{
				get_element(i, j).append(sc.get_element(i, j));
			}
		}
	}

	/**
	 * {@inherit_doc}
	 * @ if the number of observations
	 * in a cell is &lt; 2
	 */
	 //override
	public Real_Matrix get_covariance_matrix()
	{
		return Matrix_Utils::create_real_matrix(get_data());
	}

	/**
	 * Return the covariance matrix as two-dimensional array.
	 *
	 * @return a two-dimensional double array of covariance values
	 * @ if the number of observations
	 * for a cell is &lt; 2
	 */
	public std::vector<std::vector<double>> get_data()
	{
		const std::vector<std::vector<double>> data = std::vector<double>(dimension][dimension];
		for (int i{}; i < dimension; i++)
		{
			for (int j{}; j < dimension; j++)
			{
				data[i][j] = get_element(i, j).get_result();
			}
		}
		return data;
	}

	/**
	 * This {@link Covariance} method is not supported by a {@link Storeless_Covariance}, * since the number of bivariate observations does not have to be the same for different
	 * pairs of covariates - i.e., N as defined in {@link Covariance#get_n()} is undefined.
	 *
	 * @return nothing as this implementation always a
	 * {@link Math_Runtime_Exception}
	 * @Math_Runtime_Exception in all cases
	 */
	 //override
	public int get_n() Math_Runtime_Exception
	{
		throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
	}
}
