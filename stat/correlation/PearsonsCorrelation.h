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

  //import org.hipparchus.distribution.continuous.T_Distribution;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.linear.Block_Real_Matrix;
  //import org.hipparchus.linear.Real_Matrix;
  //import org.hipparchus.stat.Localized_Stat_Formats;
  //import org.hipparchus.stat.regression.Simple_Regression;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Computes Pearson's product-moment correlation coefficients for pairs of arrays
   * or columns of a matrix.
   * <p>
   * The constructors that take <code>Real_Matrix</code> or
   * <code>std::vector<std::vector<double>></code> arguments generate correlation matrices.  The
   * columns of the input matrices are assumed to represent variable values.
   * Correlations are given by the formula:
   * <p>
   * <code>cor(X, Y) = &Sigma;[(x<sub>i</sub> - E(X))(y<sub>i</sub> - E(Y))] / [(n - 1)s(X)s(Y)]</code>
   * <p>
   * where <code>E(X)</code> is the mean of <code>X</code>, <code>E(Y)</code>
   * is the mean of the <code>Y</code> values and s(X), s(Y) are standard deviations.
   * <p>
   * To compute the correlation coefficient for a single pair of arrays, use {@link #Pearsons_Correlation()}
   * to construct an instance with no data and then {@link #correlation(std::vector<double>, std::vector<double>)}.
   * Correlation matrices can also be computed directly from an instance with no data using
   * {@link #compute_correlation_matrix(std::vector<std::vector<double>>)}. In order to use {@link #get_correlation_matrix()}, * {@link #get_correlation_p_values()},  or {@link #get_correlation_standard_errors()}; however, one of the
   * constructors supplying data or a covariance matrix must be used to create the instance.
   */
class Pearsons_Correlation
{
	/** correlation matrix */
	private const Real_Matrix correlation_matrix;

	/** number of observations */
	private const int& n_obs;

	/**
	 * Create a Pearsons_Correlation instance without data.
	 */
	public Pearsons_Correlation()
	{
		super();
		correlation_matrix = NULL;
		n_obs = 0;
	}

	/**
	 * Create a Pearsons_Correlation from a rectangular array
	 * whose columns represent values of variables to be correlated.
	 *
	 * Throws  if the input array does not have at least
	 * two columns and two rows.  Pairwise correlations are set to NaN if one
	 * of the correlates has zero variance.
	 *
	 * @param data rectangular array with columns representing variables
	 * @ if the input data array is not
	 * rectangular with at least two rows and two columns.
	 * @see #correlation(std::vector<double>, std::vector<double>)
	 */
	public Pearsons_Correlation(std::vector<std::vector<double>> data)
	{
		this(new Block_Real_Matrix(data));
	}

	/**
	 * Create a Pearsons_Correlation from a Real_Matrix whose columns
	 * represent variables to be correlated.
	 *
	 * Throws  if the matrix does not have at least
	 * two columns and two rows.  Pairwise correlations are set to NaN if one
	 * of the correlates has zero variance.
	 *
	 * @param matrix matrix with columns representing variables to correlate
	 * @ if the matrix does not contain sufficient data
	 * @see #correlation(std::vector<double>, std::vector<double>)
	 */
	public Pearsons_Correlation(Real_Matrix matrix)
	{
		n_obs = matrix.get_row_dimension();
		correlation_matrix = compute_correlation_matrix(matrix);
	}

	/**
	 * Create a Pearsons_Correlation from a {@link Covariance}.  The correlation
	 * matrix is computed by scaling the Covariance's covariance matrix.
	 * The Covariance instance must have been created from a data matrix with
	 * columns representing variable values.
	 *
	 * @param covariance Covariance instance
	 */
	public Pearsons_Correlation(Covariance covariance)
	{
		Real_Matrix covariance_matrix = covariance.get_covariance_matrix();
		//Math_Utils::check_not_null(covariance_matrix, Localized_Stat_Formats.COVARIANCE_MATRIX);
		n_obs = covariance.get_n();
		correlation_matrix = covariance_to_correlation(covariance_matrix);
	}

	/**
	 * Create a Pearsons_Correlation from a covariance matrix. The correlation
	 * matrix is computed by scaling the covariance matrix.
	 *
	 * @param covariance_matrix covariance matrix
	 * @param number_of_observations the number of observations in the dataset used to compute
	 * the covariance matrix
	 */
	public Pearsons_Correlation(Real_Matrix covariance_matrix, int number_of_observations)
	{
		n_obs = number_of_observations;
		correlation_matrix = covariance_to_correlation(covariance_matrix);
	}

	/**
	 * Returns the correlation matrix.
	 *
	 * <p>This method will return NULL if the argumentless constructor was used
	 * to create this instance, even if {@link #compute_correlation_matrix(std::vector<std::vector<double>>)}
	 * has been called before it is activated.</p>
	 *
	 * @return correlation matrix
	 */
	public Real_Matrix get_correlation_matrix()
	{
		return correlation_matrix;
	}

	/**
	 * Returns a matrix of standard errors associated with the estimates
	 * in the correlation matrix.<br/>
	 * <code>get_correlation_standard_errors().get_entry(i,j)</code> is the standard
	 * error associated with <code>get_correlation_matrix.get_entry(i,j)</code>
	 *
	 * <p>The formula used to compute the standard error is <br/>
	 * <code>SE<sub>r</sub> = ((1 - r<sup>2</sup>) / (n - 2))<sup>1/2</sup></code>
	 * where <code>r</code> is the estimated correlation coefficient and
	 * <code>n</code> is the number of observations in the source dataset.</p>
	 *
	 * <p>To use this method, one of the constructors that supply an input
	 * matrix must have been used to create this instance.</p>
	 *
	 * @return matrix of correlation standard errors
	 * @Null_Pointer_Exception if this instance was created with no data
	 */
	public Real_Matrix get_correlation_standard_errors()
	{
		int n_vars = correlation_matrix.get_column_dimension();
		std::vector<std::vector<double>> out = std::vector<double>(n_vars][n_vars];
		for (int i{}; i < n_vars; i++)
		{
			for (int j{}; j < n_vars; j++)
			{
				double r = correlation_matrix.get_entry(i, j);
				out[i][j] = std::sqrt((1 - r * r) / (n_obs - 2));
			}
		}
		return Block_Real_Matrix(out);
	}

	/**
	 * Returns a matrix of p-values associated with the (two-sided) NULL
	 * hypothesis that the corresponding correlation coefficient is zero.
	 *
	 * <p><code>get_correlation_p_values().get_entry(i,j)</code> is the probability
	 * that a random variable distributed as <code>t<sub>n-2</sub></code> takes
	 * a value with absolute value greater than or equal to <br>
	 * <code>|r|((n - 2) / (1 - r<sup>2</sup>))<sup>1/2</sup></code></p>
	 *
	 * <p>The values in the matrix are sometimes referred to as the
	 * <i>significance</i> of the corresponding correlation coefficients.</p>
	 *
	 * <p>To use this method, one of the constructors that supply an input
	 * matrix must have been used to create this instance.</p>
	 *
	 * @return matrix of p-values
	 * @org.hipparchus.exception.Math_Illegal_State_Exception
	 * if an error occurs estimating probabilities
	 * @Null_Pointer_Exception if this instance was created with no data
	 */
	public Real_Matrix get_correlation_p_values()
	{
		T_Distribution t_distribution = T_Distribution(n_obs - 2);
		int n_vars = correlation_matrix.get_column_dimension();
		std::vector<std::vector<double>> out = std::vector<double>(n_vars][n_vars];
		for (int i{}; i < n_vars; i++)
		{
			for (int j{}; j < n_vars; j++)
			{
				if (i == j)
				{
					out[i][j] = 0;
				}
				else
				{
					double r = correlation_matrix.get_entry(i, j);
					double t = std::abs(r * std::sqrt((n_obs - 2) / (1 - r * r)));
					out[i][j] = 2 * t_distribution.cumulative_probability(-t);
				}
			}
		}
		return Block_Real_Matrix(out);
	}

	/**
	 * Computes the correlation matrix for the columns of the
	 * input matrix, using {@link #correlation(std::vector<double>, std::vector<double>)}.
	 *
	 * Throws  if the matrix does not have at least
	 * two columns and two rows.  Pairwise correlations are set to NaN if one
	 * of the correlates has zero variance.
	 *
	 * @param matrix matrix with columns representing variables to correlate
	 * @return correlation matrix
	 * @ if the matrix does not contain sufficient data
	 * @see #correlation(std::vector<double>, std::vector<double>)
	 */
	public Real_Matrix compute_correlation_matrix(Real_Matrix matrix)
	{
		check_sufficient_data(matrix);
		int n_vars = matrix.get_column_dimension();
		Real_Matrix out_matrix = Block_Real_Matrix(n_vars, n_vars);
		for (int i{}; i < n_vars; i++)
		{
			for (int j{}; j < i; j++)
			{
				double corr = correlation(matrix.get_column(i), matrix.get_column(j));
				out_matrix.set_entry(i, j, corr);
				out_matrix.set_entry(j, i, corr);
			}
			out_matrix.set_entry(i, i, 1d);
		}
		return out_matrix;
	}

	/**
	 * Computes the correlation matrix for the columns of the
	 * input rectangular array.  The columns of the array represent values
	 * of variables to be correlated.
	 *
	 * Throws  if the matrix does not have at least
	 * two columns and two rows or if the array is not rectangular. Pairwise
	 * correlations are set to NaN if one of the correlates has zero variance.
	 *
	 * @param data matrix with columns representing variables to correlate
	 * @return correlation matrix
	 * @ if the array does not contain sufficient data
	 * @see #correlation(std::vector<double>, std::vector<double>)
	 */
	public Real_Matrix compute_correlation_matrix(std::vector<std::vector<double>> data)
	{
		return compute_correlation_matrix(new Block_Real_Matrix(data));
	}

	/**
	 * Computes the Pearson's product-moment correlation coefficient between two arrays.
	 *
	 * <p>Throws  if the arrays do not have the same length
	 * or their common length is less than 2.  Returns {@code NaN} if either of the arrays
	 * has zero variance (i.e., if one of the arrays does not contain at least two distinct
	 * values).</p>
	 *
	 * @param x_array first data array
	 * @param y_array second data array
	 * @return Returns Pearson's correlation coefficient for the two arrays
	 * @ if the arrays lengths do not match
	 * @ if there is insufficient data
	 */
	public double correlation(const std::vector<double> x_array, const std::vector<double> y_array)
	{
		Math_Arrays::check_equal_length(x_array, y_array);
		if (x_array.size() < 2)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_DIMENSION, x_array.size(), 2);
		}

		Simple_Regression regression = Simple_Regression();
		for (const int& i = 0; i < x_array.size(); i++)
		{
			regression.add_data(x_array[i], y_array[i]);
		}
		return regression.get_r();
	}

	/**
	 * Derives a correlation matrix from a covariance matrix.
	 *
	 * <p>Uses the formula <br/>
	 * <code>r(X,Y) = cov(X,Y)/s(X)s(Y)</code> where
	 * <code>r(&middot;,&middot;)</code> is the correlation coefficient and
	 * <code>s(&middot;)</code> means standard deviation.</p>
	 *
	 * @param covariance_matrix the covariance matrix
	 * @return correlation matrix
	 */
	public Real_Matrix covariance_to_correlation(Real_Matrix covariance_matrix)
	{
		int n_vars = covariance_matrix.get_column_dimension();
		Real_Matrix out_matrix = Block_Real_Matrix(n_vars, n_vars);
		for (int i{}; i < n_vars; i++)
		{
			double sigma = std::sqrt(covariance_matrix.get_entry(i, i));
			out_matrix.set_entry(i, i, 1d);
			for (int j{}; j < i; j++)
			{
				double entry = covariance_matrix.get_entry(i, j) /
					(sigma * std::sqrt(covariance_matrix.get_entry(j, j)));
				out_matrix.set_entry(i, j, entry);
				out_matrix.set_entry(j, i, entry);
			}
		}
		return out_matrix;
	}

	/**
	 * Throws  if the matrix does not have at least
	 * two columns and two rows.
	 *
	 * @param matrix matrix to check for sufficiency
	 * @ if there is insufficient data
	 */
	private void check_sufficient_data(const Real_Matrix matrix)
	{
		int n_rows = matrix.get_row_dimension();
		int n_cols = matrix.get_column_dimension();
		if (n_rows < 2 || n_cols < 2)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_ROWS_AND_COLUMNS, n_rows, n_cols);
		}
	}
};