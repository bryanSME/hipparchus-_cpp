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
#include "../../core/linear/BlockRealMatrix.h"
#include "../descriptive/moment/Mean.h"
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.linear.Block_Real_Matrix;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.stat.descriptive.moment.Mean;
//import org.hipparchus.stat.descriptive.moment.Variance;

/**
 * Computes covariances for pairs of arrays or columns of a matrix.
 * <p>
 * The constructors that take {@code Real_Matrix} or {@code std::vector<std::vector<double>>}
 * arguments generate covariance matrices. The columns of the input
 * matrices are assumed to represent variable values.
 * <p>
 * The constructor argument {@code bias_corrected} determines whether or
 * not computed covariances are bias-corrected.
 * <p>
 * Unbiased covariances are given by the formula:
 * <p>
 * <code>cov(X, Y) = &Sigma;[(x<sub>i</sub> - E(X))(y<sub>i</sub> - E(Y))] / (n - 1)</code>
 * <p>
 * where {@code E(X)} is the mean of {@code X} and {@code E(Y)}
 * is the mean of the <code>Y</code> values.
 * <p>
 * Non-bias-corrected estimates use {@code n} in place of {@code n - 1}.
 */
class Covariance 
{
private:
    /** The covariance matrix. */
    const Real_Matrix my_covariance_matrix;

    /** Number of observations (length of covariate vectors). */
    const int my_n;

    /**
     * Throws  if the matrix does not have at least
     * one column and two rows.
     *
     * @param matrix matrix to check
     * @ if the matrix does not contain sufficient data
     * to compute covariance
     */
    void check_sufficient_data(const Real_Matrix& matrix)
    {
        int n_rows = matrix.get_row_dimension();
        int n_cols = matrix.get_column_dimension();
        if (n_rows < 2 || n_cols < 1)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_ROWS_AND_COLUMNS, n_rows, n_cols);
        }
    }

protected:
    /**
     * Compute a covariance matrix from a matrix whose columns represent covariates.
     *
     * @param matrix input matrix (must have at least one column and two rows)
     * @param bias_corrected determines whether or not covariance estimates are bias-corrected
     * @return covariance matrix
     * @ if the matrix does not contain sufficient data
     */
    Real_Matrix compute_covariance_matrix(const Real_Matrix& matrix, bool bias_corrected)
    {
        int dimension = matrix.get_column_dimension();
        Variance variance = Variance(bias_corrected);
        Real_Matrix out_matrix = Block_Real_Matrix(dimension, dimension);
        for (int i{}; i < dimension; i++)
        {
            for (int j{}; j < i; j++)
            {
                double cov = covariance(matrix.get_column(i), matrix.get_column(j), bias_corrected);
                out_matrix.set_entry(i, j, cov);
                out_matrix.set_entry(j, i, cov);
            }
            out_matrix.set_entry(i, i, variance.evaluate(matrix.get_column(i)));
        }
        return out_matrix;
    }

    /**
     * Create a covariance matrix from a matrix whose columns represent
     * covariates. Covariances are computed using the bias-corrected formula.
     *
     * @param matrix input matrix (must have at least one column and two rows)
     * @return covariance matrix
     * @ if matrix does not contain sufficient data
     * @see #Covariance
     */
    Real_Matrix compute_covariance_matrix(const Real_Matrix& matrix)
    {
        return compute_covariance_matrix(matrix, true);
    }

    /**
     * Compute a covariance matrix from a rectangular array whose columns represent covariates.
     *
     * @param data input array (must have at least one column and two rows)
     * @param bias_corrected determines whether or not covariance estimates are bias-corrected
     * @return covariance matrix
     * @ if the data array does not contain sufficient data
     * @ if the input data array is not
     * rectangular with at least one row and one column.
     */
    Real_Matrix compute_covariance_matrix(const std::vector<std::vector<double>>& data, bool bias_corrected)
    {
        return compute_covariance_matrix(Block_Real_Matrix(data), bias_corrected);
    }

    /**
     * Create a covariance matrix from a rectangular array whose columns represent
     * covariates. Covariances are computed using the bias-corrected formula.
     *
     * @param data input array (must have at least one column and two rows)
     * @return covariance matrix
     * @ if the data array does not contain sufficient data
     * @ if the input data array is not
     * rectangular with at least one row and one column.
     * @see #Covariance
     */
    Real_Matrix compute_covariance_matrix(const std::vector<std::vector<double>>& data)
    {
        return compute_covariance_matrix(data, true);
    }

public:
    /**
     * Create a Covariance with no data.
     */
    Covariance() 
        :
        my_covariance_matrix{},
        my_n{ 0 }
    {
        super();
    }

    /**
     * Create a Covariance matrix from a rectangular array
     * whose columns represent covariates.
     * <p>
     * The <code>bias_corrected</code> parameter determines whether or not
     * covariance estimates are bias-corrected.
     * <p>
     * The input array must be rectangular with at least one column
     * and two rows.
     *
     * @param data rectangular array with columns representing covariates
     * @param bias_corrected true means covariances are bias-corrected
     * @ if the input data array is not
     * rectangular with at least two rows and one column.
     * @ if the input data array is not
     * rectangular with at least one row and one column.
     */
    Covariance(const std::vector<std::vector<double>>& data, bool bias_corrected)
    {
        Covariance(Block_Real_Matrix(data), bias_corrected);
    }

    /**
     * Create a Covariance matrix from a rectangular array
     * whose columns represent covariates.
     * <p>
     * The input array must be rectangular with at least one column
     * and two rows.
     *
     * @param data rectangular array with columns representing covariates
     * @ if the input data array is not
     * rectangular with at least two rows and one column.
     * @ if the input data array is not
     * rectangular with at least one row and one column.
     */
    Covariance(const std::vector<std::vector<double>>& data)  
    {
        Covariance(data, true);
    }

    /**
     * Create a covariance matrix from a matrix whose columns
     * represent covariates.
     * <p>
     * The <code>bias_corrected</code> parameter determines whether or not
     * covariance estimates are bias-corrected.
     * <p>
     * The matrix must have at least one column and two rows.
     *
     * @param matrix matrix with columns representing covariates
     * @param bias_corrected true means covariances are bias-corrected
     * @ if the input matrix does not have
     * at least two rows and one column
     */
    Covariance(const Real_Matrix& matrix, bool bias_corrected)
    {
        check_sufficient_data(matrix);
        my_n = matrix.get_row_dimension();
        my_covariance_matrix = compute_covariance_matrix(matrix, bias_corrected);
    }

    /**
     * Create a covariance matrix from a matrix whose columns
     * represent covariates.
     * <p>
     * The matrix must have at least one column and two rows.
     *
     * @param matrix matrix with columns representing covariates
     * @ if the input matrix does not have
     * at least two rows and one column
     */
    Covariance(const Real_Matrix& matrix)  
    {
        Covariance(matrix, true);
    }

    /**
     * Returns the covariance matrix
     *
     * @return covariance matrix
     */
    Real_Matrix get_covariance_matrix() const
    {
        return my_covariance_matrix;
    }

    /**
     * Returns the number of observations (length of covariate vectors)
     *
     * @return number of observations
     */
    int get_n() const
    {
        return my_n;
    }

    /**
     * Computes the covariance between the two arrays.
     * <p>
     * Array lengths must match and the common length must be at least 2.
     *
     * @param x_array first data array
     * @param y_array second data array
     * @param bias_corrected if true, returned value will be bias-corrected
     * @return returns the covariance for the two arrays
     * @  if the arrays lengths do not match or
     * there is insufficient data
     */
    double covariance(const std::vector<double> x_array, const std::vector<double> y_array, bool bias_corrected)
    {
        Mean mean = Mean();
        double result = 0;
        int length = x_array.size();
        if (length != y_array.size()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, length, y_array.size());
        }
        if (length < 2) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_OBSERVED_POINTS_IN_SAMPLE, length, 2);
        }
        double x_mean = mean.evaluate(x_array);
        double y_mean = mean.evaluate(y_array);
        for (int i{}; i < length; i++) 
        {
            double x_dev = x_array[i] - x_mean;
            double y_dev = y_array[i] - y_mean;
            result += (x_dev * y_dev - result) / (i + 1);
        }
        return bias_corrected
            ? result * (static_cast<double>( length / static_cast<double>((length - 1))
                : result;
    }

    /**
     * Computes the covariance between the two arrays, using the bias-corrected
     * formula.
     * <p>
     * Array lengths must match and the common length must be at least 2.
     *
     * @param x_array first data array
     * @param y_array second data array
     * @return returns the covariance for the two arrays
     * @ if the arrays lengths do not match or
     * there is insufficient data
     */
    double covariance(const std::vector<double> x_array, const std::vector<double> y_array)
    {
        return covariance(x_array, y_array, true);
    }
};