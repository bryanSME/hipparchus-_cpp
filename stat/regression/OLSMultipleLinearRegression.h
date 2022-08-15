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

//import org.hipparchus.exception.;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.linear.LU_Decomposition;
//import org.hipparchus.linear.QR_Decomposition;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.stat.Stat_Utils;
//import org.hipparchus.stat.descriptive.moment.Second_Moment;
#include "AbstractMultipleLinearRegression.h"
#include "../../core/linear/QRDecomposition.h"
#include "../StatUtils.h"

/**
 * <p>Implements ordinary least squares (OLS) to estimate the parameters of a
 * multiple linear regression model.</p>
 *
 * <p>The regression coefficients, <code>b</code>, satisfy the normal equations:
 * <pre><code> X<sup>T</sup> X b = X<sup>T</sup> y </code></pre></p>
 *
 * <p>To solve the normal equations, this implementation uses QR decomposition
 * of the <code>X</code> matrix. (See {@link QR_Decomposition} for details on the
 * decomposition algorithm.) The <code>X</code> matrix, also known as the <i>design matrix,</i>
 * has rows corresponding to sample observations and columns corresponding to independent
 * variables.  When the model is estimated using an intercept term (i.e. when
 * {@link #is_no_intercept() is_no_intercept} is false as it is by default), the <code>X</code>
 * matrix includes an initial column identically equal to 1.  We solve the normal equations
 * as follows:
 * <pre><code> X<sup>T</sup>X b = X<sup>T</sup> y
 * (QR)<sup>T</sup> (QR) b = (QR)<sup>T</sup>y
 * R<sup>T</sup> (Q<sup>T</sup>Q) R b = R<sup>T</sup> Q<sup>T</sup> y
 * R<sup>T</sup> R b = R<sup>T</sup> Q<sup>T</sup> y
 * (R<sup>T</sup>)<sup>-1</sup> R<sup>T</sup> R b = (R<sup>T</sup>)<sup>-1</sup> R<sup>T</sup> Q<sup>T</sup> y
 * R b = Q<sup>T</sup> y </code></pre></p>
 *
 * <p>Given <code>Q</code> and <code>R</code>, the last equation is solved by back-substitution.</p>
 *
 */
class OLS_Multiple_Linear_Regression : public Abstract_Multiple_Linear_Regression 
{
private:
    /** Cached QR decomposition of X matrix */
    QR_Decomposition my_qr;

    /** Singularity threshold for QR decomposition */
    const double my_threshold;

public:
    /**
     * Create an empty OLS_Multiple_Linear_Regression instance.
     */
    OLS_Multiple_Linear_Regression() 
    {
        OLS_Multiple_Linear_Regression(0);
    }

    /**
     * Create an empty OLS_Multiple_Linear_Regression instance, using the given
     * singularity threshold for the QR decomposition.
     *
     * @param threshold the singularity threshold
     */
    OLS_Multiple_Linear_Regression(const double& threshold) : my_threshold{ threshold } {};


    /**
     * Loads model x and y sample data, overriding any previous sample.
     *
     * Computes and caches QR decomposition of the X matrix.
     * @param y the [n,1] array representing the y sample
     * @param x the [n,k] array representing the x sample
     * @ if the x and y array data are not
     *             compatible for the regression
     */
    void new_sample_data(std::vector<double> y, std::vector<std::vector<double>> x)  
    {
        validate_sample_data(x, y);
        new_y_sample_data(y);
        new_x_sample_data(x);
    }

    /**
     * {@inherit_doc}
     * <p>This implementation computes and caches the QR decomposition of the X matrix.</p>
     */
    //override
    void new_sample_data(std::vector<double> data, int nobs, int nvars) 
    {
        super.new_sample_data(data, nobs, nvars);
        my_qr = QR_Decomposition(get_x(), my_threshold);
    }

    /**
     * <p>Compute the "hat" matrix.
     * </p>
     * <p>The hat matrix is defined in terms of the design matrix X
     *  by X(X<sup>T</sup>X)<sup>-1</sup>X<sup>T</sup>
     * </p>
     * <p>The implementation here uses the QR decomposition to compute the
     * hat matrix as Q I<sub>p</sub>Q<sup>T</sup> where I<sub>p</sub> is the
     * p-dimensional identity matrix augmented by 0's.  This computational
     * formula is from "The Hat Matrix in Regression and ANOVA", * David C. Hoaglin and Roy E. Welsch, * <i>The American Statistician</i>, Vol. 32, No. 1 (Feb., 1978), pp. 17-22.
     * </p>
     * <p>Data for the model must have been successfully loaded using one of
     * the {@code new_sample_data} methods before invoking this method; otherwise
     * a {@code Null_Pointer_Exception} will be thrown.</p>
     *
     * @return the hat matrix
     * @Null_Pointer_Exception unless method {@code new_sample_data} has been
     * called beforehand.
     */
    Real_Matrix calculate_hat() 
    {
        // Create augmented identity matrix
        Real_Matrix Q = my_qr.get_q();
        const int p = my_qr.get_r().get_column_dimension();
        const int n = Q.get_column_dimension();
        // No try-catch or advertised  - NPE above if n < 3
        Array_2D_Row_Real_Matrix aug_i = Array_2D_Row_Real_Matrix(n, n);
        std::vector<std::vector<double>> aug_i_data = aug_i.get_data_ref();
        for (int i{}; i < n; i++) 
        {
            for (int j{}; j < n; j++) 
            {
                aug_i_data[i][j] = (i == j && i < p)
                    ? 1
                    : 0;
            }
        }

        // Compute and return Hat matrix
        // No DME advertised - args valid if we get here
        return Q.multiply(aug_i).multiply_transposed(Q);
    }

    /**
     * <p>Returns the sum of squared deviations of Y from its mean.</p>
     *
     * <p>If the model has no intercept term, <code>0</code> is used for the
     * mean of Y - i.e., what is returned is the sum of the squared Y values.</p>
     *
     * <p>The value returned by this method is the SSTO value used in
     * the {@link #calculate_r_squared() R-squared} computation.</p>
     *
     * @return SSTO - the total sum of squares
     * @Null_Pointer_Exception if the sample has not been set
     * @see #is_no_intercept()
     */
    double calculate_total_sum_of_squares() 
    {
        if (is_no_intercept()) 
        {
            return Stat_Utils.sum_sq(get_y().to_array());
        }
        return Second_Moment().evaluate(get_y().to_array());
    }

    /**
     * Returns the sum of squared residuals.
     *
     * @return residual sum of squares
     * @org.hipparchus.exception. if the design matrix is singular
     * @Null_Pointer_Exception if the data for the model have not been loaded
     */
    double calculate_residual_sum_of_squares() 
    {
        const Real_Vector residuals = calculate_residuals();
        // No advertised DME, args are valid
        return residuals.dot_product(residuals);
    }

    /**
     * Returns the R-_Squared statistic, defined by the formula <pre>
     * R<sup>2</sup> = 1 - SSR / SSTO
     * </pre>
     * where SSR is the {@link #calculate_residual_sum_of_squares() sum of squared residuals}
     * and SSTO is the {@link #calculate_total_sum_of_squares() total sum of squares}
     *
     * <p>If there is no variance in y, i.e., SSTO = 0, NaN is returned.</p>
     *
     * @return R-square statistic
     * @Null_Pointer_Exception if the sample has not been set
     * @org.hipparchus.exception. if the design matrix is singular
     */
    double calculate_r_squared() 
    {
        return 1 - calculate_residual_sum_of_squares() / calculate_total_sum_of_squares();
    }

    /**
     * <p>Returns the adjusted R-squared statistic, defined by the formula <pre>
     * R<sup>2</sup><sub>adj</sub> = 1 - [SSR (n - 1)] / [SSTO (n - p)]
     * </pre>
     * where SSR is the {@link #calculate_residual_sum_of_squares() sum of squared residuals}, * SSTO is the {@link #calculate_total_sum_of_squares() total sum of squares}, n is the number
     * of observations and p is the number of parameters estimated (including the intercept).</p>
     *
     * <p>If the regression is estimated without an intercept term, what is returned is <pre>
     * <code> 1 - (1 - {@link #calculate_r_squared()}) * (n / (n - p)) </code>
     * </pre></p>
     *
     * <p>If there is no variance in y, i.e., SSTO = 0, NaN is returned.</p>
     *
     * @return adjusted R-_Squared statistic
     * @Null_Pointer_Exception if the sample has not been set
     * @org.hipparchus.exception. if the design matrix is singular
     * @see #is_no_intercept()
     */
    double calculate_adjusted_r_squared() 
    {
        const double n = get_x().get_row_dimension();
        if (is_no_intercept()) 
        {
            return 1 - (1 - calculate_r_squared()) * (n / (n - get_x().get_column_dimension()));
        }

            return 1 - (calculate_residual_sum_of_squares() * (n - 1)) /
                (calculate_total_sum_of_squares() * (n - get_x().get_column_dimension()));
        
    }

protected:
    /**
     * {@inherit_doc}
     * <p>This implementation computes and caches the QR decomposition of the X matrix
     * once it is successfully loaded.</p>
     */
    //override
    void new_x_sample_data(std::vector<std::vector<double>> x) 
    {
        super.new_x_sample_data(x);
        my_qr = QR_Decomposition(get_x(), threshold);
    }

    /**
     * Calculates the regression coefficients using OLS.
     *
     * <p>Data for the model must have been successfully loaded using one of
     * the {@code new_sample_data} methods before invoking this method; otherwise
     * a {@code Null_Pointer_Exception} will be thrown.</p>
     *
     * @return beta
     * @org.hipparchus.exception. if the design matrix is singular
     * @Null_Pointer_Exception if the data for the model have not been loaded
     */
    //override
    Real_Vector calculate_beta() 
    {
        return my_qr.get_solver().solve(get_y());
    }

    /**
     * <p>Calculates the variance-covariance matrix of the regression parameters.
     * </p>
     * <p>Var(b) = (X<sup>T</sup>X)<sup>-1</sup>
     * </p>
     * <p>Uses QR decomposition to reduce (X<sup>T</sup>X)<sup>-1</sup>
     * to (R<sup>T</sup>R)<sup>-1</sup>, with only the top p rows of
     * R included, where p = the length of the beta vector.</p>
     *
     * <p>Data for the model must have been successfully loaded using one of
     * the {@code new_sample_data} methods before invoking this method; otherwise
     * a {@code Null_Pointer_Exception} will be thrown.</p>
     *
     * @return The beta variance-covariance matrix
     * @org.hipparchus.exception. if the design matrix is singular
     * @Null_Pointer_Exception if the data for the model have not been loaded
     */
    //override
    Real_Matrix calculate_beta_variance() 
    {
        int p = get_x().get_column_dimension();
        Real_Matrix Raug = qr.get_r().get_sub_matrix(0, p - 1 , 0, p - 1);
        Real_Matrix Rinv = LU_Decomposition(Raug).get_solver().get_inverse();
        return Rinv.multiply_transposed(Rinv);
    }

};