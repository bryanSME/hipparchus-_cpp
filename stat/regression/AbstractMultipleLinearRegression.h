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

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.linear.Array_Real_Vector;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.stat.Localized_Stat_Formats;
//import org.hipparchus.stat.descriptive.moment.Variance;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
#include <vector>
#include "../../core/util/MathUtils.h"
#include "MultipleLinearRegression.h"
#include "../../core/linear/RealMatrix.h"
#include "../../core/linear/RealVector.h"

/**
 * Abstract base class for implementations of Multiple_Linear_Regression.
 */
class Abstract_Multiple_Linear_Regression : Multiple_Linear_Regression 
{
private:
    /** X sample data. */
    Real_Matrix my_x_matrix;

    /** Y sample data. */
    Real_Vector my_y_vector;

    /** Whether or not the regression model includes an intercept.  True means no intercept. */
    bool my_no_intercept;

protected:
    /**
     * @return the X sample data.
     */
    Real_Matrix get_x() const
    {
        return my_x_matrix;
    }

    /**
     * @return the Y sample data.
     */
    Real_Vector get_y() const
    {
        return my_y_vector;
    }

public:
    /**
     * @return true if the model has no intercept term; false otherwise
     */
    bool is_no_intercept() const
    {
        return my_no_intercept;
    }

    /**
     * @param no_intercept true means the model is to be estimated without an intercept term
     */
    void set_no_intercept(bool no_intercept) 
    {
        my_no_intercept = no_intercept;
    }

    /**
     * <p>Loads model x and y sample data from a flat input array, overriding any previous sample.
     * </p>
     * <p>Assumes that rows are concatenated with y values first in each row.  For example, an input
     * <code>data</code> array containing the sequence of values (1, 2, 3, 4, 5, 6, 7, 8, 9) with
     * <code>nobs = 3</code> and <code>nvars = 2</code> creates a regression dataset with two
     * independent variables, as below:
     * <pre>
     *   y   x[0]  x[1]
     *   --------------
     *   1     2     3
     *   4     5     6
     *   7     8     9
     * </pre>
     * </p>
     * <p>Note that there is no need to add an initial unitary column (column of 1's) when
     * specifying a model including an intercept term.  If {@link #is_no_intercept()} is <code>true</code>, * the X matrix will be created without an initial column of "1"s; otherwise this column will
     * be added.
     * </p>
     * <p>Throws Illegal_Argument_Exception if any of the following preconditions fail:
     * <ul><li><code>data</code> cannot be NULL</li>
     * <li><code>data.size() = nobs * (nvars + 1)</li>
     * <li><code>nobs &gt; nvars</code></li></ul>
     * </p>
     *
     * @param data input data array
     * @param nobs number of observations (rows)
     * @param nvars number of independent variables (columns, not counting y)
     * @Null_Argument_Exception if the data array is NULL
     * @ if the length of the data array is not equal
     * to <code>nobs * (nvars + 1)</code>
     * @ if <code>nobs</code> is less than
     * <code>nvars + 1</code>
     */
    void new_sample_data(const std::vector<double>& data, const int& nobs, const int& nvars) 
    {
        //Math_Utils::check_not_null(data, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
        Math_Utils::check_dimension(data.size(), nobs * (nvars + 1));
        if (nobs <= nvars) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_OBSERVED_POINTS_IN_SAMPLE, nobs, nvars + 1);
        }
        std::vector<double> y = std::vector<double>(nobs];
        const int cols = my_no_intercept ? nvars: nvars + 1;
        std::vector<std::vector<double>> x = std::vector<double>(nobs][cols];
        int pointer = 0;
        for (int i{}; i < nobs; i++) 
        {
            y[i] = data[pointer++];
            if (!my_no_intercept) 
            {
                x[i][0] = 1.0;
            }
            for (int j = my_no_intercept ? 0 : 1; j < cols; j++)
            {
                x[i][j] = data[pointer++];
            }
        }
        my_x_matrix = Array_2D_Row_Real_Matrix(x);
        my_y_vector = Array_Real_Vector(y);
    }

    /**
     * Loads y sample data, overriding any previous data.
     *
     * @param y the array representing the y sample
     * @Null_Argument_Exception if y is NULL
     * @ if y is empty
     */
    protected void new_y_sample_data(std::vector<double> y) 
    {
        if (y == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        if (y.size() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        my_y_vector = Array_Real_Vector(y);
    }

    /**
     * <p>Loads x sample data, overriding any previous data.
     * </p>
     * The input <code>x</code> array should have one row for each sample
     * observation, with columns corresponding to independent variables.
     * For example, if <pre>
     * <code> x = std::vector<std::vector<double>> {{1, 2}, {3, 4}, {5, 6}} </code></pre>
     * then <code>set_x_sample_data(x) </code> results in a model with two independent
     * variables and 3 observations:
     * <pre>
     *   x[0]  x[1]
     *   ----------
     *     1    2
     *     3    4
     *     5    6
     * </pre>
     * </p>
     * <p>Note that there is no need to add an initial unitary column (column of 1's) when
     * specifying a model including an intercept term.
     * </p>
     * @param x the rectangular array representing the x sample
     * @Null_Argument_Exception if x is NULL
     * @ if x is empty
     * @ if x is not rectangular
     */
    protected void new_x_sample_data(std::vector<std::vector<double>> x) 
    {
        if (x == NULL) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        if (x.size() == 0) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        if (my_no_intercept)
        {
            my_x_matrix = Array_2D_Row_Real_Matrix(x, true);
        }
else { // Augment design matrix with initial unitary column
            const int& n_vars = x[0].size();
            const std::vector<std::vector<double>> x_aug = std::vector<double>(x.size()][n_vars + 1];
            for (int i{}; i < x.size(); i++) 
            {
                Math_Utils::check_dimension(x[i].size(), n_vars);
                x_aug[i][0] = 1.0;
                System.arraycopy(x[i], 0, x_aug[i], 1, n_vars);
            }
            my_x_matrix = Array_2D_Row_Real_Matrix(x_aug, false);
        }
    }

    /**
     * Validates sample data.  Checks that
     * <ul><li>Neither x nor y is NULL or empty;</li>
     * <li>The length (i.e. number of rows) of x equals the length of y</li>
     * <li>x has at least one more row than it has columns (i.e. there is
     * sufficient data to estimate regression coefficients for each of the
     * columns in x plus an intercept.</li>
     * </ul>
     *
     * @param x the [n,k] array representing the x data
     * @param y the [n,1] array representing the y data
     * @Null_Argument_Exception if {@code x} or {@code y} is NULL
     * @ if {@code x} and {@code y} do not
     * have the same length
     * @ if {@code x} or {@code y} are zero-length
     * @ if the number of rows of {@code x}
     * is not larger than the number of columns + 1 if the model has an intercept;
     * or the number of columns if there is no intercept term
     */
    protected void validate_sample_data(std::vector<std::vector<double>> x, std::vector<double> y)  
    {
        if ((x == NULL) || (y == NULL)) 
        {
            throw std::exception("not implemented");
            //throw Null_Argument_Exception();
        }
        Math_Utils::check_dimension(x.size(), y.size());
        if (x.size() == 0) {  // Must be no y data either
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        if (x[0].size() + (my_no_intercept ? 0 : 1) > x.size())
        {
            throw std::exception("not implemented");
            //throw (Localized_Stat_Formats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS, x.size(), x[0].size());
        }
    }

    /**
     * Validates that the x data and covariance matrix have the same
     * number of rows and that the covariance matrix is square.
     *
     * @param x the [n,k] array representing the x sample
     * @param covariance the [n,n] array representing the covariance matrix
     * @ if the number of rows in x is not equal
     * to the number of rows in covariance
     * @ if the covariance matrix is not square
     */
    protected void validate_covariance_data(std::vector<std::vector<double>> x, std::vector<std::vector<double>> covariance) 
    {
        Math_Utils::check_dimension(x.size(), covariance.size());
        if (covariance.size() > 0 && covariance.size() != covariance[0].size()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, covariance.size(), covariance[0].size());
        }
    }

    /**
     * {@inherit_doc}
     */
    //override
    public std::vector<double> estimate_regression_parameters() 
    {
        Real_Vector b = calculate_beta();
        return b.to_array();
    }

    /**
     * {@inherit_doc}
     */
    //override
    public std::vector<double> estimate_residuals() 
    {
        Real_Vector b = calculate_beta();
        Real_Vector e = y_vector.subtract(x_matrix.operate(b));
        return e.to_array();
    }

    /**
     * {@inherit_doc}
     */
    //override
    public std::vector<std::vector<double>> estimate_regression_parameters_variance() 
    {
        return calculate_beta_variance().get_data();
    }

    /**
     * {@inherit_doc}
     */
    //override
    public std::vector<double> estimate_regression_parameters_standard_errors() 
    {
        std::vector<std::vector<double>> beta_variance = estimate_regression_parameters_variance();
        double sigma = calculate_error_variance();
        int length = beta_variance[0].size();
        std::vector<double> result = std::vector<double>(length];
        for (int i{}; i < length; i++) 
        {
            result[i] = std::sqrt(sigma * beta_variance[i][i]);
        }
        return result;
    }

    /**
     * {@inherit_doc}
     */
    //override
    public double estimate_regressand_variance() 
    {
        return calculate_y_variance();
    }

    /**
     * Estimates the variance of the error.
     *
     * @return estimate of the error variance
     */
    public double estimate_error_variance() 
    {
        return calculate_error_variance();

    }

    /**
     * Estimates the standard error of the regression.
     *
     * @return regression standard error
     */
    public double estimate_regression_standard_error() 
    {
        return std::sqrt(estimate_error_variance());
    }

    /**
     * Calculates the beta of multiple linear regression in matrix notation.
     *
     * @return beta
     */
    protected virtual Real_Vector calculate_beta();

    /**
     * Calculates the beta variance of multiple linear regression in matrix
     * notation.
     *
     * @return beta variance
     */
    protected virtual Real_Matrix calculate_beta_variance();


    /**
     * Calculates the variance of the y values.
     *
     * @return Y variance
     */
    protected double calculate_y_variance() 
    {
        return Variance().evaluate(my_y_vector.to_array());
    }

    /**
     * <p>Calculates the variance of the error term.</p>
     * Uses the formula <pre>
     * var(u) = u &middot; u / (n - k)
     * </pre>
     * where n and k are the row and column dimensions of the design
     * matrix X.
     *
     * @return error variance estimate
     */
    protected double calculate_error_variance() 
    {
        Real_Vector residuals = calculate_residuals();
        return residuals.dot_product(residuals) /
               (my_x_matrix.get_row_dimension() - my_x_matrix.get_column_dimension());
    }

    /**
     * Calculates the residuals of multiple linear regression in matrix
     * notation.
     *
     * <pre>
     * u = y - X * b
     * </pre>
     *
     * @return The residuals [n,1] matrix
     */
    protected Real_Vector calculate_residuals() 
    {
        Real_Vector b = calculate_beta();
        return my_y_vector.subtract(my_x_matrix.operate(b));
    }

};