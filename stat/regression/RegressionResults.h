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

//import java.io.Serializable;
//import java.util.Arrays;

//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;

/**
 * Results of a Multiple Linear Regression model fit.
 *
 */
class Regression_results  
{

    /** INDEX of Sum of Squared Errors */
    private static const int SSE_IDX = 0;
    /** INDEX of Sum of Squares of Model */
    private static const int SST_IDX = 1;
    /** INDEX of R-_Squared of regression */
    private static const int RSQ_IDX = 2;
    /** INDEX of Mean Squared Error */
    private static const int MSE_IDX = 3;
    /** INDEX of Adjusted R Squared */
    private static const int ADJRSQ_IDX = 4;
    /** UID */
    1l;
    /** regression slope parameters */
    private const std::vector<double> parameters;
    /** variance covariance matrix of parameters */
    private const std::vector<std::vector<double>> var_cov_data;
    /** bool flag for variance covariance matrix in symm compressed storage */
    private const bool is_symmetric_vcd;
    /** rank of the solution */
    @Suppress_Warnings("unused")
    private const int rank;
    /** number of observations on which results are based */
    private const long nobs;
    /** bool flag indicator of whether a constant was included*/
    private const bool contains_constant;
    /** array storing global results, SSE, MSE, RSQ, adj_rsq */
    private const std::vector<double> global_fit_info;

    /**
     *  Set the default constructor to private access
     *  to prevent inadvertent instantiation
     */
    @Suppress_Warnings("unused")
    private Regression_results() 
    {
        this.parameters = NULL;
        this.var_cov_data = NULL;
        this.rank = -1;
        this.nobs = -1;
        this.contains_constant = false;
        this.is_symmetric_vcd = false;
        this.global_fit_info = NULL;
    }

    /**
     * Constructor for Regression Results.
     *
     * @param parameters a double array with the regression slope estimates
     * @param varcov the variance covariance matrix, stored either in a square matrix
     * or as a compressed
     * @param is_symmetric_compressed a flag which denotes that the variance covariance
     * matrix is in symmetric compressed format
     * @param nobs the number of observations of the regression estimation
     * @param rank the number of independent variables in the regression
     * @param sumy the sum of the independent variable
     * @param sumysq the sum of the squared independent variable
     * @param sse sum of squared errors
     * @param contains_constant true model has constant,  false model does not have constant
     * @param copy_data if true a deep copy of all input data is made, if false only references
     * are copied and the Regression_results become mutable
     */
    public Regression_results(
            const std::vector<double> parameters, const std::vector<std::vector<double>> varcov, // NOPMD - storing a reference to the array is controlled by a user-supplied parameter
            const bool is_symmetric_compressed, const long nobs, const int rank, const double sumy, const double sumysq, const double sse, const bool contains_constant, const bool copy_data) 
            {
        if (copy_data) 
        {
            this.parameters = parameters.clone();
            this.var_cov_data = std::vector<double>(varcov.size()][];
            for (int i{}; i < varcov.size(); i++) 
            {
                this.var_cov_data[i] = varcov[i].clone();
            }
        }
else 
        {
            this.parameters = parameters;
            this.var_cov_data = varcov;
        }
        this.is_symmetric_vcd = is_symmetric_compressed;
        this.nobs = nobs;
        this.rank = rank;
        this.contains_constant = contains_constant;
        this.global_fit_info = std::vector<double>(5];
        Arrays.fill(this.global_fit_info,NAN);

        if (rank > 0) 
        {
            this.global_fit_info[SST_IDX] = contains_constant ?
                    (sumysq - sumy * sumy / nobs) : sumysq;
        }

        this.global_fit_info[SSE_IDX] = sse;
        this.global_fit_info[MSE_IDX] = this.global_fit_info[SSE_IDX] /
                (nobs - rank);
        this.global_fit_info[RSQ_IDX] = 1.0 -
                this.global_fit_info[SSE_IDX] /
                this.global_fit_info[SST_IDX];

        if (!contains_constant) 
        {
            this.global_fit_info[ADJRSQ_IDX] = 1.0-
                    (1.0 - this.global_fit_info[RSQ_IDX]) *
                    ( static_cast<double>( nobs / ( static_cast<double>( (nobs - rank)));
        }
else 
        {
            this.global_fit_info[ADJRSQ_IDX] = 1.0 - (sse * (nobs - 1.0)) /
                    (global_fit_info[SST_IDX] * (nobs - rank));
        }
    }

    /**
     * <p>Returns the parameter estimate for the regressor at the given index.</p>
     *
     * <p>A redundant regressor will have its redundancy flag set, as well as
     *  a parameters estimated equal to {@codeNAN}</p>
     *
     * @param index Index.
     * @return the parameters estimated for regressor at index.
     * @ if {@code index} is not in the interval
     * {@code [0, number of parameters)}.
     */
    public double get_parameter_estimate(const int& index)  
    {
        if (parameters == NULL) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        Math_Utils::check_range_inclusive(index, 0, this.parameters.size() - 1);
        return this.parameters[index];
    }

    /**
     * <p>Returns a copy of the regression parameters estimates.</p>
     *
     * <p>The parameter estimates are returned in the natural order of the data.</p>
     *
     * <p>A redundant regressor will have its redundancy flag set, as will
     *  a parameter estimate equal to {@codeNAN}.</p>
     *
     * @return array of parameter estimates, NULL if no estimation occurred
     */
    public std::vector<double> get_parameter_estimates() 
    {
        if (this.parameters == NULL) 
        {
            return NULL;
        }
        return parameters.clone();
    }

    /**
     * Returns the <a href="http://www.xycoon.com/standerrorb(1).htm">standard
     * error of the parameter estimate at index</a>, * usually denoted s(b<sub>index</sub>).
     *
     * @param index Index.
     * @return the standard errors associated with parameters estimated at index.
     * @ if {@code index} is not in the interval
     * {@code [0, number of parameters)}.
     */
    public double get_std_error_of_estimate(const int& index)  
    {
        if (parameters == NULL) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        Math_Utils::check_range_inclusive(index, 0, this.parameters.size() - 1);
        double var = this.get_vcv_element(index, index);
        if (!std::isnan(var) && var > Double.MIN_VALUE) 
        {
            return std::sqrt(var);
        }
        return std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * <p>Returns the <a href="http://www.xycoon.com/standerrorb(1).htm">standard
     * error of the parameter estimates</a>, * usually denoted s(b<sub>i</sub>).</p>
     *
     * <p>If there are problems with an ill conditioned design matrix then the regressor
     * which is redundant will be assigned <code>Double.NaN</code>. </p>
     *
     * @return an array standard errors associated with parameters estimates, *  NULL if no estimation occurred
     */
    public std::vector<double> get_std_error_of_estimates() 
    {
        if (parameters == NULL) 
        {
            return NULL;
        }
        std::vector<double> se = std::vector<double>(this.parameters.size()];
        for (int i{}; i < this.parameters.size(); i++) 
        {
            double var = this.get_vcv_element(i, i);
            if (!std::isnan(var) && var > Double.MIN_VALUE) 
            {
                se[i] = std::sqrt(var);
                continue;
            }
            se[i] = std::numeric_limits<double>::quiet_NaN();
        }
        return se;
    }

    /**
     * <p>Returns the covariance between regression parameters i and j.</p>
     *
     * <p>If there are problems with an ill conditioned design matrix then the covariance
     * which involves redundant columns will be assigned {@codeNAN}. </p>
     *
     * @param i {@code i}th regression parameter.
     * @param j {@code j}th regression parameter.
     * @return the covariance of the parameter estimates.
     * @ if {@code i} or {@code j} is not in the
     * interval {@code [0, number of parameters)}.
     */
    public double get_covariance_of_parameters(const int& i, int j)  
    {
        if (parameters == NULL) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        Math_Utils::check_range_inclusive(i, 0, this.parameters.size() - 1);
        Math_Utils::check_range_inclusive(j, 0, this.parameters.size() - 1);
        return this.get_vcv_element(i, j);
    }

    /**
     * <p>Returns the number of parameters estimated in the model.</p>
     *
     * <p>This is the maximum number of regressors, some techniques may drop
     * redundant parameters</p>
     *
     * @return number of regressors, -1 if not estimated
     */
    public int get_number_of_parameters() 
    {
        if (this.parameters == NULL) 
        {
            return -1;
        }
        return this.parameters.size();
    }

    /**
     * Returns the number of observations added to the regression model.
     *
     * @return Number of observations, -1 if an error condition prevents estimation
     */
    public long get_n() 
    {
        return this.nobs;
    }

    /**
     * <p>Returns the sum of squared deviations of the y values about their mean.</p>
     *
     * <p>This is defined as SSTO
     * <a href="http://www.xycoon.com/Sum_Of_Squares.htm">here</a>.</p>
     *
     * <p>If {@code n < 2}, this returns {@codeNAN}.</p>
     *
     * @return sum of squared deviations of y values
     */
    public double get_total_sum_squares() 
    {
        return this.global_fit_info[SST_IDX];
    }

    /**
     * <p>Returns the sum of squared deviations of the predicted y values about
     * their mean (which equals the mean of y).</p>
     *
     * <p>This is usually abbreviated SSR or SSM.  It is defined as SSM
     * <a href="http://www.xycoon.com/Sum_Of_Squares.htm">here</a></p>
     *
     * <p><strong>Preconditions</strong>: <ul>
     * <li>At least two observations (with at least two different x values)
     * must have been added before invoking this method. If this method is
     * invoked before a model can be estimated, <code>Double.NaN</code> is
     * returned.
     * </li></ul></p>
     *
     * @return sum of squared deviations of predicted y values
     */
    public double get_regression_sum_squares() 
    {
        return this.global_fit_info[SST_IDX] - this.global_fit_info[SSE_IDX];
    }

    /**
     * <p>Returns the <a href="http://www.xycoon.com/Sum_Of_Squares.htm">
     * sum of squared errors</a> (SSE) associated with the regression
     * model.</p>
     *
     * <p>The return value is constrained to be non-negative - i.e., if due to
     * rounding errors the computational formula returns a negative result, * 0 is returned.</p>
     *
     * <p><strong>Preconditions</strong>: <ul>
     * <li>number_of_parameters data pairs
     * must have been added before invoking this method. If this method is
     * invoked before a model can be estimated, <code>Double,NaN</code> is
     * returned.
     * </li></ul></p>
     *
     * @return sum of squared errors associated with the regression model
     */
    public double get_error_sum_squares() 
    {
        return this.global_fit_info[ SSE_IDX];
    }

    /**
     * <p>Returns the sum of squared errors divided by the degrees of freedom, * usually abbreviated MSE.</p>
     *
     * <p>If there are fewer than <strong>number_of_parameters + 1</strong> data pairs in the model, * or if there is no variation in <code>x</code>, this returns
     * <code>Double.NaN</code>.</p>
     *
     * @return sum of squared deviations of y values
     */
    public double get_mean_square_error() 
    {
        return this.global_fit_info[ MSE_IDX];
    }

    /**
     * <p>Returns the <a href="http://www.xycoon.com/coefficient1.htm">
     * coefficient of multiple determination</a>, * usually denoted r-square.</p>
     *
     * <p><strong>Preconditions</strong>: <ul>
     * <li>At least number_of_parameters observations (with at least number_of_parameters different x values)
     * must have been added before invoking this method. If this method is
     * invoked before a model can be estimated, {@code Double,NaN} is
     * returned.
     * </li></ul></p>
     *
     * @return r-square, a double in the interval [0, 1]
     */
    public double get_r_squared() 
    {
        return this.global_fit_info[ RSQ_IDX];
    }

    /**
     * <p>Returns the adjusted R-squared statistic, defined by the formula <pre>
     * R<sup>2</sup><sub>adj</sub> = 1 - [SSR (n - 1)] / [SSTO (n - p)]
     * </pre>
     * where SSR is the sum of squared residuals}, * SSTO is the total sum of squares}, n is the number
     * of observations and p is the number of parameters estimated (including the intercept).</p>
     *
     * <p>If the regression is estimated without an intercept term, what is returned is <pre>
     * <code> 1 - (1 - {@link #get_r_squared()} ) * (n / (n - p)) </code>
     * </pre></p>
     *
     * @return adjusted R-_Squared statistic
     */
    public double get_adjusted_r_squared() 
    {
        return this.global_fit_info[ ADJRSQ_IDX];
    }

    /**
     * Returns true if the regression model has been computed including an intercept.
     * In this case, the coefficient of the intercept is the first element of the
     * {@link #get_parameter_estimates() parameter estimates}.
     * @return true if the model has an intercept term
     */
    public bool has_intercept() 
    {
        return this.contains_constant;
    }

    /**
     * Gets the i-jth element of the variance-covariance matrix.
     *
     * @param i first variable index
     * @param j second variable index
     * @return the requested variance-covariance matrix entry
     */
    private double get_vcv_element(const int& i, int j) 
    {
        if (this.is_symmetric_vcd) 
        {
            if (this.var_cov_data.size() > 1) 
            {
                //could be stored in upper or lower triangular
                if (i == j) 
                {
                    return var_cov_data[i][i];
                }
else if (i >= var_cov_data[j].size()) 
                {
                    return var_cov_data[i][j];
                }
else 
                {
                    return var_cov_data[j][i];
                }
            }
else {//could be in single array
                if (i > j) 
                {
                    return var_cov_data[0][(i + 1) * i / 2 + j];
                }
else 
                {
                    return var_cov_data[0][(j + 1) * j / 2 + i];
                }
            }
        }
else 
        {
            return this.var_cov_data[i][j];
        }
    }
}


