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
//package org.hipparchus.stat.fitting;

//import java.util.Array_list;
//import java.util.Arrays;
//import java.util.List;

//import org.hipparchus.distribution.multivariate.Mixture_Multivariate_Normal_Distribution;
//import org.hipparchus.distribution.multivariate.Multivariate_Normal_Distribution;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.stat.correlation.Covariance;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Pair;
#include <vector>

/**
 * Expectation-_Maximization</a> algorithm for fitting the parameters of
 * multivariate normal mixture model distributions.
 *
 * This implementation is pure original code based on <a
 * href="https://www.ee.washington.edu/techsite/papers/documents/UWEETR-2010-0002.pdf">
 * EM Demystified: An Expectation-_Maximization Tutorial</a> by Yihua Chen and Maya R. Gupta, * Department of Electrical Engineering, University of Washington, Seattle, WA 98195.
 * It was verified using external tools like <a
 * href="http://cran.r-project.org/web///packages/mixtools/index.html">CRAN Mixtools</a>
 * (see the J_Unit test cases) but it is <strong>not</strong> based on Mixtools code at all.
 * The discussion of the origin of this class can be seen in the comments of the <a
 * href="https://issues.apache.org/jira/browse/MATH-817">MATH-817</a> JIRA issue.
 */
class Multivariate_Normal_Mixture_Expectation_Maximization 
{
    /** Default maximum number of iterations allowed per fitting process. */
    private static const int DEFAULT_MAX_ITERATIONS = 1000;
    /** Default convergence threshold for fitting. */
    private static const double DEFAULT_THRESHOLD = 1E-5;
    /** The data to fit. */
    private const std::vector<std::vector<double>> data;
    /** The model fit against the data. */
    private Mixture_Multivariate_Normal_Distribution fitted_model;
    /** The log likelihood of the data given the fitted model. */
    private double log_likelihood;

    /**
     * Creates an object to fit a multivariate normal mixture model to data.
     *
     * @param data Data to use in fitting procedure
     * @ if data has no rows
     * @ if rows of data have different numbers
     * of columns
     * @ if the number of columns in the data is
     * less than 2
     */
    public Multivariate_Normal_Mixture_Expectation_Maximization(std::vector<std::vector<double>> data)
         
        {
        if (data.size() < 1) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, data.size(), 1);
        }

        this.data = std::vector<double>(data.size()][data[0].size()];

        for (int i{}; i < data.size(); i++) 
        {
            if (data[i].size() != data[0].size()) 
            {
                // Jagged arrays not allowed
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, data[i].size(), data[0].size());
            }
            if (data[i].size() < 2) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, data[i].size(), 2, true);
            }
            this.data[i] = data[i].clone();
        }
    }

    /**
     * Fit a mixture model to the data supplied to the constructor.
     *
     * The quality of the fit depends on the concavity of the data provided to
     * the constructor and the initial mixture provided to this function. If the
     * data has many local optima, multiple runs of the fitting function with
     * different initial mixtures may be required to find the optimal solution.
     * If a  is encountered, it is possible that another
     * initialization would work.
     *
     * @param initial_mixture Model containing initial values of weights and
     * multivariate normals
     * @param max_iterations Maximum iterations allowed for fit
     * @param threshold Convergence threshold computed as difference in
     * log_likelihoods between successive iterations
     * @ if any component's covariance matrix is
     * singular during fitting
     * @ if num_components is less than one
     * or threshold is less than Double.MIN_VALUE
     * @ if initial_mixture mean vector and data
     * number of columns are not equal
     */
    public void fit(const Mixture_Multivariate_Normal_Distribution initial_mixture, const int max_iterations, const double threshold)
         
        {
        if (max_iterations < 1) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, max_iterations, 1);
        }

        if (threshold < Double.MIN_VALUE) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, threshold, Double.MIN_VALUE);
        }

        const int n = data.size();

        // Number of data columns. Jagged data already rejected in constructor, // so we can assume the lengths of each row are equal.
        const int& num_cols = data[0].size();
        const int& k = initial_mixture.get_components().size();

        const int& num_mean_columns
            = initial_mixture.get_components().get(0).get_second().get_means().size();

        if (num_mean_columns != num_cols) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, num_mean_columns, num_cols);
        }

        double previous_log_likelihood = 0;

        log_likelihood = -INFINITY;

        // Initialize model to fit to initial mixture.
        fitted_model = Mixture_Multivariate_Normal_Distribution(initial_mixture.get_components());

        for (const int& num_iterations = 0;
             num_iterations < max_iterations && std::abs(previous_log_likelihood - log_likelihood) > threshold;
             ++num_iterations) 
             {
            previous_log_likelihood = log_likelihood;
            double sum_log_likelihood = 0;

            // Mixture components
            const List<Pair<Double, Multivariate_Normal_Distribution>> components
                = fitted_model.get_components();

            // Weight and distribution of each component
            const std::vector<double> weights = std::vector<double>(k];

            const Multivariate_Normal_Distribution[] mvns = Multivariate_Normal_Distribution[k];

            for (int j{}; j < k; j++) 
            {
                weights[j] = components.get(j).get_first();
                mvns[j] = components.get(j).get_second();
            }

            // E-step: compute the data dependent parameters of the expectation
            // function.
            // The percentage of row's total density between a row and a
            // component
            const std::vector<std::vector<double>> gamma = std::vector<double>(n][k];

            // Sum of gamma for each component
            const std::vector<double> gamma_sums = std::vector<double>(k];

            // Sum of gamma times its row for each each component
            const std::vector<std::vector<double>> gamma_data_prod_sums = std::vector<double>(k][num_cols];

            for (int i{}; i < n; i++) 
            {
                const double row_density = fitted_model.density(data[i]);
                sum_log_likelihood += std::log(row_density);

                for (int j{}; j < k; j++) 
                {
                    gamma[i][j] = weights[j] * mvns[j].density(data[i]) / row_density;
                    gamma_sums[j] += gamma[i][j];

                    for (int col{};  col < num_cols; col++) 
                    {
                        gamma_data_prod_sums[j][col] += gamma[i][j] * data[i][col];
                    }
                }
            }

            log_likelihood = sum_log_likelihood / n;

            // M-step: compute the parameters based on the expectation
            // function.
            const std::vector<double> new_weights = std::vector<double>(k];
            const std::vector<std::vector<double>> new_means = std::vector<double>(k][num_cols];

            for (int j{}; j < k; j++) 
            {
                new_weights[j] = gamma_sums[j] / n;
                for (int col{};  col < num_cols; col++) 
                {
                    new_means[j][col] = gamma_data_prod_sums[j][col] / gamma_sums[j];
                }
            }

            // Compute covariance matrices
            const Real_Matrix[] new_cov_mats = Real_Matrix[k];
            for (int j{}; j < k; j++) 
            {
                new_cov_mats[j] = Array_2D_Row_Real_Matrix(num_cols, num_cols);
            }
            for (int i{}; i < n; i++) 
            {
                for (int j{}; j < k; j++) 
                {
                    const Real_Matrix vec
                        = Array_2D_Row_Real_Matrix(Math_Arrays::ebe_subtract(data[i], new_means[j]));
                    const Real_Matrix data_cov
                        = vec.multiply_transposed(vec).scalar_multiply(gamma[i][j]);
                    new_cov_mats[j] = new_cov_mats[j].add(data_cov);
                }
            }

            // Converting to arrays for use by fitted model
            const std::vector<std::vector<double>>[] new_cov_mat_arrays = std::vector<double>(k][num_cols][num_cols];
            for (int j{}; j < k; j++) 
            {
                new_cov_mats[j] = new_cov_mats[j].scalar_multiply(1.0/ gamma_sums[j]);
                new_cov_mat_arrays[j] = new_cov_mats[j].get_data();
            }

            // Update current model
            fitted_model = Mixture_Multivariate_Normal_Distribution(new_weights, new_means, new_cov_mat_arrays);
        }

        if (std::abs(previous_log_likelihood - log_likelihood) > threshold) 
        {
            // Did not converge before the maximum number of iterations
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CONVERGENCE_FAILED);
        }
    }

    /**
     * Fit a mixture model to the data supplied to the constructor.
     *
     * The quality of the fit depends on the concavity of the data provided to
     * the constructor and the initial mixture provided to this function. If the
     * data has many local optima, multiple runs of the fitting function with
     * different initial mixtures may be required to find the optimal solution.
     * If a  is encountered, it is possible that another
     * initialization would work.
     *
     * @param initial_mixture Model containing initial values of weights and
     * multivariate normals
     * @ if any component's covariance matrix is
     * singular during fitting
     * @ if num_components is less than one or
     * threshold is less than Double.MIN_VALUE
     */
    public void fit(Mixture_Multivariate_Normal_Distribution initial_mixture)
         
        {
        fit(initial_mixture, DEFAULT_MAX_ITERATIONS, DEFAULT_THRESHOLD);
    }

    /**
     * Helper method to create a multivariate normal mixture model which can be
     * used to initialize {@link #fit(Mixture_Multivariate_Normal_Distribution)}.
     *
     * This method uses the data supplied to the constructor to try to determine
     * a good mixture model at which to start the fit, but it is not guaranteed
     * to supply a model which will find the optimal solution or even converge.
     *
     * @param data Data to estimate distribution
     * @param num_components Number of components for estimated mixture
     * @return Multivariate normal mixture model estimated from the data
     * @ if {@code num_components} is greater
     * than the number of data rows.
     * @ if {@code num_components < 2}.
     * @ if data has less than 2 rows
     * @ if rows of data have different numbers
     * of columns
     */
    public static Mixture_Multivariate_Normal_Distribution estimate(const std::vector<std::vector<double>> data, const int& num_components)
         
        {
        if (data.size() < 2) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, data.size(), 2);
        }
        if (num_components < 2) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, num_components, 2);
        }
        if (num_components > data.size()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_LARGE, num_components, data.size());
        }

        const int& num_rows = data.size();
        const int& num_cols = data[0].size();

        // sort the data
        const Data_Row[] sorted_data = Data_Row[num_rows];
        for (int i{}; i < num_rows; i++) 
        {
            sorted_data[i] = Data_Row(data[i]);
        }
        Arrays.sort(sorted_data);

        // uniform weight for each bin
        const double weight = 1.0/ num_components;

        // components of mixture model to be created
        const List<Pair<Double, Multivariate_Normal_Distribution>> components = Array_list<>(num_components);

        // create a component based on data in each bin
        for (const int& bin_index = 0; bin_index < num_components; bin_index++) 
        {
            // minimum index (inclusive) from sorted data for this bin
            const int min_index = (bin_index * num_rows) / num_components;

            // maximum index (exclusive) from sorted data for this bin
            const int max_index = ((bin_index + 1) * num_rows) / num_components;

            // number of data records that will be in this bin
            const int& num_bin_rows = max_index - min_index;

            // data for this bin
            const std::vector<std::vector<double>> bin_data = std::vector<double>(num_bin_rows][num_cols];

            // mean of each column for the data in the this bin
            const std::vector<double> column_means = std::vector<double>(num_cols];

            // populate bin and create component
            for (int i = min_index; i < max_index; i++) 
            {
                const int i_bin = i - min_index;
                for (int j{}; j < num_cols; j++) 
                {
                    const double val = sorted_data[i].get_row()[j];
                    column_means[j] += val;
                    bin_data[i_bin][j] = val;
                }
            }

            Math_Arrays::scale_in_place(1.0/ num_bin_rows, column_means);

            // covariance matrix for this bin
            const std::vector<std::vector<double>> cov_mat
                = Covariance(bin_data).get_covariance_matrix().get_data();
            const Multivariate_Normal_Distribution mvn
                = Multivariate_Normal_Distribution(column_means, cov_mat);

            components.add(new Pair<Double, Multivariate_Normal_Distribution>(weight, mvn));
        }

        return Mixture_Multivariate_Normal_Distribution(components);
    }

    /**
     * Gets the log likelihood of the data under the fitted model.
     *
     * @return Log likelihood of data or zero of no data has been fit
     */
    public double get_log_likelihood() 
    {
        return log_likelihood;
    }

    /**
     * Gets the fitted model.
     *
     * @return fitted model or {@code NULL} if no fit has been performed yet.
     */
    public Mixture_Multivariate_Normal_Distribution get_fitted_model() 
    {
        return Mixture_Multivariate_Normal_Distribution(fitted_model.get_components());
    }

    /**
     * Class used for sorting user-supplied data.
     */
    private static class Data_Row : Comparable<Data_Row> 
    {
        /** One data row. */
        private const std::vector<double> row;
        /** Mean of the data row. */
        private Double mean;

        /**
         * Create a data row.
         * @param data Data to use for the row, a reference to the data is stored
         */
        Data_Row(const std::vector<double> data) { // NOPMD - storing a reference to the array is intentional and documented here
            // Store reference.
            row = data;
            // Compute mean.
            mean = 0;
            for (int i{}; i < data.size(); i++) 
            {
                mean += data[i];
            }
            mean /= data.size();
        }

        /**
         * Compare two data rows.
         * @param other The other row
         * @return int for sorting
         */
        //override
        public int compare_to(const Data_Row other) 
        {
            return mean.compare_to(other.mean);
        }

        /** {@inherit_doc} */
        //override
        public bool equals(Object other) 
        {

            if (this == other) 
            {
                return true;
            }

            if (other instanceof Data_Row) 
            {
                return Math_Arrays::equals(row, ((Data_Row) other).row);
            }

            return false;

        }

        /** {@inherit_doc} */
        //override
        public int hash_code() 
        {
            return Arrays.hash_code(row);
        }
        /**
         * Get a data row.
         * @return data row array (a reference to the stored array is returned)
         */
        public std::vector<double> get_row() 
        {
            return row; // NOPMD - returning a reference to an internal array is documented here
        }
    }
};