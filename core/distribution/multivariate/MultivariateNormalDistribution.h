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

#include <cmath>

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.linear.Eigen_Decomposition;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.random.Random_Generator;
//import org.hipparchus.random.Well19937c;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;

/**
 * Implementation of the multivariate normal (Gaussian) distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">
 * Multivariate normal distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/Multivariate_Normal_Distribution.html">
 * Multivariate normal distribution (MathWorld)</a>
 */
class Multivariate_Normal_Distribution : public AbstractMultivariate_Real_Distribution 
{
    /** Default singular matrix tolerance check value **/
    private static const double DEFAULT_TOLERANCE = Precision.EPSILON;

    /** Vector of means. */
    private const std::vector<double> means;
    /** Covariance matrix. */
    private const Real_Matrix covariance_matrix;
    /** The matrix inverse of the covariance matrix. */
    private const Real_Matrix covariance_matrix_inverse;
    /** The determinant of the covariance matrix. */
    private const double covariance_matrix_determinant;
    /** Matrix used in computation of samples. */
    private const Real_Matrix sampling_matrix;
    /** Inverse singular check tolerance when testing if invertable **/
    private const double singular_matrix_check_tolerance;

    /**
     * Creates a multivariate normal distribution with the given mean vector and
     * covariance matrix.
     * <br/>
     * The number of dimensions is equal to the length of the mean vector
     * and to the number of rows and columns of the covariance matrix.
     * It is frequently written as "p" in formulae.
     * <p>
     * <b>Note:</b> this constructor will implicitly create an instance of
     * {@link Well19937c} as random generator to be used for sampling only (see
     * {@link #sample()} and {@link #samplestatic_cast<int>(}). In case no sampling is
     * needed for the created distribution, it is advised to pass {@code NULL}
     * as random generator via the appropriate constructors to avoid the
     * additional initialisation overhead.
     *
     * @param means Vector of means.
     * @param covariances Covariance matrix.
     * @ if the arrays length are
     * inconsistent.
     * @ if the eigenvalue decomposition cannot
     * be performed on the provided covariance matrix.
     * @ if any of the eigenvalues is
     * negative.
     */
    public Multivariate_Normal_Distribution(const std::vector<double>& means, const std::vector<std::vector<double>>& covariances)
    {
        this(means, covariances, DEFAULT_TOLERANCE);
    }

    /**
     * Creates a multivariate normal distribution with the given mean vector and
     * covariance matrix.
     * <br/>
     * The number of dimensions is equal to the length of the mean vector
     * and to the number of rows and columns of the covariance matrix.
     * It is frequently written as "p" in formulae.
     * <p>
     * <b>Note:</b> this constructor will implicitly create an instance of
     * {@link Well19937c} as random generator to be used for sampling only (see
     * {@link #sample()} and {@link #samplestatic_cast<int>(}). In case no sampling is
     * needed for the created distribution, it is advised to pass {@code NULL}
     * as random generator via the appropriate constructors to avoid the
     * additional initialisation overhead.
     *
     * @param means Vector of means.
     * @param covariances Covariance matrix.
     * @param singular_matrix_check_tolerance Tolerance used during the singular matrix check before inversion
     * @ if the arrays length are
     * inconsistent.
     * @ if the eigenvalue decomposition cannot
     * be performed on the provided covariance matrix.
     * @ if any of the eigenvalues is
     * negative.
     */
    public Multivariate_Normal_Distribution(const std::vector<double> means, const std::vector<std::vector<double>> covariances, const double singular_matrix_check_tolerance)
    {
        this(new Well19937c(), means, covariances, singular_matrix_check_tolerance);
    }


    /**
     * Creates a multivariate normal distribution with the given mean vector and
     * covariance matrix.
     * <br/>
     * The number of dimensions is equal to the length of the mean vector
     * and to the number of rows and columns of the covariance matrix.
     * It is frequently written as "p" in formulae.
     *
     * @param rng Random Number Generator.
     * @param means Vector of means.
     * @param covariances Covariance matrix.
     * @ if the arrays length are
     * inconsistent.
     * @ if the eigenvalue decomposition cannot
     * be performed on the provided covariance matrix.
     * @ if any of the eigenvalues is
     * negative.
     */
    public Multivariate_Normal_Distribution(Random_Generator rng, const std::vector<double> means, const std::vector<std::vector<double>> covariances) 
    {
        this(rng, means, covariances, DEFAULT_TOLERANCE);
    }

    /**
     * Creates a multivariate normal distribution with the given mean vector and
     * covariance matrix.
     * <br/>
     * The number of dimensions is equal to the length of the mean vector
     * and to the number of rows and columns of the covariance matrix.
     * It is frequently written as "p" in formulae.
     *
     * @param rng Random Number Generator.
     * @param means Vector of means.
     * @param covariances Covariance matrix.
     * @param singular_matrix_check_tolerance Tolerance used during the singular matrix check before inversion
     * @ if the arrays length are
     * inconsistent.
     * @ if the eigenvalue decomposition cannot
     * be performed on the provided covariance matrix.
     * @ if any of the eigenvalues is
     * negative.
     */
    public Multivariate_Normal_Distribution(Random_Generator rng, const std::vector<double> means, const std::vector<std::vector<double>> covariances, const double singular_matrix_check_tolerance)
    {
        super(rng, means.size());

        const int dim = means.size();

        if (covariances.size() != dim) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, covariances.size(), dim);
        }

        for (int i{}; i < dim; i++) 
        {
            if (dim != covariances[i].size()) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, covariances[i].size(), dim);
            }
        }

        this.means = means.clone();
        this.singular_matrix_check_tolerance = singular_matrix_check_tolerance;

        covariance_matrix = Array_2D_Row_Real_Matrix(covariances);

        // Covariance matrix eigen decomposition.
        const Eigen_Decomposition cov_mat_dec = Eigen_Decomposition(covariance_matrix, singular_matrix_check_tolerance);

        // Compute and store the inverse.
        covariance_matrix_inverse = cov_mat_dec.get_solver().get_inverse();
        // Compute and store the determinant.
        covariance_matrix_determinant = cov_mat_dec.get_determinant();

        // Eigenvalues of the covariance matrix.
        const std::vector<double> cov_mat_eigenvalues = cov_mat_dec.get_real_eigenvalues();

        for (int i{}; i < cov_mat_eigenvalues.size(); i++) 
        {
            if (cov_mat_eigenvalues[i] < 0) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_DEFINITE_MATRIX);
            }
        }

        // Matrix where each column is an eigenvector of the covariance matrix.
        const Array_2D_Row_Real_Matrix cov_mat_eigenvectors = Array_2D_Row_Real_Matrix(dim, dim);
        for (const int& v = 0; v < dim; v++) 
        {
            const std::vector<double> evec = cov_mat_dec.get_eigenvector(v).to_array();
            cov_mat_eigenvectors.set_column(v, evec);
        }

        const Real_Matrix tmp_matrix = cov_mat_eigenvectors.transpose();

        // Scale each eigenvector by the square root of its eigenvalue.
        for (int row{}; row < dim; row++) 
        {
            const double factor = std::sqrt(cov_mat_eigenvalues[row]);
            for (int col{};  col < dim; col++) 
            {
                tmp_matrix.multiply_entry(row, col, factor);
            }
        }

        sampling_matrix = cov_mat_eigenvectors.multiply(tmp_matrix);
    }

    /**
     * Gets the mean vector.
     *
     * @return the mean vector.
     */
    public std::vector<double> get_means() 
    {
        return means.clone();
    }

    /**
     * Gets the covariance matrix.
     *
     * @return the covariance matrix.
     */
    public Real_Matrix get_covariances() 
    {
        return covariance_matrix.copy();
    }

    /**
     * Gets the current setting for the tolerance check used during singular checks before inversion
     * @return tolerance
     */
    public double get_singular_matrix_check_tolerance() { return singular_matrix_check_tolerance; }

    /** {@inherit_doc} */
    //override
    public double density(const std::vector<double>& vals)  
    {
        const int dim = get_dimension();
        if (vals.size() != dim) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, vals.size(), dim);
        }

        return std::pow(2 * std::numbers::pi, -0.5 * dim) *
            std::pow(covariance_matrix_determinant, -0.5) *
            get_exponent_term(vals);
    }

    /**
     * Gets the square root of each element on the diagonal of the covariance
     * matrix.
     *
     * @return the standard deviations.
     */
    public std::vector<double> get_standard_deviations() 
    {
        const int dim = get_dimension();
        const std::vector<double> std = std::vector<double>(dim];
        const std::vector<std::vector<double>> s = covariance_matrix.get_data();
        for (int i{}; i < dim; i++) 
        {
            std[i] = std::sqrt(s[i][i]);
        }
        return std;
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> sample() 
    {
        const int dim = get_dimension();
        const auto normal_vals = std::vector<double>(dim);

        for (int i{}; i < dim; i++) 
        {
            normal_vals[i] = random.next_gaussian();
        }

        const std::vector<double>& vals = sampling_matrix.operate(normal_vals);

        for (int i{}; i < dim; i++) 
        {
            vals[i] += means[i];
        }

        return vals;
    }

    /**
     * Computes the term used in the exponent (see definition of the distribution).
     *
     * @param values Values at which to compute density.
     * @return the multiplication factor of density calculations.
     */
    private double get_exponent_term(const std::vector<double>& values) 
    {
        const auto centered = std::vector<double>(values.size());
        for (int i{}; i < centered.size(); i++) 
        {
            centered[i] = values[i] - get_means()[i];
        }
        const std::vector<double> pre_multiplied = covariance_matrix_inverse.pre_multiply(centered);
        double sum{};
        for (int i{}; i < pre_multiplied.size(); i++) 
        {
            sum += pre_multiplied[i] * centered[i];
        }
        return std::exp(-0.5 * sum);
    }
};