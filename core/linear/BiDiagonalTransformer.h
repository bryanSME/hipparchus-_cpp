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

//package org.hipparchus.linear;

//import org.hipparchus.util.FastMath;
#include <vector>
#include "MatrixUtils.h"

/**
 * Class transforming any matrix to bi-diagonal shape.
 * <p>Any m &times; n matrix A can be written as the product of three matrices:
 * A = U &times; B &times; V<sup>T</sup> with U an m &times; m orthogonal matrix, * B an m &times; n bi-diagonal matrix (lower diagonal if m &lt; n, upper diagonal
 * otherwise), and V an n &times; n orthogonal matrix.</p>
 * <p>Transformation to bi-diagonal shape is often not a goal by itself, but it is
 * an intermediate step in more general decomposition algorithms like {@link
 * Singular_Value_Decomposition Singular Value Decomposition}. This class is therefore
 * intended for internal use by the library and is not public. As a consequence of
 * this explicitly limited scope, many methods directly returns references to
 * internal arrays, not copies.</p>
 */
class Bi_Diagonal_Transformer 
{
private:
    /** Householder vectors. */
    const std::vector<std::vector<double>> my_householder_vectors;

    /** Main diagonal. */
    const std::vector<double> my_main;

    /** Secondary diagonal. */
    const std::vector<double> my_secondary;

    /** Cached value of U. */
    Real_Matrix my_cached_u;

    /** Cached value of B. */
    Real_Matrix my_cached_b;

    /** Cached value of V. */
    Real_Matrix my_cached_v;

    /**
     * Transform original matrix to upper bi-diagonal form.
     * <p>Transformation is done using alternate Householder transforms
     * on columns and rows.</p>
     */
    void transform_to_upper_bi_diagonal()
    {
        const int m = householder_vectors.size();
        const int n = householder_vectors[0].size();
        for (int k{}; k < n; k++)
        {

            //zero-out a column
            double x_norm_sqr{};
            for (int i = k; i < m; ++i)
            {
                const double c = householder_vectors[i][k];
                x_norm_sqr += c * c;
            }
            const auto h_k = householder_vectors[k];
            const double a = (h_k[k] > 0)
                ? -std::sqrt(x_norm_sqr)
                : std::sqrt(x_norm_sqr);
            main[k] = a;
            if (a != 0.0)
            {
                h_k[k] -= a;
                for (int j{ k + 1 }; j < n; ++j)
                {
                    double alpha = 0;
                    for (int i = k; i < m; ++i)
                    {
                        const std::vector<double> h_i = householder_vectors[i];
                        alpha -= h_i[j] * h_i[k];
                    }
                    alpha /= a * householder_vectors[k][k];
                    for (int i = k; i < m; ++i)
                    {
                        const std::vector<double> h_i = householder_vectors[i];
                        h_i[j] -= alpha * h_i[k];
                    }
                }
            }

            if (k < n - 1)
            {
                //zero-out a row
                x_norm_sqr = 0;
                for (int j{ k + 1 }; j < n; ++j)
                {
                    const double c = h_k[j];
                    x_norm_sqr += c * c;
                }
                const double b = (h_k[k + 1] > 0) ? -std::sqrt(x_norm_sqr) : std::sqrt(x_norm_sqr);
                secondary[k] = b;
                if (b != 0.0)
                {
                    h_k[k + 1] -= b;
                    for (int i = k + 1; i < m; ++i)
                    {
                        const std::vector<double> h_i = householder_vectors[i];
                        double beta = 0;
                        for (int j{ k + 1 }; j < n; ++j)
                        {
                            beta -= h_i[j] * h_k[j];
                        }
                        beta /= b * h_k[k + 1];
                        for (int j{ k + 1 }; j < n; ++j)
                        {
                            h_i[j] -= beta * h_k[j];
                        }
                    }
                }
            }

        }
    }

    /**
     * Transform original matrix to lower bi-diagonal form.
     * <p>Transformation is done using alternate Householder transforms
     * on rows and columns.</p>
     */
    void transform_to_lower_bi_diagonal()
    {
        const int m = householder_vectors.size();
        const int n = householder_vectors[0].size();
        for (int k{}; k < m; k++)
        {

            //zero-out a row
            const std::vector<double> h_k = householder_vectors[k];
            double x_norm_sqr{};
            for (int j{ k }; j < n; ++j)
            {
                const double c = h_k[j];
                x_norm_sqr += c * c;
            }
            const double& a = (h_k[k] > 0) ? -std::sqrt(x_norm_sqr) : std::sqrt(x_norm_sqr);
            main[k] = a;
            if (a != 0.0)
            {
                h_k[k] -= a;
                for (int i = k + 1; i < m; ++i)
                {
                    const std::vector<double> h_i = householder_vectors[i];
                    double alpha = 0;
                    for (int j{ k }; j < n; ++j)
                    {
                        alpha -= h_i[j] * h_k[j];
                    }
                    alpha /= a * householder_vectors[k][k];
                    for (int j{ k }; j < n; ++j)
                    {
                        h_i[j] -= alpha * h_k[j];
                    }
                }
            }

            if (k < m - 1)
            {
                //zero-out a column
                const std::vector<double> h_kp1 = householder_vectors[k + 1];
                x_norm_sqr = 0;
                for (int i = k + 1; i < m; ++i)
                {
                    const double c = householder_vectors[i][k];
                    x_norm_sqr += c * c;
                }
                const double b = (h_kp1[k] > 0) ? -std::sqrt(x_norm_sqr) : std::sqrt(x_norm_sqr);
                secondary[k] = b;
                if (b != 0.0)
                {
                    h_kp1[k] -= b;
                    for (int j{ k + 1 }; j < n; ++j)
                    {
                        double beta = 0;
                        for (int i = k + 1; i < m; ++i)
                        {
                            const std::vector<double> h_i = householder_vectors[i];
                            beta -= h_i[j] * h_i[k];
                        }
                        beta /= b * h_kp1[k];
                        for (int i = k + 1; i < m; ++i)
                        {
                            const std::vector<double> h_i = householder_vectors[i];
                            h_i[j] -= beta * h_i[k];
                        }
                    }
                }
            }

        }
    }

public:
    /**
     * Build the transformation to bi-diagonal shape of a matrix.
     * @param matrix the matrix to transform.
     */
    Bi_Diagonal_Transformer(Real_Matrix matrix) 
    {
        const int m = matrix.get_row_dimension();
        const int n = matrix.get_column_dimension();
        const int p = std::min(m, n);
        householder_vectors = matrix.get_data();
        main      = std::vector<double>(p];
        secondary = std::vector<double>(p - 1];
        cached_u   = NULL;
        cached_b   = NULL;
        cached_v   = NULL;

        // transform matrix
        if (m >= n) 
        {
            transform_to_upper_bi_diagonal();
        }
        else 
        {
            transform_to_lower_bi_diagonal();
        }

    }

    /**
     * Returns the matrix U of the transform.
     * <p>U is an orthogonal matrix, i.e. its transpose is also its inverse.</p>
     * @return the U matrix
     */
    Real_Matrix get_u() 
    {
        if (cached_u == NULL) 
        {

            const int m = householder_vectors.size();
            const int n = householder_vectors[0].size();
            const int p = main.size();
            const int diag_offset    = (m >= n) ? 0 : 1;
            const std::vector<double> diagonal = (m >= n) ? main : secondary;
            std::vector<std::vector<double>> ua = std::vector<double>(m][m];

            // fill up the part of the matrix not affected by Householder transforms
            for (int k = m - 1; k >= p; --k) 
            {
                ua[k][k] = 1;
            }

            // build up first part of the matrix by applying Householder transforms
            for (int k = p - 1; k >= diag_offset; --k) 
            {
                const std::vector<double> h_k = householder_vectors[k];
                ua[k][k] = 1;
                if (h_k[k - diag_offset] != 0.0) 
                {
                    for (int j{ k }; j < m; ++j) 
                    {
                        double alpha = 0;
                        for (int i = k; i < m; ++i) 
                        {
                            alpha -= ua[i][j] * householder_vectors[i][k - diag_offset];
                        }
                        alpha /= diagonal[k - diag_offset] * h_k[k - diag_offset];

                        for (int i = k; i < m; ++i) 
                        {
                            ua[i][j] += -alpha * householder_vectors[i][k - diag_offset];
                        }
                    }
                }
            }
            if (diag_offset > 0) 
            {
                ua[0][0] = 1;
            }
            cached_u = Matrix_Utils::create_real_matrix(ua);
        }

        // return the cached matrix
        return cached_u;
    }

    /**
     * Returns the bi-diagonal matrix B of the transform.
     * @return the B matrix
     */
    Real_Matrix get_b() 
    {

        if (cached_b == NULL) 
        {

            const int m = householder_vectors.size();
            const int n = householder_vectors[0].size();
            std::vector<std::vector<double>> ba = std::vector<double>(m][n];
            for (int i{}; i < main.size(); ++i) 
            {
                ba[i][i] = main[i];
                if (m < n) 
                {
                    if (i > 0) 
                    {
                        ba[i][i-1] = secondary[i - 1];
                    }
                }
                else 
                {
                    if (i < main.size() - 1) 
                    {
                        ba[i][i+1] = secondary[i];
                    }
                }
            }
            cached_b = Matrix_Utils::create_real_matrix(ba);
        }

        // return the cached matrix
        return cached_b;

    }

    /**
     * Returns the matrix V of the transform.
     * <p>V is an orthogonal matrix, i.e. its transpose is also its inverse.</p>
     * @return the V matrix
     */
    Real_Matrix get_v() 
    {

        if (cached_v == NULL) 
        {

            const int m = householder_vectors.size();
            const int n = householder_vectors[0].size();
            const int p = main.size();
            const int diag_offset    = (m >= n) ? 1 : 0;
            const std::vector<double> diagonal = (m >= n) ? secondary : main;
            std::vector<std::vector<double>> va = std::vector<double>(n][n];

            // fill up the part of the matrix not affected by Householder transforms
            for (int k{ n - 1 }; k >= p; --k) 
            {
                va[k][k] = 1;
            }

            // build up first part of the matrix by applying Householder transforms
            for (int k = p - 1; k >= diag_offset; --k) 
            {
                const std::vector<double> h_k = householder_vectors[k - diag_offset];
                va[k][k] = 1;
                if (h_k[k] != 0.0) 
                {
                    for (int j{ k }; j < n; ++j) 
                    {
                        double beta = 0;
                        for (int i = k; i < n; ++i) 
                        {
                            beta -= va[i][j] * h_k[i];
                        }
                        beta /= diagonal[k - diag_offset] * h_k[k];

                        for (int i = k; i < n; ++i) 
                        {
                            va[i][j] += -beta * h_k[i];
                        }
                    }
                }
            }
            if (diag_offset > 0) 
            {
                va[0][0] = 1;
            }
            cached_v = Matrix_Utils::create_real_matrix(va);
        }

        // return the cached matrix
        return cached_v;

    }

    /**
     * Get the Householder vectors of the transform.
     * <p>Note that since this class is only intended for internal use, * it returns directly a reference to its internal arrays, not a copy.</p>
     * @return the main diagonal elements of the B matrix
     */
    std::vector<std::vector<double>> get_householder_vectors_ref() const
    {
        return my_householder_vectors; // NOPMD - returning an internal array is intentional and documented here
    }

    /**
     * Get the main diagonal elements of the matrix B of the transform.
     * <p>Note that since this class is only intended for internal use, * it returns directly a reference to its internal arrays, not a copy.</p>
     * @return the main diagonal elements of the B matrix
     */
    std::vector<double> get_main_diagonal_ref() const
    {
        return my_main; // NOPMD - returning an internal array is intentional and documented here
    }

    /**
     * Get the secondary diagonal elements of the matrix B of the transform.
     * <p>Note that since this class is only intended for internal use, * it returns directly a reference to its internal arrays, not a copy.</p>
     * @return the secondary diagonal elements of the B matrix
     */
    std::vector<double> get_secondary_diagonal_ref() const
    {
        return my_secondary; // NOPMD - returning an internal array is intentional and documented here
    }

    /**
     * Check if the matrix is transformed to upper bi-diagonal.
     * @return true if the matrix is transformed to upper bi-diagonal
     */
    bool is_upper_bi_diagonal() const
    {
        return my_householder_vectors.size() >=  my_householder_vectors[0].size();
    }
};