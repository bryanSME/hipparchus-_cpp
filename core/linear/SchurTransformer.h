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
#include "MatrixUtils.h"
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;
#include "RealMatrix.h"

/**
 * Class transforming a general real matrix to Schur form.
 * <p>A m &times; m matrix A can be written as the product of three matrices: A = P
 * &times; T &times; P<sup>T</sup> with P an orthogonal matrix and T an quasi-triangular
 * matrix. Both P and T are m &times; m matrices.</p>
 * <p>Transformation to Schur form is often not a goal by itself, but it is an
 * intermediate step in more general decomposition algorithms like
 * {@link Eigen_Decomposition eigen decomposition}. This class is therefore
 * intended for internal use by the library and is not public. As a consequence
 * of this explicitly limited scope, many methods directly returns references to
 * internal arrays, not copies.</p>
 * <p>This class is based on the method hqr2 in class Eigenvalue_Decomposition
 * from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library.</p>
 *
 * @see <a href="http://mathworld.wolfram.com/SchurDecomposition.html">Schur Decomposition - MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Schur_decomposition">Schur Decomposition - Wikipedia</a>
 * @see <a href="http://en.wikipedia.org/wiki/Householder_transformation">Householder Transformations</a>
 */
class Schur_Transformer 
{
    /** Maximum allowed iterations for convergence of the transformation. */
    private static const int MAX_ITERATIONS{ 100 };

    /** P matrix. */
    private const double matrix_p[][];
    /** T matrix. */
    private const double matrixstd::vector<std::vector<T>>;
    /** Cached value of P. */
    private Real_Matrix cached_p;
    /** Cached value of T. */
    private Real_Matrix cached_t;
    /** Cached value of PT. */
    private Real_Matrix cached_pt;

    /** Epsilon criteria taken from JAMA code (originally was 2^-52). */
    private const double epsilon = Precision.EPSILON;

    /**
     * Build the transformation to Schur form of a general real matrix.
     *
     * @param matrix matrix to transform
     * @ if the matrix is not square
     */
    Schur_Transformer(const Real_Matrix matrix) 
    {
        if (!matrix.is_square()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, matrix.get_row_dimension(), matrix.get_column_dimension());
        }

        Hessenberg_Transformer transformer = Hessenberg_Transformer(matrix);
        matrix_t = transformer.get_h().get_data();
        matrix_p = transformer.get_p().get_data();
        cached_t = NULL;
        cached_p = NULL;
        cached_pt = NULL;

        // transform matrix
        transform();
    }

    /**
     * Returns the matrix P of the transform.
     * <p>P is an orthogonal matrix, i.e. its inverse is also its transpose.</p>
     *
     * @return the P matrix
     */
    public Real_Matrix get_p() 
    {
        if (cached_p == NULL) 
        {
            cached_p = Matrix_Utils::create_real_matrix(matrix_p);
        }
        return cached_p;
    }

    /**
     * Returns the transpose of the matrix P of the transform.
     * <p>P is an orthogonal matrix, i.e. its inverse is also its transpose.</p>
     *
     * @return the transpose of the P matrix
     */
    public Real_Matrix get_p_t() 
    {
        if (cached_pt == NULL) 
        {
            cached_pt = get_p().transpose();
        }

        // return the cached matrix
        return cached_pt;
    }

    /**
     * Returns the quasi-triangular Schur matrix T of the transform.
     *
     * @return the T matrix
     */
    public Real_Matrix get_t() 
    {
        if (cached_t == NULL) 
        {
            cached_t = Matrix_Utils::create_real_matrix(matrix_t);
        }

        // return the cached matrix
        return cached_t;
    }

    /**
     * Transform original matrix to Schur form.
     * @Math_Illegal_State_Exception if the transformation does not converge
     */
    private void transform() 
    {
        const int n = matrix_t.size();

        // compute matrix norm
        const double norm = get_norm();

        // shift information
        const Shift_Info shift = Shift_Info();

        // Outer loop over eigenvalue index
        int iteration = 0;
        int iu = n - 1;
        while (iu >= 0) 
        {

            // Look for single small sub-diagonal element
            const int il = find_small_sub_diagonal_element(iu, norm);

            // Check for convergence
            if (il == iu) 
            {
                // One root found
                matrix_t[iu][iu] += shift.ex_shift;
                iu--;
                iteration = 0;
            }
            else if (il == iu - 1) 
            {
                // Two roots found
                double p = (matrix_t[iu - 1][iu - 1] - matrix_t[iu][iu]) / 2.0;
                double q = p * p + matrix_t[iu][iu - 1] * matrix_t[iu - 1][iu];
                matrix_t[iu][iu] += shift.ex_shift;
                matrix_t[iu - 1][iu - 1] += shift.ex_shift;

                if (q >= 0) 
                {
                    double z = std::sqrt(std::abs(q));
                    if (p >= 0) 
                    {
                        z = p + z;
                    }
                    else 
                    {
                        z = p - z;
                    }
                    const double x = matrix_t[iu][iu - 1];
                    const double s = std::abs(x) + std::abs(z);
                    p = x / s;
                    q = z / s;
                    const double r = std::sqrt(p * p + q * q);
                    p /= r;
                    q /= r;

                    // Row modification
                    for (int j = iu - 1; j < n; j++) 
                    {
                        z = matrix_t[iu - 1][j];
                        matrix_t[iu - 1][j] = q * z + p * matrix_t[iu][j];
                        matrix_t[iu][j] = q * matrix_t[iu][j] - p * z;
                    }

                    // Column modification
                    for (int i{}; i <= iu; i++) 
                    {
                        z = matrix_t[i][iu - 1];
                        matrix_t[i][iu - 1] = q * z + p * matrix_t[i][iu];
                        matrix_t[i][iu] = q * matrix_t[i][iu] - p * z;
                    }

                    // Accumulate transformations
                    for (int i{}; i <= n - 1; i++) 
                    {
                        z = matrix_p[i][iu - 1];
                        matrix_p[i][iu - 1] = q * z + p * matrix_p[i][iu];
                        matrix_p[i][iu] = q * matrix_p[i][iu] - p * z;
                    }
                }
                iu -= 2;
                iteration = 0;
            }
            else 
            {
                // No convergence yet
                compute_shift(il, iu, iteration, shift);

                // stop transformation after too many iterations
                ++iteration;
                if (iteration > MAX_ITERATIONS) 
                {
                    throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CONVERGENCE_FAILED, MAX_ITERATIONS);
                }

                // the initial house_holder vector for the QR step
                const std::vector<double> h_vec = std::vector<double>(3];

                const int im = init_q_r_step(il, iu, shift, h_vec);
                perform_double_q_r_step(il, im, iu, shift, h_vec);
            }
        }
    }

    /**
     * Computes the L1 norm of the (quasi-)triangular matrix T.
     *
     * @return the L1 norm of matrix T
     */
    private double get_norm() 
    {
        double norm = 0.0;
        for (int i{}; i < matrix_t.size(); i++) 
        {
            // as matrix T is (quasi-)triangular, also take the sub-diagonal element into account
            for (int j = std::max(i - 1, 0); j < matrix_t.size(); j++) 
            {
                norm += std::abs(matrix_t[i][j]);
            }
        }
        return norm;
    }

    /**
     * Find the first small sub-diagonal element and returns its index.
     *
     * @param start_idx the starting index for the search
     * @param norm the L1 norm of the matrix
     * @return the index of the first small sub-diagonal element
     */
    private int find_small_sub_diagonal_element(const int start_idx, const double norm) 
    {
        int l = start_idx;
        while (l > 0) 
        {
            double s = std::abs(matrix_t[l - 1][l - 1]) + std::abs(matrix_t[l][l]);
            if (s == 0.0) 
            {
                s = norm;
            }
            if (std::abs(matrix_t[l][l - 1]) < epsilon * s) 
            {
                break;
            }
            l--;
        }
        return l;
    }

    /**
     * Compute the shift for the current iteration.
     *
     * @param l the index of the small sub-diagonal element
     * @param idx the current eigenvalue index
     * @param iteration the current iteration
     * @param shift holder for shift information
     */
    private void compute_shift(const int l, const int idx, const int iteration, const Shift_Info shift) 
    {
        // Form shift
        shift.x = matrix_t[idx][idx];
        shift.y = shift.w = 0.0;
        if (l < idx) 
        {
            shift.y = matrix_t[idx - 1][idx - 1];
            shift.w = matrix_t[idx][idx - 1] * matrix_t[idx - 1][idx];
        }

        // Wilkinson's original ad hoc shift
        if (iteration == 10) 
        {
            shift.ex_shift += shift.x;
            for (int i{}; i <= idx; i++) 
            {
                matrix_t[i][i] -= shift.x;
            }
            const double s = std::abs(matrix_t[idx][idx - 1]) + std::abs(matrix_t[idx - 1][idx - 2]);
            shift.x = 0.75 * s;
            shift.y = 0.75 * s;
            shift.w = -0.4375 * s * s;
        }

        // MATLAB's ad hoc shift
        if (iteration == 30) 
        {
            double s = (shift.y - shift.x) / 2.0;
            s = s * s + shift.w;
            if (s > 0.0) 
            {
                s = std::sqrt(s);
                if (shift.y < shift.x) 
                {
                    s = -s;
                }
                s = shift.x - shift.w / ((shift.y - shift.x) / 2.0 + s);
                for (int i{}; i <= idx; i++) 
                {
                    matrix_t[i][i] -= s;
                }
                shift.ex_shift += s;
                shift.x = shift.y = shift.w = 0.964;
            }
        }
    }

    /**
     * Initialize the householder vectors for the QR step.
     *
     * @param il the index of the small sub-diagonal element
     * @param iu the current eigenvalue index
     * @param shift shift information holder
     * @param h_vec the initial house_holder vector
     * @return the start index for the QR step
     */
    private int init_q_r_step(const int& il, const int iu, const Shift_Info shift, std::vector<double> h_vec) 
    {
        // Look for two consecutive small sub-diagonal elements
        int im = iu - 2;
        while (im >= il) 
        {
            const double z = matrix_t[im][im];
            const double r = shift.x - z;
            double s = shift.y - z;
            h_vec[0] = (r * s - shift.w) / matrix_t[im + 1][im] + matrix_t[im][im + 1];
            h_vec[1] = matrix_t[im + 1][im + 1] - z - r - s;
            h_vec[2] = matrix_t[im + 2][im + 1];

            if (im == il) 
            {
                break;
            }

            const double lhs = std::abs(matrix_t[im][im - 1]) * (std::abs(h_vec[1]) + std::abs(h_vec[2]));
            const double rhs = std::abs(h_vec[0]) * (std::abs(matrix_t[im - 1][im - 1]) +
                                                        std::abs(z) +
                                                        std::abs(matrix_t[im + 1][im + 1]));

            if (lhs < epsilon * rhs) 
            {
                break;
            }
            im--;
        }

        return im;
    }

    /**
     * Perform a double QR step involving rows l:idx and columns m:n
     *
     * @param il the index of the small sub-diagonal element
     * @param im the start index for the QR step
     * @param iu the current eigenvalue index
     * @param shift shift information holder
     * @param h_vec the initial house_holder vector
     */
    private void perform_double_q_r_step(const int il, const int im, const int iu, const Shift_Info shift, const std::vector<double> h_vec) 
    {

        const int n = matrix_t.size();
        double p = h_vec[0];
        double q = h_vec[1];
        double r = h_vec[2];

        for (int k = im; k <= iu - 1; k++) 
        {
            bool notlast = k != (iu - 1);
            if (k != im) 
            {
                p = matrix_t[k][k - 1];
                q = matrix_t[k + 1][k - 1];
                r = notlast ? matrix_t[k + 2][k - 1] : 0.0;
                shift.x = std::abs(p) + std::abs(q) + std::abs(r);
                if (Precision::equals(shift.x, 0.0, epsilon)) 
                {
                    continue;
                }
                p /= shift.x;
                q /= shift.x;
                r /= shift.x;
            }
            double s = std::sqrt(p * p + q * q + r * r);
            if (p < 0.0) 
            {
                s = -s;
            }
            if (s != 0.0) 
            {
                if (k != im) 
                {
                    matrix_t[k][k - 1] = -s * shift.x;
                }
                else if (il != im) 
                {
                    matrix_t[k][k - 1] = -matrix_t[k][k - 1];
                }
                p += s;
                shift.x = p / s;
                shift.y = q / s;
                double z = r / s;
                q /= p;
                r /= p;

                // Row modification
                for (int j{ k }; j < n; j++) 
                {
                    p = matrix_t[k][j] + q * matrix_t[k + 1][j];
                    if (notlast) 
                    {
                        p += r * matrix_t[k + 2][j];
                        matrix_t[k + 2][j] -= p * z;
                    }
                    matrix_t[k][j] -= p * shift.x;
                    matrix_t[k + 1][j] -= p * shift.y;
                }

                // Column modification
                for (int i{}; i <= std::min(iu, k + 3); i++) 
                {
                    p = shift.x * matrix_t[i][k] + shift.y * matrix_t[i][k + 1];
                    if (notlast) 
                    {
                        p += z * matrix_t[i][k + 2];
                        matrix_t[i][k + 2] -= p * r;
                    }
                    matrix_t[i][k] -= p;
                    matrix_t[i][k + 1] -= p * q;
                }

                // Accumulate transformations
                const int high = matrix_t.size() - 1;
                for (int i{}; i <= high; i++) 
                {
                    p = shift.x * matrix_p[i][k] + shift.y * matrix_p[i][k + 1];
                    if (notlast) 
                    {
                        p += z * matrix_p[i][k + 2];
                        matrix_p[i][k + 2] -= p * r;
                    }
                    matrix_p[i][k] -= p;
                    matrix_p[i][k + 1] -= p * q;
                }
            }  // (s != 0)
        }  // k loop

        // clean up pollution due to round-off errors
        for (int i = im + 2; i <= iu; i++) 
        {
            matrix_t[i][i-2] = 0.0;
            if (i > im + 2) 
            {
                matrix_t[i][i-3] = 0.0;
            }
        }
    }

    /**
     * Internal data structure holding the current shift information.
     * Contains variable names as present in the original JAMA code.
     */
    private static class Shift_Info 
    {
        // CHECKSTYLE: stop all

        /** x shift info */
        double x;
        /** y shift info */
        double y;
        /** w shift info */
        double w;
        /** Indicates an exceptional shift. */
        double ex_shift;

        // CHECKSTYLE: resume all
    }
};