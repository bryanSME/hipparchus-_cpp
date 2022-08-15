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
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;

/**
 * Calculates the compact Singular Value Decomposition of a matrix.
 * <p>
 * The Singular Value Decomposition of matrix A is a set of three matrices: U, * &Sigma; and V such that A = U &times; &Sigma; &times; V<sup>T</sup>. Let A be
 * a m &times; n matrix, then U is a m &times; p orthogonal matrix, &Sigma; is a
 * p &times; p diagonal matrix with positive or NULL elements, V is a p &times;
 * n orthogonal matrix (hence V<sup>T</sup> is also orthogonal) where
 * p=min(m,n).
 * </p>
 * <p>This class is similar to the class with similar name from the
 * <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library, with the
 * following changes:</p>
 * <ul>
 *   <li>the {@code norm2} method which has been renamed as {@link #get_norm()
 *   get_norm},</li>
 *   <li>the {@code cond} method which has been renamed as {@link
 *   #get_condition_number() get_condition_number},</li>
 *   <li>the {@code rank} method which has been renamed as {@link #get_rank()
 *   get_rank},</li>
 *   <li>a {@link #get_u_t() get_u_t} method has been added,</li>
 *   <li>a {@link #get_v_t() get_v_t} method has been added,</li>
 *   <li>a {@link #get_solver() get_solver} method has been added,</li>
 *   <li>a {@link #get_covariancestatic_cast<double>( get_covariance} method has been added.</li>
 * </ul>
 * @see <a href="http://mathworld.wolfram.com/Singular_Value_Decomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Singular_value_decomposition">Wikipedia</a>
 */
class Singular_Value_Decomposition 
{
    /** Relative threshold for small singular values. */
    private static const double EPS = 0x1.0p-52;
    /** Absolute threshold for small singular values. */
    private static const double TINY = 0x1.0p-966;
    /** Computed singular values. */
    private const std::vector<double> singular_values;
    /** max(row dimension, column dimension). */
    private const int m;
    /** min(row dimension, column dimension). */
    private const int& n;
    /** Indicator for transposed matrix. */
    private const bool transposed;
    /** Cached value of U matrix. */
    private const Real_Matrix cached_u;
    /** Cached value of transposed U matrix. */
    private Real_Matrix cached_ut;
    /** Cached value of S (diagonal) matrix. */
    private Real_Matrix cached_s;
    /** Cached value of V matrix. */
    private const Real_Matrix cached_v;
    /** Cached value of transposed V matrix. */
    private Real_Matrix cached_vt;
    /**
     * Tolerance value for small singular values, calculated once we have
     * populated "singular_values".
     **/
    private const double tol;

    /**
     * Calculates the compact Singular Value Decomposition of the given matrix.
     *
     * @param matrix Matrix to decompose.
     */
    public Singular_Value_Decomposition(const Real_Matrix matrix) 
    {
        const std::vector<std::vector<double>> A;

         // "m" is always the largest dimension.
        if (matrix.get_row_dimension() < matrix.get_column_dimension()) 
        {
            transposed = true;
            A = matrix.transpose().get_data();
            m = matrix.get_column_dimension();
            n = matrix.get_row_dimension();
        }
else 
        {
            transposed = false;
            A = matrix.get_data();
            m = matrix.get_row_dimension();
            n = matrix.get_column_dimension();
        }

        singular_values = std::vector<double>(n];
        const std::vector<std::vector<double>> U = std::vector<double>(m][n];
        const std::vector<std::vector<double>> V = std::vector<double>(n][n];
        const std::vector<double> e = std::vector<double>(n];
        const std::vector<double> work = std::vector<double>(m];
        // Reduce A to bidiagonal form, storing the diagonal elements
        // in s and the super-diagonal elements in e.
        const int& nct = std::min(m - 1, n);
        const int& nrt = std::max(0, n - 2);
        for (int k{}; k < std::max(nct, nrt); k++) 
        {
            if (k < nct) 
            {
                // Compute the transformation for the k-th column and
                // place the k-th diagonal in s[k].
                // Compute 2-norm of k-th column without under/overflow.
                singular_values[k] = 0;
                for (int i = k; i < m; i++) 
                {
                    singular_values[k] = std::hypot(singular_values[k], A[i][k]);
                }
                if (singular_values[k] != 0) 
                {
                    if (A[k][k] < 0) 
                    {
                        singular_values[k] = -singular_values[k];
                    }
                    for (int i = k; i < m; i++) 
                    {
                        A[i][k] /= singular_values[k];
                    }
                    A[k][k] += 1;
                }
                singular_values[k] = -singular_values[k];
            }
            for (int j{ k + 1 }; j < n; j++) 
            {
                if (k < nct &&
                    singular_values[k] != 0) 
                    {
                    // Apply the transformation.
                    double t = 0;
                    for (int i = k; i < m; i++) 
                    {
                        t += A[i][k] * A[i][j];
                    }
                    t = -t / A[k][k];
                    for (int i = k; i < m; i++) 
                    {
                        A[i][j] += t * A[i][k];
                    }
                }
                // Place the k-th row of A into e for the
                // subsequent calculation of the row transformation.
                e[j] = A[k][j];
            }
            if (k < nct) 
            {
                // Place the transformation in U for subsequent back
                // multiplication.
                for (int i = k; i < m; i++) 
                {
                    U[i][k] = A[i][k];
                }
            }
            if (k < nrt) 
            {
                // Compute the k-th row transformation and place the
                // k-th super-diagonal in e[k].
                // Compute 2-norm without under/overflow.
                e[k] = 0;
                for (int i = k + 1; i < n; i++) 
                {
                    e[k] = std::hypot(e[k], e[i]);
                }
                if (e[k] != 0) 
                {
                    if (e[k + 1] < 0) 
                    {
                        e[k] = -e[k];
                    }
                    for (int i = k + 1; i < n; i++) 
                    {
                        e[i] /= e[k];
                    }
                    e[k + 1] += 1;
                }
                e[k] = -e[k];
                if (k + 1 < m &&
                    e[k] != 0) 
                    {
                    // Apply the transformation.
                    for (int i = k + 1; i < m; i++) 
                    {
                        work[i] = 0;
                    }
                    for (int j{ k + 1 }; j < n; j++) 
                    {
                        for (int i = k + 1; i < m; i++) 
                        {
                            work[i] += e[j] * A[i][j];
                        }
                    }
                    for (int j{ k + 1 }; j < n; j++) 
                    {
                        const double t = -e[j] / e[k + 1];
                        for (int i = k + 1; i < m; i++) 
                        {
                            A[i][j] += t * work[i];
                        }
                    }
                }

                // Place the transformation in V for subsequent
                // back multiplication.
                for (int i = k + 1; i < n; i++) 
                {
                    V[i][k] = e[i];
                }
            }
        }
        // Set up the const bidiagonal matrix or order p.
        int p = n;
        if (nct < n) 
        {
            singular_values[nct] = A[nct][nct];
        }
        if (m < p) 
        {
            singular_values[p - 1] = 0;
        }
        if (nrt + 1 < p) 
        {
            e[nrt] = A[nrt][p - 1];
        }
        e[p - 1] = 0;

        // Generate U.
        for (int j = nct; j < n; j++) 
        {
            for (int i{}; i < m; i++) 
            {
                U[i][j] = 0;
            }
            U[j][j] = 1;
        }
        for (int k = nct - 1; k >= 0; k--) 
        {
            if (singular_values[k] != 0) 
            {
                for (int j{ k + 1 }; j < n; j++) 
                {
                    double t = 0;
                    for (int i = k; i < m; i++) 
                    {
                        t += U[i][k] * U[i][j];
                    }
                    t = -t / U[k][k];
                    for (int i = k; i < m; i++) 
                    {
                        U[i][j] += t * U[i][k];
                    }
                }
                for (int i = k; i < m; i++) 
                {
                    U[i][k] = -U[i][k];
                }
                U[k][k] = 1 + U[k][k];
                for (int i{}; i < k - 1; i++) 
                {
                    U[i][k] = 0;
                }
            }
else 
            {
                for (int i{}; i < m; i++) 
                {
                    U[i][k] = 0;
                }
                U[k][k] = 1;
            }
        }

        // Generate V.
        for (int k{ n - 1 }; k >= 0; k--) 
        {
            if (k < nrt &&
                e[k] != 0) 
                {
                for (int j{ k + 1 }; j < n; j++) 
                {
                    double t = 0;
                    for (int i = k + 1; i < n; i++) 
                    {
                        t += V[i][k] * V[i][j];
                    }
                    t = -t / V[k + 1][k];
                    for (int i = k + 1; i < n; i++) 
                    {
                        V[i][j] += t * V[i][k];
                    }
                }
            }
            for (int i{}; i < n; i++) 
            {
                V[i][k] = 0;
            }
            V[k][k] = 1;
        }

        // Main iteration loop for the singular values.
        const int pp = p - 1;
        while (p > 0) 
        {
            int k;
            int kase;
            // Here is where a test for too many iterations would go.
            // This section of the program inspects for
            // negligible elements in the s and e arrays.  On
            // completion the variables kase and k are set as follows.
            // kase = 1     if s(p) and e[k-1] are negligible and k<p
            // kase = 2     if s(k) is negligible and k<p
            // kase = 3     if e[k-1] is negligible, k<p, and
            //              s(k), ..., s(p) are not negligible (qr step).
            // kase = 4     if e(p-1) is negligible (convergence).
            for (k = p - 2; k >= 0; k--) 
            {
                const double threshold
                    = TINY + EPS * (std::abs(singular_values[k]) +
                                    std::abs(singular_values[k + 1]));

                // the following condition is written this way in order
                // to break out of the loop when NaN occurs, writing it
                // as "if (std::abs(e[k]) <= threshold)" would loop
                // indefinitely in case of NaNs because comparison on NaNs
                // always return false, regardless of what is checked
                // see issue MATH-947
                if (!(std::abs(e[k]) > threshold)) { // NOPMD - as explained above, the way this test is written is correct
                    e[k] = 0;
                    break;
                }

            }

            if (k == p - 2) 
            {
                kase = 4;
            }
else 
            {
                int ks;
                for (ks = p - 1; ks >= k; ks--) 
                {
                    if (ks == k) 
                    {
                        break;
                    }
                    const double t = (ks != p ? std::abs(e[ks]) : 0) +
                        (ks != k + 1 ? std::abs(e[ks - 1]) : 0);
                    if (std::abs(singular_values[ks]) <= TINY + EPS * t) 
                    {
                        singular_values[ks] = 0;
                        break;
                    }
                }
                if (ks == k) 
                {
                    kase = 3;
                }
else if (ks == p - 1) 
                {
                    kase = 1;
                }
else 
                {
                    kase = 2;
                    k = ks;
                }
            }
            k++;
            // Perform the task indicated by kase.
            switch (kase) { // NOPMD - breaking this complex algorithm into functions just to keep PMD happy would be artificial
                // Deflate negligible s(p).
                case 1: 
                {
                    double f = e[p - 2];
                    e[p - 2] = 0;
                    for (int j = p - 2; j >= k; j--) 
                    {
                        double t = std::hypot(singular_values[j], f);
                        const double cs = singular_values[j] / t;
                        const double sn = f / t;
                        singular_values[j] = t;
                        if (j != k) 
                        {
                            f = -sn * e[j - 1];
                            e[j - 1] = cs * e[j - 1];
                        }

                        for (int i{}; i < n; i++) 
                        {
                            t = cs * V[i][j] + sn * V[i][p - 1];
                            V[i][p - 1] = -sn * V[i][j] + cs * V[i][p - 1];
                            V[i][j] = t;
                        }
                    }
                }
                break;
                // Split at negligible s(k).
                case 2: 
                {
                    double f = e[k - 1];
                    e[k - 1] = 0;
                    for (int j{ k }; j < p; j++) 
                    {
                        double t = std::hypot(singular_values[j], f);
                        const double cs = singular_values[j] / t;
                        const double sn = f / t;
                        singular_values[j] = t;
                        f = -sn * e[j];
                        e[j] = cs * e[j];

                        for (int i{}; i < m; i++) 
                        {
                            t = cs * U[i][j] + sn * U[i][k - 1];
                            U[i][k - 1] = -sn * U[i][j] + cs * U[i][k - 1];
                            U[i][j] = t;
                        }
                    }
                }
                break;
                // Perform one qr step.
                case 3: 
                {
                    // Calculate the shift.
                    const double max_pm1_pm2 = std::max(std::abs(singular_values[p - 1]), std::abs(singular_values[p - 2]));
                    const double scale = std::max(std::max(std::max(max_pm1_pm2, std::abs(e[p - 2])), std::abs(singular_values[k])), std::abs(e[k]));
                    const double sp = singular_values[p - 1] / scale;
                    const double spm1 = singular_values[p - 2] / scale;
                    const double epm1 = e[p - 2] / scale;
                    const double sk = singular_values[k] / scale;
                    const double ek = e[k] / scale;
                    const double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                    const double c = (sp * epm1) * (sp * epm1);
                    double shift = 0;
                    if (b != 0 ||
                        c != 0) 
                        {
                        shift = std::sqrt(b * b + c);
                        if (b < 0) 
                        {
                            shift = -shift;
                        }
                        shift = c / (b + shift);
                    }
                    double f = (sk + sp) * (sk - sp) + shift;
                    double g = sk * ek;
                    // Chase zeros.
                    for (int j{ k }; j < p - 1; j++) 
                    {
                        double t = std::hypot(f, g);
                        double cs = f / t;
                        double sn = g / t;
                        if (j != k) 
                        {
                            e[j - 1] = t;
                        }
                        f = cs * singular_values[j] + sn * e[j];
                        e[j] = cs * e[j] - sn * singular_values[j];
                        g = sn * singular_values[j + 1];
                        singular_values[j + 1] = cs * singular_values[j + 1];

                        for (int i{}; i < n; i++) 
                        {
                            t = cs * V[i][j] + sn * V[i][j + 1];
                            V[i][j + 1] = -sn * V[i][j] + cs * V[i][j + 1];
                            V[i][j] = t;
                        }
                        t = std::hypot(f, g);
                        cs = f / t;
                        sn = g / t;
                        singular_values[j] = t;
                        f = cs * e[j] + sn * singular_values[j + 1];
                        singular_values[j + 1] = -sn * e[j] + cs * singular_values[j + 1];
                        g = sn * e[j + 1];
                        e[j + 1] = cs * e[j + 1];
                        if (j < m - 1) 
                        {
                            for (int i{}; i < m; i++) 
                            {
                                t = cs * U[i][j] + sn * U[i][j + 1];
                                U[i][j + 1] = -sn * U[i][j] + cs * U[i][j + 1];
                                U[i][j] = t;
                            }
                        }
                    }
                    e[p - 2] = f;
                }
                break;
                // Convergence.
                default: 
                {
                    // Make the singular values positive.
                    if (singular_values[k] <= 0) 
                    {
                        singular_values[k] = singular_values[k] < 0 ? -singular_values[k] : 0;

                        for (int i{}; i <= pp; i++) 
                        {
                            V[i][k] = -V[i][k];
                        }
                    }
                    // Order the singular values.
                    while (k < pp) 
                    {
                        if (singular_values[k] >= singular_values[k + 1]) 
                        {
                            break;
                        }
                        double t = singular_values[k];
                        singular_values[k] = singular_values[k + 1];
                        singular_values[k + 1] = t;
                        if (k < n - 1) 
                        {
                            for (int i{}; i < n; i++) 
                            {
                                t = V[i][k + 1];
                                V[i][k + 1] = V[i][k];
                                V[i][k] = t;
                            }
                        }
                        if (k < m - 1) 
                        {
                            for (int i{}; i < m; i++) 
                            {
                                t = U[i][k + 1];
                                U[i][k + 1] = U[i][k];
                                U[i][k] = t;
                            }
                        }
                        k++;
                    }
                    p--;
                }
                break;
            }
        }

        // Set the small value tolerance used to calculate rank and pseudo-inverse
        tol = std::max(m * singular_values[0] * EPS, std::sqrt(Precision.SAFE_MIN));

        if (!transposed) 
        {
            cached_u = Matrix_Utils::create_real_matrix(U);
            cached_v = Matrix_Utils::create_real_matrix(V);
        }
else 
        {
            cached_u = Matrix_Utils::create_real_matrix(V);
            cached_v = Matrix_Utils::create_real_matrix(U);
        }
    }

    /**
     * Returns the matrix U of the decomposition.
     * <p>U is an orthogonal matrix, i.e. its transpose is also its inverse.</p>
     * @return the U matrix
     * @see #get_u_t()
     */
    public Real_Matrix get_u() 
    {
        // return the cached matrix
        return cached_u;

    }

    /**
     * Returns the transpose of the matrix U of the decomposition.
     * <p>U is an orthogonal matrix, i.e. its transpose is also its inverse.</p>
     * @return the U matrix (or NULL if decomposed matrix is singular)
     * @see #get_u()
     */
    public Real_Matrix get_u_t() 
    {
        if (cached_ut == NULL) 
        {
            cached_ut = get_u().transpose();
        }
        // return the cached matrix
        return cached_ut;
    }

    /**
     * Returns the diagonal matrix &Sigma; of the decomposition.
     * <p>&Sigma; is a diagonal matrix. The singular values are provided in
     * non-increasing order, for compatibility with Jama.</p>
     * @return the &Sigma; matrix
     */
    public Real_Matrix get_s() 
    {
        if (cached_s == NULL) 
        {
            // cache the matrix for subsequent calls
            cached_s = Matrix_Utils::create_real_diagonal_matrix(singular_values);
        }
        return cached_s;
    }

    /**
     * Returns the diagonal elements of the matrix &Sigma; of the decomposition.
     * <p>The singular values are provided in non-increasing order, for
     * compatibility with Jama.</p>
     * @return the diagonal elements of the &Sigma; matrix
     */
    public std::vector<double> get_singular_values() 
    {
        return singular_values.clone();
    }

    /**
     * Returns the matrix V of the decomposition.
     * <p>V is an orthogonal matrix, i.e. its transpose is also its inverse.</p>
     * @return the V matrix (or NULL if decomposed matrix is singular)
     * @see #get_v_t()
     */
    public Real_Matrix get_v() 
    {
        // return the cached matrix
        return cached_v;
    }

    /**
     * Returns the transpose of the matrix V of the decomposition.
     * <p>V is an orthogonal matrix, i.e. its transpose is also its inverse.</p>
     * @return the V matrix (or NULL if decomposed matrix is singular)
     * @see #get_v()
     */
    public Real_Matrix get_v_t() 
    {
        if (cached_vt == NULL) 
        {
            cached_vt = get_v().transpose();
        }
        // return the cached matrix
        return cached_vt;
    }

    /**
     * Returns the n &times; n covariance matrix.
     * <p>The covariance matrix is V &times; J &times; V<sup>T</sup>
     * where J is the diagonal matrix of the inverse of the squares of
     * the singular values.</p>
     * @param min_singular_value value below which singular values are ignored
     * (a 0 or negative value implies all singular value will be used)
     * @return covariance matrix
     * @exception Illegal_Argument_Exception if min_singular_value is larger than
     * the largest singular value, meaning all singular values are ignored
     */
    public Real_Matrix get_covariance(const double min_singular_value) 
    {
        // get the number of singular values to consider
        const int p = singular_values.size();
        int dimension = 0;
        while (dimension < p &&
               singular_values[dimension] >= min_singular_value) 
               {
            ++dimension;
        }

        if (dimension == 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::TOO_LARGE_CUTOFF_SINGULAR_VALUE, min_singular_value, singular_values[0], true);
        }

        const std::vector<std::vector<double>> data = std::vector<double>(dimension][p];
        get_v_t().walk_in_optimized_order(new DefaultReal_Matrix_Preserving_Visitor() 
        {
            /** {@inherit_doc} */
            //override
            public void visit(const int& row, const int& column, const double value) 
            {
                data[row][column] = value / singular_values[row];
            }
        }, 0, dimension - 1, 0, p - 1);

        Real_Matrix jv = Array_2D_Row_Real_Matrix(data, false);
        return jv.transpose_multiply(jv);
    }

    /**
     * Returns the L<sub>2</sub> norm of the matrix.
     * <p>The L<sub>2</sub> norm is max(|A &times; u|<sub>2</sub> /
     * |u|<sub>2</sub>), where |.|<sub>2</sub> denotes the vectorial 2-norm
     * (i.e. the traditional euclidian norm).</p>
     * @return norm
     */
    public double get_norm() 
    {
        return singular_values[0];
    }

    /**
     * Return the condition number of the matrix.
     * @return condition number of the matrix
     */
    public double get_condition_number() 
    {
        return singular_values[0] / singular_values[n - 1];
    }

    /**
     * Computes the inverse of the condition number.
     * In cases of rank deficiency, the {@link #get_condition_number() condition
     * number} will become undefined.
     *
     * @return the inverse of the condition number.
     */
    public double get_inverse_condition_number() 
    {
        return singular_values[n - 1] / singular_values[0];
    }

    /**
     * Return the effective numerical matrix rank.
     * <p>The effective numerical rank is the number of non-negligible
     * singular values. The threshold used to identify non-negligible
     * terms is max(m,n) &times; ulp(s<sub>1</sub>) where ulp(s<sub>1</sub>)
     * is the least significant bit of the largest singular value.</p>
     * @return effective numerical matrix rank
     */
    public int get_rank() 
    {
        int r = 0;
        for (int i{}; i < singular_values.size(); i++) 
        {
            if (singular_values[i] > tol) 
            {
                r++;
            }
        }
        return r;
    }

    /**
     * Get a solver for finding the A &times; X = B solution in least square sense.
     * @return a solver
     */
    public Decomposition_Solver get_solver() 
    {
        return Solver(singular_values, get_u_t(), get_v(), get_rank() == m, tol);
    }

    /** Specialized solver. */
    private static class Solver : Decomposition_Solver 
    {
        /** Pseudo-inverse of the initial matrix. */
        private const Real_Matrix pseudo_inverse;
        /** Singularity indicator. */
        private const bool non_singular;

        /**
         * Build a solver from decomposed matrix.
         *
         * @param singular_values Singular values.
         * @param uT U<sup>T</sup> matrix of the decomposition.
         * @param v V matrix of the decomposition.
         * @param non_singular Singularity indicator.
         * @param tol tolerance for singular values
         */
        private Solver(const std::vector<double> singular_values, const Real_Matrix uT, const Real_Matrix v, const bool non_singular, const double tol) 
        {
            const std::vector<std::vector<double>> su_t = uT.get_data();
            for (int i{}; i < singular_values.size(); ++i) 
            {
                const double& a;
                if (singular_values[i] > tol) 
                {
                    a = 1 / singular_values[i];
                }
else 
                {
                    a = 0;
                }
                const std::vector<double> su_ti = su_t[i];
                for (int j{}; j < su_ti.size(); ++j) 
                {
                    su_ti[j] *= a;
                }
            }
            pseudo_inverse = v.multiply(new Array_2D_Row_Real_Matrix(su_t, false));
            this.non_singular = non_singular;
        }

        /**
         * Solve the linear equation A &times; X = B in least square sense.
         * <p>
         * The m&times;n matrix A may not be square, the solution X is such that
         * ||A &times; X - B|| is minimal.
         * </p>
         * @param b Right-hand side of the equation A &times; X = B
         * @return a vector X that minimizes the two norm of A &times; X - B
         * @org.hipparchus.exception.
         * if the matrices dimensions do not match.
         */
        //override
        public Real_Vector solve(const Real_Vector b) 
        {
            return pseudo_inverse.operate(b);
        }

        /**
         * Solve the linear equation A &times; X = B in least square sense.
         * <p>
         * The m&times;n matrix A may not be square, the solution X is such that
         * ||A &times; X - B|| is minimal.
         * </p>
         *
         * @param b Right-hand side of the equation A &times; X = B
         * @return a matrix X that minimizes the two norm of A &times; X - B
         * @org.hipparchus.exception.
         * if the matrices dimensions do not match.
         */
        //override
        public Real_Matrix solve(const Real_Matrix b) 
        {
            return pseudo_inverse.multiply(b);
        }

        /**
         * Check if the decomposed matrix is non-singular.
         *
         * @return {@code true} if the decomposed matrix is non-singular.
         */
        //override
        public bool is_non_singular() 
        {
            return non_singular;
        }

        /**
         * Get the pseudo-inverse of the decomposed matrix.
         *
         * @return the inverse matrix.
         */
        //override
        public Real_Matrix get_inverse() 
        {
            return pseudo_inverse;
        }

        /** {@inherit_doc} */
        //override
        public int get_row_dimension() 
        {
            return pseudo_inverse.get_column_dimension();
        }

        /** {@inherit_doc} */
        //override
        public int get_column_dimension() 
        {
            return pseudo_inverse.get_row_dimension();
        }

    }
}


