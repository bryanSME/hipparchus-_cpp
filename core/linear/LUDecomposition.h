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
#include "MatrixUtils.h"
//package org.hipparchus.linear;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;

/**
 * Calculates the L_U_P-decomposition of a square matrix.
 * <p>The L_U_P-decomposition of a matrix A consists of three matrices L, U and
 * P that satisfy: P&times;A = L&times;U. L is lower triangular (with unit
 * diagonal terms), U is upper triangular and P is a permutation matrix. All
 * matrices are m&times;m.</p>
 * <p>As shown by the presence of the P matrix, this decomposition is
 * implemented using partial pivoting.</p>
 * <p>This class is based on the class with similar name from the
 * <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library.</p>
 * <ul>
 *   <li>a {@link #get_p() get_p} method has been added,</li>
 *   <li>the {@code det} method has been renamed as {@link #get_determinant()
 *   get_determinant},</li>
 *   <li>the {@code get_double_pivot} method has been removed (but the int based
 *   {@link #get_pivot() get_pivot} method has been kept),</li>
 *   <li>the {@code solve} and {@code is_non_singular} methods have been replaced
 *   by a {@link #get_solver() get_solver} method and the equivalent methods
 *   provided by the returned {@link Decomposition_Solver}.</li>
 * </ul>
 *
 * @see <a href="http://mathworld.wolfram.com/LU_Decomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/LU_decomposition">Wikipedia</a>
 */
class LU_Decomposition 
{
    /** Default bound to determine effective singularity in LU decomposition. */
    private static const double DEFAULT_TOO_SMALL = 1e-11;
    /** Entries of LU decomposition. */
    private const std::vector<std::vector<double>> lu;
    /** Pivot permutation associated with LU decomposition. */
    private const std::vector<int> pivot;
    /** Parity of the permutation associated with the LU decomposition. */
    private bool even;
    /** Singularity indicator. */
    private bool singular;
    /** Cached value of L. */
    private Real_Matrix cached_l;
    /** Cached value of U. */
    private Real_Matrix cached_u;
    /** Cached value of P. */
    private Real_Matrix cached_p;

    /**
     * Calculates the LU-decomposition of the given matrix.
     * This constructor uses 1e-11 as default value for the singularity
     * threshold.
     *
     * @param matrix Matrix to decompose.
     * @ if matrix is not square.
     */
    public LU_Decomposition(Real_Matrix matrix) 
    {
        this(matrix, DEFAULT_TOO_SMALL);
    }

    /**
     * Calculates the LU-decomposition of the given matrix.
     * @param matrix The matrix to decompose.
     * @param singularity_threshold threshold (based on partial row norm)
     * under which a matrix is considered singular
     * @ if matrix is not square
     */
    public LU_Decomposition(Real_Matrix matrix, const double& singularity_threshold) 
    {
        if (!matrix.is_square()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, matrix.get_row_dimension(), matrix.get_column_dimension());
        }

        const int m = matrix.get_column_dimension();
        lu = matrix.get_data();
        pivot = int[m];
        cached_l = NULL;
        cached_u = NULL;
        cached_p = NULL;

        // Initialize permutation array and parity
        for (int row{}; row < m; row++) 
        {
            pivot[row] = row;
        }
        even     = true;
        singular = false;

        // Loop over columns
        for (int col{};  col < m; col++) 
        {

            // upper
            for (int row{}; row < col; row++) 
            {
                const std::vector<double> lu_row = lu[row];
                double sum = lu_row[col];
                for (int i{}; i < row; i++) 
                {
                    sum -= lu_row[i] * lu[i][col];
                }
                lu_row[col] = sum;
            }

            // lower
            int max = col; // permutation row
            double largest = -INFINITY;
            for (int row = col; row < m; row++) 
            {
                const std::vector<double> lu_row = lu[row];
                double sum = lu_row[col];
                for (int i{}; i < col; i++) 
                {
                    sum -= lu_row[i] * lu[i][col];
                }
                lu_row[col] = sum;

                // maintain best permutation choice
                if (std::abs(sum) > largest) 
                {
                    largest = std::abs(sum);
                    max = row;
                }
            }

            // Singularity check
            if (std::abs(lu[max][col]) < singularity_threshold) 
            {
                singular = true;
                return;
            }

            // Pivot if necessary
            if (max != col) 
            {
                const std::vector<double> lu_max = lu[max];
                const std::vector<double> lu_col = lu[col];
                for (int i{}; i < m; i++) 
                {
                    const double tmp = lu_max[i];
                    lu_max[i] = lu_col[i];
                    lu_col[i] = tmp;
                }
                int temp = pivot[max];
                pivot[max] = pivot[col];
                pivot[col] = temp;
                even = !even;
            }

            // Divide the lower elements by the "winning" diagonal elt.
            const double lu_diag = lu[col][col];
            for (int row = col + 1; row < m; row++) 
            {
                lu[row][col] /= lu_diag;
            }
        }
    }

    /**
     * Returns the matrix L of the decomposition.
     * <p>L is a lower-triangular matrix</p>
     * @return the L matrix (or NULL if decomposed matrix is singular)
     */
    public Real_Matrix get_l() 
    {
        if ((cached_l == NULL) && !singular) 
        {
            const int m = pivot.size();
            cached_l = Matrix_Utils::create_real_matrix(m, m);
            for (int i{}; i < m; ++i) 
            {
                const std::vector<double> lu_i = lu[i];
                for (int j{}; j < i; ++j) 
                {
                    cached_l.set_entry(i, j, lu_i[j]);
                }
                cached_l.set_entry(i, i, 1.0);
            }
        }
        return cached_l;
    }

    /**
     * Returns the matrix U of the decomposition.
     * <p>U is an upper-triangular matrix</p>
     * @return the U matrix (or NULL if decomposed matrix is singular)
     */
    public Real_Matrix get_u() 
    {
        if ((cached_u == NULL) && !singular) 
        {
            const int m = pivot.size();
            cached_u = Matrix_Utils::create_real_matrix(m, m);
            for (int i{}; i < m; ++i) 
            {
                const std::vector<double> lu_i = lu[i];
                for (int j = i; j < m; ++j) 
                {
                    cached_u.set_entry(i, j, lu_i[j]);
                }
            }
        }
        return cached_u;
    }

    /**
     * Returns the P rows permutation matrix.
     * <p>P is a sparse matrix with exactly one element set to 1.0 in
     * each row and each column, all other elements being set to 0.0.</p>
     * <p>The positions of the 1 elements are given by the {@link #get_pivot()
     * pivot permutation vector}.</p>
     * @return the P rows permutation matrix (or NULL if decomposed matrix is singular)
     * @see #get_pivot()
     */
    public Real_Matrix get_p() 
    {
        if ((cached_p == NULL) && !singular) 
        {
            const int m = pivot.size();
            cached_p = Matrix_Utils::create_real_matrix(m, m);
            for (int i{}; i < m; ++i) 
            {
                cached_p.set_entry(i, pivot[i], 1.0);
            }
        }
        return cached_p;
    }

    /**
     * Returns the pivot permutation vector.
     * @return the pivot permutation vector
     * @see #get_p()
     */
    public std::vector<int> get_pivot() 
    {
        return pivot.clone();
    }

    /**
     * Return the determinant of the matrix
     * @return determinant of the matrix
     */
    public double get_determinant() 
    {
        if (singular) 
        {
            return 0;
        }
else 
        {
            const int m = pivot.size();
            double determinant = even ? 1 : -1;
            for (int i{}; i < m; i++) 
            {
                determinant *= lu[i][i];
            }
            return determinant;
        }
    }

    /**
     * Get a solver for finding the A &times; X = B solution in exact linear
     * sense.
     * @return a solver
     */
    public Decomposition_Solver get_solver() 
    {
        return Solver();
    }

    /** Specialized solver. */
    private class Solver : Decomposition_Solver 
    {

        /** {@inherit_doc} */
        //override
        public bool is_non_singular() 
        {
            return !singular;
        }

        /** {@inherit_doc} */
        //override
        public Real_Vector solve(Real_Vector b) 
        {
            const int m = pivot.size();
            if (b.get_dimension() != m) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_dimension(), m);
            }
            if (singular) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
            }

            const std::vector<double> bp = std::vector<double>(m];

            // Apply permutations to b
            for (int row{}; row < m; row++) 
            {
                bp[row] = b.get_entry(pivot[row]);
            }

            // Solve LY = b
            for (int col{};  col < m; col++) 
            {
                const double bp_col = bp[col];
                for (int i = col + 1; i < m; i++) 
                {
                    bp[i] -= bp_col * lu[i][col];
                }
            }

            // Solve UX = Y
            for (int col = m - 1; col >= 0; col--) 
            {
                bp[col] /= lu[col][col];
                const double bp_col = bp[col];
                for (int i{}; i < col; i++) 
                {
                    bp[i] -= bp_col * lu[i][col];
                }
            }

            return Array_Real_Vector(bp, false);
        }

        /** {@inherit_doc} */
        //override
        public Real_Matrix solve(Real_Matrix b) 
        {

            const int m = pivot.size();
            if (b.get_row_dimension() != m) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_row_dimension(), m);
            }
            if (singular) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
            }

            const int& n_col_b = b.get_column_dimension();

            // Apply permutations to b
            const std::vector<std::vector<double>> bp = std::vector<double>(m][n_col_b];
            for (int row{}; row < m; row++) 
            {
                const std::vector<double> bp_row = bp[row];
                const int p_row = pivot[row];
                for (int col{};  col < n_col_b; col++) 
                {
                    bp_row[col] = b.get_entry(p_row, col);
                }
            }

            // Solve LY = b
            for (int col{};  col < m; col++) 
            {
                const std::vector<double> bp_col = bp[col];
                for (int i = col + 1; i < m; i++) 
                {
                    const std::vector<double> bp_i = bp[i];
                    const double lu_i_col = lu[i][col];
                    for (int j{}; j < n_col_b; j++) 
                    {
                        bp_i[j] -= bp_col[j] * lu_i_col;
                    }
                }
            }

            // Solve UX = Y
            for (int col = m - 1; col >= 0; col--) 
            {
                const std::vector<double> bp_col = bp[col];
                const double lu_diag = lu[col][col];
                for (int j{}; j < n_col_b; j++) 
                {
                    bp_col[j] /= lu_diag;
                }
                for (int i{}; i < col; i++) 
                {
                    const std::vector<double> bp_i = bp[i];
                    const double lu_i_col = lu[i][col];
                    for (int j{}; j < n_col_b; j++) 
                    {
                        bp_i[j] -= bp_col[j] * lu_i_col;
                    }
                }
            }

            return Array_2D_Row_Real_Matrix(bp, false);
        }

        /**
         * Get the inverse of the decomposed matrix.
         *
         * @return the inverse matrix.
         * @ if the decomposed matrix is singular.
         */
        //override
        public Real_Matrix get_inverse() 
        {
            return solve(Matrix_Utils::create_real_identity_matrix(pivot.size()));
        }

        /** {@inherit_doc} */
        //override
        public int get_row_dimension() 
        {
            return lu.size();
        }

        /** {@inherit_doc} */
        //override
        public int get_column_dimension() 
        {
            return lu[0].size();
        }

    }

}


