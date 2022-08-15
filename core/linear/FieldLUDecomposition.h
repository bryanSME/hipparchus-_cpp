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

//import java.util.function.Predicate;
#include "MatrixUtils.h"
//import org.hipparchus.Field;
//import org.hipparchus.Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;

/**
 * Calculates the L_U_P-decomposition of a square matrix.
 * <p>The L_U_P-decomposition of a matrix A consists of three matrices
 * L, U and P that satisfy: PA = LU, L is lower triangular, and U is
 * upper triangular and P is a permutation matrix. All matrices are
 * m&times;m.</p>
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
 * @param <T> the type of the field elements
 * @see <a href="http://mathworld.wolfram.com/LU_Decomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/LU_decomposition">Wikipedia</a>
 */
template<typename T, typename std::enable_if<std::is_base_of<Field_Element<T>, T>::value>::type* = nullptr>
class FieldLU_Decomposition
{

    /** Field to which the elements belong. */
    private const Field<T> field;

    /** Entries of LU decomposition. */
    private std::vector<std::vector<T>> lu;

    /** Pivot permutation associated with LU decomposition. */
    private std::vector<int> pivot;

    /** Parity of the permutation associated with the LU decomposition. */
    private bool even;

    /** Singularity indicator. */
    private bool singular;

    /** Cached value of L. */
    private Field_Matrix<T> cached_l;

    /** Cached value of U. */
    private Field_Matrix<T> cached_u;

    /** Cached value of P. */
    private Field_Matrix<T> cached_p;

    /**
     * Calculates the LU-decomposition of the given matrix.
     * <p>
     * By default, <code>numeric_permutation_choice</code> is set to <code>true</code>.
     * </p>
     * @param matrix The matrix to decompose.
     * @ if matrix is not square
     * @see #FieldLU_Decomposition(Field_Matrix, Predicate)
     * @see #FieldLU_Decomposition(Field_Matrix, Predicate, bool)
     */
    public FieldLU_Decomposition(Field_Matrix<T> matrix) 
    {
        this(matrix, e -> e.is_zero());
    }

    /**
     * Calculates the LU-decomposition of the given matrix.
     * <p>
     * By default, <code>numeric_permutation_choice</code> is set to <code>true</code>.
     * </p>
     * @param matrix The matrix to decompose.
     * @param zero_checker checker for zero elements
     * @ if matrix is not square
     * @see #FieldLU_Decomposition(Field_Matrix, Predicate, bool)
     */
    public FieldLU_Decomposition(Field_Matrix<T> matrix, const Predicate<T> zero_checker ) 
    {
        this(matrix, zero_checker, true);
    }

    /**
     * Calculates the LU-decomposition of the given matrix.
     * @param matrix The matrix to decompose.
     * @param zero_checker checker for zero elements
     * @param numeric_permutation_choice if <code>true</code> choose permutation index with numeric calculations, otherwise choose with <code>zero_checker</code>
     * @ if matrix is not square
     */
    public FieldLU_Decomposition(Field_Matrix<T> matrix, const Predicate<T> zero_checker, bool numeric_permutation_choice) 
    {
        if (!matrix.is_square()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, matrix.get_row_dimension(), matrix.get_column_dimension());
        }

        const int m = matrix.get_column_dimension();
        field = matrix.get_field();
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
                const std::vector<T> lu_row = lu[row];
                T sum = lu_row[col];
                for (int i{}; i < row; i++) 
                {
                    sum = sum.subtract(lu_row[i].multiply(lu[i][col]));
                }
                lu_row[col] = sum;
            }

            int max = col; // permutation row
            if (numeric_permutation_choice) 
            {

                // lower
                double largest = -INFINITY;

                for (int row = col; row < m; row++) 
                {
                    const std::vector<T> lu_row = lu[row];
                    T sum = lu_row[col];
                    for (int i{}; i < col; i++) 
                    {
                        sum = sum.subtract(lu_row[i].multiply(lu[i][col]));
                    }
                    lu_row[col] = sum;

                    // maintain best permutation choice
                    double abs_sum = std::abs(sum.get_real());
                    if (abs_sum > largest) 
                    {
                        largest = abs_sum;
                        max = row;
                    }
                }

            }
else 
            {

                // lower
                int non_zero = col; // permutation row
                for (int row = col; row < m; row++) 
                {
                    const std::vector<T> lu_row = lu[row];
                    T sum = lu_row[col];
                    for (int i{}; i < col; i++) 
                    {
                        sum = sum.subtract(lu_row[i].multiply(lu[i][col]));
                    }
                    lu_row[col] = sum;

                    if (zero_checker.test(lu[non_zero][col])) 
                    {
                        // try to select a better permutation choice
                        ++non_zero;
                    }
                }
                max = std::min(m - 1, non_zero);

            }

            // Singularity check
            if (zero_checker.test(lu[max][col])) 
            {
                singular = true;
                return;
            }

            // Pivot if necessary
            if (max != col) 
            {
                const std::vector<T> lu_max = lu[max];
                const std::vector<T> lu_col = lu[col];
                for (int i{}; i < m; i++) 
                {
                    const T tmp = lu_max[i];
                    lu_max[i] = lu_col[i];
                    lu_col[i] = tmp;
                }
                int temp = pivot[max];
                pivot[max] = pivot[col];
                pivot[col] = temp;
                even = !even;
            }

            // Divide the lower elements by the "winning" diagonal elt.
            const T lu_diag = lu[col][col];
            for (int row = col + 1; row < m; row++) 
            {
                lu[row][col] = lu[row][col].divide(lu_diag);
            }
        }

    }

    /**
     * Returns the matrix L of the decomposition.
     * <p>L is a lower-triangular matrix</p>
     * @return the L matrix (or NULL if decomposed matrix is singular)
     */
    public Field_Matrix<T> get_l() 
    {
        if ((cached_l == NULL) && !singular) 
        {
            const int m = pivot.size();
            cached_l = Array2DRowField_Matrix<>(field, m, m);
            for (int i{}; i < m; ++i) 
            {
                const std::vector<T> lu_i = lu[i];
                for (int j{}; j < i; ++j) 
                {
                    cached_l.set_entry(i, j, lu_i[j]);
                }
                cached_l.set_entry(i, i, field.get_one());
            }
        }
        return cached_l;
    }

    /**
     * Returns the matrix U of the decomposition.
     * <p>U is an upper-triangular matrix</p>
     * @return the U matrix (or NULL if decomposed matrix is singular)
     */
    public Field_Matrix<T> get_u() 
    {
        if ((cached_u == NULL) && !singular) 
        {
            const int m = pivot.size();
            cached_u = Array2DRowField_Matrix<>(field, m, m);
            for (int i{}; i < m; ++i) 
            {
                const std::vector<T> lu_i = lu[i];
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
    public Field_Matrix<T> get_p() 
    {
        if ((cached_p == NULL) && !singular) 
        {
            const int m = pivot.size();
            cached_p = Array2DRowField_Matrix<>(field, m, m);
            for (int i{}; i < m; ++i) 
            {
                cached_p.set_entry(i, pivot[i], field.get_one());
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
     * Return the determinant of the matrix.
     * @return determinant of the matrix
     */
    public T get_determinant() 
    {
        if (singular) 
        {
            return field.get_zero();
        }
else 
        {
            const int m = pivot.size();
            T determinant = even ? field.get_one() : field.get_zero().subtract(field.get_one());
            for (int i{}; i < m; i++) 
            {
                determinant = determinant.multiply(lu[i][i]);
            }
            return determinant;
        }
    }

    /**
     * Get a solver for finding the A &times; X = B solution in exact linear sense.
     * @return a solver
     */
    public FieldDecomposition_Solver<T> get_solver() 
    {
        return Solver();
    }

    /** Specialized solver.
     */
    private class Solver : FieldDecomposition_Solver<T> 
    {

        /** {@inherit_doc} */
        //override
        public bool is_non_singular() 
        {
            return !singular;
        }

        /** {@inherit_doc} */
        //override
        public Field_Vector<T> solve(Field_Vector<T> b) 
        {
            if (b instanceof ArrayField_Vector) 
            {
                return solve((ArrayField_Vector<T>) b);
            }
else 
            {

                const int m = pivot.size();
                if (b.get_dimension() != m) 
                {
                    throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_dimension(), m);
                }
                if (singular) 
                {
                    throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
                }

                // Apply permutations to b
                const std::vector<T> bp = Math_Arrays::build_array(field, m);
                for (int row{}; row < m; row++) 
                {
                    bp[row] = b.get_entry(pivot[row]);
                }

                // Solve LY = b
                for (int col{};  col < m; col++) 
                {
                    const T& bp_col = bp[col];
                    for (int i = col + 1; i < m; i++) 
                    {
                        bp[i] = bp[i].subtract(bp_col.multiply(lu[i][col]));
                    }
                }

                // Solve UX = Y
                for (int col = m - 1; col >= 0; col--) 
                {
                    bp[col] = bp[col].divide(lu[col][col]);
                    const T& bp_col = bp[col];
                    for (int i{}; i < col; i++) 
                    {
                        bp[i] = bp[i].subtract(bp_col.multiply(lu[i][col]));
                    }
                }

                return ArrayField_Vector<T>(field, bp, false);

            }
        }

        /** Solve the linear equation A &times; X = B.
         * <p>The A matrix is implicit here. It is </p>
         * @param b right-hand side of the equation A &times; X = B
         * @return a vector X such that A &times; X = B
         * @ if the matrices dimensions do not match.
         * @ if the decomposed matrix is singular.
         */
        public ArrayField_Vector<T> solve(ArrayField_Vector<T> b) 
        {
            const int m = pivot.size();
            const int length = b.get_dimension();
            if (length != m) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, length, m);
            }
            if (singular) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
            }

            // Apply permutations to b
            const std::vector<T> bp = Math_Arrays::build_array(field, m);
            for (int row{}; row < m; row++) 
            {
                bp[row] = b.get_entry(pivot[row]);
            }

            // Solve LY = b
            for (int col{};  col < m; col++) 
            {
                const T& bp_col = bp[col];
                for (int i = col + 1; i < m; i++) 
                {
                    bp[i] = bp[i].subtract(bp_col.multiply(lu[i][col]));
                }
            }

            // Solve UX = Y
            for (int col = m - 1; col >= 0; col--) 
            {
                bp[col] = bp[col].divide(lu[col][col]);
                const T& bp_col = bp[col];
                for (int i{}; i < col; i++) 
                {
                    bp[i] = bp[i].subtract(bp_col.multiply(lu[i][col]));
                }
            }

            return ArrayField_Vector<T>(bp, false);
        }

        /** {@inherit_doc} */
        //override
        public Field_Matrix<T> solve(Field_Matrix<T> b) 
        {
            const int m = pivot.size();
            if (b.get_row_dimension() != m) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_row_dimension(), m);
            }
            if (singular) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
            }

            const int& n_col_b = b.get_column_dimension();

            // Apply permutations to b
            const std::vector<std::vector<T>> bp = Math_Arrays::build_array(field, m, n_col_b);
            for (int row{}; row < m; row++) 
            {
                const std::vector<T> bp_row = bp[row];
                const int p_row = pivot[row];
                for (int col{};  col < n_col_b; col++) 
                {
                    bp_row[col] = b.get_entry(p_row, col);
                }
            }

            // Solve LY = b
            for (int col{};  col < m; col++) 
            {
                const std::vector<T> bp_col = bp[col];
                for (int i = col + 1; i < m; i++) 
                {
                    const std::vector<T> bp_i = bp[i];
                    const T lu_i_col = lu[i][col];
                    for (int j{}; j < n_col_b; j++) 
                    {
                        bp_i[j] = bp_i[j].subtract(bp_col[j].multiply(lu_i_col));
                    }
                }
            }

            // Solve UX = Y
            for (int col = m - 1; col >= 0; col--) 
            {
                const std::vector<T> bp_col = bp[col];
                const T lu_diag = lu[col][col];
                for (int j{}; j < n_col_b; j++) 
                {
                    bp_col[j] = bp_col[j].divide(lu_diag);
                }
                for (int i{}; i < col; i++) 
                {
                    const std::vector<T> bp_i = bp[i];
                    const T lu_i_col = lu[i][col];
                    for (int j{}; j < n_col_b; j++) 
                    {
                        bp_i[j] = bp_i[j].subtract(bp_col[j].multiply(lu_i_col));
                    }
                }
            }

            return Array2DRowField_Matrix<T>(field, bp, false);

        }

        /** {@inherit_doc} */
        //override
        public Field_Matrix<T> get_inverse() 
        {
            return solve(Matrix_Utils::create_field_identity_matrix(field, pivot.size()));
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


