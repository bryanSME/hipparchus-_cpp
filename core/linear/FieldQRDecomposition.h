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
#include <type_traits>
#include "../CalculusFieldElement.hpp"
//import java.util.Arrays;
//import java.util.function.Predicate;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;


/**
 * Calculates the QR-decomposition of a field matrix.
 * <p>The QR-decomposition of a matrix A consists of two matrices Q and R
 * that satisfy: A = QR, Q is orthogonal (Q<sup>T</sup>Q = I), and R is
 * upper triangular. If A is m&times;n, Q is m&times;m and R m&times;n.</p>
 * <p>This class compute the decomposition using Householder reflectors.</p>
 * <p>For efficiency purposes, the decomposition in packed form is transposed.
 * This allows inner loop to iterate inside rows, which is much more cache-efficient
 * in Java.</p>
 * <p>This class is based on the class {@link QR_Decomposition}.</p>
 *
 * @param <T> type of the underlying field elements
 * @see <a href="http://mathworld.wolfram.com/QR_Decomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/QR_decomposition">Wikipedia</a>
 *
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldQR_Decomposition 
{
    /**
     * A packed TRANSPOSED representation of the QR decomposition.
     * <p>The elements BELOW the diagonal are the elements of the UPPER triangular
     * matrix R, and the rows ABOVE the diagonal are the Householder reflector vectors
     * from which an explicit form of Q can be recomputed if desired.</p>
     */
    private std::vector<std::vector<T>> qrt;
    /** The diagonal elements of R. */
    private std::vector<T> r_diag;
    /** Cached value of Q. */
    private Field_Matrix<T> cached_q;
    /** Cached value of QT. */
    private Field_Matrix<T> cached_q_t;
    /** Cached value of R. */
    private Field_Matrix<T> cached_r;
    /** Cached value of H. */
    private Field_Matrix<T> cached_h;
    /** Singularity threshold. */
    private const T threshold;
    /** checker for zero. */
    private const Predicate<T> zero_checker;

    /**
     * Calculates the QR-decomposition of the given matrix.
     * The singularity threshold defaults to zero.
     *
     * @param matrix The matrix to decompose.
     *
     * @see #FieldQR_Decomposition(Field_Matrix, Calculus_Field_Element)
     */
    public FieldQR_Decomposition(Field_Matrix<T> matrix) 
    {
        this(matrix, matrix.get_field().get_zero());
    }

    /**
     * Calculates the QR-decomposition of the given matrix.
     *
     * @param matrix The matrix to decompose.
     * @param threshold Singularity threshold.
     */
    public FieldQR_Decomposition(Field_Matrix<T> matrix, T threshold) 
    {
        this(matrix, threshold, e -> e.is_zero());
    }

    /**
     * Calculates the QR-decomposition of the given matrix.
     *
     * @param matrix The matrix to decompose.
     * @param threshold Singularity threshold.
     * @param zero_checker checker for zero
     */
    public FieldQR_Decomposition(Field_Matrix<T> matrix, T threshold, Predicate<T> zero_checker) 
    {
        this.threshold   = threshold;
        this.zero_checker = zero_checker;

        const int m = matrix.get_row_dimension();
        const int n = matrix.get_column_dimension();
        qrt = matrix.transpose().get_data();
        r_diag = Math_Arrays::build_array(threshold.get_field(),std::min(m, n));
        cached_q  = NULL;
        cached_q_t = NULL;
        cached_r  = NULL;
        cached_h  = NULL;

        decompose(qrt);

    }

    /** Decompose matrix.
     * @param matrix transposed matrix
     */
    protected void decompose(std::vector<std::vector<T>> matrix) 
    {
        for (const int& minor = 0; minor < std::min(matrix.size(), matrix[0].size()); minor++) 
        {
            perform_householder_reflection(minor, matrix);
        }
    }

    /** Perform Householder reflection for a minor A(minor, minor) of A.
     * @param minor minor index
     * @param matrix transposed matrix
     */
    protected void perform_householder_reflection(const int& minor, std::vector<std::vector<T>> matrix) 
    {

        const std::vector<T> qrt_minor = matrix[minor];
        const T zero = threshold.get_field().get_zero();
        /*
         * Let x be the first column of the minor, and a^2 = |x|^2.
         * x will be in the positions qr[minor][minor] through qr[m][minor].
         * The first column of the transformed minor will be (a,0,0,..)'
         * The sign of a is chosen to be opposite to the sign of the first
         * component of x. Let's find a:
         */
        T x_norm_sqr = zero;
        for (int row = minor; row < qrt_minor.size(); row++) 
        {
            const T c = qrt_minor[row];
            x_norm_sqr = x_norm_sqr.add(c.multiply(c));
        }
        const T a = (qrt_minor[minor].get_real() > 0) ? x_norm_sqr.sqrt().negate() : x_norm_sqr.sqrt();
        r_diag[minor] = a;

        if (!zero_checker.test(a)) 
        {

            /*
             * Calculate the normalized reflection vector v and transform
             * the first column. We know the norm of v beforehand: v = x-ae
             * so |v|^2 = <x-ae,x-ae> = <x,x>-2a<x,e>+a^2<e,e> =
             * a^2+a^2-2a<x,e> = 2a*(a - <x,e>).
             * Here <x, e> is now qr[minor][minor].
             * v = x-ae is stored in the column at qr:
             */
            qrt_minor[minor] = qrt_minor[minor].subtract(a); // now |v|^2 = -2a*(qr[minor][minor])

            /*
             * Transform the rest of the columns of the minor:
             * They will be transformed by the matrix H = I-2vv'/|v|^2.
             * If x is a column vector of the minor, then
             * Hx = (I-2vv'/|v|^2)x = x-2vv'x/|v|^2 = x - 2<x,v>/|v|^2 v.
             * Therefore the transformation is easily calculated by
             * subtracting the column vector (2<x,v>/|v|^2)v from x.
             *
             * Let 2<x,v>/|v|^2 = alpha. From above we have
             * |v|^2 = -2a*(qr[minor][minor]), so
             * alpha = -<x,v>/(a*qr[minor][minor])
             */
            for (int col = minor+1; col < matrix.size(); col++) 
            {
                const std::vector<T> qrt_col = matrix[col];
                T alpha = zero;
                for (int row = minor; row < qrt_col.size(); row++) 
                {
                    alpha = alpha.subtract(qrt_col[row].multiply(qrt_minor[row]));
                }
                alpha = alpha.divide(a.multiply(qrt_minor[minor]));

                // Subtract the column vector alpha*v from x.
                for (int row = minor; row < qrt_col.size(); row++) 
                {
                    qrt_col[row] = qrt_col[row].subtract(alpha.multiply(qrt_minor[row]));
                }
            }
        }
    }


    /**
     * Returns the matrix R of the decomposition.
     * <p>R is an upper-triangular matrix</p>
     * @return the R matrix
     */
    public Field_Matrix<T> get_r() 
    {

        if (cached_r == NULL) 
        {

            // R is supposed to be m x n
            const int n = qrt.size();
            const int m = qrt[0].size();
            std::vector<std::vector<T>> ra = Math_Arrays::build_array(threshold.get_field(), m, n);
            // copy the diagonal from r_diag and the upper triangle of qr
            for (int row = std::min(m, n) - 1; row >= 0; row--) 
            {
                ra[row][row] = r_diag[row];
                for (int col = row + 1; col < n; col++) 
                {
                    ra[row][col] = qrt[col][row];
                }
            }
            cached_r = Matrix_Utils::create_field_matrix(ra);
        }

        // return the cached matrix
        return cached_r;
    }

    /**
     * Returns the matrix Q of the decomposition.
     * <p>Q is an orthogonal matrix</p>
     * @return the Q matrix
     */
    public Field_Matrix<T> get_q() 
    {
        if (cached_q == NULL) 
        {
            cached_q = get_q_t().transpose();
        }
        return cached_q;
    }

    /**
     * Returns the transpose of the matrix Q of the decomposition.
     * <p>Q is an orthogonal matrix</p>
     * @return the transpose of the Q matrix, Q<sup>T</sup>
     */
    public Field_Matrix<T> get_q_t() 
    {
        if (cached_q_t == NULL) 
        {

            // QT is supposed to be m x m
            const int n = qrt.size();
            const int m = qrt[0].size();
            std::vector<std::vector<T>> qta = Math_Arrays::build_array(threshold.get_field(), m, m);

            /*
             * Q = Q1 Q2 ... Q_m, so Q is formed by first constructing Q_m and then
             * applying the Householder transformations Q_(m-1),Q_(m-2),...,Q1 in
             * succession to the result
             */
            for (const int& minor = m - 1; minor >= std::min(m, n); minor--) 
            {
                qta[minor][minor] = threshold.get_field().get_one();
            }

            for (const int& minor = std::min(m, n)-1; minor >= 0; minor--)
            {
                const std::vector<T> qrt_minor = qrt[minor];
                qta[minor][minor] = threshold.get_field().get_one();
                if (!qrt_minor[minor].is_zero()) 
                {
                    for (int col = minor; col < m; col++) 
                    {
                        T alpha = threshold.get_field().get_zero();
                        for (int row = minor; row < m; row++) 
                        {
                            alpha = alpha.subtract(qta[col][row].multiply(qrt_minor[row]));
                        }
                        alpha = alpha.divide(r_diag[minor].multiply(qrt_minor[minor]));

                        for (int row = minor; row < m; row++) 
                        {
                            qta[col][row] = qta[col][row].add(alpha.negate().multiply(qrt_minor[row]));
                        }
                    }
                }
            }
            cached_q_t = Matrix_Utils::create_field_matrix(qta);
        }

        // return the cached matrix
        return cached_q_t;
    }

    /**
     * Returns the Householder reflector vectors.
     * <p>H is a lower trapezoidal matrix whose columns represent
     * each successive Householder reflector vector. This matrix is used
     * to compute Q.</p>
     * @return a matrix containing the Householder reflector vectors
     */
    public Field_Matrix<T> get_h() 
    {
        if (cached_h == NULL) 
        {

            const int n = qrt.size();
            const int m = qrt[0].size();
            std::vector<std::vector<T>> ha = Math_Arrays::build_array(threshold.get_field(), m, n);
            for (int i{}; i < m; ++i) 
            {
                for (int j{}; j < std::min(i + 1, n); ++j) 
                {
                    ha[i][j] = qrt[j][i].divide(r_diag[j].negate());
                }
            }
            cached_h = Matrix_Utils::create_field_matrix(ha);
        }

        // return the cached matrix
        return cached_h;
    }

    /**
     * Get a solver for finding the A &times; X = B solution in least square sense.
     * <p>
     * Least Square sense means a solver can be computed for an overdetermined system, * (i.e. a system with more equations than unknowns, which corresponds to a tall A
     * matrix with more rows than columns). In any case, if the matrix is singular
     * within the tolerance set at {@link #FieldQR_Decomposition(Field_Matrix, * Calculus_Field_Element) construction}, an error will be triggered when
     * the {@link Decomposition_Solver#solve(Real_Vector) solve} method will be called.
     * </p>
     * @return a solver
     */
    public FieldDecomposition_Solver<T> get_solver() 
    {
        return Field_Solver();
    }

    /**
     * Specialized solver.
     */
    private class Field_Solver : FieldDecomposition_Solver<T>
    {

        /** {@inherit_doc} */
        //override
        public bool is_non_singular() 
        {
            return !check_singular(r_diag, threshold, false);
        }

        /** {@inherit_doc} */
        //override
        public Field_Vector<T> solve(Field_Vector<T> b) 
        {
            const int n = qrt.size();
            const int m = qrt[0].size();
            if (b.get_dimension() != m) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_dimension(), m);
            }
            check_singular(r_diag, threshold, true);

            const std::vector<T> x =Math_Arrays::build_array(threshold.get_field(),n);
            const std::vector<T> y = b.to_array();

            // apply Householder transforms to solve Q.y = b
            for (const int& minor = 0; minor < std::min(m, n); minor++) 
            {

                const std::vector<T> qrt_minor = qrt[minor];
                T dot_product = threshold.get_field().get_zero();
                for (int row = minor; row < m; row++) 
                {
                    dot_product = dot_product.add(y[row].multiply(qrt_minor[row]));
                }
                dot_product =  dot_product.divide(r_diag[minor].multiply(qrt_minor[minor]));

                for (int row = minor; row < m; row++) 
                {
                    y[row] = y[row].add(dot_product.multiply(qrt_minor[row]));
                }
            }

            // solve triangular system R.x = y
            for (int row = r_diag.size() - 1; row >= 0; --row) 
            {
                y[row] = y[row].divide(r_diag[row]);
                const T y_row = y[row];
                const std::vector<T> qrt_row = qrt[row];
                x[row] = y_row;
                for (int i{}; i < row; i++) 
                {
                    y[i] = y[i].subtract(y_row.multiply(qrt_row[i]));
                }
            }

            return ArrayField_Vector<T>(x, false);
        }

        /** {@inherit_doc} */
        //override
        public Field_Matrix<T> solve(Field_Matrix<T> b) 
        {
            const int n = qrt.size();
            const int m = qrt[0].size();
            if (b.get_row_dimension() != m) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_row_dimension(), m);
            }
            check_singular(r_diag, threshold, true);

            const int columns        = b.get_column_dimension();
            const int block_size      = BlockField_Matrix.BLOCK_SIZE;
            const int c_blocks        = (columns + block_size - 1) / block_size;
            const std::vector<std::vector<T>> x_blocks = BlockField_Matrix.create_blocks_layout(threshold.get_field(),n, columns);
            const std::vector<std::vector<T>> y       = Math_Arrays::build_array(threshold.get_field(), b.get_row_dimension(), block_size);
            const std::vector<T>   alpha   = Math_Arrays::build_array(threshold.get_field(), block_size);

            for (const int& k_block = 0; k_block < c_blocks; ++k_block) 
            {
                const int& k_start = k_block * block_size;
                const int& k_end   = std::min(k_start + block_size, columns);
                const int& k_width = k_end - k_start;

                // get the right hand side vector
                b.copy_sub_matrix(0, m - 1, k_start, k_end - 1, y);

                // apply Householder transforms to solve Q.y = b
                for (const int& minor = 0; minor < std::min(m, n); minor++) 
                {
                    const std::vector<T> qrt_minor = qrt[minor];
                    const T factor     = r_diag[minor].multiply(qrt_minor[minor]).reciprocal();

                    Arrays.fill(alpha, 0, k_width, threshold.get_field().get_zero());
                    for (int row = minor; row < m; ++row) 
                    {
                        const T   d    = qrt_minor[row];
                        const std::vector<T> y_row = y[row];
                        for (int k{}; k < k_width; ++k) 
                        {
                            alpha[k] = alpha[k].add(d.multiply(y_row[k]));
                        }
                    }

                    for (int k{}; k < k_width; ++k) 
                    {
                        alpha[k] = alpha[k].multiply(factor);
                    }

                    for (int row = minor; row < m; ++row) 
                    {
                        const T   d    = qrt_minor[row];
                        const std::vector<T> y_row = y[row];
                        for (int k{}; k < k_width; ++k) 
                        {
                            y_row[k] = y_row[k].add(alpha[k].multiply(d));
                        }
                    }
                }

                // solve triangular system R.x = y
                for (int j = r_diag.size() - 1; j >= 0; --j) 
                {
                    const int      j_block = j / block_size;
                    const int      j_start = j_block * block_size;
                    const T   factor = r_diag[j].reciprocal();
                    const std::vector<T> yJ     = y[j];
                    const std::vector<T> x_block = x_blocks[j_block * c_blocks + k_block];
                    int index = (j - j_start) * k_width;
                    for (int k{}; k < k_width; ++k) 
                    {
                        yJ[k]           =yJ[k].multiply(factor);
                        x_block[index++] = yJ[k];
                    }

                    const std::vector<T> qrtJ = qrt[j];
                    for (int i{}; i < j; ++i) 
                    {
                        const T rIJ  = qrtJ[i];
                        const std::vector<T> y_i = y[i];
                        for (int k{}; k < k_width; ++k) 
                        {
                            y_i[k] = y_i[k].subtract(yJ[k].multiply(rIJ));
                        }
                    }
                }
            }

            return BlockField_Matrix<T>(n, columns, x_blocks, false);
        }

        /**
         * {@inherit_doc}
         * @ if the decomposed matrix is singular.
         */
        //override
        public Field_Matrix<T> get_inverse() 
        {
            return solve(Matrix_Utils::create_field_identity_matrix(threshold.get_field(), qrt[0].size()));
        }

        /**
         * Check singularity.
         *
         * @param diag Diagonal elements of the R matrix.
         * @param min Singularity threshold.
         * @param raise Whether to raise a {@link }
         * if any element of the diagonal fails the check.
         * @return {@code true} if any element of the diagonal is smaller
         * or equal to {@code min}.
         * @ if the matrix is singular and
         * {@code raise} is {@code true}.
         */
        private bool check_singular(std::vector<T> diag, T min, bool raise) 
        {
            const int len = diag.size();
            for (int i{}; i < len; i++) 
            {
                const T d = diag[i];
                if (std::abs(d.get_real()) <= min.get_real()) 
                {
                    if (raise) 
                    {
                        throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
                    }
else 
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        /** {@inherit_doc} */
        //override
        public int get_row_dimension() 
        {
            return qrt[0].size();
        }

        /** {@inherit_doc} */
        //override
        public int get_column_dimension() 
        {
            return qrt.size();
        }

    }
}


