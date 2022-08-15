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

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
#include "MatrixUtils.h"

/**
 * Calculates the Cholesky decomposition of a matrix.
 * <p>The Cholesky decomposition of a real symmetric positive-definite
 * matrix A consists of a lower triangular matrix L with same size such
 * that: A = LL<sup>T</sup>. In a sense, this is the square root of A.</p>
 * <p>This class is based on the class with similar name from the
 * <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library, with the
 * following changes:</p>
 * <ul>
 *   <li>a {@link #get_l_t() get_l_t} method has been added,</li>
 *   <li>the {@code isspd} method has been removed, since the constructor of
 *   this class a {@link } when a
 *   matrix cannot be decomposed,</li>
 *   <li>a {@link #get_determinant() get_determinant} method has been added,</li>
 *   <li>the {@code solve} method has been replaced by a {@link #get_solver()
 *   get_solver} method and the equivalent method provided by the returned
 *   {@link Decomposition_Solver}.</li>
 * </ul>
 *
 * @see <a href="http://mathworld.wolfram.com/Cholesky_Decomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Cholesky_decomposition">Wikipedia</a>
 */
class Cholesky_Decomposition 
{
public:
    /**
     * Default threshold above which off-diagonal elements are considered too different
     * and matrix not symmetric.
     */
    public static constexpr double DEFAULT_RELATIVE_SYMMETRY_THRESHOLD{ 1.0e-15 };
    /**
     * Default threshold below which diagonal elements are considered NULL
     * and matrix not positive definite.
     */
    public static constexpr double DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD{ 1.0e-10 };
private:
    /** Row-oriented storage for L<sup>T</sup> matrix data. */
    const std::vector<std::vector<double>> mu_l_t_data;
    /** Cached value of L. */
    Real_Matrix my_cached_l;
    /** Cached value of LT. */
    Real_Matrix my_cached_l_t;

public:
    /**
     * Calculates the Cholesky decomposition of the given matrix.
     * <p>
     * Calling this constructor is equivalent to call {@link
     * #Cholesky_Decomposition(Real_Matrix, double, double)} with the
     * thresholds set to the default values {@link
     * #DEFAULT_RELATIVE_SYMMETRY_THRESHOLD} and {@link
     * #DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD}
     * </p>
     * @param matrix the matrix to decompose
     * @ if the matrix is not square.
     * @ if the matrix is not symmetric.
     * @ if the matrix is not
     * strictly positive definite.
     * @see #Cholesky_Decomposition(Real_Matrix, double, double)
     * @see #DEFAULT_RELATIVE_SYMMETRY_THRESHOLD
     * @see #DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD
     */
    Cholesky_Decomposition(const Real_Matrix& matrix) 
    {
        this(matrix, DEFAULT_RELATIVE_SYMMETRY_THRESHOLD, DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD);
    }

    /**
     * Calculates the Cholesky decomposition of the given matrix.
     * @param matrix the matrix to decompose
     * @param relative_symmetry_threshold threshold above which off-diagonal
     * elements are considered too different and matrix not symmetric
     * @param absolute_positivity_threshold threshold below which diagonal
     * elements are considered NULL and matrix not positive definite
     * @ if the matrix is not square.
     * @ if the matrix is not symmetric.
     * @ if the matrix is not
     * strictly positive definite.
     * @see #Cholesky_Decomposition(Real_Matrix)
     * @see #DEFAULT_RELATIVE_SYMMETRY_THRESHOLD
     * @see #DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD
     */
    Cholesky_Decomposition(const Real_Matrix& matrix, const double& relative_symmetry_threshold, const double& absolute_positivity_threshold) 
    {
        if (!matrix.is_square()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, matrix.get_row_dimension(), matrix.get_column_dimension());
        }

        const int order = matrix.get_row_dimension();
        my_l_t_data   = matrix.get_data();
        my_cached_l  = NULL;
        my_cached_l_t = NULL;

        // check the matrix before transformation
        for (int i{}; i < order; ++i) 
        {
            const auto lI = my_l_t_data[i];

            // check off-diagonal elements (and reset them to 0)
            for (int j = i + 1; j < order; ++j) 
            {
                const auto lJ = my_l_t_data[j];
                const double lij = lI[j];
                const double lji = lJ[i];
                const double max_delta =
                    relative_symmetry_threshold * std::max(std::abs(lij), std::abs(lji));
                if (std::abs(lij - lji) > max_delta) 
                {
                    throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SYMMETRIC_MATRIX, i, j, relative_symmetry_threshold);
                }
                lJ[i] = 0;
           }
        }

        // transform the matrix
        for (int i{}; i < order; ++i) 
        {
            auto lt_i = my_l_t_data[i];

            // check diagonal element
            if (lt_i[i] <= absolute_positivity_threshold) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_DEFINITE_MATRIX);
            }

            lt_i[i] = std::sqrt(lt_i[i]);
            const double inverse = 1.0 / lt_i[i];

            for (int q = order - 1; q > i; --q) 
            {
                lt_i[q] *= inverse;
                const auto lt_q = my_l_t_data[q];
                for (const int& p = q; p < order; ++p) 
                {
                    lt_q[p] -= lt_i[q] * lt_i[p];
                }
            }
        }
    }

    /**
     * Returns the matrix L of the decomposition.
     * <p>L is an lower-triangular matrix</p>
     * @return the L matrix
     */
    Real_Matrix get_l() 
    {
        if (cached_l == NULL) 
        {
            cached_l = get_l_t().transpose();
        }
        return cached_l;
    }

    /**
     * Returns the transpose of the matrix L of the decomposition.
     * <p>L<sup>T</sup> is an upper-triangular matrix</p>
     * @return the transpose of the matrix L of the decomposition
     */
    Real_Matrix get_l_t() 
    {

        if (cached_l_t == NULL) 
        {
            cached_l_t = Matrix_Utils::create_real_matrix(l_t_data);
        }

        // return the cached matrix
        return cached_l_t;
    }

    /**
     * Return the determinant of the matrix
     * @return determinant of the matrix
     */
    double get_determinant() 
    {
        double determinant = 1.0;
        for (int i{}; i < l_t_data.size(); ++i) 
        {
            double l_tii = l_t_data[i][i];
            determinant *= l_tii * l_tii;
        }
        return determinant;
    }

    /**
     * Get a solver for finding the A &times; X = B solution in least square sense.
     * @return a solver
     */
    Decomposition_Solver get_solver() 
    {
        return Solver();
    }

    /** Specialized solver. */
    class Solver : public Decomposition_Solver 
    {
    public:
        /** {@inherit_doc} */
        //override
        bool is_non_singular() 
        {
            // if we get this far, the matrix was positive definite, hence non-singular
            return true;
        }

        /** {@inherit_doc} */
        //override
        Real_Vector solve(const Real_Vector b) 
        {
            const int m = l_t_data.size();
            if (b.get_dimension() != m) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_dimension(), m);
            }

            const std::vector<double> x = b.to_array();

            // Solve LY = b
            for (int j{}; j < m; j++) 
            {
                const std::vector<double> lJ = l_t_data[j];
                x[j] /= lJ[j];
                const double xJ = x[j];
                for (int i = j + 1; i < m; i++) 
                {
                    x[i] -= xJ * lJ[i];
                }
            }

            // Solve LTX = Y
            for (int j = m - 1; j >= 0; j--) 
            {
                x[j] /= l_t_data[j][j];
                const double xJ = x[j];
                for (int i{}; i < j; i++) 
                {
                    x[i] -= xJ * l_t_data[i][j];
                }
            }

            return Array_Real_Vector(x, false);
        }

        /** {@inherit_doc} */
        //override
        Real_Matrix solve(Real_Matrix b) 
        {
            const int m = l_t_data.size();
            if (b.get_row_dimension() != m) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_row_dimension(), m);
            }

            const int& n_col_b = b.get_column_dimension();
            const std::vector<std::vector<double>> x = b.get_data();

            // Solve LY = b
            for (int j{}; j < m; j++) 
            {
                const std::vector<double> lJ = l_t_data[j];
                const double lJJ = lJ[j];
                const std::vector<double> xJ = x[j];
                for (int k{}; k < n_col_b; ++k) 
                {
                    xJ[k] /= lJJ;
                }
                for (int i = j + 1; i < m; i++) 
                {
                    const std::vector<double> xI = x[i];
                    const double lji = lJ[i];
                    for (int k{}; k < n_col_b; ++k) 
                    {
                        xI[k] -= xJ[k] * lji;
                    }
                }
            }

            // Solve LTX = Y
            for (int j = m - 1; j >= 0; j--) 
            {
                const double lJJ = l_t_data[j][j];
                const std::vector<double> xJ = x[j];
                for (int k{}; k < n_col_b; ++k) 
                {
                    xJ[k] /= lJJ;
                }
                for (int i{}; i < j; i++) 
                {
                    const std::vector<double> xI = x[i];
                    const double lij = l_t_data[i][j];
                    for (int k{}; k < n_col_b; ++k) 
                    {
                        xI[k] -= xJ[k] * lij;
                    }
                }
            }

            return Array_2D_Row_Real_Matrix(x);
        }

        /**
         * Get the inverse of the decomposed matrix.
         *
         * @return the inverse matrix.
         */
        //override
        Real_Matrix get_inverse() 
        {
            return solve(Matrix_Utils::create_real_identity_matrix(l_t_data.size()));
        }

        /** {@inherit_doc} */
        //override
        int get_row_dimension() 
        {
            return l_t_data.size();
        }

        /** {@inherit_doc} */
        //override
        int get_column_dimension() 
        {
            return l_t_data[0].size();
        }

    }

}


