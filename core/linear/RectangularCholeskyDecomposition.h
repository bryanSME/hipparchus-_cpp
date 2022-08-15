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
#include <vector>
#include "MatrixUtils.h"
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;

/**
 * Calculates the rectangular Cholesky decomposition of a matrix.
 * <p>The rectangular Cholesky decomposition of a real symmetric positive
 * semidefinite matrix A consists of a rectangular matrix B with the same
 * number of rows such that: A is almost equal to BB<sup>T</sup>, depending
 * on a user-defined tolerance. In a sense, this is the square root of A.</p>
 * <p>The difference with respect to the regular {@link Cholesky_Decomposition}
 * is that rows/columns may be permuted (hence the rectangular shape instead
 * of the traditional triangular shape) and there is a threshold to ignore
 * small diagonal elements. This is used for example to generate {@link
 * org.hipparchus.random.CorrelatedRandom_Vector_Generator correlated
 * random n-dimensions vectors} in a p-dimension subspace (p &lt; n).
 * In other words, it allows generating random vectors from a covariance
 * matrix that is only positive semidefinite, and not positive definite.</p>
 * <p>Rectangular Cholesky decomposition is <em>not</em> suited for solving
 * linear systems, so it does not provide any {@link Decomposition_Solver
 * decomposition solver}.</p>
 *
 * @see <a href="http://mathworld.wolfram.com/Cholesky_Decomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Cholesky_decomposition">Wikipedia</a>
 */
class RectangularCholesky_Decomposition 
{
private:
    /** Permutated Cholesky root of the symmetric positive semidefinite matrix. */
    Real_Matrix my_root;

    /** Rank of the symmetric positive semidefinite matrix. */
    int my_rank;

public:
    /**
     * Decompose a symmetric positive semidefinite matrix.
     * <p>
     * <b>Note:</b> this constructor follows the linpack method to detect dependent
     * columns by proceeding with the Cholesky algorithm until a nonpositive diagonal
     * element is encountered.
     *
     * @see <a href="http://eprints.ma.man.ac.uk/1193/01/covered/MIMS_ep2008_56.pdf">
     * Analysis of the Cholesky Decomposition of a Semi-definite Matrix</a>
     *
     * @param matrix Symmetric positive semidefinite matrix.
     * @exception  if the matrix is not
     * positive semidefinite.
     */
    RectangularCholesky_Decomposition(const Real_Matrix& matrix)
    {
        RectangularCholesky_Decomposition(matrix, 0);
    }

    /**
     * Decompose a symmetric positive semidefinite matrix.
     *
     * @param matrix Symmetric positive semidefinite matrix.
     * @param small Diagonal elements threshold under which columns are
     * considered to be dependent on previous ones and are discarded.
     * @exception  if the matrix is not
     * positive semidefinite.
     */
    RectangularCholesky_Decomposition(const Real_Matrix& matrix, const std::vector<double>& sma)
    {
        const auto order = matrix.get_row_dimension();
        auto c = matrix.get_data();
        auto b = std::vector<double>(order][order];

        auto index = std::vector<int>(order);
        for (int i{}; i < order; ++i) 
        {
            index[i] = i;
        }

        int r{};
        for (bool loop{ true }; loop;)
        {

            // find maximal diagonal element
            int swap_r{ r };
            for (int i{ r + 1 }; i < order; ++i)
            {
                int ii{ index[i] };
                int isr{ index[swap_r] };
                if (c[ii][ii] > c[isr][isr]) 
                {
                    swap_r = i;
                }
            }


            // swap elements
            if (swap_r != r) 
            {
                const int tmp_index = index[r];
                index[r] = index[swap_r];
                index[swap_r] = tmp_index;
                const std::vector<double> tmp_row = b[r];
                b[r] = b[swap_r];
                b[swap_r] = tmp_row;
            }

            // check diagonal element
            int ir{ index[r] };
            if (c[ir][ir] <= small) 
            {

                if (r == 0) 
                {
                    throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_DEFINITE_MATRIX);
                }

                // check remaining diagonal elements
                for (int i = r; i < order; ++i) 
                {
                    if (c[index[i]][index[i]] < -small) 
                    {
                        // there is at least one sufficiently negative diagonal element, // the symmetric positive semidefinite matrix is wrong
                        throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_DEFINITE_MATRIX);
                    }
                }

                // all remaining diagonal elements are close to zero, we consider we have
                // found the rank of the symmetric positive semidefinite matrix
                loop = false;

            }
            else 
            {

                // transform the matrix
                const double sqrt = std::sqrt(c[ir][ir]);
                b[r][r] = sqrt;
                const double inverse  = 1 / sqrt;
                const double inverse2 = 1 / c[ir][ir];
                for (int i = r + 1; i < order; ++i) 
                {
                    const int ii = index[i];
                    const double e = inverse * c[ii][ir];
                    b[i][r] = e;
                    c[ii][ii] -= c[ii][ir] * c[ii][ir] * inverse2;
                    for (int j = r + 1; j < i; ++j) 
                    {
                        const int ij = index[j];
                        const double f = c[ii][ij] - e * b[j][r];
                        c[ii][ij] = f;
                        c[ij][ii] = f;
                    }
                }

                // prepare next iteration
                loop = ++r < order;
            }
        }

        // build the root matrix
        my)rank = r;
        my_root = Matrix_Utils::create_real_matrix(order, r);
        for (int i{}; i < order; ++i) 
        {
            for (int j{}; j < r; ++j) 
            {
                my_root.set_entry(index[i], j, b[i][j]);
            }
        }

    }

    /** Get the root of the covariance matrix.
     * The root is the rectangular matrix <code>B</code> such that
     * the covariance matrix is equal to <code>B.B<sup>T</sup></code>
     * @return root of the square matrix
     * @see #get_rank()
     */
    Real_Matrix get_root_matrix() const 
    {
        return my_root;
    }

    /** Get the rank of the symmetric positive semidefinite matrix.
     * The r is the number of independent rows in the symmetric positive semidefinite
     * matrix, it is also the number of columns of the rectangular
     * matrix of the decomposition.
     * @return r of the square matrix.
     * @see #get_root_matrix()
     */
    int get_rank() const 
    {
        return my_rank;
    }
};