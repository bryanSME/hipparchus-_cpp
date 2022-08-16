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
//import org.hipparchus.util.Precision;

/**
 * Class transforming a general real matrix to Hessenberg form.
 * <p>A m &times; m matrix A can be written as the product of three matrices: A = P
 * &times; H &times; P<sup>T</sup> with P an orthogonal matrix and H a Hessenberg
 * matrix. Both P and H are m &times; m matrices.</p>
 * <p>Transformation to Hessenberg form is often not a goal by itself, but it is an
 * intermediate step in more general decomposition algorithms like
 * {@link Eigen_Decomposition eigen decomposition}. This class is therefore
 * intended for internal use by the library and is not public. As a consequence
 * of this explicitly limited scope, many methods directly returns references to
 * internal arrays, not copies.</p>
 * <p>This class is based on the method orthes in class Eigenvalue_Decomposition
 * from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library.</p>
 *
 * @see <a href="http://mathworld.wolfram.com/HessenbergDecomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Householder_transformation">Householder Transformations</a>
 */
class Hessenberg_Transformer 
{
private:
    /** Householder vectors. */
    const std::vector<std::vector<double>> my_householder_vectors;
    /** Temporary storage vector. */
    const std::vector<double> my_ort;
    /** Cached value of P. */
    Real_Matrix my_cached_p;
    /** Cached value of Pt. */
    Real_Matrix my_cached_pt;
    /** Cached value of H. */
    Real_Matrix my_cached_h;

public:
    /**
     * Build the transformation to Hessenberg form of a general matrix.
     *
     * @param matrix matrix to transform
     * @ if the matrix is not square
     */
    Hessenberg_Transformer(const Real_Matrix matrix) 
    {
        if (!matrix.is_square()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, matrix.get_row_dimension(), matrix.get_column_dimension());
        }

        const int m = matrix.get_row_dimension();
        householder_vectors = matrix.get_data();
        ort = std::vector<double>(m];
        cached_p = NULL;
        cached_pt = NULL;
        cached_h = NULL;

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
            const int n = householder_vectors.size();
            const int high = n - 1;
            const std::vector<std::vector<double>> pa = std::vector<double>(n][n];

            for (int i{}; i < n; i++) 
            {
                for (int j{}; j < n; j++) 
                {
                    pa[i][j] = (i == j) ? 1 : 0;
                }
            }

            for (const int& m = high - 1; m >= 1; m--) 
            {
                if (householder_vectors[m][m - 1] != 0.0) 
                {
                    for (int i = m + 1; i <= high; i++) 
                    {
                        ort[i] = householder_vectors[i][m - 1];
                    }

                    for (int j = m; j <= high; j++) 
                    {
                        double g = 0.0;

                        for (int i = m; i <= high; i++) 
                        {
                            g += ort[i] * pa[i][j];
                        }

                        // Double division avoids possible underflow
                        g = (g / ort[m]) / householder_vectors[m][m - 1];

                        for (int i = m; i <= high; i++) 
                        {
                            pa[i][j] += g * ort[i];
                        }
                    }
                }
            }

            cached_p = Matrix_Utils::create_real_matrix(pa);
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
     * Returns the Hessenberg matrix H of the transform.
     *
     * @return the H matrix
     */
    public Real_Matrix get_h() 
    {
        if (cached_h == NULL) 
        {
            const int m = householder_vectors.size();
            const std::vector<std::vector<double>> h = std::vector<double>(m][m];
            for (int i{}; i < m; ++i) 
            {
                if (i > 0) 
                {
                    // copy the entry of the lower sub-diagonal
                    h[i][i - 1] = householder_vectors[i][i - 1];
                }

                // copy upper triangular part of the matrix
                for (int j = i; j < m; ++j) 
                {
                    h[i][j] = householder_vectors[i][j];
                }
            }
            cached_h = Matrix_Utils::create_real_matrix(h);
        }

        // return the cached matrix
        return cached_h;
    }

    /**
     * Get the Householder vectors of the transform.
     * <p>Note that since this class is only intended for internal use, it returns
     * directly a reference to its internal arrays, not a copy.</p>
     *
     * @return the main diagonal elements of the B matrix
     */
    std::vector<std::vector<double>> get_householder_vectors_ref() 
    {
        return householder_vectors; // NOPMD - returning an internal array is intentional and documented here
    }

    /**
     * Transform original matrix to Hessenberg form.
     * <p>Transformation is done using Householder transforms.</p>
     */
    private void transform() 
    {
        const int n = householder_vectors.size();
        const int high = n - 1;

        for (const int& m = 1; m <= high - 1; m++) 
        {
            // Scale column.
            double scale = 0;
            for (int i = m; i <= high; i++) 
            {
                scale += std::abs(householder_vectors[i][m - 1]);
            }

            if (!Precision.equals(scale, 0)) 
            {
                // Compute Householder transformation.
                double h = 0;
                for (int i = high; i >= m; i--) 
                {
                    ort[i] = householder_vectors[i][m - 1] / scale;
                    h += ort[i] * ort[i];
                }
                const double g = (ort[m] > 0) ? -std::sqrt(h) : std::sqrt(h);

                h -= ort[m] * g;
                ort[m] -= g;

                // Apply Householder similarity transformation
                // H = (I - u*u' / h) * H * (I - u*u' / h)

                for (int j = m; j < n; j++) 
                {
                    double f = 0;
                    for (int i = high; i >= m; i--) 
                    {
                        f += ort[i] * householder_vectors[i][j];
                    }
                    f /= h;
                    for (int i = m; i <= high; i++) 
                    {
                        householder_vectors[i][j] -= f * ort[i];
                    }
                }

                for (int i{}; i <= high; i++) 
                {
                    double f = 0;
                    for (int j = high; j >= m; j--) 
                    {
                        f += ort[j] * householder_vectors[i][j];
                    }
                    f /= h;
                    for (int j = m; j <= high; j++) 
                    {
                        householder_vectors[i][j] -= f * ort[j];
                    }
                }

                ort[m] = scale * ort[m];
                householder_vectors[m][m - 1] = scale * g;
            }
        }
    }
}


