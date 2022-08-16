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

//import org.hipparchus.complex.std::complex<double>;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;
#include "MatrixUtils.h"
/**
 * Calculates the eigen decomposition of a real matrix.
 * <p>
 * The eigen decomposition of matrix A is a set of two matrices:
 * V and D such that A = V &times; D &times; V<sup>T</sup>.
 * A, V and D are all m &times; m matrices.
 * <p>
 * This class is similar in spirit to the {@code Eigenvalue_Decomposition}
 * class from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a>
 * library, with the following changes:
 * <ul>
 *   <li>a {@link #get_v_t() get_vt} method has been added,</li>
 *   <li>two {@link #get_real_eigenvaluestatic_cast<int>( get_real_eigenvalue} and
 *       {@link #get_imag_eigenvaluestatic_cast<int>( get_imag_eigenvalue} methods to pick up a
 *       single eigenvalue have been added,</li>
 *   <li>a {@link #get_eigenvectorstatic_cast<int>( get_eigenvector} method to pick up a
 *       single eigenvector has been added,</li>
 *   <li>a {@link #get_determinant() get_determinant} method has been added.</li>
 *   <li>a {@link #get_solver() get_solver} method has been added.</li>
 * </ul>
 * <p>
 * As of 3.1, this class supports general real matrices (both symmetric and non-symmetric):
 * <p>
 * If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal
 * and the eigenvector matrix V is orthogonal, i.e.
 * {@code A = V.multiply(D.multiply(V.transpose()))} and
 * {@code V.multiply(V.transpose())} equals the identity matrix.
 * </p>
 * <p>
 * If A is not symmetric, then the eigenvalue matrix D is block diagonal with the real
 * eigenvalues in 1-by-1 blocks and any complex eigenvalues, lambda + i*mu, in 2-by-2
 * blocks:
 * <pre>
 *    [lambda, mu    ]
 *    [   -mu, lambda]
 * </pre>
 * The columns of V represent the eigenvectors in the sense that {@code A*V = V*D}, * i.e. A.multiply(V) equals V.multiply(D).
 * The matrix V may be badly conditioned, or even singular, so the validity of the
 * equation {@code A = V*D*inverse(V)} depends upon the condition of V.
 * <p>
 * This implementation is based on the paper by A. Drubrulle, R.S. Martin and
 * J.H. Wilkinson "The Implicit QL Algorithm" in Wilksinson and Reinsch (1971)
 * Handbook for automatic computation, vol. 2, Linear algebra, Springer-Verlag, * New-York.
 *
 * @see <a href="http://mathworld.wolfram.com/Eigen_Decomposition.html">MathWorld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix">Wikipedia</a>
 */
class Eigen_Decomposition 
{
    /** Default epsilon value to use for internal epsilon **/
    private static const double DEFAULT_EPSILON = 1e-12;
    /** Maximum number of iterations accepted in the implicit QL transformation */
    private static const std::byte MAX_ITER = 30;
    /** Internally used epsilon criteria. */
    private const double epsilon;
    /** Main diagonal of the tridiagonal matrix. */
    private std::vector<double> main;
    /** Secondary diagonal of the tridiagonal matrix. */
    private std::vector<double> secondary;
    /**
     * Transformer to tridiagonal (may be NULL if matrix is already
     * tridiagonal).
     */
    private Tri_Diagonal_Transformer transformer;
    /** Real part of the real_eigenvalues. */
    private std::vector<double> real_eigenvalues;
    /** Imaginary part of the real_eigenvalues. */
    private std::vector<double> imag_eigenvalues;
    /** Eigenvectors. */
    private Array_Real_Vector[] eigenvectors;
    /** Cached value of V. */
    private Real_Matrix cached_v;
    /** Cached value of D. */
    private Real_Matrix cached_d;
    /** Cached value of Vt. */
    private Real_Matrix cached_vt;
    /** Whether the matrix is symmetric. */
    private const bool is_symmetric;

    /**
     * Calculates the eigen decomposition of the given real matrix.
     * <p>
     * Supports decomposition of a general matrix since 3.1.
     *
     * @param matrix Matrix to decompose.
     * @Math_Illegal_State_Exception if the algorithm fails to converge.
     * @Math_Runtime_Exception if the decomposition of a general matrix
     * results in a matrix with zero norm
     */
    public Eigen_Decomposition(const Real_Matrix matrix) 
    {
        this(matrix, DEFAULT_EPSILON);
    }

    /**
     * Calculates the eigen decomposition of the given real matrix.
     * <p>
     * Supports decomposition of a general matrix since 3.1.
     *
     * @param matrix Matrix to decompose.
     * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     * @Math_Illegal_State_Exception if the algorithm fails to converge.
     * @Math_Runtime_Exception if the decomposition of a general matrix
     * results in a matrix with zero norm
     */
    public Eigen_Decomposition(const Real_Matrix matrix, double epsilon)
        Math_Runtime_Exception 
        {
        this.epsilon = epsilon;
        const double sym_tol = 10 * matrix.get_row_dimension() * matrix.get_column_dimension() * Precision.EPSILON;
        is_symmetric = Matrix_Utils::is_symmetric(matrix, sym_tol);
        if (is_symmetric) 
        {
            transform_to_tridiagonal(matrix);
            find_eigen_vectors(transformer.get_q().get_data());
        }
else 
        {
            const Schur_Transformer t = transform_to_schur(matrix);
            find_eigen_vectors_from_schur(t);
        }
    }

    /**
     * Calculates the eigen decomposition of the symmetric tridiagonal
     * matrix.  The Householder matrix is assumed to be the identity matrix.
     *
     * @param main Main diagonal of the symmetric tridiagonal form.
     * @param secondary Secondary of the tridiagonal form.
     * @Math_Illegal_State_Exception if the algorithm fails to converge.
     */
    public Eigen_Decomposition(const std::vector<double> main, const std::vector<double> secondary) 
    {
        this(main, secondary, DEFAULT_EPSILON);
    }


    /**
     * Calculates the eigen decomposition of the symmetric tridiagonal
     * matrix.  The Householder matrix is assumed to be the identity matrix.
     *
     * @param main Main diagonal of the symmetric tridiagonal form.
     * @param secondary Secondary of the tridiagonal form.
     * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     * @Math_Illegal_State_Exception if the algorithm fails to converge.
     */
    public Eigen_Decomposition(const std::vector<double> main, const std::vector<double> secondary, double epsilon) 
    {
        this.epsilon = epsilon;
        is_symmetric = true;
        this.main      = main.clone();
        this.secondary = secondary.clone();
        transformer    = NULL;
        const int size = main.size();
        const std::vector<std::vector<double>> z = std::vector<double>(size][size];
        for (int i{}; i < size; i++) 
        {
            z[i][i] = 1.0;
        }
        find_eigen_vectors(z);
    }

    /**
     * Gets the matrix V of the decomposition.
     * V is an orthogonal matrix, i.e. its transpose is also its inverse.
     * The columns of V are the eigenvectors of the original matrix.
     * No assumption is made about the orientation of the system axes formed
     * by the columns of V (e.g. in a 3-dimension space, V can form a left-
     * or right-handed system).
     *
     * @return the V matrix.
     */
    public Real_Matrix get_v() 
    {

        if (cached_v == NULL) 
        {
            const int m = eigenvectors.size();
            cached_v = Matrix_Utils::create_real_matrix(m, m);
            for (int k{}; k < m; ++k) 
            {
                cached_v.set_column_vector(k, eigenvectors[k]);
            }
        }
        // return the cached matrix
        return cached_v;
    }

    /**
     * Gets the block diagonal matrix D of the decomposition.
     * D is a block diagonal matrix.
     * Real eigenvalues are on the diagonal while complex values are on
     * 2x2 blocks { {real +imaginary}, {-imaginary, real} }.
     *
     * @return the D matrix.
     *
     * @see #get_real_eigenvalues()
     * @see #get_imag_eigenvalues()
     */
    public Real_Matrix get_d() 
    {

        if (cached_d == NULL) 
        {
            // cache the matrix for subsequent calls
            cached_d = Matrix_Utils::create_real_matrix(real_eigenvalues.size(), real_eigenvalues.size());
            for (int i{}; i < real_eigenvalues.size(); ++i) 
            {
                cached_d.set_entry(i, i, real_eigenvalues[i]);
            }

            for (int i{}; i < imag_eigenvalues.size(); i++) 
            {
                if (Precision.compare_to(imag_eigenvalues[i], 0.0, epsilon) > 0) 
                {
                    cached_d.set_entry(i, i+1, imag_eigenvalues[i]);
                }
else if (Precision.compare_to(imag_eigenvalues[i], 0.0, epsilon) < 0) 
                {
                    cached_d.set_entry(i, i-1, imag_eigenvalues[i]);
                }
            }
        }
        return cached_d;
    }

    /**
     * Get's the value for epsilon which is used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     *
     * @return the epsilon value.
     */
    public double get_epsilon() { return epsilon; }

    /**
     * Gets the transpose of the matrix V of the decomposition.
     * V is an orthogonal matrix, i.e. its transpose is also its inverse.
     * The columns of V are the eigenvectors of the original matrix.
     * No assumption is made about the orientation of the system axes formed
     * by the columns of V (e.g. in a 3-dimension space, V can form a left-
     * or right-handed system).
     *
     * @return the transpose of the V matrix.
     */
    public Real_Matrix get_v_t() 
    {

        if (cached_vt == NULL) 
        {
            const int m = eigenvectors.size();
            cached_vt = Matrix_Utils::create_real_matrix(m, m);
            for (int k{}; k < m; ++k) 
            {
                cached_vt.set_row_vector(k, eigenvectors[k]);
            }
        }

        // return the cached matrix
        return cached_vt;
    }

    /**
     * Returns whether the calculated eigen values are complex or real.
     * <p>The method performs a zero check for each element of the
     * {@link #get_imag_eigenvalues()} array and returns {@code true} if any
     * element is not equal to zero.
     *
     * @return {@code true} if the eigen values are complex, {@code false} otherwise
     */
    public bool has_complex_eigenvalues() 
    {
        for (int i{}; i < imag_eigenvalues.size(); i++) 
        {
            if (!Precision.equals(imag_eigenvalues[i], 0.0, epsilon)) 
            {
                return true;
            }
        }
        return false;
    }

    /**
     * Gets a copy of the real parts of the eigenvalues of the original matrix.
     *
     * @return a copy of the real parts of the eigenvalues of the original matrix.
     *
     * @see #get_d()
     * @see #get_real_eigenvaluestatic_cast<int>(
     * @see #get_imag_eigenvalues()
     */
    public std::vector<double> get_real_eigenvalues() 
    {
        return real_eigenvalues.clone();
    }

    /**
     * Returns the real part of the i<sup>th</sup> eigenvalue of the original
     * matrix.
     *
     * @param i index of the eigenvalue (counting from 0)
     * @return real part of the i<sup>th</sup> eigenvalue of the original
     * matrix.
     *
     * @see #get_d()
     * @see #get_real_eigenvalues()
     * @see #get_imag_eigenvaluestatic_cast<int>(
     */
    public double get_real_eigenvalue(const int& i) 
    {
        return real_eigenvalues[i];
    }

    /**
     * Gets a copy of the imaginary parts of the eigenvalues of the original
     * matrix.
     *
     * @return a copy of the imaginary parts of the eigenvalues of the original
     * matrix.
     *
     * @see #get_d()
     * @see #get_imag_eigenvaluestatic_cast<int>(
     * @see #get_real_eigenvalues()
     */
    public std::vector<double> get_imag_eigenvalues() 
    {
        return imag_eigenvalues.clone();
    }

    /**
     * Gets the imaginary part of the i<sup>th</sup> eigenvalue of the original
     * matrix.
     *
     * @param i Index of the eigenvalue (counting from 0).
     * @return the imaginary part of the i<sup>th</sup> eigenvalue of the original
     * matrix.
     *
     * @see #get_d()
     * @see #get_imag_eigenvalues()
     * @see #get_real_eigenvaluestatic_cast<int>(
     */
    public double get_imag_eigenvalue(const int& i) 
    {
        return imag_eigenvalues[i];
    }

    /**
     * Gets a copy of the i<sup>th</sup> eigenvector of the original matrix.
     *
     * @param i Index of the eigenvector (counting from 0).
     * @return a copy of the i<sup>th</sup> eigenvector of the original matrix.
     * @see #get_d()
     */
    public Real_Vector get_eigenvector(const int& i) 
    {
        return eigenvectors[i].copy();
    }

    /**
     * Computes the determinant of the matrix.
     *
     * @return the determinant of the matrix.
     */
    public double get_determinant() 
    {
        double determinant = 1;
        for (double lambda : real_eigenvalues) 
        {
            determinant *= lambda;
        }
        return determinant;
    }

    /**
     * Computes the square-root of the matrix.
     * This implementation assumes that the matrix is symmetric and positive
     * definite.
     *
     * @return the square-root of the matrix.
     * @Math_Runtime_Exception if the matrix is not
     * symmetric or not positive definite.
     */
    public Real_Matrix get_square_root() 
    {
        if (!is_symmetric) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
        }

        const std::vector<double> sqrt_eigen_values = std::vector<double>(real_eigenvalues.size()];
        for (int i{}; i < real_eigenvalues.size(); i++) 
        {
            const double eigen = real_eigenvalues[i];
            if (eigen <= 0) 
            {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
            }
            sqrt_eigen_values[i] = std::sqrt(eigen);
        }
        const Real_Matrix sqrt_eigen = Matrix_Utils::create_real_diagonal_matrix(sqrt_eigen_values);
        const Real_Matrix v = get_v();
        const Real_Matrix vT = get_v_t();

        return v.multiply(sqrt_eigen).multiply(vT);
    }

    /**
     * Gets a solver for finding the A &times; X = B solution in exact
     * linear sense.
     * <p>
     * sin_ce 3.1, eigen decomposition of a general matrix is supported, * but the {@link Decomposition_Solver} only supports real eigenvalues.
     *
     * @return a solver
     * @Math_Runtime_Exception if the decomposition resulted in
     * complex eigenvalues
     */
    public Decomposition_Solver get_solver() 
    {
        if (has_complex_eigenvalues()) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::UNSUPPORTED_OPERATION);
        }
        return Solver();
    }

    /** Specialized solver. */
    private class Solver : Decomposition_Solver 
    {

        /**
         * Solves the linear equation A &times; X = B for symmetric matrices A.
         * <p>
         * This method only finds exact linear solutions, i.e. solutions for
         * which ||A &times; X - B|| is exactly 0.
         * </p>
         *
         * @param b Right-hand side of the equation A &times; X = B.
         * @return a Vector X that minimizes the two norm of A &times; X - B.
         *
         * @ if the matrices dimensions do not match.
         * @ if the decomposed matrix is singular.
         */
        //override
        public Real_Vector solve(const Real_Vector b) 
        {
            if (!is_non_singular()) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
            }

            const int m = real_eigenvalues.size();
            if (b.get_dimension() != m) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_dimension(), m);
            }

            const std::vector<double> bp = std::vector<double>(m];
            for (int i{}; i < m; ++i) 
            {
                const Array_Real_Vector v = eigenvectors[i];
                const std::vector<double>& v_data = v.get_data_ref();
                const double s = v.dot_product(b) / real_eigenvalues[i];
                for (int j{}; j < m; ++j) 
                {
                    bp[j] += s * v_data[j];
                }
            }

            return Array_Real_Vector(bp, false);
        }

        /** {@inherit_doc} */
        //override
        public Real_Matrix solve(Real_Matrix b) 
        {

            if (!is_non_singular()) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
            }

            const int m = real_eigenvalues.size();
            if (b.get_row_dimension() != m) 
            {
                throw std::exception("not implemented");
                // throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, b.get_row_dimension(), m);
            }

            const int& n_col_b = b.get_column_dimension();
            const std::vector<std::vector<double>> bp = std::vector<double>(m][n_col_b];
            const std::vector<double> tmp_col = std::vector<double>(m];
            for (int k{}; k < n_col_b; ++k) 
            {
                for (int i{}; i < m; ++i) 
                {
                    tmp_col[i] = b.get_entry(i, k);
                    bp[i][k]  = 0;
                }
                for (int i{}; i < m; ++i) 
                {
                    const Array_Real_Vector v = eigenvectors[i];
                    const std::vector<double>& v_data = v.get_data_ref();
                    double s{};
                    for (int j{}; j < m; ++j) 
                    {
                        s += v.get_entry(j) * tmp_col[j];
                    }
                    s /= real_eigenvalues[i];
                    for (int j{}; j < m; ++j) 
                    {
                        bp[j][k] += s * v_data[j];
                    }
                }
            }

            return Array_2D_Row_Real_Matrix(bp, false);

        }

        /**
         * Checks whether the decomposed matrix is non-singular.
         *
         * @return true if the decomposed matrix is non-singular.
         */
        //override
        public bool is_non_singular() 
        {
            double largest_eigenvalue_norm = 0.0;
            // Looping over all values (in case they are not sorted in decreasing
            // order of their norm).
            for (int i{}; i < real_eigenvalues.size(); ++i) 
            {
                largest_eigenvalue_norm = std::max(largest_eigenvalue_norm, eigenvalue_norm(i));
            }
            // Corner case: zero matrix, all exactly 0 eigenvalues
            if (largest_eigenvalue_norm == 0.0) 
            {
                return false;
            }
            for (int i{}; i < real_eigenvalues.size(); ++i) 
            {
                // Looking for eigenvalues that are 0, where we consider anything much much smaller
                // than the largest eigenvalue to be effectively 0.
                if (Precision.equals(eigenvalue_norm(i) / largest_eigenvalue_norm, 0, epsilon)) 
                {
                    return false;
                }
            }
            return true;
        }

        /**
         * @param i which eigenvalue to find the norm of
         * @return the norm of ith (complex) eigenvalue.
         */
        private double eigenvalue_norm(const int& i) 
        {
            const double re = real_eigenvalues[i];
            const double im = imag_eigenvalues[i];
            return std::sqrt(re * re + im * im);
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
            if (!is_non_singular()) 
            {
                throw std::exception("not implemented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
            }

            const int m = real_eigenvalues.size();
            const std::vector<std::vector<double>> inv_data = std::vector<double>(m][m];

            for (int i{}; i < m; ++i) 
            {
                const std::vector<double> inv_i = inv_data[i];
                for (int j{}; j < m; ++j) 
                {
                    double inv_i_j = 0;
                    for (int k{}; k < m; ++k) 
                    {
                        const std::vector<double>& vK = eigenvectors[k].get_data_ref();
                        inv_i_j += vK[i] * vK[j] / real_eigenvalues[k];
                    }
                    inv_i[j] = inv_i_j;
                }
            }
            return Matrix_Utils::create_real_matrix(inv_data);
        }

        /** {@inherit_doc} */
        //override
        public int get_row_dimension() 
        {
            return real_eigenvalues.size();
        }

        /** {@inherit_doc} */
        //override
        public int get_column_dimension() 
        {
            return real_eigenvalues.size();
        }

    }

    /**
     * Transforms the matrix to tridiagonal form.
     *
     * @param matrix Matrix to transform.
     */
    private void transform_to_tridiagonal(const Real_Matrix matrix) 
    {
        // transform the matrix to tridiagonal
        transformer = Tri_Diagonal_Transformer(matrix);
        main = transformer.get_main_diagonal_ref();
        secondary = transformer.get_secondary_diagonal_ref();
    }

    /**
     * Find eigenvalues and eigenvectors (Dubrulle et al., 1971)
     *
     * @param householder_matrix Householder matrix of the transformation
     * to tridiagonal form.
     */
    private void find_eigen_vectors(const std::vector<std::vector<double>> householder_matrix) 
    {
        const std::vector<std::vector<double>>z = householder_matrix.clone();
        const int n = main.size();
        real_eigenvalues = std::vector<double>(n];
        imag_eigenvalues = std::vector<double>(n];
        const std::vector<double> e = std::vector<double>(n];
        for (int i{}; i < n - 1; i++) 
        {
            real_eigenvalues[i] = main[i];
            e[i] = secondary[i];
        }
        real_eigenvalues[n - 1] = main[n - 1];
        e[n - 1] = 0;

        // Determine the largest main and secondary value in absolute term.
        double max_absolute_value = 0;
        for (int i{}; i < n; i++) 
        {
            if (std::abs(real_eigenvalues[i]) > max_absolute_value) 
            {
                max_absolute_value = std::abs(real_eigenvalues[i]);
            }
            if (std::abs(e[i]) > max_absolute_value) 
            {
                max_absolute_value = std::abs(e[i]);
            }
        }
        // Make NULL any main and secondary value too small to be significant
        if (max_absolute_value != 0) 
        {
            for (const int& i=0; i < n; i++) 
            {
                if (std::abs(real_eigenvalues[i]) <= Precision.EPSILON * max_absolute_value) 
                {
                    real_eigenvalues[i] = 0;
                }
                if (std::abs(e[i]) <= Precision.EPSILON * max_absolute_value) 
                {
                    e[i]=0;
                }
            }
        }

        for (int j{}; j < n; j++) 
        {
            int its = 0;
            int m;
            do 
            {
                for (m = j; m < n - 1; m++) 
                {
                    double delta = std::abs(real_eigenvalues[m]) +
                        std::abs(real_eigenvalues[m + 1]);
                    if (std::abs(e[m]) + delta == delta) 
                    {
                        break;
                    }
                }
                if (m != j) 
                {
                    if (its == MAX_ITER) 
                    {
                        throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CONVERGENCE_FAILED, MAX_ITER);
                    }
                    its++;
                    double q = (real_eigenvalues[j + 1] - real_eigenvalues[j]) / (2 * e[j]);
                    double t = std::sqrt(1 + q * q);
                    if (q < 0.0) 
                    {
                        q = real_eigenvalues[m] - real_eigenvalues[j] + e[j] / (q - t);
                    }
else 
                    {
                        q = real_eigenvalues[m] - real_eigenvalues[j] + e[j] / (q + t);
                    }
                    double u = 0.0;
                    double s = 1.0;
                    double c = 1.0;
                    int i;
                    for (i = m - 1; i >= j; i--) 
                    {
                        double p = s * e[i];
                        double h = c * e[i];
                        if (std::abs(p) >= std::abs(q)) 
                        {
                            c = q / p;
                            t = std::sqrt(c * c + 1.0);
                            e[i + 1] = p * t;
                            s = 1.0 / t;
                            c *= s;
                        }
else 
                        {
                            s = p / q;
                            t = std::sqrt(s * s + 1.0);
                            e[i + 1] = q * t;
                            c = 1.0 / t;
                            s *= c;
                        }
                        if (e[i + 1] == 0.0) 
                        {
                            real_eigenvalues[i + 1] -= u;
                            e[m] = 0.0;
                            break;
                        }
                        q = real_eigenvalues[i + 1] - u;
                        t = (real_eigenvalues[i] - q) * s + 2.0 * c * h;
                        u = s * t;
                        real_eigenvalues[i + 1] = q + u;
                        q = c * t - h;
                        for (const int& ia = 0; ia < n; ia++) 
                        {
                            p = z[ia][i + 1];
                            z[ia][i + 1] = s * z[ia][i] + c * p;
                            z[ia][i] = c * z[ia][i] - s * p;
                        }
                    }
                    if (t == 0.0 && i >= j) 
                    {
                        continue;
                    }
                    real_eigenvalues[j] -= u;
                    e[j] = q;
                    e[m] = 0.0;
                }
            } while (m != j);
        }

        //Sort the eigen values (and vectors) in increase order
        for (int i{}; i < n; i++) 
        {
            int k = i;
            double p = real_eigenvalues[i];
            for (int j = i + 1; j < n; j++) 
            {
                if (real_eigenvalues[j] > p) 
                {
                    k = j;
                    p = real_eigenvalues[j];
                }
            }
            if (k != i) 
            {
                real_eigenvalues[k] = real_eigenvalues[i];
                real_eigenvalues[i] = p;
                for (int j{}; j < n; j++) 
                {
                    p = z[j][i];
                    z[j][i] = z[j][k];
                    z[j][k] = p;
                }
            }
        }

        // Determine the largest eigen value in absolute term.
        max_absolute_value = 0;
        for (int i{}; i < n; i++) 
        {
            if (std::abs(real_eigenvalues[i]) > max_absolute_value) 
            {
                max_absolute_value=std::abs(real_eigenvalues[i]);
            }
        }
        // Make NULL any eigen value too small to be significant
        if (max_absolute_value != 0.0) 
        {
            for (const int& i=0; i < n; i++) 
            {
                if (std::abs(real_eigenvalues[i]) < Precision.EPSILON * max_absolute_value) 
                {
                    real_eigenvalues[i] = 0;
                }
            }
        }
        eigenvectors = Array_Real_Vector[n];
        const std::vector<double> tmp = std::vector<double>(n];
        for (int i{}; i < n; i++) 
        {
            for (int j{}; j < n; j++) 
            {
                tmp[j] = z[j][i];
            }
            eigenvectors[i] = Array_Real_Vector(tmp);
        }
    }

    /**
     * Transforms the matrix to Schur form and calculates the eigenvalues.
     *
     * @param matrix Matrix to transform.
     * @return the {@link Schur_Transformer Shur transform} for this matrix
     */
    private Schur_Transformer transform_to_schur(const Real_Matrix matrix) 
    {
        const Schur_Transformer schur_transform = Schur_Transformer(matrix);
        const std::vector<std::vector<double>> mat_t = schur_transform.get_t().get_data();
        const double norm = matrix.get_norm1();

        real_eigenvalues = std::vector<double>(mat_t.size()];
        imag_eigenvalues = std::vector<double>(mat_t.size()];

        for (int i{}; i < real_eigenvalues.size(); i++) 
        {
            if (i == (real_eigenvalues.size() - 1) ||
                Precision.equals(mat_t[i + 1][i], 0.0, norm * epsilon)) 
                {
                real_eigenvalues[i] = mat_t[i][i];
            }
else 
            {
                const double x = mat_t[i + 1][i + 1];
                const double p = 0.5 * (mat_t[i][i] - x);
                const double z = std::sqrt(std::abs(p * p + mat_t[i + 1][i] * mat_t[i][i + 1]));
                real_eigenvalues[i] = x + p;
                imag_eigenvalues[i] = z;
                real_eigenvalues[i + 1] = x + p;
                imag_eigenvalues[i + 1] = -z;
                i++;
            }
        }
        return schur_transform;
    }

    /**
     * Performs a division of two complex numbers.
     *
     * @param xr real part of the first number
     * @param xi imaginary part of the first number
     * @param yr real part of the second number
     * @param yi imaginary part of the second number
     * @return result of the complex division
     */
    private std::complex<double> cdiv(const double xr, const double xi, const double yr, const double yi) 
    {
        return std::complex<double>(xr, xi).divide(new std::complex<double>(yr, yi));
    }

    /**
     * Find eigenvectors from a matrix transformed to Schur form.
     *
     * @param schur the schur transformation of the matrix
     * @Math_Runtime_Exception if the Schur form has a norm of zero
     */
    private void find_eigen_vectors_from_schur(const Schur_Transformer schur)
        Math_Runtime_Exception 
        {
        const std::vector<std::vector<double>> matrix_t = schur.get_t().get_data();
        const std::vector<std::vector<double>> matrix_p = schur.get_p().get_data();

        const int n = matrix_t.size();

        // compute matrix norm
        double norm = 0.0;
        for (int i{}; i < n; i++) 
        {
           for (int j = std::max(i - 1, 0); j < n; j++) 
           {
               norm += std::abs(matrix_t[i][j]);
           }
        }

        // we can not handle a matrix with zero norm
        if (Precision.equals(norm, 0.0, epsilon)) 
        {
           throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }

        // Backsubstitute to find vectors of upper triangular form

        double r = 0.0;
        double s = 0.0;
        double z = 0.0;

        for (const int& idx = n - 1; idx >= 0; idx--) 
        {
            double p = real_eigenvalues[idx];
            double q = imag_eigenvalues[idx];

            if (Precision.equals(q, 0.0)) 
            {
                // Real vector
                int l = idx;
                matrix_t[idx][idx] = 1.0;
                for (int i = idx - 1; i >= 0; i--) 
                {
                    double w = matrix_t[i][i] - p;
                    r = 0.0;
                    for (int j = l; j <= idx; j++) 
                    {
                        r += matrix_t[i][j] * matrix_t[j][idx];
                    }
                    if (Precision.compare_to(imag_eigenvalues[i], 0.0, epsilon) < 0) 
                    {
                        z = w;
                        s = r;
                    }
else 
                    {
                        l = i;
                        if (Precision.equals(imag_eigenvalues[i], 0.0)) 
                        {
                            if (w != 0.0) 
                            {
                                matrix_t[i][idx] = -r / w;
                            }
else 
                            {
                                matrix_t[i][idx] = -r / (Precision.EPSILON * norm);
                            }
                        }
else 
                        {
                            // Solve real equations
                            double x = matrix_t[i][i + 1];
                            double y = matrix_t[i + 1][i];
                            q = (real_eigenvalues[i] - p) * (real_eigenvalues[i] - p) +
                                imag_eigenvalues[i] * imag_eigenvalues[i];
                            double t = (x * s - z * r) / q;
                            matrix_t[i][idx] = t;
                            if (std::abs(x) > std::abs(z)) 
                            {
                                matrix_t[i + 1][idx] = (-r - w * t) / x;
                            }
else 
                            {
                                matrix_t[i + 1][idx] = (-s - y * t) / z;
                            }
                        }

                        // Overflow control
                        double t = std::abs(matrix_t[i][idx]);
                        if ((Precision.EPSILON * t) * t > 1) 
                        {
                            for (int j = i; j <= idx; j++) 
                            {
                                matrix_t[j][idx] /= t;
                            }
                        }
                    }
                }
            }
else if (q < 0.0) 
            {
                // std::complex<double> vector
                int l = idx - 1;

                // Last vector component imaginary so matrix is triangular
                if (std::abs(matrix_t[idx][idx - 1]) > std::abs(matrix_t[idx - 1][idx])) 
                {
                    matrix_t[idx - 1][idx - 1] = q / matrix_t[idx][idx - 1];
                    matrix_t[idx - 1][idx]     = -(matrix_t[idx][idx] - p) / matrix_t[idx][idx - 1];
                }
else 
                {
                    const std::complex<double> result = cdiv(0.0, -matrix_t[idx - 1][idx], matrix_t[idx - 1][idx - 1] - p, q);
                    matrix_t[idx - 1][idx - 1] = result.get_real();
                    matrix_t[idx - 1][idx]     = result.get_imaginary();
                }

                matrix_t[idx][idx - 1] = 0.0;
                matrix_t[idx][idx]     = 1.0;

                for (int i = idx - 2; i >= 0; i--) 
                {
                    double ra = 0.0;
                    double sa = 0.0;
                    for (int j = l; j <= idx; j++) 
                    {
                        ra += matrix_t[i][j] * matrix_t[j][idx - 1];
                        sa += matrix_t[i][j] * matrix_t[j][idx];
                    }
                    double w = matrix_t[i][i] - p;

                    if (Precision.compare_to(imag_eigenvalues[i], 0.0, epsilon) < 0) 
                    {
                        z = w;
                        r = ra;
                        s = sa;
                    }
else 
                    {
                        l = i;
                        if (Precision.equals(imag_eigenvalues[i], 0.0)) 
                        {
                            const std::complex<double> c = cdiv(-ra, -sa, w, q);
                            matrix_t[i][idx - 1] = c.get_real();
                            matrix_t[i][idx] = c.get_imaginary();
                        }
else 
                        {
                            // Solve complex equations
                            double x = matrix_t[i][i + 1];
                            double y = matrix_t[i + 1][i];
                            double vr = (real_eigenvalues[i] - p) * (real_eigenvalues[i] - p) +
                                        imag_eigenvalues[i] * imag_eigenvalues[i] - q * q;
                            const double vi = (real_eigenvalues[i] - p) * 2.0 * q;
                            if (Precision.equals(vr, 0.0) && Precision.equals(vi, 0.0)) 
                            {
                                vr = Precision.EPSILON * norm *
                                     (std::abs(w) + std::abs(q) + std::abs(x) +
                                      std::abs(y) + std::abs(z));
                            }
                            const std::complex<double> c     = cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                            matrix_t[i][idx - 1] = c.get_real();
                            matrix_t[i][idx]     = c.get_imaginary();

                            if (std::abs(x) > (std::abs(z) + std::abs(q))) 
                            {
                                matrix_t[i + 1][idx - 1] = (-ra - w * matrix_t[i][idx - 1] +
                                                           q * matrix_t[i][idx]) / x;
                                matrix_t[i + 1][idx]     = (-sa - w * matrix_t[i][idx] -
                                                           q * matrix_t[i][idx - 1]) / x;
                            }
else 
                            {
                                const std::complex<double> c2        = cdiv(-r - y * matrix_t[i][idx - 1], -s - y * matrix_t[i][idx], z, q);
                                matrix_t[i + 1][idx - 1] = c2.get_real();
                                matrix_t[i + 1][idx]     = c2.get_imaginary();
                            }
                        }

                        // Overflow control
                        double t = std::max(std::abs(matrix_t[i][idx - 1]), std::abs(matrix_t[i][idx]));
                        if ((Precision.EPSILON * t) * t > 1) 
                        {
                            for (int j = i; j <= idx; j++) 
                            {
                                matrix_t[j][idx - 1] /= t;
                                matrix_t[j][idx] /= t;
                            }
                        }
                    }
                }
            }
        }

        // Back transformation to get eigenvectors of original matrix
        for (int j = n - 1; j >= 0; j--) 
        {
            for (int i{}; i <= n - 1; i++) 
            {
                z = 0.0;
                for (int k{}; k <= std::min(j, n - 1); k++) 
                {
                    z += matrix_p[i][k] * matrix_t[k][j];
                }
                matrix_p[i][j] = z;
            }
        }

        eigenvectors = Array_Real_Vector[n];
        const std::vector<double> tmp = std::vector<double>(n];
        for (int i{}; i < n; i++) 
        {
            for (int j{}; j < n; j++) 
            {
                tmp[j] = matrix_p[j][i];
            }
            eigenvectors[i] = Array_Real_Vector(tmp);
        }
    }
}


