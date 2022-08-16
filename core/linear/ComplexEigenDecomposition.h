#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
//package org.hipparchus.linear;

//import java.lang.reflect.Array;

//import org.hipparchus.complex.std::complex<double>;
//import org.hipparchus.complex.std::complex<double>_Field;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;
#include "MatrixUtils.h"
#include <complex>

/**
 * Given a matrix A, it computes a complex eigen decomposition AV = VD.
 *
 * <p>
 * std::complex<double> Eigen Decomposition differs from the {@link Eigen_Decomposition} since it
 * computes the eigen vectors as complex eigen vectors (if applicable).
 * </p>
 *
 * <p>
 * Beware that in the complex case, you do not always have \(V 	imes V^{T} = I\) or even a
 * diagonal matrix, even if the eigenvectors that form the columns of the V
 * matrix are independent. On example is the square matrix
 * \[
 * A = \left(\begin{matrix}
 * 3 & -2\\
 * 4 & -1
 * \end{matrix}\right)
 * \]
 * which has two conjugate eigenvalues \(\lambda_1=1+2i\) and \(\lambda_2=1-2i\)
 * with associated eigenvectors \(v_1^T = (1, 1-i)\) and \(v_2^T = (1, 1+i)\).
 * \[
 * V	imesV^T = \left(\begin{matrix}
 * 2 & 2\\
 * 2 & 0
 * \end{matrix}\right)
 * \]
 * which is not the identity matrix. Therefore, despite \(A 	imes V = V 	imes D\), * \(A \ne V 	imes D 	ime V^T\), which would hold for real eigendecomposition.
 * </p>
 *
 * Compute complex eigen values from the Schur transform. Compute complex eigen
 * vectors based on eigen values and the inverse iteration method.
 *
 * see: https://en.wikipedia.org/wiki/Inverse_iteration
 * https://en.wikiversity.org/wiki/Shifted_inverse_iteration
 * http://www.robots.ox.ac.uk/~sjrob/Teaching/EngComp/ecl4.pdf
 * http://www.math.ohiou.edu/courses/math3600/lecture16.pdf
 *
 */
class Complex_Eigen_Decomposition 
{
    /** Default threshold below which eigenvectors are considered equal. */
    public static const double DEFAULT_EIGENVECTORS_EQUALITY{ 1.0e-5 };
    /** Default value to use for internal epsilon. */
    public static const double DEFAULT_EPSILON = 1e-12;
    /** Internally used epsilon criteria for const AV=VD check. */
    public static const double DEFAULT_EPSILON_AV_VD_CHECK = 1e-6;
    /** Maximum number of inverse iterations. */
    private static const int MAX_ITER = 10;
    /** complex eigenvalues. */
    private std::vector<std::complex<double>>eigenvalues;
    /** Eigenvectors. */
    private Field_Vector<std::complex<double>>[] eigenvectors;
    /** Cached value of V. */
    private Field_Matrix<std::complex<double>> V;
    /** Cached value of D. */
    private Field_Matrix<std::complex<double>> D;
    /** Internally used threshold below which eigenvectors are considered equal. */
    private const double eigen_vectors_equality;
    /** Internally used epsilon criteria. */
    private const double epsilon;
    /** Internally used epsilon criteria for const AV=VD check. */
    private const double epsilon_a_v_v_d_check;

    /**
     * Constructor for decomposition.
     * <p>
     * This constructor uses the default values {@link #DEFAULT_EIGENVECTORS_EQUALITY}, * {@link #DEFAULT_EPSILON} and {@link #DEFAULT_EPSILON_AV_VD_CHECK}
     * </p>
     * @param matrix
     *            real matrix.
     */
    public Complex_Eigen_Decomposition(const Real_Matrix matrix) 
    {
        this(matrix, DEFAULT_EIGENVECTORS_EQUALITY, DEFAULT_EPSILON, DEFAULT_EPSILON_AV_VD_CHECK);
    }

    /**
     * Constructor for decomposition.
     * <p>
     * The {@code eigen_vectors_equality} threshold is used to ensure the L∞-normalized
     * eigenvectors found using inverse iteration are different from each other.
     * if \(min(|e_i-e_j|,|e_i+e_j|)\) is smaller than this threshold, the algorithm
     * considers it has found again an already known vector, so it drops it and attempts
     * a inverse iteration with a different start vector. This value should be
     * much larger than {@code epsilon} which is used for convergence
     * </p>
     * @param matrix real matrix.
     * @param eigen_vectors_equality threshold below which eigenvectors are considered equal
     * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     * @param epsilon_a_v_v_d_check Epsilon criteria for const AV=VD check
     * @since 1.8
     */
    public Complex_Eigen_Decomposition(const Real_Matrix matrix, const double eigen_vectors_equality, const double epsilon, const double epsilon_a_v_v_d_check) 
    {

        if (!matrix.is_square()) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, matrix.get_row_dimension(), matrix.get_column_dimension());
        }
        this.eigen_vectors_equality = eigen_vectors_equality;
        this.epsilon              = epsilon;
        this.epsilon_a_v_v_d_check     = epsilon_a_v_v_d_check;

        // computing the eigen values
        find_eigen_values(matrix);
        // computing the eigen vectors
        find_eigen_vectors(convertToField_Complex<double>(matrix));

        // V
        const int m = eigenvectors.size();
        V = Matrix_Utils::create_field_matrix(std::complex<double>_Field.get_instance(), m, m);
        for (int k{}; k < m; ++k) 
        {
            V.set_column_vector(k, eigenvectors[k]);
        }

        // D
        D = Matrix_Utils::createFieldDiagonal_Matrix(eigenvalues);

        check_definition(matrix);
    }

    /**
     * Getter of the eigen values.
     *
     * @return eigen values.
     */
    public std::vector<std::complex<double>>get_eigenvalues() 
    {
        return eigenvalues.clone();
    }

    /**
     * Getter of the eigen vectors.
     *
     * @param i
     *            which eigen vector.
     * @return eigen vector.
     */
    public Field_Vector<std::complex<double>> get_eigenvector(const int& i) 
    {
        return eigenvectors[i].copy();
    }

    /** Reset eigenvalues and eigen vectors from matrices.
     * <p>
     * This method is intended to be called by sub-classes (mainly {@link OrderedComplex_Eigen_Decomposition})
     * that reorder the matrices elements. It rebuild the eigenvalues and eigen vectors arrays
     * from the D and V matrices.
     * </p>
     * @since 2.1
     */
    protected void matrices_to_eigen_arrays() 
    {
        for (int i{}; i < eigenvalues.size(); ++i) 
        {
            eigenvalues[i] = D.get_entry(i, i);
        }
        for (int i{}; i < eigenvectors.size(); ++i) 
        {
            for (int j{}; j < eigenvectors[i].get_dimension(); ++j) 
            {
                eigenvectors[i].set_entry(j, V.get_entry(j, i));
            }
        }
    }

    /**
     * Confirm if there are complex eigen values.
     *
     * @return true if there are complex eigen values.
     */
    public bool has_complex_eigenvalues() 
    {
        for (int i{}; i < eigenvalues.size(); i++) 
        {
            if (!Precision.equals(eigenvalues[i].get_imaginary(), 0.0, epsilon)) 
            {
                return true;
            }
        }
        return false;
    }

    /**
     * Computes the determinant.
     *
     * @return the determinant.
     */
    public double get_determinant() 
    {
        std::complex<double> determinant = std::complex<double>(1, 0);
        for (std::complex<double> lambda : eigenvalues) 
        {
            determinant = determinant.multiply(lambda);
        }
        return determinant.get_real();
    }

    /**
     * Getter V.
     *
     * @return V.
     */
    public Field_Matrix<std::complex<double>> get_v() 
    {
        return V;
    }

    /**
     * Getter D.
     *
     * @return D.
     */
    public Field_Matrix<std::complex<double>> get_d() 
    {
        return D;
    }

    /**
     * Getter VT.
     *
     * @return VT.
     */
    public Field_Matrix<std::complex<double>> get_v_t() 
    {
        return V.transpose();
    }

    /**
     * Compute eigen values using the Schur transform.
     *
     * @param matrix
     *            real matrix to compute eigen values.
     */
    protected void find_eigen_values(const Real_Matrix matrix) 
    {
        const Schur_Transformer schur_transform = Schur_Transformer(matrix);
        const std::vector<std::vector<double>> mat_t = schur_transform.get_t().get_data();

        eigenvalues = std::complex<double>[mat_t.size()];

        for (int i{}; i < eigenvalues.size(); i++) 
        {
            if (i == (eigenvalues.size() - 1) || Precision.equals(mat_t[i + 1][i], 0.0, epsilon)) 
            {
                eigenvalues[i] = std::complex<double>(mat_t[i][i]);
            }
else 
            {
                const double x = mat_t[i + 1][i + 1];
                const double p = 0.5 * (mat_t[i][i] - x);
                const double z = std::sqrt(std::abs(p * p + mat_t[i + 1][i] * mat_t[i][i + 1]));
                eigenvalues[i] = std::complex<double>(x + p, z);
                eigenvalues[i + 1] = std::complex<double>(x + p, -z);
                i++;
            }
        }

    }

    /**
     * Compute the eigen vectors using the inverse power method.
     *
     * @param matrix
     *            real matrix to compute eigen vectors.
     */
    //@Suppress_Warnings("unchecked")
    protected void find_eigen_vectors(const Field_Matrix<std::complex<double>> matrix) 
    {
        // number of eigen values/vectors
        int n = eigenvalues.size();

        // eigen vectors
        eigenvectors = (Field_Vector<std::complex<double>>[]) Array.new_instance(Field_Vector.class, n);

        // computing eigen vector based on eigen values and inverse iteration
        for (int i{}; i < eigenvalues.size(); i++) 
        {

            // shifted non-singular matrix matrix A-(λ+ε)I that is close to the singular matrix A-λI
            std::complex<double> mu = eigenvalues[i].add(epsilon);
            const Field_Matrix<std::complex<double>> shifted = matrix.copy();
            for (int k{}; k < matrix.get_column_dimension(); ++k) 
            {
                shifted.set_entry(k, k, shifted.get_entry(k, k).subtract(mu));
            }

            // solver for linear system (A - (λ+ε)I) Bₖ₊₁ = Bₖ
            FieldDecomposition_Solver<std::complex<double>> solver = FieldQR_Decomposition<>(shifted).get_solver();

            // loop over possible start vectors
            for (const int& p = 0; eigenvectors[i] == NULL && p < matrix.get_column_dimension(); ++p) 
            {

                // find a vector to start iterations
                Field_Vector<std::complex<double>> b = findStart(p);

                if (get_norm(b).norm() > Precision.SAFE_MIN) 
                {
                    // start vector is a good candidate for inverse iteration

                    // perform inverse iteration
                    double delta = INFINITY;
                    for (int k{}; delta > epsilon && k < MAX_ITER; k++) 
                    {

                        // solve (A - (λ+ε)) Bₖ₊₁ = Bₖ
                        const Field_Vector<std::complex<double>> bNext = solver.solve(b);

                        // normalize according to L∞ norm
                        normalize(bNext);

                        // compute convergence criterion, comparing Bₖ and both ±Bₖ₊₁
                        // as iterations sometimes flip between two opposite vectors
                        delta = separation(b, bNext);

                        // prepare next iteration
                        b = bNext;

                    }

                    // check we have not found again an already known vector
                    for (int j = 0; b != NULL && j < i; ++j) 
                    {
                        if (separation(eigenvectors[j], b) <= eigen_vectors_equality) 
                        {
                            // the selected start vector leads us to found a known vector again, // we must try another start
                            b = NULL;
                        }
                    }
                    eigenvectors[i] = b;

                }
            }

        }
    }

    /** Find a start vector orthogonal to all already found normalized eigenvectors.
     * @param index index of the vector
     * @return start vector
     */
    private Field_Vector<std::complex<double>> findStart(const int index) 
    {

        // create vector
        const Field_Vector<std::complex<double>> start =
                        Matrix_Utils::create_field_vector(std::complex<double>_Field.get_instance(), eigenvalues.size());

        // initialize with a canonical vector
        start.set_entry(index, std::complex<double>.ONE);

        return start;

    }

    /**
     * Compute the L∞ norm of the a given vector.
     *
     * @param vector
     *            vector.
     * @return L∞ norm.
     */
    private std::complex<double> get_norm(Field_Vector<std::complex<double>> vector) 
    {
        double  normR = 0;
        std::complex<double> norm_c = std::complex<double>.ZERO;
        for (int i{}; i < vector.get_dimension(); i++) 
        {
            const std::complex<double> ci = vector.get_entry(i);
            const double  ni = std::hypot(ci.get_real(), ci.get_imaginary());
            if (ni > normR) 
            {
                normR = ni;
                norm_c = ci;
            }
        }
        return norm_c;
    }

    /** Normalize a vector with respect to L∞ norm.
     * @param v vector to normalized
     */
    private void normalize(const Field_Vector<std::complex<double>> v) 
    {
        const std::complex<double> inv_norm = get_norm(v).reciprocal();
        for (int j{}; j < v.get_dimension(); ++j) 
        {
            v.set_entry(j, v.get_entry(j).multiply(inv_norm));
        }
    }

    /** Compute the separation between two normalized vectors (which may be in opposite directions).
     * @param v1 first normalized vector
     * @param v2 second normalized vector
     * @return min (|v1 - v2|, |v1+v2|)
     */
    private double separation(const Field_Vector<std::complex<double>> v1, const Field_Vector<std::complex<double>> v2) 
    {
        double delta_plus  = 0;
        double delta_minus = 0;
        for (int j{}; j < v1.get_dimension(); ++j) 
        {
            const std::complex<double> bCurrj = v1.get_entry(j);
            const std::complex<double> bNextj = v2.get_entry(j);
            delta_plus  = std::max(delta_plus, std::hypot(bNextj.get_real()      + bCurrj.get_real(), bNextj.get_imaginary() + bCurrj.get_imaginary()));
            delta_minus = std::max(delta_minus, std::hypot(bNextj.get_real()      - bCurrj.get_real(), bNextj.get_imaginary() - bCurrj.get_imaginary()));
        }
        return std::min(delta_plus, delta_minus);
    }

    /**
     * Check definition of the decomposition in runtime.
     *
     * @param matrix
     *            matrix to be decomposed.
     */
    protected void check_definition(const Real_Matrix matrix) 
    {
        Field_Matrix<std::complex<double>> matrix_c = convertToField_Complex<double>(matrix);

        // checking definition of the decomposition
        // testing A*V = V*D
        Field_Matrix<std::complex<double>> AV = matrix_c.multiply(get_v());
        Field_Matrix<std::complex<double>> VD = get_v().multiply(get_d());
        if (!equalsWith_precision(AV, VD, epsilon_a_v_v_d_check)) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::FAILED_DECOMPOSITION, matrix.get_row_dimension(), matrix.get_column_dimension());

        }

    }

    /**
     * Helper method that checks with two matrix is equals taking into account a
     * given precision.
     *
     * @param matrix1 first matrix to compare
     * @param matrix2 second matrix to compare
     * @param tolerance tolerance on matrices entries
     * @return true is matrices entries are equal within tolerance, * false otherwise
     */
    private bool equalsWith_precision(const Field_Matrix<std::complex<double>> matrix1, const Field_Matrix<std::complex<double>> matrix2, const double& tolerance) 
    {
        bool to_ret = true;
        for (int i{}; i < matrix1.get_row_dimension(); i++) 
        {
            for (int j{}; j < matrix1.get_column_dimension(); j++) 
            {
                std::complex<double> c1 = matrix1.get_entry(i, j);
                std::complex<double> c2 = matrix2.get_entry(i, j);
                if (c1.add(c2.negate()).norm() > tolerance) 
                {
                    to_ret = false;
                    break;
                }
            }
        }
        return to_ret;
    }

    /**
     * It converts a real matrix into a complex field matrix.
     *
     * @param matrix
     *            real matrix.
     * @return complex matrix.
     */
    private Field_Matrix<std::complex<double>> convertToField_Complex<double>(Real_Matrix matrix) 
    {
        const Field_Matrix<std::complex<double>> to_ret =
                        Matrix_Utils::create_field_identity_matrix(std::complex<double>_Field.get_instance(), matrix.get_row_dimension());
        for (int i{}; i < to_ret.get_row_dimension(); i++) 
        {
            for (int j{}; j < to_ret.get_column_dimension(); j++) 
            {
                to_ret.set_entry(i, j, std::complex<double>(matrix.get_entry(i, j)));
            }
        }
        return to_ret;
    }
}


