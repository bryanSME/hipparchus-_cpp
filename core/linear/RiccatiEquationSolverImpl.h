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

//import org.hipparchus.complex.std::complex<double>;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.util.FastMath;
#include "MatrixUtils.h"
#include "RiccatiEquationSolver.h"
/**
 *
 * This solver computes the solution using the following approach:
 *
 * 1. Compute the Hamiltonian matrix 2. Extract its complex eigen vectors (not
 * the best solution, a better solution would be ordered Schur transformation)
 * 3. Approximate the initial solution given by 2 using the Kleinman algorithm
 * (an iterative method)
 */
class Riccati_Equation_SolverImpl : public Riccati_Equation_Solver 
{
private:
    /** Internally used maximum iterations. */
    static constexpr int MAX_ITERATIONS{ 100 };

    /** Internally used epsilon criteria. */
    static constexpr double EPSILON{ 1e-8 };

    /** The solution of the algebraic Riccati equation. */
    const Real_Matrix P;

    /** The computed K. */
    const Real_Matrix K;

public:
    /**
     * Constructor of the solver. A and B should be compatible. B and R must be
     * multiplicative compatible. A and Q must be multiplicative compatible. R
     * must be invertible.
     *
     * @param A state transition matrix
     * @param B control multipliers matrix
     * @param Q state cost matrix
     * @param R control cost matrix
     */
    Riccati_Equation_SolverImpl(const Real_Matrix A, const Real_Matrix B, const Real_Matrix Q, const Real_Matrix R) 
    {

        // checking A
        if (!A.is_square()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NON_SQUARE_MATRIX, A.get_row_dimension(), A.get_column_dimension());
        }
        if (A.get_column_dimension() != B.get_row_dimension()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, A.get_row_dimension(), B.get_row_dimension());
        }
        Matrix_Utils::check_multiplication_compatible(B, R);
        Matrix_Utils::check_multiplication_compatible(A, Q);

        // checking R
        const Singular_Value_Decomposition svd = Singular_Value_Decomposition(R);
        if (!svd.get_solver().is_non_singular()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
        }

        const Real_Matrix R_inv = svd.get_solver().get_inverse();

        P = compute_p(A, B, Q, R, R_inv, MAX_ITERATIONS, EPSILON);

        K = R_inv.multiply_transposed(B).multiply(P);
    }

    /**
     * Compute an initial stable solution and then applies the Kleinman
     * algorithm to approximate it using an EPSILON.
     *
     * @param A state transition matrix
     * @param B control multipliers matrix
     * @param Q state cost matrix
     * @param R control cost matrix
     * @param R_inv inverse of matrix R
     * @param max_iterations maximum number of iterations
     * @param epsilon epsilon to be used
     * @return matrix P, solution of the algebraic Riccati equation
     */
    private Real_Matrix compute_p(const Real_Matrix A, const Real_Matrix B, const Real_Matrix Q, const Real_Matrix R, const Real_Matrix R_inv, const int max_iterations, const double epsilon) 
    {
        const Real_Matrix P_ = compute_initial_p(A, B, Q, R_inv);
        return approximate_p(A, B, Q, R, R_inv, P_, max_iterations, epsilon);
    }

    /**
     * Compute initial P using the Hamiltonian and the ordered eigen values
     * decomposition.
     *
     * @param A state transition matrix
     * @param B control multipliers matrix
     * @param Q state cost matrix
     * @param R_inv inverse of matrix R
     * @return initial solution
     */
    private Real_Matrix compute_initial_p(const Real_Matrix A, const Real_Matrix B, const Real_Matrix Q, const Real_Matrix R_inv) 
    {
        const Real_Matrix B_tran = B.transpose();

        // computing the Hamiltonian Matrix
        const Real_Matrix m11 = A;
        const Real_Matrix m12 = B.multiply(R_inv).multiply(B_tran).scalar_multiply(-1).scalar_add(0);
        const Real_Matrix m21 = Q.scalar_multiply(-1).scalar_add(0);
        const Real_Matrix m22 = A.transpose().scalar_multiply(-1).scalar_add(0);
        if (m11.get_row_dimension() != m12.get_row_dimension()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, m11.get_row_dimension(), m12.get_row_dimension());
        }
        if (m21.get_row_dimension() != m22.get_row_dimension()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, m21.get_row_dimension(), m22.get_row_dimension());
        }
        if (m11.get_column_dimension() != m21.get_column_dimension()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, m11.get_column_dimension(), m21.get_column_dimension());
        }
        if (m21.get_column_dimension() != m22.get_column_dimension()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, m21.get_column_dimension(), m22.get_column_dimension());
        }

        // defining M
        const Real_Matrix m = Matrix_Utils::create_real_matrix(m11.get_row_dimension() + m21.get_row_dimension(), m11.get_column_dimension() + m12.get_column_dimension());
        // defining submatrixes
        m.set_sub_matrix(m11.get_data(), 0, 0);
        m.set_sub_matrix(m12.get_data(), 0, m11.get_column_dimension());
        m.set_sub_matrix(m21.get_data(), m11.get_row_dimension(), 0);
        m.set_sub_matrix(m22.get_data(), m11.get_row_dimension(), m11.get_column_dimension());

        // eigen decomposition
        // numerically bad, but it is used for the initial stable solution for
        // the
        // Kleinman Algorithm
        // it must be ordered in order to work with submatrices
        const OrderedComplex_Eigen_Decomposition eigen_decomposition =
                        OrderedComplex_Eigen_Decomposition(m, EPSILON, Complex_Eigen_Decomposition.DEFAULT_EPSILON, Complex_Eigen_Decomposition.DEFAULT_EPSILON_AV_VD_CHECK);
        const Field_Matrix<std::complex<double>> u = eigen_decomposition.get_v();

        // solving linear system
        // P = U_21*U_11^-1
        const Field_Matrix<std::complex<double>> u11 = u.get_sub_matrix(0, m11.get_row_dimension() - 1, 0, m11.get_column_dimension() - 1);
        const Field_Matrix<std::complex<double>> u21 = u.get_sub_matrix(m11.get_row_dimension(), 2 * m11.get_row_dimension() - 1, 0, m11.get_column_dimension() - 1);

        const FieldDecomposition_Solver<std::complex<double>> solver = FieldLU_Decomposition<std::complex<double>>(u11).get_solver();

        if (!solver.is_non_singular()) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
        }

        // solving U_11^{-1}
        Field_Matrix<std::complex<double>> u11_inv = solver.get_inverse();

        // P = U_21*U_11^-1
        Field_Matrix<std::complex<double>> p = u21.multiply(u11_inv);

        // converting to realmatrix - ignoring precision errors in imaginary
        // components
        return convert_to_real__matrix(p, Double.MAX_VALUE);

    }

    /**
     * Applies the Kleinman's algorithm.
     *
     * @param A state transition matrix
     * @param B control multipliers matrix
     * @param Q state cost matrix
     * @param R control cost matrix
     * @param R_inv inverse of matrix R
     * @param initial_p initial solution
     * @param max_iterations maximum number of iterations allowed
     * @param epsilon convergence threshold
     * @return improved solution
     */
    private Real_Matrix approximate_p(const Real_Matrix A, const Real_Matrix B, const Real_Matrix Q, const Real_Matrix R, const Real_Matrix R_inv, const Real_Matrix initial_p, const int max_iterations, const double epsilon) 
    {
        Real_Matrix P_ = initial_p;

        double error = 1;
        int i = 1;
        while (error > epsilon) 
        {
            const Real_Matrix K_ = P_.multiply(B).multiply(R_inv).scalar_multiply(-1);

            // X = AA+BB*K1';
            const Real_Matrix X = A.add(B.multiply_transposed(K_));
            // Y = -K1*RR*K1' - QQ;
            const Real_Matrix Y = K_.multiply(R).multiply_transposed(K_).scalar_multiply(-1).subtract(Q);

            const Array_2D_Row_Real_Matrix X_ = (Array_2D_Row_Real_Matrix) X.transpose();
            const Array_2D_Row_Real_Matrix Y_ = (Array_2D_Row_Real_Matrix) Y;
            const Array_2D_Row_Real_Matrix eye_x =
                            (Array_2D_Row_Real_Matrix) Matrix_Utils::create_real_identity_matrix(X_.get_row_dimension());

            // X1=kron(X',eye(size(X))) + kron(eye(size(X)),X');
            const Real_Matrix X__ = X_.kronecker_product(eye_x).add(eye_x.kronecker_product(X_));
            // Y1=reshape(Y,prod(size(Y)),1); %%stack
            const Real_Matrix Y__ = Y_.stack();

            // PX = inv(X1)*Y1;
            // sensitive to numerical erros
            // const Real_Matrix PX = Matrix_Utils::inverse(X__).multiply(Y__);
            Decomposition_Solver solver = LU_Decomposition(X__).get_solver();
            if (!solver.is_non_singular()) 
            {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::SINGULAR_MATRIX);
            }
            const Real_Matrix PX = solver.solve(Y__);

            // P = reshape(PX,sqrt(length(PX)),sqrt(length(PX))); %%unstack
            const Real_Matrix P__ = ((Array_2D_Row_Real_Matrix) PX).unstack_square();

            // aerror = norm(P - P1);
            const Real_Matrix diff = P__.subtract(P_);
            // calculationg l2 norm
            const Singular_Value_Decomposition svd = Singular_Value_Decomposition(diff);
            error = svd.get_norm();

            P_ = P__;
            i++;
            if (i > max_iterations) 
            {
                throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::CONVERGENCE_FAILED);
            }
        }

        return P_;
    }

    /** {inherit_doc} */
    //override
    public Real_Matrix get_p() 
    {
        return P;
    }

    /** {inherit_doc} */
    //override
    public Real_Matrix get_k() 
    {
        return K;
    }

    /**
     * Converts a given complex matrix into a real matrix taking into account a
     * precision for the imaginary components.
     *
     * @param matrix complex field matrix
     * @param tolerance tolerance on the imaginary part
     * @return real matrix.
     */
    private Real_Matrix convert_to_real__matrix(Field_Matrix<std::complex<double>> matrix, double tolerance) 
    {
        const Real_Matrix to_ret = Matrix_Utils::create_real_matrix(matrix.get_row_dimension(), matrix.get_row_dimension());
        for (int i{}; i < to_ret.get_row_dimension(); i++) 
        {
            for (int j{}; j < to_ret.get_column_dimension(); j++) 
            {
                std::complex<double> c = matrix.get_entry(i, j);
                if (c.get_imaginary() != 0 && std::abs(c.get_imaginary()) > tolerance) 
                {
                    throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::COMPLEX_CANNOT_BE_CONSIDERED_A_REAL_NUMBER, c.get_real(), c.get_imaginary());
                }
                to_ret.set_entry(i, j, c.get_real());
            }
        }
        return to_ret;
    }

}


