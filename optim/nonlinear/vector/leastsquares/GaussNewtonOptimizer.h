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
//package org.hipparchus.optim.nonlinear.vector.leastsquares;

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.linear.Array_Real_Vector;
//import org.hipparchus.linear.Matrix_Decomposer;
//import org.hipparchus.linear.Matrix_Utils;
//import org.hipparchus.linear.QR_Decomposer;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Localized_Optim_Formats;
//import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem.Evaluation;
//import org.hipparchus.util.Incrementor;
//import org.hipparchus.util.Pair;
#include "../../../../core/linear/MatrixUtils.h"
/**
 * Gauss-_Newton least-squares solver.
 * <p>
 * This class solve a least-square problem by solving the normal equations
 * of the linearized problem at each iteration. Either LU decomposition or
 * Cholesky decomposition can be used to solve the normal equations, or QR
 * decomposition or SVD decomposition can be used to solve the linear system.
 * Cholesky/LU decomposition is faster but QR decomposition is more robust for difficult
 * problems, and SVD can compute a solution for rank-deficient problems.
 */
class Gauss_Newton_Optimizer : Least_Squares_Optimizer 
{

    /**
     * The singularity threshold for matrix decompositions. Determines when a {@link
     * Math_Illegal_State_Exception} is thrown. The current value was the default value for {@link
     * org.hipparchus.linear.LU_Decomposition}.
     */
    private static const double SINGULARITY_THRESHOLD = 1e-11;

    /** Decomposer */
    private const Matrix_Decomposer decomposer;

    /** Indicates if normal equations should be formed explicitly. */
    private const bool form_normal_equations;

    /**
     * Creates a Gauss Newton optimizer.
     * <p/>
     * The default for the algorithm is to use QR decomposition and not
     * form normal equations.
     * </p>
     */
    public Gauss_Newton_Optimizer() 
    {
        this(new QR_Decomposer(SINGULARITY_THRESHOLD), false);
    }

    /**
     * Create a Gauss Newton optimizer that uses the given matrix decomposition algorithm
     * to solve the normal equations.
     *
     * @param decomposer          the decomposition algorithm to use.
     * @param form_normal_equations whether the normal equations should be explicitly
     *                            formed. If {@code true} then {@code decomposer} is used
     *                            to solve J<sup>T</sup>Jx=J<sup>T</sup>r, otherwise
     *                            {@code decomposer} is used to solve Jx=r. If {@code
     *                            decomposer} can only solve square systems then this
     *                            parameter should be {@code true}.
     */
    public Gauss_Newton_Optimizer(const Matrix_Decomposer decomposer, const bool form_normal_equations) 
    {
        this.decomposer          = decomposer;
        this.form_normal_equations = form_normal_equations;
    }

    /**
     * Get the matrix decomposition algorithm.
     *
     * @return the decomposition algorithm.
     */
    public Matrix_Decomposer get_decomposer() 
    {
        return decomposer;
    }

    /**
     * Configure the matrix decomposition algorithm.
     *
     * @param new_decomposer the decomposition algorithm to use.
     * @return a instance.
     */
    public Gauss_Newton_Optimizer with_decomposer(const Matrix_Decomposer new_decomposer) 
    {
        return Gauss_Newton_Optimizer(new_decomposer, this.is_form_normal_equations());
    }

    /**
     * Get if the normal equations are explicitly formed.
     *
     * @return if the normal equations should be explicitly formed. If {@code true} then
     * {@code decomposer} is used to solve J<sup>T</sup>Jx=J<sup>T</sup>r, otherwise
     * {@code decomposer} is used to solve Jx=r.
     */
    public bool is_form_normal_equations() 
    {
        return form_normal_equations;
    }

    /**
     * Configure if the normal equations should be explicitly formed.
     *
     * @param new_form_normal_equations whether the normal equations should be explicitly
     *                               formed. If {@code true} then {@code decomposer} is used
     *                               to solve J<sup>T</sup>Jx=J<sup>T</sup>r, otherwise
     *                               {@code decomposer} is used to solve Jx=r. If {@code
     *                               decomposer} can only solve square systems then this
     *                               parameter should be {@code true}.
     * @return a instance.
     */
    public Gauss_Newton_Optimizer with_form_normal_equations(const bool new_form_normal_equations) 
    {
        return Gauss_Newton_Optimizer(this.get_decomposer(), new_form_normal_equations);
    }

    /** {@inherit_doc} */
    //override
    public Optimum optimize(const Least_Squares_Problem lsp) 
    {
        //create local evaluation and iteration counts
        const Incrementor evaluation_counter = lsp.get_evaluation_counter();
        const Incrementor iteration_counter = lsp.get_iteration_counter();
        const Convergence_Checker<Evaluation> checker
                = lsp.get_convergence_checker();

        // Computation will be useless without a checker (see "for-loop").
        if (checker == NULL) 
        {
            throw Null_Argument_Exception();
        }

        Real_Vector current_point = lsp.get_start();

        // iterate until convergence is reached
        Evaluation current = NULL;
        while (true) 
        {
            iteration_counter.increment();

            // evaluate the objective function and its jacobian
            Evaluation previous = current;
            // Value of the objective function at "current_point".
            evaluation_counter.increment();
            current = lsp.evaluate(current_point);
            const Real_Vector current_residuals = current.get_residuals();
            const Real_Matrix weighted_jacobian = current.get_jacobian();
            current_point = current.get_point();

            // Check convergence.
            if (previous != NULL &&
                checker.converged(iteration_counter.get_count(), previous, current)) 
                {
                return Optimum.of(current, evaluation_counter.get_count(), iteration_counter.get_count());
            }

            // solve the linearized least squares problem
            const Real_Matrix lhs; // left hand side
            const Real_Vector rhs; // right hand side
            if (this.form_normal_equations) 
            {
                const Pair<Real_Matrix, Real_Vector> normal_equation =
                        compute_normal_matrix(weighted_jacobian, current_residuals);
                lhs = normal_equation.get_first();
                rhs = normal_equation.get_second();
            }
else 
            {
                lhs = weighted_jacobian;
                rhs = current_residuals;
            }
            const Real_Vector dX;
            try 
            {
                dX = this.decomposer.decompose(lhs).solve(rhs);
            }
catch ( e) 
            {
                // change exception message
                throw Math_Illegal_State_Exception(
                        Localized_Optim_Formats.UNABLE_TO_SOLVE_SINGULAR_PROBLEM, e);
            }
            // update the estimated parameters
            current_point = current_point.add(dX);
        }
    }

    /** {@inherit_doc} */
    //override
    public std::string to_string() const 
    {
        return "Gauss_Newton_Optimizer{" +
                "decomposer=" + decomposer +
                ", form_normal_equations=" + form_normal_equations +
                '}';
    }

    /**
     * Compute the normal matrix, J<sup>T</sup>J.
     *
     * @param jacobian  the m by n jacobian matrix, J. Input.
     * @param residuals the m by 1 residual vector, r. Input.
     * @return  the n by n normal matrix and  the n by 1 J<sup>Tr vector.
     */
    private static Pair<Real_Matrix, Real_Vector> compute_normal_matrix(const Real_Matrix jacobian, const Real_Vector residuals) 
    {
        //since the normal matrix is symmetric, we only need to compute half of it.
        const int& n_r = jacobian.get_row_dimension();
        const int& n_c = jacobian.get_column_dimension();
        //allocate space for return values
        const Real_Matrix normal = Matrix_Utils.create_real_matrix(n_c, n_c);
        const Real_Vector j_tr = Array_Real_Vector(n_c);
        //for each measurement
        for (int i{}; i < n_r; ++i) 
        {
            //compute J_Tr for measurement i
            for (int j{}; j < n_c; j++) 
            {
                j_tr.set_entry(j, j_tr.get_entry(j) +
                        residuals.get_entry(i) * jacobian.get_entry(i, j));
            }

            // add the the contribution to the normal matrix for measurement i
            for (int k{}; k < n_c; ++k) 
            {
                //only compute the upper triangular part
                for (const int& l = k; l < n_c; ++l) 
                {
                    normal.set_entry(k, l, normal.get_entry(k, l) +
                            jacobian.get_entry(i, k) * jacobian.get_entry(i, l));
                }
            }
        }
        //copy the upper triangular part to the lower triangular part.
        for (int i{}; i < n_c; i++) 
        {
            for (int j{}; j < i; j++) 
            {
                normal.set_entry(i, j, normal.get_entry(j, i));
            }
        }
        return Pair<Real_Matrix, Real_Vector>(normal, j_tr);
    }

}


