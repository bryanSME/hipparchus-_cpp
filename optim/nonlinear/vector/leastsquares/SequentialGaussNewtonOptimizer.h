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
//package org.hipparchus.optim.nonlinear.vector.leastsquares;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.linear.Array_Real_Vector;
//import org.hipparchus.linear.Cholesky_Decomposition;
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
 * Sequential Gauss-_Newton least-squares solver.
 * <p>
 * This class solve a least-square problem by solving the normal equations of
 * the linearized problem at each iteration.
 * </p>
 *
 */
class Sequential_Gauss_Newton_Optimizer : Least_Squares_Optimizer 
{

    /**
     * The singularity threshold for matrix decompositions. Determines when a
     * {@link Math_Illegal_State_Exception} is thrown. The current value was the
     * default value for {@link LU_Decomposition}.
     */
    private static const double SINGULARITY_THRESHOLD = 1e-11;

    /** Decomposer. */
    private const Matrix_Decomposer decomposer;

    /** Indicates if normal equations should be formed explicitly. */
    private const bool form_normal_equations;

    /** Old evaluation previously computed. */
    private const Evaluation old_evaluation;

    /** Old jacobian previously computed. */
    private const Real_Matrix old_lhs;

    /** Old residuals previously computed. */
    private const Real_Vector old_rhs;

    /**
     * Create a sequential Gauss Newton optimizer.
     * <p/>
     * The default for the algorithm is to use QR decomposition, not
     * form normal equations and have no previous evaluation
     * </p>
     *
     */
    public Sequential_Gauss_Newton_Optimizer() 
    {
        this(new QR_Decomposer(SINGULARITY_THRESHOLD), false, NULL);
    }

    /**
     * Create a sequential Gauss Newton optimizer that uses the given matrix
     * decomposition algorithm to solve the normal equations.
     * <p>
     * The {@code decomposer} is used to solve J<sup>T</sup>Jx=J<sup>T</sup>r.
     * </p>
     *
     * @param decomposer the decomposition algorithm to use.
     * @param form_normal_equations whether the normal equations should be explicitly
     *                            formed. If {@code true} then {@code decomposer} is used
     *                            to solve J<sup>T</sup>Jx=J<sup>T</sup>r, otherwise
     *                            {@code decomposer} is used to solve Jx=r. If {@code
     *                            decomposer} can only solve square systems then this
     *                            parameter should be {@code true}.
     * @param evaluation old evaluation previously computed, NULL if there are no previous evaluations.
     */
    public Sequential_Gauss_Newton_Optimizer(const Matrix_Decomposer decomposer, const bool form_normal_equations, const Evaluation evaluation) 
    {
        this.decomposer          = decomposer;
        this.form_normal_equations = form_normal_equations;
        this.old_evaluation       = evaluation;
        if (evaluation == NULL) 
        {
            this.old_lhs = NULL;
            this.old_rhs = NULL;
        }
else 
        {
            if (form_normal_equations) 
            {
                const Pair<Real_Matrix, Real_Vector> normal_equation =
                                compute_normal_matrix(evaluation.get_jacobian(), evaluation.get_residuals());
                // solve the linearized least squares problem
                this.old_lhs = normal_equation.get_first();
                this.old_rhs = normal_equation.get_second();
            }
else 
            {
                this.old_lhs = evaluation.get_jacobian();
                this.old_rhs = evaluation.get_residuals();
            }
        }
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
    public Sequential_Gauss_Newton_Optimizer with_decomposer(const Matrix_Decomposer new_decomposer) 
    {
        return Sequential_Gauss_Newton_Optimizer(new_decomposer, this.is_form_normal_equations(), this.get_old_evaluation());
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
    public Sequential_Gauss_Newton_Optimizer with_form_normal_equations(const bool new_form_normal_equations) 
    {
        return Sequential_Gauss_Newton_Optimizer(this.get_decomposer(), new_form_normal_equations, this.get_old_evaluation());
    }

    /**
     * Get the previous evaluation used by the optimizer.
     *
     * @return the previous evaluation.
     */
    public Evaluation get_old_evaluation() 
    {
        return old_evaluation;
    }

    /**
     * Configure the previous evaluation used by the optimizer.
     * <p>
     * This building method uses a complete evaluation to retrieve
     * a priori data. Note that as {@link #with_a_priori_data(Real_Vector, Real_Matrix)}
     * generates a fake evaluation and calls this method, either
     * {@link #with_a_priori_data(Real_Vector, Real_Matrix)} or {@link #with_evaluation(Least_Squares_Problem.Evaluation)}
     * should be called, but not both as the last one called will //override the previous one.
     * </p>
     * @param previous_evaluation the previous evaluation used by the optimizer.
     * @return a instance.
     */
    public Sequential_Gauss_Newton_Optimizer with_evaluation(const Evaluation previous_evaluation) 
    {
        return Sequential_Gauss_Newton_Optimizer(this.get_decomposer(), this.is_form_normal_equations(), previous_evaluation);
    }

    /**
     * Configure from a priori state and covariance.
     * <p>
     * This building method generates a fake evaluation and calls
     * {@link #with_evaluation(Least_Squares_Problem.Evaluation)}, so either
     * {@link #with_a_priori_data(Real_Vector, Real_Matrix)} or {@link #with_evaluation(Least_Squares_Problem.Evaluation)}
     * should be called, but not both as the last one called will //override the previous one.
     * </p>
     * @param a_priori_state a priori state to use
     * @param a_priori_covariance a priori covariance to use
     * @return a instance.
     */
    public Sequential_Gauss_Newton_Optimizer with_a_priori_data(const Real_Vector a_priori_state, const Real_Matrix& a_priori_covariance) 
    {

        // we consider the a priori state and covariance come from a
        // previous estimation with exactly one observation of each state
        // component, so partials are the identity matrix, weight is the
        // square root of inverse of covariance, and residuals are zero

        // create a fake weighted Jacobian
        const Real_Matrix jTj              = get_decomposer().decompose(a_priori_covariance).get_inverse();
        const Real_Matrix weighted_jacobian = Cholesky_Decomposition(jTj).get_l_t();

        // create fake zero residuals
        const Real_Vector residuals        = Matrix_Utils.create_real__vector(a_priori_state.get_dimension());

        // combine everything as an evaluation
        const Evaluation fake_evaluation   = Abstract_Evaluation(a_priori_state.get_dimension()) 
        {

            /** {@inherit_doc} */
            //override
            public Real_Vector get_residuals() 
            {
                return residuals;
            }

            /** {@inherit_doc} */
            //override
            public Real_Vector get_point() 
            {
                return a_priori_state;
            }

            /** {@inherit_doc} */
            //override
            public Real_Matrix get_jacobian() 
            {
                return weighted_jacobian;
            }
        };

        return with_evaluation(fake_evaluation);

    }

    /** {@inherit_doc} */
    //override
    public Optimum optimize(const Least_Squares_Problem lsp) 
    {
        // create local evaluation and iteration counts
        const Incrementor evaluation_counter = lsp.get_evaluation_counter();
        const Incrementor iteration_counter = lsp.get_iteration_counter();
        const Convergence_Checker<Evaluation> checker =
            lsp.get_convergence_checker();

        // Computation will be useless without a checker (see "for-loop").
        if (checker == NULL) 
        {
            throw Null_Argument_Exception();
        }

        Real_Vector current_point = lsp.get_start();

        if (old_evaluation != NULL &&
            current_point.get_dimension() != old_evaluation.get_point().get_dimension()) 
            {
            throw Math_Illegal_State_Exception(Localized_Core_Formats.DIMENSIONS_MISMATCH, current_point.get_dimension(), old_evaluation.get_point().get_dimension());
        }

        // iterate until convergence is reached
        Evaluation current = NULL;
        while (true) 
        {
            iteration_counter.increment();

            // evaluate the objective function and its jacobian
            const Evaluation previous = current;

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
                // combine old and evaluations
                const Evaluation combined_evaluation = old_evaluation == NULL ?
                                                      current :
                                                      Combined_Evaluation(old_evaluation, current);
                return Optimum.of(combined_evaluation, evaluation_counter.get_count(), iteration_counter.get_count());
            }

           // solve the linearized least squares problem
            const Real_Matrix lhs; // left hand side
            const Real_Vector rhs; // right hand side
            if (this.form_normal_equations) 
            {
                const Pair<Real_Matrix, Real_Vector> normal_equation =
                                compute_normal_matrix(weighted_jacobian, current_residuals);

                lhs = old_lhs == NULL ?
                      normal_equation.get_first() :
                      normal_equation.get_first().add(old_lhs); // left hand side
                rhs = old_rhs == NULL ?
                      normal_equation.get_second() :
                      normal_equation.get_second().add(old_rhs); // right hand side
            }
else 
            {
                lhs = old_lhs == NULL ?
                      weighted_jacobian :
                      combine_jacobians(old_lhs, weighted_jacobian);
                rhs = old_rhs == NULL ?
                      current_residuals :
                      combine_residuals(old_rhs, current_residuals);
            }

            const Real_Vector dX;
            try 
            {
                dX = this.decomposer.decompose(lhs).solve(rhs);
            }
catch ( e) 
            {
                // change exception message
                throw Math_Illegal_State_Exception(Localized_Optim_Formats.UNABLE_TO_SOLVE_SINGULAR_PROBLEM, e);
            }
            // update the estimated parameters
            current_point = current_point.add(dX);

        }
    }

    /** {@inherit_doc} */
    //override
    public std::string to_string() const 
    {
        return "Sequential_Gauss_Newton_Optimizer{" +
               "decomposer=" + decomposer + '}';
    }

    /**
     * Compute the normal matrix, J<sup>T</sup>J.
     *
     * @param jacobian  the m by n jacobian matrix, J. Input.
     * @param residuals the m by 1 residual vector, r. Input.
     * @return  the n by n normal matrix and the n by 1 J<sup>Tr</sup> vector.
     */
    private static Pair<Real_Matrix, Real_Vector>
        compute_normal_matrix(const Real_Matrix jacobian, const Real_Vector residuals) 
        {
        // since the normal matrix is symmetric, we only need to compute half of
        // it.
        const int& n_r = jacobian.get_row_dimension();
        const int& n_c = jacobian.get_column_dimension();
        // allocate space for return values
        const Real_Matrix normal = Matrix_Utils.create_real_matrix(n_c, n_c);
        const Real_Vector j_tr = Array_Real_Vector(n_c);
        // for each measurement
        for (int i{}; i < n_r; ++i) 
        {
            // compute J_Tr for measurement i
            for (int j{}; j < n_c; j++) 
            {
                j_tr.set_entry(j, j_tr.get_entry(j) +
                                residuals.get_entry(i) *
                                               jacobian.get_entry(i, j));
            }

            // add the the contribution to the normal matrix for measurement i
            for (int k{}; k < n_c; ++k) 
            {
                // only compute the upper triangular part
                for (const int& l = k; l < n_c; ++l) 
                {
                    normal
                        .set_entry(k, l, normal.get_entry(k, l) +
                                        jacobian.get_entry(i, k) *
                                                       jacobian.get_entry(i, l));
                }
            }
        }
        // copy the upper triangular part to the lower triangular part.
        for (int i{}; i < n_c; i++) 
        {
            for (int j{}; j < i; j++) 
            {
                normal.set_entry(i, j, normal.get_entry(j, i));
            }
        }
        return Pair<Real_Matrix, Real_Vector>(normal, j_tr);
    }

    /** Combine Jacobian matrices
     * @param old_jacobian old Jacobian matrix
     * @param new_jacobian Jacobian matrix
     * @return combined Jacobian matrix
     */
    private static Real_Matrix combine_jacobians(const Real_Matrix old_jacobian, const Real_Matrix new_jacobian) 
    {
        const int old_row_dimension    = old_jacobian.get_row_dimension();
        const int old_column_dimension = old_jacobian.get_column_dimension();
        const Real_Matrix jacobian =
                        Matrix_Utils.create_real_matrix(old_row_dimension + new_jacobian.get_row_dimension(), old_column_dimension);
        jacobian.set_sub_matrix(old_jacobian.get_data(), 0,               0);
        jacobian.set_sub_matrix(new_jacobian.get_data(), old_row_dimension, 0);
        return jacobian;
    }

    /** Combine residuals vectors
     * @param old_residuals old residuals vector
     * @param new_residuals residuals vector
     * @return combined residuals vector
     */
    private static Real_Vector combine_residuals(const Real_Vector old_residuals, const Real_Vector new_residuals) 
    {
        return old_residuals.append(new_residuals);
    }

    /**
     * Container with an old and a evaluation and combine both of them
     */
    private static class Combined_Evaluation extends Abstract_Evaluation 
    {

        /** Point of evaluation. */
        private const Real_Vector point;

        /** Derivative at point. */
        private const Real_Matrix jacobian;

        /** Computed residuals. */
        private const Real_Vector residuals;

        /**
         * Create an {@link Evaluation} with no weights.
         *
         * @param old_evaluation the old evaluation.
         * @param new_evaluation the evaluation
         */
        private Combined_Evaluation(const Evaluation old_evaluation, const Evaluation new_evaluation) 
        {

            super(old_evaluation.get_residuals().get_dimension() +
                  new_evaluation.get_residuals().get_dimension());

            this.point    = new_evaluation.get_point();
            this.jacobian = combine_jacobians(old_evaluation.get_jacobian(), new_evaluation.get_jacobian());
            this.residuals = combine_residuals(old_evaluation.get_residuals(), new_evaluation.get_residuals());
        }

        /** {@inherit_doc} */
        //override
        public Real_Matrix get_jacobian() 
        {
            return jacobian;
        }

        /** {@inherit_doc} */
        //override
        public Real_Vector get_point() 
        {
            return point;
        }

        /** {@inherit_doc} */
        //override
        public Real_Vector get_residuals() 
        {
            return residuals;
        }

    }

}


