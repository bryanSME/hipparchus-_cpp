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

//import java.util.Arrays;

//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.linear.Array_Real_Vector;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.optim.Localized_Optim_Formats;
//import org.hipparchus.optim.nonlinear.vector.leastsquares.Least_Squares_Problem.Evaluation;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Incrementor;
//import org.hipparchus.util.Precision;


/**
 * This class solves a least-squares problem using the Levenberg-_Marquardt
 * algorithm.
 *
 * <p>This implementation <em>should</em> work even for over-determined systems
 * (i.e. systems having more point than equations). Over-determined systems
 * are solved by ignoring the point which have the smallest impact according
 * to their jacobian column norm. Only the rank of the matrix and some loop bounds
 * are changed to implement this.</p>
 *
 * <p>The resolution engine is a simple translation of the MINPACK <a
 * href="http://www.netlib.org/minpack/lmder.f">lmder</a> routine with minor
 * changes. The changes include the over-determined resolution, the use of
 * inherited convergence checker and the Q.R. decomposition which has been
 * rewritten following the algorithm described in the
 * P. Lascaux and R. Theodor book <i>Analyse num&eacute;rique matricielle
 * appliqu&eacute;e &agrave; l'art de l'ing&eacute;nieur</i>, Masson 1986.</p>
 * <p>The authors of the original fortran version are:
 * <ul>
 * <li>Argonne National Laboratory. MINPACK project. March 1980</li>
 * <li>Burton S. Garbow</li>
 * <li>Kenneth E. Hillstrom</li>
 * <li>Jorge J. More</li>
 * </ul>
 * The redistribution policy for MINPACK is available <a
 * href="http://www.netlib.org/minpack/disclaimer">here</a>, for convenience, it
 * is reproduced below.</p>
 *
 * <table border="0" width="80%" cellpadding="10" align="center" bgcolor="#E0E0E0">
 * <tr><td>
 *    Minpack Copyright Notice (1999) University of Chicago.
 *    All rights reserved
 * </td></tr>
 * <tr><td>
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * <ol>
 *  <li>Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.</li>
 * <li>Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.</li>
 * <li>The end-user documentation included with the redistribution, if any, *     must include the following acknowledgment:
 *     <code>This product includes software developed by the University of
 *           Chicago, as Operator of Argonne National Laboratory.</code>
 *     Alternately, this acknowledgment may appear in the software itself, *     if and wherever such third-party acknowledgments normally appear.</li>
 * <li><strong>WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
 *     WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
 *     UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
 *     THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
 *     IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
 *     OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
 *     OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
 *     USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
 *     THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
 *     DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
 *     UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
 *     BE CORRECTED.</strong></li>
 * <li><strong>LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
 *     HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
 *     ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT, *     INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
 *     ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
 *     PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
 *     SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
 *     (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE, *     EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
 *     POSSIBILITY OF SUCH LOSS OR DAMAGES.</strong></li>
 * <ol></td></tr>
 * </table>
 *
 */
class Levenberg_Marquardt_Optimizer : Least_Squares_Optimizer 
{

    /** Twice the "epsilon machine". */
    private static const double TWO_EPS = 2 * Precision.EPSILON;

    /* configuration parameters */
    /** Positive input variable used in determining the initial step bound. */
    private const double initial_step_bound_factor;
    /** Desired relative error in the sum of squares. */
    private const double cost_relative_tolerance;
    /**  Desired relative error in the approximate solution parameters. */
    private const double par_relative_tolerance;
    /** Desired max cosine on the orthogonality between the function vector
     * and the columns of the jacobian. */
    private const double ortho_tolerance;
    /** Threshold for QR ranking. */
    private const double qr_ranking_threshold;

    /** Default constructor.
     * <p>
     * The default values for the algorithm settings are:
     * <ul>
     *  <li>Initial step bound factor: 100</li>
     *  <li>Cost relative tolerance: 1e-10</li>
     *  <li>Parameters relative tolerance: 1e-10</li>
     *  <li>Orthogonality tolerance: 1e-10</li>
     *  <li>QR ranking threshold: {@link Precision#SAFE_MIN}</li>
     * </ul>
     **/
    public Levenberg_Marquardt_Optimizer() 
    {
        this(100, 1e-10, 1e-10, 1e-10, Precision.SAFE_MIN);
    }

    /**
     * Construct an instance with all parameters specified.
     *
     * @param initial_step_bound_factor initial step bound factor
     * @param cost_relative_tolerance  cost relative tolerance
     * @param par_relative_tolerance   parameters relative tolerance
     * @param ortho_tolerance         orthogonality tolerance
     * @param qr_ranking_threshold     threshold in the QR decomposition. Columns with a 2
     *                               norm less than this threshold are considered to be
     *                               all 0s.
     */
    public Levenberg_Marquardt_Optimizer(
            const double initial_step_bound_factor, const double cost_relative_tolerance, const double par_relative_tolerance, const double ortho_tolerance, const double qr_ranking_threshold) 
            {
        this.initial_step_bound_factor = initial_step_bound_factor;
        this.cost_relative_tolerance = cost_relative_tolerance;
        this.par_relative_tolerance = par_relative_tolerance;
        this.ortho_tolerance = ortho_tolerance;
        this.qr_ranking_threshold = qr_ranking_threshold;
    }

    /**
     * @param new_initial_step_bound_factor Positive input variable used in
     * determining the initial step bound. This bound is set to the
     * product of initial_step_bound_factor and the euclidean norm of
     * {@code diag * x} if non-zero, or else to {@code new_initial_step_bound_factor}
     * itself. In most cases factor should lie in the interval
     * {@code (0.1, 100.0)}. {@code 100} is a generally recommended value.
     * of the matrix is reduced.
     * @return a instance.
     */
    public Levenberg_Marquardt_Optimizer with_initial_step_bound_factor(double new_initial_step_bound_factor) 
    {
        return Levenberg_Marquardt_Optimizer(
                new_initial_step_bound_factor, cost_relative_tolerance, par_relative_tolerance, ortho_tolerance, qr_ranking_threshold);
    }

    /**
     * @param new_cost_relative_tolerance Desired relative error in the sum of squares.
     * @return a instance.
     */
    public Levenberg_Marquardt_Optimizer with_cost_relative_tolerance(double new_cost_relative_tolerance) 
    {
        return Levenberg_Marquardt_Optimizer(
                initial_step_bound_factor, new_cost_relative_tolerance, par_relative_tolerance, ortho_tolerance, qr_ranking_threshold);
    }

    /**
     * @param new_par_relative_tolerance Desired relative error in the approximate solution
     * parameters.
     * @return a instance.
     */
    public Levenberg_Marquardt_Optimizer with_parameter_relative_tolerance(double new_par_relative_tolerance) 
    {
        return Levenberg_Marquardt_Optimizer(
                initial_step_bound_factor, cost_relative_tolerance, new_par_relative_tolerance, ortho_tolerance, qr_ranking_threshold);
    }

    /**
     * Modifies the given parameter.
     *
     * @param new_ortho_tolerance Desired max cosine on the orthogonality between
     * the function vector and the columns of the Jacobian.
     * @return a instance.
     */
    public Levenberg_Marquardt_Optimizer with_ortho_tolerance(double new_ortho_tolerance) 
    {
        return Levenberg_Marquardt_Optimizer(
                initial_step_bound_factor, cost_relative_tolerance, par_relative_tolerance, new_ortho_tolerance, qr_ranking_threshold);
    }

    /**
     * @param new_q_r_ranking_threshold Desired threshold for QR ranking.
     * If the squared norm of a column vector is smaller or equal to this
     * threshold during QR decomposition, it is considered to be a zero vector
     * and hence the rank of the matrix is reduced.
     * @return a instance.
     */
    public Levenberg_Marquardt_Optimizer with_ranking_threshold(double new_q_r_ranking_threshold) 
    {
        return Levenberg_Marquardt_Optimizer(
                initial_step_bound_factor, cost_relative_tolerance, par_relative_tolerance, ortho_tolerance, new_q_r_ranking_threshold);
    }

    /**
     * Gets the value of a tuning parameter.
     * @see #with_initial_step_bound_factorstatic_cast<double>(
     *
     * @return the parameter's value.
     */
    public double get_initial_step_bound_factor() 
    {
        return initial_step_bound_factor;
    }

    /**
     * Gets the value of a tuning parameter.
     * @see #with_cost_relative_tolerancestatic_cast<double>(
     *
     * @return the parameter's value.
     */
    public double get_cost_relative_tolerance() 
    {
        return cost_relative_tolerance;
    }

    /**
     * Gets the value of a tuning parameter.
     * @see #with_parameter_relative_tolerancestatic_cast<double>(
     *
     * @return the parameter's value.
     */
    public double get_parameter_relative_tolerance() 
    {
        return par_relative_tolerance;
    }

    /**
     * Gets the value of a tuning parameter.
     * @see #with_ortho_tolerancestatic_cast<double>(
     *
     * @return the parameter's value.
     */
    public double get_ortho_tolerance() 
    {
        return ortho_tolerance;
    }

    /**
     * Gets the value of a tuning parameter.
     * @see #with_ranking_thresholdstatic_cast<double>(
     *
     * @return the parameter's value.
     */
    public double get_ranking_threshold() 
    {
        return qr_ranking_threshold;
    }

    /** {@inherit_doc} */
    //override
    public Optimum optimize(const Least_Squares_Problem problem) 
    {
        // Pull in relevant data from the problem as locals.
        const int& n_r = problem.get_observation_size(); // Number of observed data.
        const int& n_c = problem.get_parameter_size(); // Number of parameters.
        // Counters.
        const Incrementor iteration_counter = problem.get_iteration_counter();
        const Incrementor evaluation_counter = problem.get_evaluation_counter();
        // Convergence criterion.
        const Convergence_Checker<Evaluation> checker = problem.get_convergence_checker();

        // arrays shared with the other private methods
        const int solved_cols  = std::min(n_r, n_c);
        /* Parameters evolution direction associated with lm_par. */
        std::vector<double> lm_dir = std::vector<double>(n_c];
        /* Levenberg-_Marquardt parameter. */
        double lm_par = 0;

        // local point
        double   delta   = 0;
        double   x_norm   = 0;
        std::vector<double> diag    = std::vector<double>(n_c];
        std::vector<double> old_x    = std::vector<double>(n_c];
        std::vector<double> old_res  = std::vector<double>(n_r];
        std::vector<double> qtf     = std::vector<double>(n_r];
        std::vector<double> work1   = std::vector<double>(n_c];
        std::vector<double> work2   = std::vector<double>(n_c];
        std::vector<double> work3   = std::vector<double>(n_c];


        // Evaluate the function at the starting point and calculate its norm.
        evaluation_counter.increment();
        //value will be reassigned in the loop
        Evaluation current = problem.evaluate(problem.get_start());
        std::vector<double> current_residuals = current.get_residuals().to_array();
        double current_cost = current.get_cost();
        std::vector<double> current_point = current.get_point().to_array();

        // Outer loop.
        bool first_iteration = true;
        while (true) 
        {
            iteration_counter.increment();

            const Evaluation previous = current;

            // QR decomposition of the jacobian matrix
            const Internal_Data internal_data = qr_decomposition(current.get_jacobian(), solved_cols);
            const std::vector<std::vector<double>> weighted_jacobian = internal_data.weighted_jacobian;
            const std::vector<int> permutation = internal_data.permutation;
            const std::vector<double> diag_r = internal_data.diag_r;
            const std::vector<double> jac_norm = internal_data.jac_norm;

            //residuals already have weights applied
            std::vector<double> weighted_residual = current_residuals;
            for (int i{}; i < n_r; i++) 
            {
                qtf[i] = weighted_residual[i];
            }

            // compute Qt.res
            q_ty(qtf, internal_data);

            // now we don't need Q anymore, // so let jacobian contain the R matrix with its diagonal elements
            for (int k{}; k < solved_cols; ++k) 
            {
                int pk = permutation[k];
                weighted_jacobian[k][pk] = diag_r[pk];
            }

            if (first_iteration) 
            {
                // scale the point according to the norms of the columns
                // of the initial jacobian
                x_norm = 0;
                for (int k{}; k < n_c; ++k) 
                {
                    double dk = jac_norm[k];
                    if (dk == 0) 
                    {
                        dk = 1.0;
                    }
                    double xk = dk * current_point[k];
                    x_norm  += xk * xk;
                    diag[k] = dk;
                }
                x_norm = std::sqrt(x_norm);

                // initialize the step bound delta
                delta = (x_norm == 0) ? initial_step_bound_factor : (initial_step_bound_factor * x_norm);
            }

            // check orthogonality between function vector and jacobian columns
            double max_cosine = 0;
            if (current_cost != 0) 
            {
                for (int j{}; j < solved_cols; ++j) 
                {
                    int    pj = permutation[j];
                    double s  = jac_norm[pj];
                    if (s != 0) 
                    {
                        double sum{};
                        for (int i{}; i <= j; ++i) 
                        {
                            sum += weighted_jacobian[i][pj] * qtf[i];
                        }
                        max_cosine = std::max(max_cosine, std::abs(sum) / (s * current_cost));
                    }
                }
            }
            if (max_cosine <= ortho_tolerance) 
            {
                // Convergence has been reached.
                return Optimum.of(
                        current, evaluation_counter.get_count(), iteration_counter.get_count());
            }

            // rescale if necessary
            for (int j{}; j < n_c; ++j) 
            {
                diag[j] = std::max(diag[j], jac_norm[j]);
            }

            // Inner loop.
            for (double ratio = 0; ratio < 1.0e-4;) 
            {

                // save the state
                for (int j{}; j < solved_cols; ++j) 
                {
                    int pj = permutation[j];
                    old_x[pj] = current_point[pj];
                }
                const double previous_cost = current_cost;
                std::vector<double> tmp_vec = weighted_residual;
                weighted_residual = old_res;
                old_res    = tmp_vec;

                // determine the Levenberg-_Marquardt parameter
                lm_par = determine_l_m_parameter(qtf, delta, diag, internal_data, solved_cols, work1, work2, work3, lm_dir, lm_par);

                // compute the point and the norm of the evolution direction
                double lm_norm = 0;
                for (int j{}; j < solved_cols; ++j) 
                {
                    int pj = permutation[j];
                    lm_dir[pj] = -lm_dir[pj];
                    current_point[pj] = old_x[pj] + lm_dir[pj];
                    double s = diag[pj] * lm_dir[pj];
                    lm_norm  += s * s;
                }
                lm_norm = std::sqrt(lm_norm);
                // on the first iteration, adjust the initial step bound.
                if (first_iteration) 
                {
                    delta = std::min(delta, lm_norm);
                }

                // Evaluate the function at x + p and calculate its norm.
                evaluation_counter.increment();
                current = problem.evaluate(new Array_Real_Vector(current_point));
                current_residuals = current.get_residuals().to_array();
                current_cost = current.get_cost();
                current_point = current.get_point().to_array();

                // compute the scaled actual reduction
                double act_red = -1.0;
                if (0.1 * current_cost < previous_cost) 
                {
                    double r = current_cost / previous_cost;
                    act_red = 1.0 - r * r;
                }

                // compute the scaled predicted reduction
                // and the scaled directional derivative
                for (int j{}; j < solved_cols; ++j) 
                {
                    int pj = permutation[j];
                    double dir_j = lm_dir[pj];
                    work1[j] = 0;
                    for (int i{}; i <= j; ++i) 
                    {
                        work1[i] += weighted_jacobian[i][pj] * dir_j;
                    }
                }
                double coeff1 = 0;
                for (int j{}; j < solved_cols; ++j) 
                {
                    coeff1 += work1[j] * work1[j];
                }
                double pc2 = previous_cost * previous_cost;
                coeff1 /= pc2;
                double coeff2 = lm_par * lm_norm * lm_norm / pc2;
                double pre_red = coeff1 + 2 * coeff2;
                double dir_der = -(coeff1 + coeff2);

                // ratio of the actual to the predicted reduction
                ratio = (pre_red == 0) ? 0 : (act_red / pre_red);

                // update the step bound
                if (ratio <= 0.25) 
                {
                    double tmp =
                        (act_red < 0) ? (0.5 * dir_der / (dir_der + 0.5 * act_red)) : 0.5;
                        if ((0.1 * current_cost >= previous_cost) || (tmp < 0.1)) 
                        {
                            tmp = 0.1;
                        }
                        delta = tmp * std::min(delta, 10.0 * lm_norm);
                        lm_par /= tmp;
                }
else if ((lm_par == 0) || (ratio >= 0.75)) 
                {
                    delta = 2 * lm_norm;
                    lm_par *= 0.5;
                }

                // test for successful iteration.
                if (ratio >= 1.0e-4) 
                {
                    // successful iteration, update the norm
                    first_iteration = false;
                    x_norm = 0;
                    for (int k{}; k < n_c; ++k) 
                    {
                        double xK = diag[k] * current_point[k];
                        x_norm += xK * xK;
                    }
                    x_norm = std::sqrt(x_norm);

                    // tests for convergence.
                    if (checker != NULL && checker.converged(iteration_counter.get_count(), previous, current)) 
                    {
                        return Optimum.of(current, evaluation_counter.get_count(), iteration_counter.get_count());
                    }
                }
else 
                {
                    // failed iteration, reset the previous values
                    current_cost = previous_cost;
                    for (int j{}; j < solved_cols; ++j) 
                    {
                        int pj = permutation[j];
                        current_point[pj] = old_x[pj];
                    }
                    tmp_vec    = weighted_residual;
                    weighted_residual = old_res;
                    old_res    = tmp_vec;
                    // Reset "current" to previous values.
                    current = previous;
                }

                // Default convergence criteria.
                if ((std::abs(act_red) <= cost_relative_tolerance &&
                     pre_red <= cost_relative_tolerance &&
                     ratio <= 2.0) ||
                    delta <= par_relative_tolerance * x_norm) 
                    {
                    return Optimum.of(current, evaluation_counter.get_count(), iteration_counter.get_count());
                }

                // tests for termination and stringent tolerances
                if (std::abs(act_red) <= TWO_EPS &&
                    pre_red <= TWO_EPS &&
                    ratio <= 2.0) 
                    {
                    throw Math_Illegal_State_Exception(Localized_Optim_Formats.TOO_SMALL_COST_RELATIVE_TOLERANCE, cost_relative_tolerance);
                }
else if (delta <= TWO_EPS * x_norm) 
                {
                    throw Math_Illegal_State_Exception(Localized_Optim_Formats.TOO_SMALL_PARAMETERS_RELATIVE_TOLERANCE, par_relative_tolerance);
                }
else if (max_cosine <= TWO_EPS) 
                {
                    throw Math_Illegal_State_Exception(Localized_Optim_Formats.TOO_SMALL_ORTHOGONALITY_TOLERANCE, ortho_tolerance);
                }
            }
        }
    }

    /**
     * Holds internal data.
     * This structure was created so that all optimizer fields can be "const".
     * Code should be further refactored in order to not pass around arguments
     * that will modified in-place (cf. "work" arrays).
     */
    private static class Internal_Data 
    {
        /** Weighted Jacobian. */
        private const std::vector<std::vector<double>> weighted_jacobian;
        /** Columns permutation array. */
        private const std::vector<int> permutation;
        /** Rank of the Jacobian matrix. */
        private const int rank;
        /** Diagonal elements of the R matrix in the QR decomposition. */
        private const std::vector<double> diag_r;
        /** Norms of the columns of the jacobian matrix. */
        private const std::vector<double> jac_norm;
        /** Coefficients of the Householder transforms vectors. */
        private const std::vector<double> beta;

        /**
         * <p>
         * All arrays are stored by reference
         * </p>
         * @param weighted_jacobian Weighted Jacobian.
         * @param permutation Columns permutation array.
         * @param rank Rank of the Jacobian matrix.
         * @param diag_r Diagonal elements of the R matrix in the QR decomposition.
         * @param jac_norm Norms of the columns of the jacobian matrix.
         * @param beta Coefficients of the Householder transforms vectors.
         */
        Internal_Data(std::vector<std::vector<double>> weighted_jacobian, // NOPMD - staring array references is intentional and documented here
                     std::vector<int> permutation,           // NOPMD - staring array references is intentional and documented here
                     int rank, std::vector<double> diag_r,              // NOPMD - staring array references is intentional and documented here
                     std::vector<double> jac_norm,            // NOPMD - staring array references is intentional and documented here
                     const std::vector<double>& beta) {             // NOPMD - staring array references is intentional and documented here
            this.weighted_jacobian = weighted_jacobian;
            this.permutation = permutation;
            this.rank = rank;
            this.diag_r = diag_r;
            this.jac_norm = jac_norm;
            this.beta = beta;
        }
    }

    /**
     * Determines the Levenberg-_Marquardt parameter.
     *
     * <p>This implementation is a translation in Java of the MINPACK
     * <a href="http://www.netlib.org/minpack/lmpar.f">lmpar</a>
     * routine.</p>
     * <p>This method sets the lm_par and lm_dir attributes.</p>
     * <p>The authors of the original fortran function are:</p>
     * <ul>
     *   <li>Argonne National Laboratory. MINPACK project. March 1980</li>
     *   <li>Burton  S. Garbow</li>
     *   <li>Kenneth E. Hillstrom</li>
     *   <li>Jorge   J. More</li>
     * </ul>
     * <p>Luc Maisonobe did the Java translation.</p>
     *
     * @param qy Array containing q_ty.
     * @param delta Upper bound on the euclidean norm of diag_r * lm_dir.
     * @param diag Diagonal matrix.
     * @param internal_data Data (modified in-place in this method).
     * @param solved_cols Number of solved point.
     * @param work1 work array
     * @param work2 work array
     * @param work3 work array
     * @param lm_dir the "returned" LM direction will be stored in this array.
     * @param lm_par the value of the LM parameter from the previous iteration.
     * @return the LM parameter
     */
    private double determine_l_m_parameter(std::vector<double> qy, double delta, std::vector<double> diag, Internal_Data internal_data, int solved_cols, std::vector<double> work1, std::vector<double> work2, std::vector<double> work3, std::vector<double> lm_dir, double lm_par) 
    {
        const std::vector<std::vector<double>> weighted_jacobian = internal_data.weighted_jacobian;
        const std::vector<int> permutation = internal_data.permutation;
        const int rank = internal_data.rank;
        const std::vector<double> diag_r = internal_data.diag_r;

        const int& n_c = weighted_jacobian[0].size();

        // compute and store in x the gauss-newton direction, if the
        // jacobian is rank-deficient, obtain a least squares solution
        for (int j{}; j < rank; ++j) 
        {
            lm_dir[permutation[j]] = qy[j];
        }
        for (int j = rank; j < n_c; ++j) 
        {
            lm_dir[permutation[j]] = 0;
        }
        for (int k = rank - 1; k >= 0; --k) 
        {
            int pk = permutation[k];
            double ypk = lm_dir[pk] / diag_r[pk];
            for (int i{}; i < k; ++i) 
            {
                lm_dir[permutation[i]] -= ypk * weighted_jacobian[i][pk];
            }
            lm_dir[pk] = ypk;
        }

        // evaluate the function at the origin, and test
        // for acceptance of the Gauss-_Newton direction
        double dx_norm = 0;
        for (int j{}; j < solved_cols; ++j) 
        {
            int pj = permutation[j];
            double s = diag[pj] * lm_dir[pj];
            work1[pj] = s;
            dx_norm += s * s;
        }
        dx_norm = std::sqrt(dx_norm);
        double fp = dx_norm - delta;
        if (fp <= 0.1 * delta) 
        {
            lm_par = 0;
            return lm_par;
        }

        // if the jacobian is not rank deficient, the Newton step provides
        // a lower bound, parl, for the zero of the function, // otherwise set this bound to zero
        double sum2;
        double parl = 0;
        if (rank == solved_cols) 
        {
            for (int j{}; j < solved_cols; ++j) 
            {
                int pj = permutation[j];
                work1[pj] *= diag[pj] / dx_norm;
            }
            sum2 = 0;
            for (int j{}; j < solved_cols; ++j) 
            {
                int pj = permutation[j];
                double sum{};
                for (int i{}; i < j; ++i) 
                {
                    sum += weighted_jacobian[i][pj] * work1[permutation[i]];
                }
                double s = (work1[pj] - sum) / diag_r[pj];
                work1[pj] = s;
                sum2 += s * s;
            }
            parl = fp / (delta * sum2);
        }

        // calculate an upper bound, paru, for the zero of the function
        sum2 = 0;
        for (int j{}; j < solved_cols; ++j) 
        {
            int pj = permutation[j];
            double sum{};
            for (int i{}; i <= j; ++i) 
            {
                sum += weighted_jacobian[i][pj] * qy[i];
            }
            sum /= diag[pj];
            sum2 += sum * sum;
        }
        double g_norm = std::sqrt(sum2);
        double paru = g_norm / delta;
        if (paru == 0) 
        {
            paru = Precision.SAFE_MIN / std::min(delta, 0.1);
        }

        // if the input par lies outside of the interval (parl,paru), // set par to the closer endpoint
        lm_par = std::min(paru, std::max(lm_par, parl));
        if (lm_par == 0) 
        {
            lm_par = g_norm / dx_norm;
        }

        for (const int& countdown = 10; countdown >= 0; --countdown) 
        {

            // evaluate the function at the current value of lm_par
            if (lm_par == 0) 
            {
                lm_par = std::max(Precision.SAFE_MIN, 0.001 * paru);
            }
            double s_par = std::sqrt(lm_par);
            for (int j{}; j < solved_cols; ++j) 
            {
                int pj = permutation[j];
                work1[pj] = s_par * diag[pj];
            }
            determine_l_m_direction(qy, work1, work2, internal_data, solved_cols, work3, lm_dir);

            dx_norm = 0;
            for (int j{}; j < solved_cols; ++j) 
            {
                int pj = permutation[j];
                double s = diag[pj] * lm_dir[pj];
                work3[pj] = s;
                dx_norm += s * s;
            }
            dx_norm = std::sqrt(dx_norm);
            double previous_f_p = fp;
            fp = dx_norm - delta;

            // if the function is small enough, accept the current value
            // of lm_par, also test for the exceptional cases where parl is zero
            if (std::abs(fp) <= 0.1 * delta ||
                (parl == 0 &&
                 fp <= previous_f_p &&
                 previous_f_p < 0)) 
                 {
                return lm_par;
            }

            // compute the Newton correction
            for (int j{}; j < solved_cols; ++j) 
            {
                int pj = permutation[j];
                work1[pj] = work3[pj] * diag[pj] / dx_norm;
            }
            for (int j{}; j < solved_cols; ++j) 
            {
                int pj = permutation[j];
                work1[pj] /= work2[j];
                double tmp = work1[pj];
                for (int i = j + 1; i < solved_cols; ++i) 
                {
                    work1[permutation[i]] -= weighted_jacobian[i][pj] * tmp;
                }
            }
            sum2 = 0;
            for (int j{}; j < solved_cols; ++j) 
            {
                double s = work1[permutation[j]];
                sum2 += s * s;
            }
            double correction = fp / (delta * sum2);

            // depending on the sign of the function, update parl or paru.
            if (fp > 0) 
            {
                parl = std::max(parl, lm_par);
            }
else if (fp < 0) 
            {
                paru = std::min(paru, lm_par);
            }

            // compute an improved estimate for lm_par
            lm_par = std::max(parl, lm_par + correction);
        }

        return lm_par;
    }

    /**
     * Solve a*x = b and d*x = 0 in the least squares sense.
     * <p>This implementation is a translation in Java of the MINPACK
     * <a href="http://www.netlib.org/minpack/qrsolv.f">qrsolv</a>
     * routine.</p>
     * <p>This method sets the lm_dir and lm_diag attributes.</p>
     * <p>The authors of the original fortran function are:</p>
     * <ul>
     *   <li>Argonne National Laboratory. MINPACK project. March 1980</li>
     *   <li>Burton  S. Garbow</li>
     *   <li>Kenneth E. Hillstrom</li>
     *   <li>Jorge   J. More</li>
     * </ul>
     * <p>Luc Maisonobe did the Java translation.</p>
     *
     * @param qy array containing q_ty
     * @param diag diagonal matrix
     * @param lm_diag diagonal elements associated with lm_dir
     * @param internal_data Data (modified in-place in this method).
     * @param solved_cols Number of sloved point.
     * @param work work array
     * @param lm_dir the "returned" LM direction is stored in this array
     */
    private void determine_l_m_direction(std::vector<double> qy, std::vector<double> diag, std::vector<double> lm_diag, Internal_Data internal_data, int solved_cols, std::vector<double> work, std::vector<double> lm_dir) 
    {
        const std::vector<int> permutation = internal_data.permutation;
        const std::vector<std::vector<double>> weighted_jacobian = internal_data.weighted_jacobian;
        const std::vector<double> diag_r = internal_data.diag_r;

        // copy R and Qty to preserve input and initialize s
        //  in particular, save the diagonal elements of R in lm_dir
        for (int j{}; j < solved_cols; ++j) 
        {
            int pj = permutation[j];
            for (int i = j + 1; i < solved_cols; ++i) 
            {
                weighted_jacobian[i][pj] = weighted_jacobian[j][permutation[i]];
            }
            lm_dir[j] = diag_r[pj];
            work[j]  = qy[j];
        }

        // eliminate the diagonal matrix d using a Givens rotation
        for (int j{}; j < solved_cols; ++j) 
        {

            // prepare the row of d to be eliminated, locating the
            // diagonal element using p from the Q.R. factorization
            int pj = permutation[j];
            double dpj = diag[pj];
            if (dpj != 0) 
            {
                Arrays.fill(lm_diag, j + 1, lm_diag.size(), 0);
            }
            lm_diag[j] = dpj;

            //  the transformations to eliminate the row of d
            // modify only a single element of Qty
            // beyond the first n, which is initially zero.
            double qtbpj = 0;
            for (int k = j; k < solved_cols; ++k) 
            {
                int pk = permutation[k];

                // determine a Givens rotation which eliminates the
                // appropriate element in the current row of d
                if (lm_diag[k] != 0) 
                {

                    const double sin;
                    const double cos;
                    double rkk = weighted_jacobian[k][pk];
                    if (std::abs(rkk) < std::abs(lm_diag[k])) 
                    {
                        const double cotan = rkk / lm_diag[k];
                        sin   = 1.0 / std::sqrt(1.0 + cotan * cotan);
                        cos   = sin * cotan;
                    }
else 
                    {
                        const double tan = lm_diag[k] / rkk;
                        cos = 1.0 / std::sqrt(1.0 + tan * tan);
                        sin = cos * tan;
                    }

                    // compute the modified diagonal element of R and
                    // the modified element of (Qty,0)
                    weighted_jacobian[k][pk] = cos * rkk + sin * lm_diag[k];
                    const double temp = cos * work[k] + sin * qtbpj;
                    qtbpj = -sin * work[k] + cos * qtbpj;
                    work[k] = temp;

                    // accumulate the tranformation in the row of s
                    for (int i = k + 1; i < solved_cols; ++i) 
                    {
                        double rik = weighted_jacobian[i][pk];
                        const double temp2 = cos * rik + sin * lm_diag[i];
                        lm_diag[i] = -sin * rik + cos * lm_diag[i];
                        weighted_jacobian[i][pk] = temp2;
                    }
                }
            }

            // store the diagonal element of s and restore
            // the corresponding diagonal element of R
            lm_diag[j] = weighted_jacobian[j][permutation[j]];
            weighted_jacobian[j][permutation[j]] = lm_dir[j];
        }

        // solve the triangular system for z, if the system is
        // singular, then obtain a least squares solution
        int n_sing = solved_cols;
        for (int j{}; j < solved_cols; ++j) 
        {
            if ((lm_diag[j] == 0) && (n_sing == solved_cols)) 
            {
                n_sing = j;
            }
            if (n_sing < solved_cols) 
            {
                work[j] = 0;
            }
        }
        if (n_sing > 0) 
        {
            for (int j = n_sing - 1; j >= 0; --j) 
            {
                int pj = permutation[j];
                double sum{};
                for (int i = j + 1; i < n_sing; ++i) 
                {
                    sum += weighted_jacobian[i][pj] * work[i];
                }
                work[j] = (work[j] - sum) / lm_diag[j];
            }
        }

        // permute the components of z back to components of lm_dir
        for (int j{}; j < lm_dir.size(); ++j) 
        {
            lm_dir[permutation[j]] = work[j];
        }
    }

    /**
     * Decompose a matrix A as A.P = Q.R using Householder transforms.
     * <p>As suggested in the P. Lascaux and R. Theodor book
     * <i>Analyse num&eacute;rique matricielle appliqu&eacute;e &agrave;
     * l'art de l'ing&eacute;nieur</i> (Masson, 1986), instead of representing
     * the Householder transforms with u<sub>k</sub> unit vectors such that:
     * <pre>
     * H<sub>k</sub> = I - 2u<sub>k</sub>.u<sub>k</sub><sup>t</sup>
     * </pre>
     * we use <sub>k</sub> non-unit vectors such that:
     * <pre>
     * H<sub>k</sub> = I - beta<sub>k</sub>v<sub>k</sub>.v<sub>k</sub><sup>t</sup>
     * </pre>
     * where v<sub>k</sub> = a<sub>k</sub> - alpha<sub>k</sub> e<sub>k</sub>.
     * The beta<sub>k</sub> coefficients are provided upon exit as recomputing
     * them from the v<sub>k</sub> vectors would be costly.</p>
     * <p>This decomposition handles rank deficient cases since the tranformations
     * are performed in non-increasing columns norms order thanks to columns
     * pivoting. The diagonal elements of the R matrix are therefore also in
     * non-increasing absolute values order.</p>
     *
     * @param jacobian Weighted Jacobian matrix at the current point.
     * @param solved_cols Number of solved point.
     * @return data used in other methods of this class.
     * @Math_Illegal_State_Exception if the decomposition cannot be performed.
     */
    private Internal_Data qr_decomposition(Real_Matrix jacobian, int solved_cols)
        Math_Illegal_State_Exception 
        {
        // Code in this class assumes that the weighted Jacobian is -(W^(1/2) J), // hence the multiplication by -1.
        const std::vector<std::vector<double>> weighted_jacobian = jacobian.scalar_multiply(-1).get_data();

        const int& n_r = weighted_jacobian.size();
        const int& n_c = weighted_jacobian[0].size();

        const std::vector<int> permutation = int[n_c];
        const std::vector<double> diag_r = std::vector<double>(n_c];
        const std::vector<double> jac_norm = std::vector<double>(n_c];
        const auto beta = std::vector<double>std::vector<double>(n_c];

        // initializations
        for (int k{}; k < n_c; ++k) 
        {
            permutation[k] = k;
            double norm2 = 0;
            for (int i{}; i < n_r; ++i) 
            {
                double akk = weighted_jacobian[i][k];
                norm2 += akk * akk;
            }
            jac_norm[k] = std::sqrt(norm2);
        }

        // transform the matrix column after column
        for (int k{}; k < n_c; ++k) 
        {

            // select the column with the greatest norm on active components
            int next_column = -1;
            double ak2 = -INFINITY;
            for (int i = k; i < n_c; ++i) 
            {
                double norm2 = 0;
                for (int j{ k }; j < n_r; ++j) 
                {
                    double aki = weighted_jacobian[j][permutation[i]];
                    norm2 += aki * aki;
                }
                if (Double.std::isinfinite(norm2) || std::isnan(norm2)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Optim_Formats.UNABLE_TO_PERFORM_QR_DECOMPOSITION_ON_JACOBIAN, n_r, n_c);
                }
                if (norm2 > ak2) 
                {
                    next_column = i;
                    ak2        = norm2;
                }
            }
            if (ak2 <= qr_ranking_threshold) 
            {
                return Internal_Data(weighted_jacobian, permutation, k, diag_r, jac_norm, beta);
            }
            int pk = permutation[next_column];
            permutation[next_column] = permutation[k];
            permutation[k] = pk;

            // choose alpha such that Hk.u = alpha ek
            double akk = weighted_jacobian[k][pk];
            double alpha = (akk > 0) ? -std::sqrt(ak2) : std::sqrt(ak2);
            double betak = 1.0 / (ak2 - akk * alpha);
            beta[pk] = betak;

            // transform the current column
            diag_r[pk] = alpha;
            weighted_jacobian[k][pk] -= alpha;

            // transform the remaining columns
            for (const int& dk = n_c - 1 - k; dk > 0; --dk) 
            {
                double gamma = 0;
                for (int j{ k }; j < n_r; ++j) 
                {
                    gamma += weighted_jacobian[j][pk] * weighted_jacobian[j][permutation[k + dk]];
                }
                gamma *= betak;
                for (int j{ k }; j < n_r; ++j) 
                {
                    weighted_jacobian[j][permutation[k + dk]] -= gamma * weighted_jacobian[j][pk];
                }
            }
        }

        return Internal_Data(weighted_jacobian, permutation, solved_cols, diag_r, jac_norm, beta);
    }

    /**
     * Compute the product Qt.y for some Q.R. decomposition.
     *
     * @param y vector to multiply (will be overwritten with the result)
     * @param internal_data Data.
     */
    private void q_ty(std::vector<double> y, Internal_Data internal_data) 
    {
        const std::vector<std::vector<double>> weighted_jacobian = internal_data.weighted_jacobian;
        const std::vector<int> permutation = internal_data.permutation;
        const auto beta = std::vector<double>internal_data.beta;

        const int& n_r = weighted_jacobian.size();
        const int& n_c = weighted_jacobian[0].size();

        for (int k{}; k < n_c; ++k) 
        {
            int pk = permutation[k];
            double gamma = 0;
            for (int i = k; i < n_r; ++i) 
            {
                gamma += weighted_jacobian[i][pk] * y[i];
            }
            gamma *= beta[pk];
            for (int i = k; i < n_r; ++i) 
            {
                y[i] -= gamma * weighted_jacobian[i][pk];
            }
        }
    }
}


