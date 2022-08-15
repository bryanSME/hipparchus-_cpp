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

//package org.hipparchus.ode.nonstiff;

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.ode.Expandable_ODE;
//import org.hipparchus.ode.Localized_ODE_Formats;
//import org.hipparchus.ode.ODE_State;
//import org.hipparchus.ode.ODE_State_And_Derivative;
//import org.hipparchus.util.FastMath;
#include "StepsizeHelper.h"
#include <string>
#include <vector>

/**
 * This class : a Gragg-_Bulirsch-_Stoer integrator for
 * Ordinary Differential Equations.
 *
 * <p>The Gragg-_Bulirsch-_Stoer algorithm is one of the most efficient
 * ones currently available for smooth problems. It uses Richardson
 * extrapolation to estimate what would be the solution if the step
 * size could be decreased down to zero.</p>
 *
 * <p>
 * This method changes both the step size and the order during
 * integration, in order to minimize computation cost. It is
 * particularly well suited when a very high precision is needed. The
 * limit where this method becomes more efficient than high-order
 * embedded Runge-Kutta methods like {@link Dormand_Prince853_Integrator
 * Dormand-Prince 8(5,3)} depends on the problem. Results given in the
 * Hairer, Norsett and Wanner book show for example that this limit
 * occurs for accuracy around 1e-6 when integrating Saltzam-Lorenz
 * equations (the authors note this problem is <i>extremely sensitive
 * to the errors in the first integration steps</i>), and around 1e-11
 * for a two dimensional celestial mechanics problems with seven
 * bodies (pleiades problem, involving quasi-collisions for which
 * <i>automatic step size control is essential</i>).
 * </p>
 *
 * <p>
 * This implementation is basically a reimplementation in Java of the
 * <a
 * href="http://www.unige.ch/math/folks/hairer/prog/nonstiff/odex.f">odex</a>
 * fortran code by E. Hairer and G. Wanner. The redistribution policy
 * for this code is available <a
 * href="http://www.unige.ch/~hairer/prog/licence.txt">here</a>, for
 * convenience, it is reproduced below.</p>
 * </p>
 *
 * <table border="0" width="80%" cellpadding="10" align="center" bgcolor="#E0E0E0">
 * <tr><td>Copyright (c) 2004, Ernst Hairer</td></tr>
 *
 * <tr><td>Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the following
 * conditions are met:
 * <ul>
 *  <li>Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.</li>
 *  <li>Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.</li>
 * </ul></td></tr>
 *
 * <tr><td><strong>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, * BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</strong></td></tr>
 * </table>
 *
 */

class Gragg_Bulirsch_Stoer_Integrator extends Adaptive_Stepsize_Integrator 
{

    /** Integrator method name. */
    private static const std::string METHOD_NAME = "Gragg-_Bulirsch-_Stoer";

    /** maximal order. */
    private int max_order;

    /** step size sequence. */
    private std::vector<int> sequence;

    /** overall cost of applying step reduction up to iteration k + 1, in number of calls. */
    private std::vector<int> cost_per_step;

    /** cost per unit step. */
    private std::vector<double> cost_per_time_unit;

    /** optimal steps for each order. */
    private std::vector<double> optimal_step;

    /** extrapolation coefficients. */
    private std::vector<std::vector<double>> coeff;

    /** stability check enabling parameter. */
    private bool perform_test;

    /** maximal number of checks for each iteration. */
    private int max_checks;

    /** maximal number of iterations for which checks are performed. */
    private int max_iter;

    /** stepsize reduction factor in case of stability check failure. */
    private double stability_reduction;

    /** first stepsize control factor. */
    private double step_control1;

    /** second stepsize control factor. */
    private double step_control2;

    /** third stepsize control factor. */
    private double step_control3;

    /** fourth stepsize control factor. */
    private double step_control4;

    /** first order control factor. */
    private double order_control_1;

    /** second order control factor. */
    private double order_control2;

    /** use interpolation error in stepsize control. */
    private bool use_interpolation_error;

    /** interpolation order control parameter. */
    private int mudif;

    /** Simple constructor.
     * Build a Gragg-_Bulirsch-_Stoer integrator with the given step
     * bounds. All tuning parameters are set to their default
     * values. The default step handler does nothing.
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scal_absolute_tolerance allowed absolute error
     * @param scal_relative_tolerance allowed relative error
     */
    public Gragg_Bulirsch_Stoer_Integrator(const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance) 
    {
        super(METHOD_NAME, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
        set_stability_check(true, -1, -1, -1);
        set_control_factors(-1, -1, -1, -1);
        set_order_control(-1, -1, -1);
        set_interpolation_control(true, -1);
    }

    /** Simple constructor.
     * Build a Gragg-_Bulirsch-_Stoer integrator with the given step
     * bounds. All tuning parameters are set to their default
     * values. The default step handler does nothing.
     * @param min_step minimal step (must be positive even for backward
     * integration), the last step can be smaller than this
     * @param max_step maximal step (must be positive even for backward
     * integration)
     * @param vec_absolute_tolerance allowed absolute error
     * @param vec_relative_tolerance allowed relative error
     */
    public Gragg_Bulirsch_Stoer_Integrator(const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance) 
    {
        super(METHOD_NAME, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
        set_stability_check(true, -1, -1, -1);
        set_control_factors(-1, -1, -1, -1);
        set_order_control(-1, -1, -1);
        set_interpolation_control(true, -1);
    }

    /** Set the stability check controls.
     * <p>The stability check is performed on the first few iterations of
     * the extrapolation scheme. If this test fails, the step is rejected
     * and the stepsize is reduced.</p>
     * <p>By default, the test is performed, at most during two
     * iterations at each step, and at most once for each of these
     * iterations. The default stepsize reduction factor is 0.5.</p>
     * @param perform_stability_check if true, stability check will be performed, if false, the check will be skipped
     * @param max_num_iter maximal number of iterations for which checks are
     * performed (the number of iterations is reset to default if negative
     * or NULL)
     * @param max_num_checks maximal number of checks for each iteration
     * (the number of checks is reset to default if negative or NULL)
     * @param stepsize_reduction_factor stepsize reduction factor in case of
     * failure (the factor is reset to default if lower than 0.0001 or
     * greater than 0.9999)
     */
    public void set_stability_check(const bool perform_stability_check, const int max_num_iter, const int max_num_checks, const double stepsize_reduction_factor) 
    {

        this.perform_test = perform_stability_check;
        this.max_iter     = (max_num_iter   <= 0) ? 2 : max_num_iter;
        this.max_checks   = (max_num_checks <= 0) ? 1 : max_num_checks;

        if ((stepsize_reduction_factor < 0.0001) || (stepsize_reduction_factor > 0.9999)) 
        {
            this.stability_reduction = 0.5;
        }
else 
        {
            this.stability_reduction = stepsize_reduction_factor;
        }

    }

    /** Set the step size control factors.

     * <p>The step size h_new is computed from the old one h by:
     * <pre>
     * h_new = h * step_control2 / (err/step_control1)^(1/(2k + 1))
     * </pre>
     * where err is the scaled error and k the iteration number of the
     * extrapolation scheme (counting from 0). The default values are
     * 0.65 for step_control1 and 0.94 for step_control2.</p>
     * <p>The step size is subject to the restriction:
     * <pre>
     * step_control3^(1/(2k + 1))/step_control4 &lt;= h_new/h &lt;= 1/step_control3^(1/(2k + 1))
     * </pre>
     * The default values are 0.02 for step_control3 and 4.0 for
     * step_control4.</p>
     * @param control1 first stepsize control factor (the factor is
     * reset to default if lower than 0.0001 or greater than 0.9999)
     * @param control2 second stepsize control factor (the factor
     * is reset to default if lower than 0.0001 or greater than 0.9999)
     * @param control3 third stepsize control factor (the factor is
     * reset to default if lower than 0.0001 or greater than 0.9999)
     * @param control4 fourth stepsize control factor (the factor
     * is reset to default if lower than 1.0001 or greater than 999.9)
     */
    public void set_control_factors(const double control1, const double control2, const double control3, const double control4) 
    {

        if ((control1 < 0.0001) || (control1 > 0.9999)) 
        {
            this.step_control1 = 0.65;
        }
else 
        {
            this.step_control1 = control1;
        }

        if ((control2 < 0.0001) || (control2 > 0.9999)) 
        {
            this.step_control2 = 0.94;
        }
else 
        {
            this.step_control2 = control2;
        }

        if ((control3 < 0.0001) || (control3 > 0.9999)) 
        {
            this.step_control3 = 0.02;
        }
else 
        {
            this.step_control3 = control3;
        }

        if ((control4 < 1.0001) || (control4 > 999.9)) 
        {
            this.step_control4 = 4.0;
        }
else 
        {
            this.step_control4 = control4;
        }

    }

    /** Set the order control parameters.
     * <p>The Gragg-_Bulirsch-_Stoer method changes both the step size and
     * the order during integration, in order to minimize computation
     * cost. Each extrapolation step increases the order by 2, so the
     * maximal order that will be used is always even, it is twice the
     * maximal number of columns in the extrapolation table.</p>
     * <pre>
     * order is decreased if w(k - 1) &lt;= w(k)     * order_control_1
     * order is increased if w(k)     &lt;= w(k - 1) * order_control2
     * </pre>
     * <p>where w is the table of work per unit step for each order
     * (number of function calls divided by the step length), and k is
     * the current order.</p>
     * <p>The default maximal order after construction is 18 (i.e. the
     * maximal number of columns is 9). The default values are 0.8 for
     * order_control_1 and 0.9 for order_control2.</p>
     * @param maximal_order maximal order in the extrapolation table (the
     * maximal order is reset to default if order &lt;= 6 or odd)
     * @param control1 first order control factor (the factor is
     * reset to default if lower than 0.0001 or greater than 0.9999)
     * @param control2 second order control factor (the factor
     * is reset to default if lower than 0.0001 or greater than 0.9999)
     */
    public void set_order_control(const int maximal_order, const double control1, const double control2) 
    {

        if (maximal_order > 6 && maximal_order % 2 == 0) 
        {
            this.max_order = maximal_order;
        }
else 
        {
            this.max_order = 18;
        }

        if ((control1 < 0.0001) || (control1 > 0.9999)) 
        {
            this.order_control_1 = 0.8;
        }
else 
        {
            this.order_control_1 = control1;
        }

        if ((control2 < 0.0001) || (control2 > 0.9999)) 
        {
            this.order_control2 = 0.9;
        }
else 
        {
            this.order_control2 = control2;
        }

        // reinitialize the arrays
        initialize_arrays();

    }

    /** Initialize the integrator internal arrays. */
    private void initialize_arrays() 
    {

        const int size = max_order / 2;

        if ((sequence == NULL) || (sequence.size() != size)) 
        {
            // all arrays should be reallocated with the right size
            sequence        = int[size];
            cost_per_step     = int[size];
            coeff           = std::vector<double>(size][];
            cost_per_time_unit = std::vector<double>(size];
            optimal_step     = std::vector<double>(size];
        }

        // step size sequence: 2, 6, 10, 14, ...
        for (int k{}; k < size; ++k) 
        {
            sequence[k] = 4 * k + 2;
        }

        // initialize the order selection cost array
        // (number of function calls for each column of the extrapolation table)
        cost_per_step[0] = sequence[0] + 1;
        for (int k{ 1 }; k < size; ++k) 
        {
            cost_per_step[k] = cost_per_step[k - 1] + sequence[k];
        }

        // initialize the extrapolation tables
        for (int k{}; k < size; ++k) 
        {
            coeff[k] = (k > 0) ? std::vector<double>(k] : NULL;
            for (const int& l = 0; l < k; ++l) 
            {
                const double ratio = (static_cast<double>( sequence[k]) / sequence[k - l - 1];
                coeff[k][l] = 1.0 / (ratio * ratio - 1.0);
            }
        }

    }

    /** Set the interpolation order control parameter.
     * The interpolation order for dense output is 2k - mudif + 1. The
     * default value for mudif is 4 and the interpolation error is used
     * in stepsize control by default.

     * @param use_interpolation_error_for_control if true, interpolation error is used
     * for stepsize control
     * @param mudif_control_parameter interpolation order control parameter (the parameter
     * is reset to default if &lt;= 0 or &gt;= 7)
     */
    public void set_interpolation_control(const bool use_interpolation_error_for_control, const int mudif_control_parameter) 
    {

        this.use_interpolation_error = use_interpolation_error_for_control;

        if ((mudif_control_parameter <= 0) || (mudif_control_parameter >= 7)) 
        {
            this.mudif = 4;
        }
else 
        {
            this.mudif = mudif_control_parameter;
        }

    }

    /** Update scaling array.
     * @param y1 first state vector to use for scaling
     * @param y2 second state vector to use for scaling
     * @param scale scaling array to update (can be shorter than state)
     */
    private void rescale(const std::vector<double> y1, const std::vector<double> y2, const std::vector<double> scale) 
    {
        const Stepsize_Helper helper = get_step_size_helper();
        for (int i{}; i < scale.size(); ++i) 
        {
            scale[i] = helper.get_tolerance(i, std::max(std::abs(y1[i]), std::abs(y2[i])));
        }
    }

    /** Perform integration over one step using substeps of a modified
     * midpoint method.
     * @param t0 initial time
     * @param y0 initial value of the state vector at t0
     * @param step global step
     * @param k iteration number (from 0 to sequence.size() - 1)
     * @param scale scaling array (can be shorter than state)
     * @param f placeholder where to put the state vector derivatives at each substep
     *          (element 0 already contains initial derivative)
     * @param y_middle placeholder where to put the state vector at the middle of the step
     * @param y_end placeholder where to put the state vector at the end
     * @return true if computation was done properly, *         false if stability check failed before end of computation
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if arrays dimensions do not match equations settings
     */
    private bool try_step(const double t0, const std::vector<double> y0, const double step, const int& k, const std::vector<double> scale, const std::vector<std::vector<double>> f, const std::vector<double> y_middle, const std::vector<double> y_end)
        , Math_Illegal_State_Exception 
        {

        const int    n        = sequence[k];
        const double sub_step  = step / n;
        const double sub_step2 = 2 * sub_step;

        // first substep
        double t = t0 + sub_step;
        for (int i{}; i < y0.size(); ++i) 
        {
            y_end[i] = y0[i] + sub_step * f[0][i];
        }
        f[1] = compute_derivatives(t, y_end);

        // other substeps
        const std::vector<double> y_tmp = y0.clone();
        for (int j{ 1 }; j < n; ++j) 
        {

            if (2 * j == n) 
            {
                // save the point at the middle of the step
                System.arraycopy(y_end, 0, y_middle, 0, y0.size());
            }

            t += sub_step;
            for (int i{}; i < y0.size(); ++i) 
            {
                const double middle = y_end[i];
                y_end[i]       = y_tmp[i] + sub_step2 * f[j][i];
                y_tmp[i]       = middle;
            }

            f[j + 1] = compute_derivatives(t, y_end);

            // stability check
            if (perform_test && (j <= max_checks) && (k < max_iter)) 
            {
                double initial_norm = 0.0;
                for (const int& l = 0; l < scale.size(); ++l) 
                {
                    const double ratio = f[0][l] / scale[l];
                    initial_norm += ratio * ratio;
                }
                double delta_norm = 0.0;
                for (const int& l = 0; l < scale.size(); ++l) 
                {
                    const double ratio = (f[j + 1][l] - f[0][l]) / scale[l];
                    delta_norm += ratio * ratio;
                }
                if (delta_norm > 4 * std::max(1.0e-15, initial_norm)) 
                {
                    return false;
                }
            }

        }

        // correction of the last substep (at t0 + step)
        for (int i{}; i < y0.size(); ++i) 
        {
            y_end[i] = 0.5 * (y_tmp[i] + y_end[i] + sub_step * f[n][i]);
        }

        return true;

    }

    /** Extrapolate a vector.
     * @param offset offset to use in the coefficients table
     * @param k index of the last updated point
     * @param diag working diagonal of the Aitken-Neville's
     * triangle, without the last element
     * @param last last element
     */
    private void extrapolate(const int offset, const int& k, const std::vector<std::vector<double>> diag, const std::vector<double> last) 
    {

        // update the diagonal
        for (int j{ 1 }; j < k; ++j) 
        {
            for (int i{}; i < last.size(); ++i) 
            {
                // Aitken-Neville's recursive formula
                diag[k - j - 1][i] = diag[k - j][i] +
                                coeff[k + offset][j - 1] * (diag[k - j][i] - diag[k - j - 1][i]);
            }
        }

        // update the last element
        for (int i{}; i < last.size(); ++i) 
        {
            // Aitken-Neville's recursive formula
            last[i] = diag[0][i] + coeff[k + offset][k - 1] * (diag[0][i] - last[i]);
        }
    }

    /** {@inherit_doc} */
    //override
    public ODE_State_And_Derivative integrate(const Expandable_ODE equations, const ODE_State initial_state, const double const_time)
        , Math_Illegal_State_Exception 
        {

        sanity_checks(initial_state, const_time);
        set_step_start(init_integration(equations, initial_state, const_time));
        const bool forward = const_time > initial_state.get_time();

        // create some internal working arrays
        std::vector<double>         y        = get_step_start().get_complete_state();
        const std::vector<double>   y1       = std::vector<double>(y.size()];
        const std::vector<std::vector<double>> diagonal = std::vector<double>(sequence.size() - 1][];
        const std::vector<std::vector<double>> y1_diag   = std::vector<double>(sequence.size() - 1][];
        for (int k{}; k < sequence.size() - 1; ++k) 
        {
            diagonal[k] = std::vector<double>(y.size()];
            y1_diag[k]   = std::vector<double>(y.size()];
        }

        const std::vector<std::vector<double>>[] fk = std::vector<double>(sequence.size()][][];
        for (int k{}; k < sequence.size(); ++k) 
        {
            fk[k] = std::vector<double>(sequence[k] + 1][];
        }

        // scaled derivatives at the middle of the step $	au$
        // (element k is $h^{k} d^{k}y(	au)/dt^{k}$ where h is step size...)
        const std::vector<std::vector<double>> y_mid_dots = std::vector<double>(1 + 2 * sequence.size()][y.size()];

        // initial scaling
        const int main_set_dimension = get_step_size_helper().get_main_set_dimension();
        const std::vector<double> scale = std::vector<double>(main_set_dimension];
        rescale(y, y, scale);

        // initial order selection
        const double tol    = get_step_size_helper().get_relative_tolerance(0);
        const double log10_r = std::log10(std::max(1.0e-10, tol));
        int target_iter = std::max(1, std::min(sequence.size() - 2, static_cast<int>( std::floor(0.5 - 0.6 * log10_r)));

        double  h_new                     = 0;
        double  max_error                 = Double.MAX_VALUE;
        bool previous_rejected         = false;
        bool first_time                = true;
        bool new_step                  = true;
        cost_per_time_unit[0] = 0;
        set_is_last_step(false);
        do 
        {

            double error;
            bool reject = false;

            if (new_step) 
            {

                // first evaluation, at the beginning of the step
                const std::vector<double> y_dot_0 = get_step_start().get_complete_derivative();
                for (int k{}; k < sequence.size(); ++k) 
                {
                    // all sequences start from the same point, so we share the derivatives
                    fk[k][0] = y_dot_0;
                }

                if (first_time) 
                {
                    h_new = initialize_step(forward, 2 * target_iter + 1, scale, get_step_start(), equations.get_mapper());
                }

                new_step = false;

            }

            set_step_size(h_new);

            // step adjustment near bounds
            if (forward) 
            {
                if (get_step_start().get_time() + get_step_size() >= const_time) 
                {
                    set_step_size(const_time - get_step_start().get_time());
                }
            }
else 
            {
                if (get_step_start().get_time() + get_step_size() <= const_time) 
                {
                    set_step_size(const_time - get_step_start().get_time());
                }
            }
            const double next_t = get_step_start().get_time() + get_step_size();
            set_is_last_step(forward ? (next_t >= const_time) : (next_t <= const_time));

            // iterate over several substep sizes
            int k = -1;
            for (bool loop = true; loop; ) 
            {

                ++k;

                // modified midpoint integration with the current substep
                if ( ! try_step(get_step_start().get_time(), y, get_step_size(), k, scale, fk[k], (k == 0) ? y_mid_dots[0] : diagonal[k - 1], (k == 0) ? y1 : y1_diag[k - 1])) 
                {

                    // the stability check failed, we reduce the global step
                    h_new   = std::abs(get_step_size_helper().filter_step(get_step_size() * stability_reduction, forward, false));
                    reject = true;
                    loop   = false;

                }
else 
                {

                    // the substep was computed successfully
                    if (k > 0) 
                    {

                        // extrapolate the state at the end of the step
                        // using last iteration data
                        extrapolate(0, k, y1_diag, y1);
                        rescale(y, y1, scale);

                        // estimate the error at the end of the step.
                        error = 0;
                        for (int j{}; j < main_set_dimension; ++j) 
                        {
                            const double e = std::abs(y1[j] - y1_diag[0][j]) / scale[j];
                            error += e * e;
                        }
                        error = std::sqrt(error / main_set_dimension);
                        if (std::isnan(error)) 
                        {
                            throw Math_Illegal_State_Exception(Localized_ODE_Formats.NAN_APPEARING_DURING_INTEGRATION, next_t);
                        }

                        if ((error > 1.0e15) || ((k > 1) && (error > max_error))) 
                        {
                            // error is too big, we reduce the global step
                            h_new   = std::abs(get_step_size_helper().filter_step(get_step_size() * stability_reduction, forward, false));
                            reject = true;
                            loop   = false;
                        }
else 
                        {

                            max_error = std::max(4 * error, 1.0);

                            // compute optimal stepsize for this order
                            const double exp = 1.0 / (2 * k + 1);
                            double fac = step_control2 / std::pow(error / step_control1, exp);
                            const double pow = std::pow(step_control3, exp);
                            fac = std::max(pow / step_control4, std::min(1 / pow, fac));
                            const bool accept_small = k < target_iter;
                            optimal_step[k]     = std::abs(get_step_size_helper().filter_step(get_step_size() * fac, forward, accept_small));
                            cost_per_time_unit[k] = cost_per_step[k] / optimal_step[k];

                            // check convergence
                            switch (k - target_iter) 
                            {

                                case -1 :
                                    if ((target_iter > 1) && ! previous_rejected) 
                                    {

                                        // check if we can stop iterations now
                                        if (error <= 1.0) 
                                        {
                                            // convergence have been reached just before target_iter
                                            loop = false;
                                        }
else 
                                        {
                                            // estimate if there is a chance convergence will
                                            // be reached on next iteration, using the
                                            // asymptotic evolution of error
                                            const double ratio = (static_cast<double>( sequence [target_iter] * sequence[target_iter + 1]) /
                                                            (sequence[0] * sequence[0]);
                                            if (error > ratio * ratio) 
                                            {
                                                // we don't expect to converge on next iteration
                                                // we reject the step immediately and reduce order
                                                reject = true;
                                                loop   = false;
                                                target_iter = k;
                                                if ((target_iter > 1) &&
                                                    (cost_per_time_unit[target_iter - 1] <
                                                                    order_control_1 * cost_per_time_unit[target_iter])) 
                                                                    {
                                                    --target_iter;
                                                }
                                                h_new = get_step_size_helper().filter_step(optimal_step[target_iter], forward, false);
                                            }
                                        }
                                    }
                                    break;

                                case 0:
                                    if (error <= 1.0) 
                                    {
                                        // convergence has been reached exactly at target_iter
                                        loop = false;
                                    }
else 
                                    {
                                        // estimate if there is a chance convergence will
                                        // be reached on next iteration, using the
                                        // asymptotic evolution of error
                                        const double ratio = (static_cast<double>( sequence[k + 1]) / sequence[0];
                                        if (error > ratio * ratio) 
                                        {
                                            // we don't expect to converge on next iteration
                                            // we reject the step immediately
                                            reject = true;
                                            loop = false;
                                            if ((target_iter > 1) &&
                                                 (cost_per_time_unit[target_iter - 1] <
                                                                 order_control_1 * cost_per_time_unit[target_iter])) 
                                                                 {
                                                --target_iter;
                                            }
                                            h_new = get_step_size_helper().filter_step(optimal_step[target_iter], forward, false);
                                        }
                                    }
                                    break;

                                case 1 :
                                    if (error > 1.0) 
                                    {
                                        reject = true;
                                        if ((target_iter > 1) &&
                                            (cost_per_time_unit[target_iter - 1] <
                                                            order_control_1 * cost_per_time_unit[target_iter])) 
                                                            {
                                            --target_iter;
                                        }
                                        h_new = get_step_size_helper().filter_step(optimal_step[target_iter], forward, false);
                                    }
                                    loop = false;
                                    break;

                                default :
                                    if ((first_time || is_last_step()) && (error <= 1.0)) 
                                    {
                                        loop = false;
                                    }
                                    break;

                            }

                        }
                    }
                }
            }

            // dense output handling
            double h_int = get_max_step();
            const Gragg_Bulirsch_Stoer_State_Interpolator interpolator;
            if (! reject) 
            {

                // extrapolate state at middle point of the step
                for (int j{ 1 }; j <= k; ++j) 
                {
                    extrapolate(0, j, diagonal, y_mid_dots[0]);
                }

                const int mu = 2 * k - mudif + 3;

                for (const int& l = 0; l < mu; ++l) 
                {

                    // derivative at middle point of the step
                    const int l2 = l / 2;
                    double factor = std::pow(0.5 * sequence[l2], l);
                    int middle_index = fk[l2].size() / 2;
                    for (int i{}; i < y.size(); ++i) 
                    {
                        y_mid_dots[l + 1][i] = factor * fk[l2][middle_index + l][i];
                    }
                    for (int j{ 1 }; j <= k - l2; ++j) 
                    {
                        factor = std::pow(0.5 * sequence[j + l2], l);
                        middle_index = fk[l2 + j].size() / 2;
                        for (int i{}; i < y.size(); ++i) 
                        {
                            diagonal[j - 1][i] = factor * fk[l2 + j][middle_index + l][i];
                        }
                        extrapolate(l2, j, diagonal, y_mid_dots[l + 1]);
                    }
                    for (int i{}; i < y.size(); ++i) 
                    {
                        y_mid_dots[l + 1][i] *= get_step_size();
                    }

                    // compute centered differences to evaluate next derivatives
                    for (int j = (l + 1) / 2; j <= k; ++j) 
                    {
                        for (const int& m = fk[j].size() - 1; m >= 2 * (l + 1); --m) 
                        {
                            for (int i{}; i < y.size(); ++i) 
                            {
                                fk[j][m][i] -= fk[j][m - 2][i];
                            }
                        }
                    }

                }

                // state at end of step
                const ODE_State_And_Derivative step_end =
                    equations.get_mapper().map_state_and_derivative(next_t, y1, compute_derivatives(next_t, y1));

                // set up interpolator covering the full step
                interpolator = Gragg_Bulirsch_Stoer_State_Interpolator(forward, get_step_start(), step_end, get_step_start(), step_end, equations.get_mapper(), y_mid_dots, mu);

                if (mu >= 0 && use_interpolation_error) 
                {
                    // use the interpolation error to limit stepsize
                    const double interp_error = interpolator.estimate_error(scale);
                    h_int = std::abs(get_step_size() /
                                        std::max(std::pow(interp_error, 1.0 / (mu + 4)), 0.01));
                    if (interp_error > 10.0) 
                    {
                        h_new   = get_step_size_helper().filter_step(h_int, forward, false);
                        reject = true;
                    }
                }

            }
else 
            {
                interpolator = NULL;
            }

            if (! reject) 
            {

                // Discrete events handling
                set_step_start(accept_step(interpolator, const_time));

                // prepare next step
                // beware that y1 is not always valid anymore here, // as some event may have triggered a reset
                // so we need to copy the step start set previously
                y = get_step_start().get_complete_state();

                int optimal_iter;
                if (k == 1) 
                {
                    optimal_iter = 2;
                    if (previous_rejected) 
                    {
                        optimal_iter = 1;
                    }
                }
else if (k <= target_iter) 
                {
                    optimal_iter = k;
                    if (cost_per_time_unit[k - 1] < order_control_1 * cost_per_time_unit[k]) 
                    {
                        optimal_iter = k - 1;
                    }
else if (cost_per_time_unit[k] < order_control2 * cost_per_time_unit[k - 1]) 
                    {
                        optimal_iter = std::min(k + 1, sequence.size() - 2);
                    }
                }
else 
                {
                    optimal_iter = k - 1;
                    if ((k > 2) && (cost_per_time_unit[k - 2] < order_control_1 * cost_per_time_unit[k - 1])) 
                    {
                        optimal_iter = k - 2;
                    }
                    if (cost_per_time_unit[k] < order_control2 * cost_per_time_unit[optimal_iter]) 
                    {
                        optimal_iter = std::min(k, sequence.size() - 2);
                    }
                }

                if (previous_rejected) 
                {
                    // after a rejected step neither order nor stepsize
                    // should increase
                    target_iter = std::min(optimal_iter, k);
                    h_new = std::min(std::abs(get_step_size()), optimal_step[target_iter]);
                }
else 
                {
                    // stepsize control
                    if (optimal_iter <= k) 
                    {
                        h_new = get_step_size_helper().filter_step(optimal_step[optimal_iter], forward, false);
                    }
else 
                    {
                        if ((k < target_iter) &&
                                        (cost_per_time_unit[k] < order_control2 * cost_per_time_unit[k - 1])) 
                                        {
                            h_new = get_step_size_helper().
                                   filter_step(optimal_step[k] * cost_per_step[optimal_iter + 1] / cost_per_step[k], forward, false);
                        }
else 
                        {
                            h_new = get_step_size_helper().
                                   filter_step(optimal_step[k] * cost_per_step[optimal_iter] / cost_per_step[k], forward, false);
                        }
                    }

                    target_iter = optimal_iter;

                }

                new_step = true;

            }

            h_new = std::min(h_new, h_int);
            if (! forward) 
            {
                h_new = -h_new;
            }

            first_time = false;

            if (reject) 
            {
                set_is_last_step(false);
                previous_rejected = true;
            }
else 
            {
                previous_rejected = false;
            }

        } while (!is_last_step());

        const ODE_State_And_Derivative const_state = get_step_start();
        reset_internal_state();
        return const_state;

    }

}


