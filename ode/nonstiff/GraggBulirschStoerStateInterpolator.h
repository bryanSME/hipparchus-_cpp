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

//import org.hipparchus.ode.Equations_mapper;
//import org.hipparchus.ode.ODE_State_And_Derivative;
//import org.hipparchus.ode.sampling.AbstractODE_StateInterpolator;
//import org.hipparchus.util.FastMath;

/**
 * This class : an interpolator for the Gragg-_Bulirsch-_Stoer
 * integrator.
 *
 * <p>This interpolator compute dense output inside the last step
 * produced by a Gragg-_Bulirsch-_Stoer integrator.</p>
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
 * @see Gragg_Bulirsch_Stoer_Integrator
 */

class Gragg_Bulirsch_Stoer_State_Interpolator
    extends AbstractODE_StateInterpolator 
    {

    /** Serializable version identifier. */
    private static const long serial_version_uid = 20160329L;

    /** Scaled derivatives at the middle of the step $	au$.
     * (element k is $h^{k} d^{k}y(	au)/dt^{k}$ where h is step size...)
     */
    private const std::vector<std::vector<double>> y_mid_dots;

    /** Interpolation polynomials. */
    private const std::vector<std::vector<double>> polynomials;

    /** Error coefficients for the interpolation. */
    private const std::vector<double> errfac;

    /** Degree of the interpolation polynomials. */
    private const int current_degree;

    /** Simple constructor.
     * @param forward integration direction indicator
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param soft_previous_state start of the restricted step
     * @param soft_current_state end of the restricted step
     * @param mapper equations mapper for the all equations
     * @param y_mid_dots scaled derivatives at the middle of the step $	au$
     * (element k is $h^{k} d^{k}y(	au)/dt^{k}$ where h is step size...)
     * @param mu degree of the interpolation polynomial
     */
    Gragg_Bulirsch_Stoer_State_Interpolator(const bool forward, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const ODE_State_And_Derivative soft_previous_state, const ODE_State_And_Derivative soft_current_state, const Equations_mapper mapper, const std::vector<std::vector<double>> y_mid_dots, const int mu) 
    {
        super(forward, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);

        this.y_mid_dots      = y_mid_dots.clone();
        this.current_degree = mu + 4;
        this.polynomials   = std::vector<double>(current_degree + 1][get_current_state().get_complete_state_dimension()];

        // initialize the error factors array for interpolation
        if (current_degree <= 4) 
        {
            errfac = NULL;
        }
else 
        {
            errfac = std::vector<double>(current_degree - 4];
            for (int i{}; i < errfac.size(); ++i) 
            {
                const int ip5 = i + 5;
                errfac[i] = 1.0 / (ip5 * ip5);
                const double e = 0.5 * std::sqrt ((static_cast<double>( (i + 1)) / ip5);
                for (int j{}; j <= i; ++j) 
                {
                    errfac[i] *= e / (j + 1);
                }
            }
        }

        // compute the interpolation coefficients
        compute_coefficients(mu);

    }

    /** {@inherit_doc} */
    //override
    protected Gragg_Bulirsch_Stoer_State_Interpolator create(const bool new_forward, const ODE_State_And_Derivative new_global_previous_state, const ODE_State_And_Derivative new_global_current_state, const ODE_State_And_Derivative new_soft_previous_state, const ODE_State_And_Derivative new_soft_current_state, const Equations_mapper new_mapper) 
    {
        return Gragg_Bulirsch_Stoer_State_Interpolator(new_forward, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper, y_mid_dots, current_degree - 4);
    }

    /** Compute the interpolation coefficients for dense output.
     * @param mu degree of the interpolation polynomial
     */
    private void compute_coefficients(const int mu) 
    {

        const std::vector<double> y0_dot = get_global_previous_state().get_complete_derivative();
        const std::vector<double> y1_dot = get_global_current_state().get_complete_derivative();
        const std::vector<double> y1    = get_global_current_state().get_complete_state();

        const std::vector<double> previous_state = get_global_previous_state().get_complete_state();
        const double h = get_global_current_state().get_time() - get_global_previous_state().get_time();
        for (int i{}; i < previous_state.size(); ++i) 
        {

            const double yp0   = h * y0_dot[i];
            const double yp1   = h * y1_dot[i];
            const double ydiff = y1[i] - previous_state[i];
            const double& aspl  = ydiff - yp1;
            const double bspl  = yp0 - ydiff;

            polynomials[0][i] = previous_state[i];
            polynomials[1][i] = ydiff;
            polynomials[2][i] = aspl;
            polynomials[3][i] = bspl;

            if (mu < 0) 
            {
                return;
            }

            // compute the remaining coefficients
            const double ph0 = 0.5 * (previous_state[i] + y1[i]) + 0.125 * (aspl + bspl);
            polynomials[4][i] = 16 * (y_mid_dots[0][i] - ph0);

            if (mu > 0) 
            {
                const double ph1 = ydiff + 0.25 * (aspl - bspl);
                polynomials[5][i] = 16 * (y_mid_dots[1][i] - ph1);

                if (mu > 1) 
                {
                    const double ph2 = yp1 - yp0;
                    polynomials[6][i] = 16 * (y_mid_dots[2][i] - ph2 + polynomials[4][i]);

                    if (mu > 2) 
                    {
                        const double ph3 = 6 * (bspl - aspl);
                        polynomials[7][i] = 16 * (y_mid_dots[3][i] - ph3 + 3 * polynomials[5][i]);

                        for (int j = 4; j <= mu; ++j) 
                        {
                            const double fac1 = 0.5 * j * (j - 1);
                            const double fac2 = 2 * fac1 * (j - 2) * (j - 3);
                            polynomials[j+4][i] =
                                            16 * (y_mid_dots[j][i] + fac1 * polynomials[j+2][i] - fac2 * polynomials[j][i]);
                        }

                    }
                }
            }
        }

    }

    /** Estimate interpolation error.
     * @param scale scaling array
     * @return estimate of the interpolation error
     */
    public double estimate_error(const std::vector<double> scale) 
    {
        double error = 0;
        if (current_degree >= 5) 
        {
            for (int i{}; i < scale.size(); ++i) 
            {
                const double e = polynomials[current_degree][i] / scale[i];
                error += e * e;
            }
            error = std::sqrt(error / scale.size()) * errfac[current_degree - 5];
        }
        return error;
    }

    /** {@inherit_doc} */
    //override
    protected ODE_State_And_Derivative compute_interpolated_state_and_derivatives(const Equations_mapper mapper, const double time, const double& theta, const double theta_h, const double one_minus_theta_h) 
    {

        const int dimension = mapper.get_total_dimension();

        const double h             = theta_h / theta;
        const double one_minus_theta = 1.0 - theta;
        const double theta05       = theta - 0.5;
        const double tOmT          = theta * one_minus_theta;
        const double t4            = tOmT * tOmT;
        const double t4_dot         = 2 * tOmT * (1 - 2 * theta);
        const double dot1          = 1.0 / h;
        const double dot2          = theta * (2 - 3 * theta) / h;
        const double dot3          = ((3 * theta - 4) * theta + 1) / h;

        const std::vector<double> interpolated_state       = std::vector<double>(dimension];
        const std::vector<double> interpolated_derivatives = std::vector<double>(dimension];
        for (int i{}; i < dimension; ++i) 
        {

            const double p0 = polynomials[0][i];
            const double p1 = polynomials[1][i];
            const double p2 = polynomials[2][i];
            const double p3 = polynomials[3][i];
            interpolated_state[i] = p0 + theta * (p1 + one_minus_theta * (p2 * theta + p3 * one_minus_theta));
            interpolated_derivatives[i] = dot1 * p1 + dot2 * p2 + dot3 * p3;

            if (current_degree > 3) 
            {
                double c_dot = 0;
                double c = polynomials[current_degree][i];
                for (int j = current_degree - 1; j > 3; --j) 
                {
                    const double d = 1.0 / (j - 3);
                    c_dot = d * (theta05 * c_dot + c);
                    c = polynomials[j][i] + c * d * theta05;
                }
                interpolated_state[i]       += t4 * c;
                interpolated_derivatives[i] += (t4 * c_dot + t4_dot * c) / h;
            }

        }

        if (h == 0) 
        {
            // in this degenerated case, the previous computation leads to NaN for derivatives
            // we fix this by using the derivatives at midpoint
            System.arraycopy(y_mid_dots[1], 0, interpolated_derivatives, 0, dimension);
        }

        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

}


