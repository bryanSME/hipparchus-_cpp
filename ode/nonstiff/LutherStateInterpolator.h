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
//import org.hipparchus.util.FastMath;

/**
 * This class represents an interpolator over the last step during an
 * ODE integration for the 6th order Luther integrator.
 *
 * <p>This interpolator computes dense output inside the last
 * step computed. The interpolation equation is consistent with the
 * integration scheme.</p>
 *
 * @see Luther_Integrator
 */

class Luther_State_Interpolator extends Runge_Kutta_State_Interpolator 
{

    /** Serializable version identifier */
    private static const long serial_version_uid = 20160328;

    /** Square root. */
    private static const double Q = std::sqrt(21);

    /** Simple constructor.
     * @param forward integration direction indicator
     * @param y_dot_k slopes at the intermediate points
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param soft_previous_state start of the restricted step
     * @param soft_current_state end of the restricted step
     * @param mapper equations mapper for the all equations
     */
    Luther_State_Interpolator(const bool forward, const std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const ODE_State_And_Derivative soft_previous_state, const ODE_State_And_Derivative soft_current_state, const Equations_mapper mapper) 
    {
        super(forward, y_dot_k, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
    }

    /** {@inherit_doc} */
    //override
    protected Luther_State_Interpolator create(const bool new_forward, const std::vector<std::vector<double>> new_y_dot_k, const ODE_State_And_Derivative new_global_previous_state, const ODE_State_And_Derivative new_global_current_state, const ODE_State_And_Derivative new_soft_previous_state, const ODE_State_And_Derivative new_soft_current_state, const Equations_mapper new_mapper) 
    {
        return Luther_State_Interpolator(new_forward, new_y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** {@inherit_doc} */
    //override
    protected ODE_State_And_Derivative compute_interpolated_state_and_derivatives(const Equations_mapper mapper, const double time, const double& theta, const double theta_h, const double one_minus_theta_h) 
    {

        // the coefficients below have been computed by solving the
        // order conditions from a theorem from Butcher (1963), using
        // the method explained in Folkmar Bornemann paper "Runge-Kutta
        // Methods, Trees, and Maple", Center of Mathematical Sciences, Munich
        // University of Technology, February 9, 2001
        //<http://wwwzenger.informatik.tu-muenchen.de/selcuk/sjam012101.html>

        // the method is implemented in the rkcheck tool
        // <https://www.spaceroots.org/software/rkcheck/index.html>.
        // Running it for order 5 gives the following order conditions
        // for an interpolator:
        // order 1 conditions
        // \sum_{i=1}^{i=s}\left(b_{i} \right) =1
        // order 2 conditions
        // \sum_{i=1}^{i=s}\left(b_{i} c_{i}\right) = \frac{	heta}{2}
        // order 3 conditions
        // \sum_{i=2}^{i=s}\left(b_{i} \sum_{j=1}^{j=i-1}{\left(a_{i,j} c_{j} \right)}\right) = \frac{	heta^{2}}{6}
        // \sum_{i=1}^{i=s}\left(b_{i} c_{i}^{2}\right) = \frac{	heta^{2}}{3}
        // order 4 conditions
        // \sum_{i=3}^{i=s}\left(b_{i} \sum_{j=2}^{j=i-1}{\left(a_{i,j} \sum_{k=1}^{k=j-1}{\left(a_{j,k} c_{k} \right)} \right)}\right) = \frac{	heta^{3}}{24}
        // \sum_{i=2}^{i=s}\left(b_{i} \sum_{j=1}^{j=i-1}{\left(a_{i,j} c_{j}^{2} \right)}\right) = \frac{	heta^{3}}{12}
        // \sum_{i=2}^{i=s}\left(b_{i} c_{i}\sum_{j=1}^{j=i-1}{\left(a_{i,j} c_{j} \right)}\right) = \frac{	heta^{3}}{8}
        // \sum_{i=1}^{i=s}\left(b_{i} c_{i}^{3}\right) = \frac{	heta^{3}}{4}
        // order 5 conditions
        // \sum_{i=4}^{i=s}\left(b_{i} \sum_{j=3}^{j=i-1}{\left(a_{i,j} \sum_{k=2}^{k=j-1}{\left(a_{j,k} \sum_{l=1}^{l=k-1}{\left(a_{k,l} c_{l} \right)} \right)} \right)}\right) = \frac{	heta^{4}}{120}
        // \sum_{i=3}^{i=s}\left(b_{i} \sum_{j=2}^{j=i-1}{\left(a_{i,j} \sum_{k=1}^{k=j-1}{\left(a_{j,k} c_{k}^{2} \right)} \right)}\right) = \frac{	heta^{4}}{60}
        // \sum_{i=3}^{i=s}\left(b_{i} \sum_{j=2}^{j=i-1}{\left(a_{i,j} c_{j}\sum_{k=1}^{k=j-1}{\left(a_{j,k} c_{k} \right)} \right)}\right) = \frac{	heta^{4}}{40}
        // \sum_{i=2}^{i=s}\left(b_{i} \sum_{j=1}^{j=i-1}{\left(a_{i,j} c_{j}^{3} \right)}\right) = \frac{	heta^{4}}{20}
        // \sum_{i=3}^{i=s}\left(b_{i} c_{i}\sum_{j=2}^{j=i-1}{\left(a_{i,j} \sum_{k=1}^{k=j-1}{\left(a_{j,k} c_{k} \right)} \right)}\right) = \frac{	heta^{4}}{30}
        // \sum_{i=2}^{i=s}\left(b_{i} c_{i}\sum_{j=1}^{j=i-1}{\left(a_{i,j} c_{j}^{2} \right)}\right) = \frac{	heta^{4}}{15}
        // \sum_{i=2}^{i=s}\left(b_{i} \left(\sum_{j=1}^{j=i-1}{\left(a_{i,j} c_{j} \right)} \right)^{2}\right) = \frac{	heta^{4}}{20}
        // \sum_{i=2}^{i=s}\left(b_{i} c_{i}^{2}\sum_{j=1}^{j=i-1}{\left(a_{i,j} c_{j} \right)}\right) = \frac{	heta^{4}}{10}
        // \sum_{i=1}^{i=s}\left(b_{i} c_{i}^{4}\right) = \frac{	heta^{4}}{5}

        // The a_{j,k} and c_{k} are given by the integrator Butcher arrays. What remains to solve
        // are the b_i for the interpolator. They are found by solving the above equations.
        // For a given interpolator, some equations are redundant, so in our case when we select
        // all equations from order 1 to 4, we still don't have enough independent equations
        // to solve from b_1 to b_7. We need to also select one equation from order 5. Here, // we selected the last equation. It appears this choice implied at least the last 3 equations
        // are fulfilled, but some of the former ones are not, so the resulting interpolator is order 5.
        // At the end, we get the b_i as polynomials in theta.

        const std::vector<double> interpolated_state;
        const std::vector<double> interpolated_derivatives;

        const double coeff_dot_1 =  1 + theta * ( -54            /   5.0 + theta * (   36                   + theta * ( -47                   + theta *   21)));
        const double coeff_dot_2 =  0;
        const double coeff_dot_3 =      theta * (-208            /  15.0 + theta * (  320            / 3.0  + theta * (-608            /  3.0 + theta *  112)));
        const double coeff_dot_4 =      theta * ( 324            /  25.0 + theta * ( -486            / 5.0  + theta * ( 972            /  5.0 + theta * -567           /  5.0)));
        const double coeff_dot5 =      theta * ((833 + 343 * Q) / 150.0 + theta * ((-637 - 357 * Q) / 30.0 + theta * ((392 + 287 * Q) / 15.0 + theta * (-49 - 49 * Q) /  5.0)));
        const double coeff_dot6 =      theta * ((833 - 343 * Q) / 150.0 + theta * ((-637 + 357 * Q) / 30.0 + theta * ((392 - 287 * Q) / 15.0 + theta * (-49 + 49 * Q) /  5.0)));
        const double coeff_dot7 =      theta * (   3            /   5.0 + theta * (   -3                   + theta *     3));

        if (get_global_previous_state() != NULL && theta <= 0.5) 
        {

            const double coeff1    =  1 + theta * ( -27            /   5.0 + theta * (   12                   + theta * ( -47            /  4.0 + theta *   21           /  5.0)));
            const double coeff2    =  0;
            const double coeff3    =      theta * (-104            /  15.0 + theta * (  320            / 9.0  + theta * (-152            /  3.0 + theta *  112           /  5.0)));
            const double coeff4    =      theta * ( 162            /  25.0 + theta * ( -162            / 5.0  + theta * ( 243            /  5.0 + theta * -567           / 25.0)));
            const double coeff5    =      theta * ((833 + 343 * Q) / 300.0 + theta * ((-637 - 357 * Q) / 90.0 + theta * ((392 + 287 * Q) / 60.0 + theta * (-49 - 49 * Q) / 25.0)));
            const double coeff6    =      theta * ((833 - 343 * Q) / 300.0 + theta * ((-637 + 357 * Q) / 90.0 + theta * ((392 - 287 * Q) / 60.0 + theta * (-49 + 49 * Q) / 25.0)));
            const double coeff7    =      theta * (   3            /  10.0 + theta * (   -1                   + theta * (   3            /  4.0)));
            interpolated_state       = previous_state_linear_combination(theta_h * coeff1, theta_h * coeff2, theta_h * coeff3, theta_h * coeff4, theta_h * coeff5, theta_h * coeff6, theta_h * coeff7);
            interpolated_derivatives = derivative_linear_combination(coeff_dot_1, coeff_dot_2, coeff_dot_3, coeff_dot_4, coeff_dot5, coeff_dot6, coeff_dot7);
        }
else 
        {

            const double coeff1    =  -1 /  20.0 + theta * (  19            /  20.0 + theta * (  -89             /  20.0  + theta * (   151            /  20.0 + theta *  -21           /   5.0)));
            const double coeff2    =  0;
            const double coeff3    = -16 /  45.0 + theta * ( -16            /  45.0 + theta * ( -328             /  45.0  + theta * (   424            /  15.0 + theta * -112           /   5.0)));
            const double coeff4    =               theta * (                          theta * (  162             /  25.0  + theta * (  -648            /  25.0 + theta *  567           /  25.0)));
            const double coeff5    = -49 / 180.0 + theta * ( -49            / 180.0 + theta * ((2254 + 1029 * Q) / 900.0  + theta * ((-1372 - 847 * Q) / 300.0 + theta * ( 49 + 49 * Q) /  25.0)));
            const double coeff6    = -49 / 180.0 + theta * ( -49            / 180.0 + theta * ((2254 - 1029 * Q) / 900.0  + theta * ((-1372 + 847 * Q) / 300.0 + theta * ( 49 - 49 * Q) /  25.0)));
            const double coeff7    =  -1 /  20.0 + theta * (  -1            /  20.0 + theta * (    1             /   4.0  + theta * (    -3            /   4.0)));
            interpolated_state       = current_state_linear_combination(one_minus_theta_h * coeff1, one_minus_theta_h * coeff2, one_minus_theta_h * coeff3, one_minus_theta_h * coeff4, one_minus_theta_h * coeff5, one_minus_theta_h * coeff6, one_minus_theta_h * coeff7);
            interpolated_derivatives = derivative_linear_combination(coeff_dot_1, coeff_dot_2, coeff_dot_3, coeff_dot_4, coeff_dot5, coeff_dot6, coeff_dot7);
        }

        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

}


