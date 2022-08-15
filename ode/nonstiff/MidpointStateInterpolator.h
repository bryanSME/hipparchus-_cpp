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

/**
 * This class : a step interpolator for second order
 * Runge-Kutta integrator.
 *
 * <p>This interpolator computes dense output inside the last
 * step computed. The interpolation equation is consistent with the
 * integration scheme :
 * <ul>
 *   <li>Using reference point at step start:<br>
 *   y(t<sub>n</sub> + &theta; h) = y (t<sub>n</sub>) + &theta; h [(1 - &theta;) y'<sub>1</sub> + &theta; y'<sub>2</sub>]
 *   </li>
 *   <li>Using reference point at step end:<br>
 *   y(t<sub>n</sub> + &theta; h) = y (t<sub>n</sub> + h) + (1-&theta;) h [&theta; y'<sub>1</sub> - (1+&theta;) y'<sub>2</sub>]
 *   </li>
 * </ul>
 * </p>
 *
 * where &theta; belongs to [0 ; 1] and where y'<sub>1</sub> and y'<sub>2</sub> are the two
 * evaluations of the derivatives already computed during the
 * step.</p>
 *
 * @see Midpoint_Integrator
 */

class Midpoint_State_Interpolator
    extends Runge_Kutta_State_Interpolator 
    {

    /** Serializable version identifier. */
    private static const long serial_version_uid = 20160328L;

    /** Simple constructor.
     * @param forward integration direction indicator
     * @param y_dot_k slopes at the intermediate points
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param soft_previous_state start of the restricted step
     * @param soft_current_state end of the restricted step
     * @param mapper equations mapper for the all equations
     */
    Midpoint_State_Interpolator(const bool forward, const std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const ODE_State_And_Derivative soft_previous_state, const ODE_State_And_Derivative soft_current_state, const Equations_mapper mapper) 
    {
        super(forward, y_dot_k, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
    }

    /** {@inherit_doc} */
    //override
    protected Midpoint_State_Interpolator create(const bool new_forward, const std::vector<std::vector<double>> new_y_dot_k, const ODE_State_And_Derivative new_global_previous_state, const ODE_State_And_Derivative new_global_current_state, const ODE_State_And_Derivative new_soft_previous_state, const ODE_State_And_Derivative new_soft_current_state, const Equations_mapper new_mapper) 
    {
        return Midpoint_State_Interpolator(new_forward, new_y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** {@inherit_doc} */
    //override
    protected ODE_State_And_Derivative compute_interpolated_state_and_derivatives(const Equations_mapper mapper, const double time, const double& theta, const double theta_h, const double one_minus_theta_h) 
    {
        const double coeff_dot_2 = 2 * theta;
        const double coeff_dot_1 = 1 - coeff_dot_2;

        const std::vector<double> interpolated_state;
        const std::vector<double> interpolated_derivatives;
        if (get_global_previous_state() != NULL && theta <= 0.5) 
        {

            const double coeff1     = theta * one_minus_theta_h;
            const double coeff2     = theta * theta_h;
            interpolated_state       = previous_state_linear_combination(coeff1, coeff2);
            interpolated_derivatives = derivative_linear_combination(coeff_dot_1, coeff_dot_2);
        }
else 
        {
            const double coeff1     =  one_minus_theta_h * theta;
            const double coeff2     = -one_minus_theta_h * (1.0 + theta);
            interpolated_state       = current_state_linear_combination(coeff1, coeff2);
            interpolated_derivatives = derivative_linear_combination(coeff_dot_1, coeff_dot_2);
        }

        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

}


