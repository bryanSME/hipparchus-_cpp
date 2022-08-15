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
 * This class represents an interpolator over the last step during an
 * ODE integration for the 5(4) Higham and Hall integrator.
 *
 * @see Higham_Hall54_Integrator
 *
 */

class Higham_Hall54_State_Interpolator
    extends Runge_Kutta_State_Interpolator 
    {

    /** Serializable version identifier */
    private static const long serial_version_uid = 20111120L;

    /** Simple constructor.
     * @param forward integration direction indicator
     * @param y_dot_k slopes at the intermediate points
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param soft_previous_state start of the restricted step
     * @param soft_current_state end of the restricted step
     * @param mapper equations mapper for the all equations
     */
    Higham_Hall54_State_Interpolator(const bool forward, const std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const ODE_State_And_Derivative soft_previous_state, const ODE_State_And_Derivative soft_current_state, const Equations_mapper mapper) 
    {
        super(forward, y_dot_k, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
    }

    /** {@inherit_doc} */
    //override
    protected Higham_Hall54_State_Interpolator create(const bool new_forward, const std::vector<std::vector<double>> new_y_dot_k, const ODE_State_And_Derivative new_global_previous_state, const ODE_State_And_Derivative new_global_current_state, const ODE_State_And_Derivative new_soft_previous_state, const ODE_State_And_Derivative new_soft_current_state, const Equations_mapper new_mapper) 
    {
        return Higham_Hall54_State_Interpolator(new_forward, new_y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** {@inherit_doc} */
    //override
    protected ODE_State_And_Derivative compute_interpolated_state_and_derivatives(const Equations_mapper mapper, const double time, const double& theta, const double theta_h, const double one_minus_theta_h) 
    {

        const double b_dot0 = 1 + theta * (-15.0/2.0 + theta * (16.0 - 10.0 * theta));
        const double b_dot1 = 0;
        const double b_dot2 = theta * (459.0/16.0 + theta * (-729.0/8.0 + 135.0/2.0 * theta));
        const double b_dot3 = theta * (-44.0 + theta * (152.0 - 120.0 * theta));
        const double b_dot4 = theta * (375.0/16.0 + theta * (-625.0/8.0 + 125.0/2.0 * theta));
        const double b_dot5 = theta * 5.0/8.0 * (2 * theta - 1);

        const std::vector<double> interpolated_state;
        const std::vector<double> interpolated_derivatives;
        if (get_global_previous_state() != NULL && theta <= 0.5) 
        {
            const double b0 = theta_h * (1.0 + theta * (-15.0/4.0  + theta * (16.0/3.0 - 5.0/2.0 * theta)));
            const double b1 = 0;
            const double b2 = theta_h * (      theta * (459.0/32.0 + theta * (-243.0/8.0 + theta * 135.0/8.0)));
            const double b3 = theta_h * (      theta * (-22.0      + theta * (152.0/3.0  + theta * -30.0)));
            const double b4 = theta_h * (      theta * (375.0/32.0 + theta * (-625.0/24.0 + theta * 125.0/8.0)));
            const double b5 = theta_h * (      theta * (-5.0/16.0  + theta *  5.0/12.0));
            interpolated_state       = previous_state_linear_combination(b0 , b1, b2, b3, b4, b5);
            interpolated_derivatives = derivative_linear_combination(b_dot0 , b_dot1, b_dot2, b_dot3, b_dot4, b_dot5);
        }
else 
        {
            const double theta2 = theta * theta;
            const double h      = theta_h / theta;
            const double b0 = h * (-1.0/12.0 + theta * (1.0 + theta * (-15.0/4.0 + theta * (16.0/3.0 + theta * -5.0/2.0))));
            const double b1 = 0;
            const double b2 = h * (-27.0/32.0 + theta2 * (459.0/32.0 + theta * (-243.0/8.0 + theta * 135.0/8.0)));
            const double b3 = h * (4.0/3.0 + theta2 * (-22.0 + theta * (152.0/3.0  + theta * -30.0)));
            const double b4 = h * (-125.0/96.0 + theta2 * (375.0/32.0 + theta * (-625.0/24.0 + theta * 125.0/8.0)));
            const double b5 = h * (-5.0/48.0 + theta2 * (-5.0/16.0 + theta * 5.0/12.0));
            interpolated_state       = current_state_linear_combination(b0 , b1, b2, b3, b4, b5);
            interpolated_derivatives = derivative_linear_combination(b_dot0 , b_dot1, b_dot2, b_dot3, b_dot4, b_dot5);
        }

        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

}


