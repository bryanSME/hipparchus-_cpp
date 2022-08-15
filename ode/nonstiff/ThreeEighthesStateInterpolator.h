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
 * This class : a step interpolator for the 3/8 fourth
 * order Runge-Kutta integrator.
 *
 * <p>This interpolator allows to compute dense output inside the last
 * step computed. The interpolation equation is consistent with the
 * integration scheme :
 * <ul>
 *   <li>Using reference point at step start:<br>
 *     y(t<sub>n</sub> + &theta; h) = y (t<sub>n</sub>)
 *                      + &theta; (h/8) [ (8 - 15 &theta; +  8 &theta;<sup>2</sup>) y'<sub>1</sub>
 *                                     +  3 * (15 &theta; - 12 &theta;<sup>2</sup>) y'<sub>2</sub>
 *                                     +        3 &theta;                           y'<sub>3</sub>
 *                                     +      (-3 &theta; +  4 &theta;<sup>2</sup>) y'<sub>4</sub>
 *                                    ]
 *   </li>
 *   <li>Using reference point at step end:<br>
 *     y(t<sub>n</sub> + &theta; h) = y (t<sub>n</sub> + h)
 *                      - (1 - &theta;) (h/8) [(1 - 7 &theta; + 8 &theta;<sup>2</sup>) y'<sub>1</sub>
 *                                         + 3 (1 +   &theta; - 4 &theta;<sup>2</sup>) y'<sub>2</sub>
 *                                         + 3 (1 +   &theta;)                         y'<sub>3</sub>
 *                                         +   (1 +   &theta; + 4 &theta;<sup>2</sup>) y'<sub>4</sub>
 *                                          ]
 *   </li>
 * </ul>
 * </p>
 *
 * where &theta; belongs to [0 ; 1] and where y'<sub>1</sub> to y'<sub>4</sub> are the four
 * evaluations of the derivatives already computed during the
 * step.</p>
 *
 * @see Three_Eighthes_Integrator
 */

class Three_Eighthes_State_Interpolator
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
    Three_Eighthes_State_Interpolator(const bool forward, const std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const ODE_State_And_Derivative soft_previous_state, const ODE_State_And_Derivative soft_current_state, const Equations_mapper mapper) 
    {
        super(forward, y_dot_k, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
    }

    /** {@inherit_doc} */
    //override
    protected Three_Eighthes_State_Interpolator create(const bool new_forward, const std::vector<std::vector<double>> new_y_dot_k, const ODE_State_And_Derivative new_global_previous_state, const ODE_State_And_Derivative new_global_current_state, const ODE_State_And_Derivative new_soft_previous_state, const ODE_State_And_Derivative new_soft_current_state, const Equations_mapper new_mapper) 
    {
        return Three_Eighthes_State_Interpolator(new_forward, new_y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** {@inherit_doc} */
    //override
    protected ODE_State_And_Derivative compute_interpolated_state_and_derivatives(const Equations_mapper mapper, const double time, const double& theta, const double theta_h, const double one_minus_theta_h) 
    {

        const double coeff_dot_3  = 0.75 * theta;
        const double coeff_dot_1  = coeff_dot_3 * (4 * theta - 5) + 1;
        const double coeff_dot_2  = coeff_dot_3 * (5 - 6 * theta);
        const double coeff_dot_4  = coeff_dot_3 * (2 * theta - 1);
        const std::vector<double> interpolated_state;
        const std::vector<double> interpolated_derivatives;

        if (get_global_previous_state() != NULL && theta <= 0.5) 
        {
            const double s          = theta_h / 8.0;
            const double four_theta_2 = 4 * theta * theta;
            const double coeff1     = s * (8 - 15 * theta + 2 * four_theta_2);
            const double coeff2     = 3 * s * (5 * theta - four_theta_2);
            const double coeff3     = 3 * s * theta;
            const double coeff4     = s * (-3 * theta + four_theta_2);
            interpolated_state       = previous_state_linear_combination(coeff1, coeff2, coeff3, coeff4);
            interpolated_derivatives = derivative_linear_combination(coeff_dot_1, coeff_dot_2, coeff_dot_3, coeff_dot_4);
        }
else 
        {
            const double s          = one_minus_theta_h / -8.0;
            const double four_theta_2 = 4 * theta * theta;
            const double coeff1     = s * (1 - 7 * theta + 2 * four_theta_2);
            const double coeff2     = 3 * s * (1 + theta - four_theta_2);
            const double coeff3     = 3 * s * (1 + theta);
            const double coeff4     = s * (1 + theta + four_theta_2);
            interpolated_state       = current_state_linear_combination(coeff1, coeff2, coeff3, coeff4);
            interpolated_derivatives = derivative_linear_combination(coeff_dot_1, coeff_dot_2, coeff_dot_3, coeff_dot_4);
        }

        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

}


