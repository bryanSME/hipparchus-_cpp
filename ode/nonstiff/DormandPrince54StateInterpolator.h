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
 * ODE integration for the 5(4) Dormand-Prince integrator.
 *
 * @see Dormand_Prince54_Integrator
 *
 */

class Dormand_Prince54_State_Interpolator
    extends Runge_Kutta_State_Interpolator 
    {

    /** Last row of the Butcher-array internal weights, element 0. */
    private static const double A70 =    35.0 /  384.0;

    // element 1 is zero, so it is neither stored nor used

    /** Last row of the Butcher-array internal weights, element 2. */
    private static const double A72 =   500.0 / 1113.0;

    /** Last row of the Butcher-array internal weights, element 3. */
    private static const double A73 =   125.0 /  192.0;

    /** Last row of the Butcher-array internal weights, element 4. */
    private static const double A74 = -2187.0 / 6784.0;

    /** Last row of the Butcher-array internal weights, element 5. */
    private static const double A75 =    11.0 /   84.0;

    /** Shampine (1986) Dense output, element 0. */
    private static const double D0 =  -12715105075.0 /  11282082432.0;

    // element 1 is zero, so it is neither stored nor used

    /** Shampine (1986) Dense output, element 2. */
    private static const double D2 =   87487479700.0 /  32700410799.0;

    /** Shampine (1986) Dense output, element 3. */
    private static const double D3 =  -10690763975.0 /   1880347072.0;

    /** Shampine (1986) Dense output, element 4. */
    private static const double D4 =  701980252875.0 / 199316789632.0;

    /** Shampine (1986) Dense output, element 5. */
    private static const double D5 =   -1453857185.0 /    822651844.0;

    /** Shampine (1986) Dense output, element 6. */
    private static const double D6 =      69997945.0 /     29380423.0;

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
    Dormand_Prince54_State_Interpolator(const bool forward, const std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const ODE_State_And_Derivative soft_previous_state, const ODE_State_And_Derivative soft_current_state, const Equations_mapper mapper) 
    {
        super(forward, y_dot_k, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
    }

    /** {@inherit_doc} */
    //override
    protected Dormand_Prince54_State_Interpolator create(const bool new_forward, const std::vector<std::vector<double>> new_y_dot_k, const ODE_State_And_Derivative new_global_previous_state, const ODE_State_And_Derivative new_global_current_state, const ODE_State_And_Derivative new_soft_previous_state, const ODE_State_And_Derivative new_soft_current_state, const Equations_mapper new_mapper) 
    {
        return Dormand_Prince54_State_Interpolator(new_forward, new_y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** {@inherit_doc} */
    //override
    protected ODE_State_And_Derivative compute_interpolated_state_and_derivatives(const Equations_mapper mapper, const double time, const double& theta, const double theta_h, const double one_minus_theta_h) 
    {

        // interpolate
        const double eta = 1 - theta;
        const double two_theta = 2 * theta;
        const double dot2 = 1 - two_theta;
        const double dot3 = theta * (2 - 3 * theta);
        const double dot4 = two_theta * (1 + theta * (two_theta - 3));

        const std::vector<double> interpolated_state;
        const std::vector<double> interpolated_derivatives;
        if (get_global_previous_state() != NULL && theta <= 0.5) 
        {
            const double f1        = theta_h;
            const double f2        = f1 * eta;
            const double f3        = f2 * theta;
            const double f4        = f3 * eta;
            const double coeff0    = f1 * A70 - f2   * (A70 - 1) + f3   * (2 * A70 - 1) + f4   * D0;
            const double coeff1    = 0;
            const double coeff2    = f1 * A72 - f2   * A72       + f3   * (2 * A72)     + f4   * D2;
            const double coeff3    = f1 * A73 - f2   * A73       + f3   * (2 * A73)     + f4   * D3;
            const double coeff4    = f1 * A74 - f2   * A74       + f3   * (2 * A74)     + f4   * D4;
            const double coeff5    = f1 * A75 - f2   * A75       + f3   * (2 * A75)     + f4   * D5;
            const double coeff6    = f4 * D6 - f3;
            const double coeff_dot0 =      A70 - dot2 * (A70 - 1) + dot3 * (2 * A70 - 1) + dot4 * D0;
            const double coeff_dot_1 = 0;
            const double coeff_dot_2 =      A72 - dot2 * A72       + dot3 * (2 * A72)     + dot4 * D2;
            const double coeff_dot_3 =      A73 - dot2 * A73       + dot3 * (2 * A73)     + dot4 * D3;
            const double coeff_dot_4 =      A74 - dot2 * A74       + dot3 * (2 * A74)     + dot4 * D4;
            const double coeff_dot5 =      A75 - dot2 * A75       + dot3 * (2 * A75)     + dot4 * D5;
            const double coeff_dot6 = dot4 * D6 - dot3;
            interpolated_state       = previous_state_linear_combination(coeff0, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6);
            interpolated_derivatives = derivative_linear_combination(coeff_dot0, coeff_dot_1, coeff_dot_2, coeff_dot_3, coeff_dot_4, coeff_dot5, coeff_dot6);
        }
else 
        {
            const double f1        = -one_minus_theta_h;
            const double f2        = one_minus_theta_h * theta;
            const double f3        = f2 * theta;
            const double f4        = f3 * eta;
            const double coeff0    = f1 * A70 - f2   * (A70 - 1) + f3   * (2 * A70 - 1) + f4   * D0;
            const double coeff1    = 0;
            const double coeff2    = f1 * A72 - f2   * A72       + f3   * (2 * A72)     + f4   * D2;
            const double coeff3    = f1 * A73 - f2   * A73       + f3   * (2 * A73)     + f4   * D3;
            const double coeff4    = f1 * A74 - f2   * A74       + f3   * (2 * A74)     + f4   * D4;
            const double coeff5    = f1 * A75 - f2   * A75       + f3   * (2 * A75)     + f4   * D5;
            const double coeff6    = f4 * D6 - f3;
            const double coeff_dot0 =      A70 - dot2 * (A70 - 1) + dot3 * (2 * A70 - 1) + dot4 * D0;
            const double coeff_dot_1 = 0;
            const double coeff_dot_2 =      A72 - dot2 * A72       + dot3 * (2 * A72)     + dot4 * D2;
            const double coeff_dot_3 =      A73 - dot2 * A73       + dot3 * (2 * A73)     + dot4 * D3;
            const double coeff_dot_4 =      A74 - dot2 * A74       + dot3 * (2 * A74)     + dot4 * D4;
            const double coeff_dot5 =      A75 - dot2 * A75       + dot3 * (2 * A75)     + dot4 * D5;
            const double coeff_dot6 = dot4 * D6 - dot3;
            interpolated_state       = current_state_linear_combination(coeff0, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6);
            interpolated_derivatives = derivative_linear_combination(coeff_dot0, coeff_dot_1, coeff_dot_2, coeff_dot_3, coeff_dot_4, coeff_dot5, coeff_dot6);
        }

        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

}


