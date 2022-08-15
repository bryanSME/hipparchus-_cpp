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

//package org.hipparchus.ode.nonstiff;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.ode.FieldEquations_mapper;
//import org.hipparchus.ode.Field_ODE_State_And_Derivative;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"

/**
 * This class represents an interpolator over the last step during an
 * ODE integration for the 5(4) Higham and Hall integrator.
 *
 * @see Higham_Hall54_Field_Integrator
 *
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Higham_Hall54_Field_State_Interpolator : public Runge_Kutta_Field_State_Interpolator<T>
{

    /** Simple constructor.
     * @param field field to which the time and state vector elements belong
     * @param forward integration direction indicator
     * @param y_dot_k slopes at the intermediate points
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param soft_previous_state start of the restricted step
     * @param soft_current_state end of the restricted step
     * @param mapper equations mapper for the all equations
     */
    Higham_Hall54_Field_State_Interpolator(const Field<T> field, const bool forward, const std::vector<std::vector<T>> y_dot_k, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const Field_ODE_State_And_Derivative<T> soft_previous_state, const Field_ODE_State_And_Derivative<T> soft_current_state, const FieldEquations_mapper<T> mapper)
    {
        super(field, forward, y_dot_k, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
    }

    /** {@inherit_doc} */
    //override
    protected Higham_Hall54_Field_State_Interpolator<T> create(const Field<T> new_field, const bool new_forward, const std::vector<std::vector<T>> new_y_dot_k, const Field_ODE_State_And_Derivative<T> new_global_previous_state, const Field_ODE_State_And_Derivative<T> new_global_current_state, const Field_ODE_State_And_Derivative<T> new_soft_previous_state, const Field_ODE_State_And_Derivative<T> new_soft_current_state, const FieldEquations_mapper<T> new_mapper)
    {
        return Higham_Hall54_Field_State_Interpolator<T>(new_field, new_forward, new_y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** {@inherit_doc} */
    //@Suppress_Warnings("unchecked")
    //override
    protected Field_ODE_State_And_Derivative<T> compute_interpolated_state_and_derivatives(const FieldEquations_mapper<T> mapper, const T time, const T theta, const T theta_h, const T one_minus_theta_h)
    {

        const T& b_dot0 = theta.multiply(theta.multiply(theta.multiply(-10.0).add(16.0)).add(-15.0 / 2.0)).add(1);
        const T& b_dot1 = time.get_field().get_zero();
        const T& b_dot2 = theta.multiply(theta.multiply(theta.multiply(135.0 / 2.0).add(-729.0 / 8.0)).add(459.0 / 16.0));
        const T& b_dot3 = theta.multiply(theta.multiply(theta.multiply(-120.0).add(152.0)).add(-44.0));
        const T& b_dot4 = theta.multiply(theta.multiply(theta.multiply(125.0 / 2.0).add(-625.0 / 8.0)).add(375.0 / 16.0));
        const T& b_dot5 = theta.multiply(5.0 / 8.0).multiply(theta.multiply(2).subtract(1));
        const std::vector<T> interpolated_state;
        const std::vector<T> interpolated_derivatives;

        if (get_global_previous_state() != NULL && theta.get_real() <= 0.5)
        {
            const T& b0 = theta_h.multiply(theta.multiply(theta.multiply(theta.multiply(-5.0 / 2.0).add(16.0 / 3.0)).add(-15.0 / 4.0)).add(1));
            const T& b1 = time.get_field().get_zero();
            const T& b2 = theta_h.multiply(theta.multiply(theta.multiply(theta.multiply(135.0 / 8.0).add(-243.0 / 8.0)).add(459.0 / 32.0)));
            const T& b3 = theta_h.multiply(theta.multiply(theta.multiply(theta.multiply(-30.0).add(152.0 / 3.0)).add(-22.0)));
            const T& b4 = theta_h.multiply(theta.multiply(theta.multiply(theta.multiply(125.0 / 8.0).add(-625.0 / 24.0)).add(375.0 / 32.0)));
            const T& b5 = theta_h.multiply(theta.multiply(theta.multiply(5.0 / 12.0).add(-5.0 / 16.0)));
            interpolated_state = previous_state_linear_combination(b0, b1, b2, b3, b4, b5);
            interpolated_derivatives = derivative_linear_combination(b_dot0, b_dot1, b_dot2, b_dot3, b_dot4, b_dot5);
        }
        else
        {
            const T theta2 = theta.multiply(theta);
            const T h = theta_h.divide(theta);
            const T& b0 = h.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply(-5.0 / 2.0).add(16.0 / 3.0)).add(-15.0 / 4.0)).add(1.0)).add(-1.0 / 12.0));
            const T& b1 = time.get_field().get_zero();
            const T& b2 = h.multiply(theta2.multiply(theta.multiply(theta.multiply(135.0 / 8.0).add(-243.0 / 8.0)).add(459.0 / 32.0)).add(-27.0 / 32.0));
            const T& b3 = h.multiply(theta2.multiply(theta.multiply(theta.multiply(-30.0).add(152.0 / 3.0)).add(-22.0)).add(4.0 / 3.0));
            const T& b4 = h.multiply(theta2.multiply(theta.multiply(theta.multiply(125.0 / 8.0).add(-625.0 / 24.0)).add(375.0 / 32.0)).add(-125.0 / 96.0));
            const T& b5 = h.multiply(theta2.multiply(theta.multiply(5.0 / 12.0).add(-5.0 / 16.0)).add(-5.0 / 48.0));
            interpolated_state = current_state_linear_combination(b0, b1, b2, b3, b4, b5);
            interpolated_derivatives = derivative_linear_combination(b_dot0, b_dot1, b_dot2, b_dot3, b_dot4, b_dot5);
        }

        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

};