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
 * This class : a linear interpolator for step.
 *
 * <p>This interpolator computes dense output inside the last
 * step computed. The interpolation equation is consistent with the
 * integration scheme :
 * <ul>
 *   <li>Using reference point at step start:<br>
 *     y(t<sub>n</sub> + &theta; h) = y (t<sub>n</sub>) + &theta; h y'
 *   </li>
 *   <li>Using reference point at step end:<br>
 *     y(t<sub>n</sub> + &theta; h) = y (t<sub>n</sub> + h) - (1-&theta;) h y'
 *   </li>
 * </ul>
 * </p>
 *
 * where &theta; belongs to [0 ; 1] and where y' is the evaluation of
 * the derivatives already computed during the step.</p>
 *
 * @see Euler_fieldIntegrator
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Euler_fieldStateInterpolator : public Runge_Kutta_Field_State_Interpolator<T>
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
    Euler_fieldStateInterpolator(const Field<T> field, const bool forward, const std::vector<std::vector<T>> y_dot_k, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const Field_ODE_State_And_Derivative<T> soft_previous_state, const Field_ODE_State_And_Derivative<T> soft_current_state, const FieldEquations_mapper<T> mapper)
    {
        super(field, forward, y_dot_k, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
    }

    /** {@inherit_doc} */
    //override
    protected Euler_fieldStateInterpolator<T> create(const Field<T> new_field, const bool new_forward, const std::vector<std::vector<T>> new_y_dot_k, const Field_ODE_State_And_Derivative<T> new_global_previous_state, const Field_ODE_State_And_Derivative<T> new_global_current_state, const Field_ODE_State_And_Derivative<T> new_soft_previous_state, const Field_ODE_State_And_Derivative<T> new_soft_current_state, const FieldEquations_mapper<T> new_mapper)
    {
        return Euler_fieldStateInterpolator<T>(new_field, new_forward, new_y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** {@inherit_doc} */
    //@Suppress_Warnings("unchecked")
    //override
    protected Field_ODE_State_And_Derivative<T> compute_interpolated_state_and_derivatives(const FieldEquations_mapper<T> mapper, const T time, const T theta, const T theta_h, const T one_minus_theta_h)
    {
        const std::vector<T> interpolated_state;
        const std::vector<T> interpolated_derivatives;
        if (get_global_previous_state() != NULL && theta.get_real() <= 0.5)
        {
            interpolated_state = previous_state_linear_combination(theta_h);
            interpolated_derivatives = derivative_linear_combination(time.get_field().get_one());
        }
        else
        {
            interpolated_state = current_state_linear_combination(one_minus_theta_h.negate());
            interpolated_derivatives = derivative_linear_combination(time.get_field().get_one());
        }

        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

};