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
 * This class : a step interpolator for the classical fourth
 * order Runge-Kutta integrator.
 *
 * <p>This interpolator allows to compute dense output inside the last
 * step computed. The interpolation equation is consistent with the
 * integration scheme :
 * <ul>
 *   <li>Using reference point at step start:<br>
 *   y(t<sub>n</sub> + &theta; h) = y (t<sub>n</sub>)
 *                    + &theta; (h/6) [  (6 - 9 &theta; + 4 &theta;<sup>2</sup>) y'<sub>1</sub>
 *                                     + (    6 &theta; - 4 &theta;<sup>2</sup>) (y'<sub>2</sub> + y'<sub>3</sub>)
 *                                     + (   -3 &theta; + 4 &theta;<sup>2</sup>) y'<sub>4</sub>
 *                                    ]
 *   </li>
 *   <li>Using reference point at step end:<br>
 *   y(t<sub>n</sub> + &theta; h) = y (t<sub>n</sub> + h)
 *                    + (1 - &theta;) (h/6) [ (-4 &theta;^2 + 5 &theta; - 1) y'<sub>1</sub>
 *                                          +(4 &theta;^2 - 2 &theta; - 2) (y'<sub>2</sub> + y'<sub>3</sub>)
 *                                          -(4 &theta;^2 +   &theta; + 1) y'<sub>4</sub>
 *                                        ]
 *   </li>
 * </ul>
 * </p>
 *
 * where &theta; belongs to [0 ; 1] and where y'<sub>1</sub> to y'<sub>4</sub> are the four
 * evaluations of the derivatives already computed during the
 * step.</p>
 *
 * @see ClassicalRunge_Kutta_Field_Integrator
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class ClassicalRunge_Kutta_Field_State_Interpolator : Runge_Kutta_Field_State_Interpolator<T>
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
    ClassicalRunge_Kutta_Field_State_Interpolator(const Field<T> field, const bool forward, const std::vector<std::vector<T>> y_dot_k, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const Field_ODE_State_And_Derivative<T> soft_previous_state, const Field_ODE_State_And_Derivative<T> soft_current_state, const FieldEquations_mapper<T> mapper)
    {
        super(field, forward, y_dot_k, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
    }

    /** {@inherit_doc} */
    //override
    protected ClassicalRunge_Kutta_Field_State_Interpolator<T> create(const Field<T> new_field, const bool new_forward, const std::vector<std::vector<T>> new_y_dot_k, const Field_ODE_State_And_Derivative<T> new_global_previous_state, const Field_ODE_State_And_Derivative<T> new_global_current_state, const Field_ODE_State_And_Derivative<T> new_soft_previous_state, const Field_ODE_State_And_Derivative<T> new_soft_current_state, const FieldEquations_mapper<T> new_mapper)
    {
        return ClassicalRunge_Kutta_Field_State_Interpolator<T>(new_field, new_forward, new_y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** {@inherit_doc} */
    //@Suppress_Warnings("unchecked")
    //override
    protected Field_ODE_State_And_Derivative<T> compute_interpolated_state_and_derivatives(const FieldEquations_mapper<T> mapper, const T time, const T theta, const T theta_h, const T one_minus_theta_h)
    {

        const T one = time.get_field().get_one();
        const T one_minus_theta = one.subtract(theta);
        const T one_minus2_theta = one.subtract(theta.multiply(2));
        const T coeff_dot_1 = one_minus_theta.multiply(one_minus2_theta);
        const T coeff_dot_23 = theta.multiply(one_minus_theta).multiply(2);
        const T coeff_dot_4 = theta.multiply(one_minus2_theta).negate();
        const std::vector<T> interpolated_state;
        const std::vector<T> interpolated_derivatives;

        if (get_global_previous_state() != NULL && theta.get_real() <= 0.5)
        {
            const T four_theta_2 = theta.multiply(theta).multiply(4);
            const T s = theta_h.divide(6.0);
            const T coeff1 = s.multiply(four_theta_2.subtract(theta.multiply(9)).add(6));
            const T coeff23 = s.multiply(theta.multiply(6).subtract(four_theta_2));
            const T coeff4 = s.multiply(four_theta_2.subtract(theta.multiply(3)));
            interpolated_state = previous_state_linear_combination(coeff1, coeff23, coeff23, coeff4);
            interpolated_derivatives = derivative_linear_combination(coeff_dot_1, coeff_dot_23, coeff_dot_23, coeff_dot_4);
        }
        else
        {
            const T four_theta = theta.multiply(4);
            const T s = one_minus_theta_h.divide(6);
            const T coeff1 = s.multiply(theta.multiply(four_theta.negate().add(5)).subtract(1));
            const T coeff23 = s.multiply(theta.multiply(four_theta.subtract(2)).subtract(2));
            const T coeff4 = s.multiply(theta.multiply(four_theta.negate().subtract(1)).subtract(1));
            interpolated_state = current_state_linear_combination(coeff1, coeff23, coeff23, coeff4);
            interpolated_derivatives = derivative_linear_combination(coeff_dot_1, coeff_dot_23, coeff_dot_23, coeff_dot_4);
        }

        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

};