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
//import org.hipparchus.ode.sampling.AbstractFieldODE_StateInterpolator;
//import org.hipparchus.util.Math_Arrays;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"
#include "../../core/util/MathArrays.h"

/** This class represents an interpolator over the last step during an
 * ODE integration for Runge-Kutta and embedded Runge-Kutta integrators.
 *
 * @see Runge_Kutta_Field_Integrator
 * @see EmbeddedRunge_Kutta_Field_Integrator
 *
 * @param <T> the type of the field elements
 */

template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Runge_Kutta_Field_State_Interpolator : public AbstractFieldODE_StateInterpolator<T>
{

    /** Field to which the time and state vector elements belong. */
    private const Field<T> field;

    /** Slopes at the intermediate points. */
    private const std::vector<std::vector<T>> y_dot_k;

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
    protected Runge_Kutta_Field_State_Interpolator(const Field<T> field, const bool forward, const std::vector<std::vector<T>> y_dot_k, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const Field_ODE_State_And_Derivative<T> soft_previous_state, const Field_ODE_State_And_Derivative<T> soft_current_state, const FieldEquations_mapper<T> mapper)
    {
        super(forward, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
        this.field = field;
        this.y_dot_k = Math_Arrays::build_array(field, y_dot_k.size(), -1);
        for (int i{}; i < y_dot_k.size(); ++i)
        {
            this.y_dot_k[i] = y_dot_k[i].clone();
        }
    }

    /** {@inherit_doc} */
    //override
    protected Runge_Kutta_Field_State_Interpolator<T> create(bool new_forward, Field_ODE_State_And_Derivative<T> new_global_previous_state, Field_ODE_State_And_Derivative<T> new_global_current_state, Field_ODE_State_And_Derivative<T> new_soft_previous_state, Field_ODE_State_And_Derivative<T> new_soft_current_state, FieldEquations_mapper<T> new_mapper)
    {
        return create(field, new_forward, y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** Create a instance.
     * @param new_field field to which the time and state vector elements belong
     * @param new_forward integration direction indicator
     * @param new_y_dot_k slopes at the intermediate points
     * @param new_global_previous_state start of the global step
     * @param new_global_current_state end of the global step
     * @param new_soft_previous_state start of the restricted step
     * @param new_soft_current_state end of the restricted step
     * @param new_mapper equations mapper for the all equations
     * @return a instance
     */
    protected virtual Runge_Kutta_Field_State_Interpolator<T> create(Field<T> new_field, bool new_forward, std::vector<std::vector<T>> new_y_dot_k, Field_ODE_State_And_Derivative<T> new_global_previous_state, Field_ODE_State_And_Derivative<T> new_global_current_state, Field_ODE_State_And_Derivative<T> new_soft_previous_state, Field_ODE_State_And_Derivative<T> new_soft_current_state, FieldEquations_mapper<T> new_mapper);

    /** Compute a state by linear combination added to previous state.
     * @param coefficients coefficients to apply to the method staged derivatives
     * @return combined state
     */
     //@Safe_Varargs
    protected const std::vector<T> previous_state_linear_combination(const T ... coefficients)
    {
        return combine(get_global_previous_state().get_complete_state(), coefficients);
    }

    /** Compute a state by linear combination added to current state.
     * @param coefficients coefficients to apply to the method staged derivatives
     * @return combined state
     */
     //@Suppress_Warnings("unchecked")
    protected std::vector<T> current_state_linear_combination(const T ... coefficients)
    {
        return combine(get_global_current_state().get_complete_state(), coefficients);
    }

    /** Compute a state derivative by linear combination.
     * @param coefficients coefficients to apply to the method staged derivatives
     * @return combined state
     */
     //@Suppress_Warnings("unchecked")
    protected std::vector<T> derivative_linear_combination(const T ... coefficients)
    {
        return combine(Math_Arrays::build_array(field, y_dot_k[0].size()), coefficients);
    }

    /** Linearly combine arrays.
     * @param a array to add to
     * @param coefficients coefficients to apply to the method staged derivatives
     * @return a itself, as a conveniency for fluent API
     */
     //@Suppress_Warnings("unchecked")
    private std::vector<T> combine(const std::vector<T> a, const T ... coefficients)
    {
        for (int i{}; i < a.size(); ++i)
        {
            for (int k{}; k < coefficients.size(); ++k)
            {
                a[i] = a[i].add(coefficients[k].multiply(y_dot_k[k][i]));
            }
        }
        return a;
    }

};