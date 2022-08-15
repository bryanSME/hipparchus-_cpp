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
//import org.hipparchus.util.Math_Arrays;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"

/**
 * This class : a simple Euler integrator for Ordinary
 * Differential Equations.
 *
 * <p>The Euler algorithm is the simplest one that can be used to
 * integrate ordinary differential equations. It is a simple inversion
 * of the forward difference expression :
 * <code>f'=(f(t+h)-f(t))/h</code> which leads to
 * <code>f(t+h)=f(t)+hf'</code>. The interpolation scheme used for
 * dense output is the linear scheme already used for integration.</p>
 *
 * <p>This algorithm looks cheap because it needs only one function
 * evaluation per step. However, as it uses linear estimates, it needs
 * very small steps to achieve high accuracy, and small steps lead to
 * numerical errors and instabilities.</p>
 *
 * <p>This algorithm is almost never used and has been included in
 * this //package only as a comparison reference for more useful
 * integrators.</p>
 *
 * @see Midpoint_Field_Integrator
 * @see ClassicalRunge_Kutta_Field_Integrator
 * @see Gill_Field_Integrator
 * @see Three_Eighthes_Field_Integrator
 * @see Luther_fieldIntegrator
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Euler_fieldIntegrator : public Runge_Kutta_Field_Integrator<T>
{

    /** Simple constructor.
     * Build an Euler integrator with the given step.
     * @param field field to which the time and state vector elements belong
     * @param step integration step
     */
    public Euler_fieldIntegrator(const Field<T> field, const T step)
    {
        super(field, "Euler", step);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> get_c()
    {
        return Math_Arrays::build_array(get_field(), 0);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<std::vector<T>> get_a()
    {
        return Math_Arrays::build_array(get_field(), 0, 0);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> get_b()
    {
        const std::vector<T> b = Math_Arrays::build_array(get_field(), 1);
        b[0] = get_field().get_one();
        return b;
    }

    /** {@inherit_doc} */
    //override
    protected Euler_fieldStateInterpolator<T>
        create_interpolator(const bool forward, std::vector<std::vector<T>> y_dot_k, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const FieldEquations_mapper<T> mapper)
    {
        return Euler_fieldStateInterpolator<T>(get_field(), forward, y_dot_k, global_previous_state, global_current_state, global_previous_state, global_current_state, mapper);
    }

};