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
 * This class : the classical fourth order Runge-Kutta
 * integrator for Ordinary Differential Equations (it is the most
 * often used Runge-Kutta method).
 *
 * <p>This method is an explicit Runge-Kutta method, its Butcher-array
 * is the following one :
 * <pre>
 *    0  |  0    0    0    0
 *   1/2 | 1/2   0    0    0
 *   1/2 |  0   1/2   0    0
 *    1  |  0    0    1    0
 *       |--------------------
 *       | 1/6  1/3  1/3  1/6
 * </pre>
 * </p>
 *
 * @see Euler_fieldIntegrator
 * @see Gill_Field_Integrator
 * @see Midpoint_Field_Integrator
 * @see Three_Eighthes_Field_Integrator
 * @see Luther_fieldIntegrator
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class ClassicalRunge_Kutta_Field_Integrator : public Runge_Kutta_Field_Integrator<T>
{

    /** Simple constructor.
     * Build a fourth-order Runge-Kutta integrator with the given step.
     * @param field field to which the time and state vector elements belong
     * @param step integration step
     */
    public ClassicalRunge_Kutta_Field_Integrator(const Field<T> field, const T step)
    {
        super(field, "classical Runge-Kutta", step);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> get_c()
    {
        const std::vector<T> c = Math_Arrays::build_array(get_field(), 3);
        c[0] = get_field().get_one().multiply(0.5);
        c[1] = c[0];
        c[2] = get_field().get_one();
        return c;
    }

    /** {@inherit_doc} */
    //override
    public std::vector<std::vector<T>> get_a()
    {
        const std::vector<std::vector<T>> a = Math_Arrays::build_array(get_field(), 3, -1);
        for (int i{}; i < a.size(); ++i)
        {
            a[i] = Math_Arrays::build_array(get_field(), i + 1);
        }
        a[0][0] = fraction(1, 2);
        a[1][0] = get_field().get_zero();
        a[1][1] = a[0][0];
        a[2][0] = get_field().get_zero();
        a[2][1] = get_field().get_zero();
        a[2][2] = get_field().get_one();
        return a;
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> get_b()
    {
        const std::vector<T> b = Math_Arrays::build_array(get_field(), 4);
        b[0] = fraction(1, 6);
        b[1] = fraction(1, 3);
        b[2] = b[1];
        b[3] = b[0];
        return b;
    }

    /** {@inherit_doc} */
    //override
    protected ClassicalRunge_Kutta_Field_State_Interpolator<T>
        create_interpolator(const bool forward, std::vector<std::vector<T>> y_dot_k, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const FieldEquations_mapper<T> mapper)
    {
        return ClassicalRunge_Kutta_Field_State_Interpolator<T>(get_field(), forward, y_dot_k, global_previous_state, global_current_state, global_previous_state, global_current_state, mapper);
    }

};