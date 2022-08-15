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
 * ODE integration for the 5(4) Dormand-Prince integrator.
 *
 * @see Dormand_Prince54_Integrator
 *
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Dormand_Prince54_Field_State_Interpolator : public Runge_Kutta_Field_State_Interpolator<T> 
{

    /** Last row of the Butcher-array internal weights, element 0. */
    private const T a70;

    // element 1 is zero, so it is neither stored nor used

    /** Last row of the Butcher-array internal weights, element 2. */
    private const T a72;

    /** Last row of the Butcher-array internal weights, element 3. */
    private const T a73;

    /** Last row of the Butcher-array internal weights, element 4. */
    private const T a74;

    /** Last row of the Butcher-array internal weights, element 5. */
    private const T a75;

    /** Shampine (1986) Dense output, element 0. */
    private const T d0;

    // element 1 is zero, so it is neither stored nor used

    /** Shampine (1986) Dense output, element 2. */
    private const T d2;

    /** Shampine (1986) Dense output, element 3. */
    private const T d3;

    /** Shampine (1986) Dense output, element 4. */
    private const T d4;

    /** Shampine (1986) Dense output, element 5. */
    private const T d5;

    /** Shampine (1986) Dense output, element 6. */
    private const T d6;

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
    Dormand_Prince54_Field_State_Interpolator(const Field<T> field, const bool forward, const std::vector<std::vector<T>> y_dot_k, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const Field_ODE_State_And_Derivative<T> soft_previous_state, const Field_ODE_State_And_Derivative<T> soft_current_state, const FieldEquations_mapper<T> mapper) 
    {
        super(field, forward, y_dot_k, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
        const T one = field.get_one();
        a70 = one.multiply(   35.0).divide( 384.0);
        a72 = one.multiply(  500.0).divide(1113.0);
        a73 = one.multiply(  125.0).divide( 192.0);
        a74 = one.multiply(-2187.0).divide(6784.0);
        a75 = one.multiply(   11.0).divide(  84.0);
        d0  = one.multiply(-12715105075.0).divide( 11282082432.0);
        d2  = one.multiply( 87487479700.0).divide( 32700410799.0);
        d3  = one.multiply(-10690763975.0).divide(  1880347072.0);
        d4  = one.multiply(701980252875.0).divide(199316789632.0);
        d5  = one.multiply( -1453857185.0).divide(   822651844.0);
        d6  = one.multiply(    69997945.0).divide(    29380423.0);
    }

    /** {@inherit_doc} */
    //override
    protected Dormand_Prince54_Field_State_Interpolator<T> create(const Field<T> new_field, const bool new_forward, const std::vector<std::vector<T>> new_y_dot_k, const Field_ODE_State_And_Derivative<T> new_global_previous_state, const Field_ODE_State_And_Derivative<T> new_global_current_state, const Field_ODE_State_And_Derivative<T> new_soft_previous_state, const Field_ODE_State_And_Derivative<T> new_soft_current_state, const FieldEquations_mapper<T> new_mapper) 
    {
        return Dormand_Prince54_Field_State_Interpolator<T>(new_field, new_forward, new_y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }
    /** {@inherit_doc} */
    //@Suppress_Warnings("unchecked")
    //override
    protected Field_ODE_State_And_Derivative<T> compute_interpolated_state_and_derivatives(const FieldEquations_mapper<T> mapper, const T time, const T theta, const T theta_h, const T one_minus_theta_h) 
    {

        // interpolate
        const T one      = time.get_field().get_one();
        const T eta      = one.subtract(theta);
        const T two_theta = theta.multiply(2);
        const T dot2     = one.subtract(two_theta);
        const T dot3     = theta.multiply(theta.multiply(-3).add(2));
        const T dot4     = two_theta.multiply(theta.multiply(two_theta.subtract(3)).add(1));
        const std::vector<T> interpolated_state;
        const std::vector<T> interpolated_derivatives;
        if (get_global_previous_state() != NULL && theta.get_real() <= 0.5) 
        {
            const T f1        = theta_h;
            const T f2        = f1.multiply(eta);
            const T f3        = f2.multiply(theta);
            const T f4        = f3.multiply(eta);
            const T coeff0    = f1.multiply(a70).
                                subtract(f2.multiply(a70.subtract(1))).
                                add(f3.multiply(a70.multiply(2).subtract(1))).
                                add(f4.multiply(d0));
            const T coeff1    = time.get_field().get_zero();
            const T coeff2    = f1.multiply(a72).
                                subtract(f2.multiply(a72)).
                                add(f3.multiply(a72.multiply(2))).
                                add(f4.multiply(d2));
            const T coeff3    = f1.multiply(a73).
                                subtract(f2.multiply(a73)).
                                add(f3.multiply(a73.multiply(2))).
                                add(f4.multiply(d3));
            const T coeff4    = f1.multiply(a74).
                                subtract(f2.multiply(a74)).
                                add(f3.multiply(a74.multiply(2))).
                                add(f4.multiply(d4));
            const T coeff5    = f1.multiply(a75).
                                subtract(f2.multiply(a75)).
                                add(f3.multiply(a75.multiply(2))).
                                add(f4.multiply(d5));
            const T coeff6    = f4.multiply(d6).subtract(f3);
            const T coeff_dot0 = a70.
                                subtract(dot2.multiply(a70.subtract(1))).
                                add(dot3.multiply(a70.multiply(2).subtract(1))).
                                add(dot4.multiply(d0));
            const T coeff_dot_1 = time.get_field().get_zero();
            const T coeff_dot_2 = a72.
                                subtract(dot2.multiply(a72)).
                                add(dot3.multiply(a72.multiply(2))).
                                add(dot4.multiply(d2));
            const T coeff_dot_3 = a73.
                                subtract(dot2.multiply(a73)).
                                add(dot3.multiply(a73.multiply(2))).
                                add(dot4.multiply(d3));
            const T coeff_dot_4 = a74.
                                subtract(dot2.multiply(a74)).
                                add(dot3.multiply(a74.multiply(2))).
                                add(dot4.multiply(d4));
            const T coeff_dot5 = a75.
                                subtract(dot2.multiply(a75)).
                                add(dot3.multiply(a75.multiply(2))).
                                add(dot4.multiply(d5));
            const T coeff_dot6 = dot4.multiply(d6).subtract(dot3);
            interpolated_state       = previous_state_linear_combination(coeff0, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6);
            interpolated_derivatives = derivative_linear_combination(coeff_dot0, coeff_dot_1, coeff_dot_2, coeff_dot_3, coeff_dot_4, coeff_dot5, coeff_dot6);
        }
else 
        {
            const T f1        = one_minus_theta_h.negate();
            const T f2        = one_minus_theta_h.multiply(theta);
            const T f3        = f2.multiply(theta);
            const T f4        = f3.multiply(eta);
            const T coeff0    = f1.multiply(a70).
                                subtract(f2.multiply(a70.subtract(1))).
                                add(f3.multiply(a70.multiply(2).subtract(1))).
                                add(f4.multiply(d0));
            const T coeff1    = time.get_field().get_zero();
            const T coeff2    = f1.multiply(a72).
                                subtract(f2.multiply(a72)).
                                add(f3.multiply(a72.multiply(2))).
                                add(f4.multiply(d2));
            const T coeff3    = f1.multiply(a73).
                                subtract(f2.multiply(a73)).
                                add(f3.multiply(a73.multiply(2))).
                                add(f4.multiply(d3));
            const T coeff4    = f1.multiply(a74).
                                subtract(f2.multiply(a74)).
                                add(f3.multiply(a74.multiply(2))).
                                add(f4.multiply(d4));
            const T coeff5    = f1.multiply(a75).
                                subtract(f2.multiply(a75)).
                                add(f3.multiply(a75.multiply(2))).
                                add(f4.multiply(d5));
            const T coeff6    = f4.multiply(d6).subtract(f3);
            const T coeff_dot0 = a70.
                                subtract(dot2.multiply(a70.subtract(1))).
                                add(dot3.multiply(a70.multiply(2).subtract(1))).
                                add(dot4.multiply(d0));
            const T coeff_dot_1 = time.get_field().get_zero();
            const T coeff_dot_2 = a72.
                                subtract(dot2.multiply(a72)).
                                add(dot3.multiply(a72.multiply(2))).
                                add(dot4.multiply(d2));
            const T coeff_dot_3 = a73.
                                subtract(dot2.multiply(a73)).
                                add(dot3.multiply(a73.multiply(2))).
                                add(dot4.multiply(d3));
            const T coeff_dot_4 = a74.
                                subtract(dot2.multiply(a74)).
                                add(dot3.multiply(a74.multiply(2))).
                                add(dot4.multiply(d4));
            const T coeff_dot5 = a75.
                                subtract(dot2.multiply(a75)).
                                add(dot3.multiply(a75.multiply(2))).
                                add(dot4.multiply(d5));
            const T coeff_dot6 = dot4.multiply(d6).subtract(dot3);
            interpolated_state       = current_state_linear_combination(coeff0, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6);
            interpolated_derivatives = derivative_linear_combination(coeff_dot0, coeff_dot_1, coeff_dot_2, coeff_dot_3, coeff_dot_4, coeff_dot5, coeff_dot6);
        }
        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

};