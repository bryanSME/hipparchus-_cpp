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

//import java.util.Arrays;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.linear.Array2DRowField_Matrix;
//import org.hipparchus.ode.FieldEquations_mapper;
//import org.hipparchus.ode.Field_ODE_State_And_Derivative;
//import org.hipparchus.ode.sampling.AbstractFieldODE_StateInterpolator;
//import org.hipparchus.util.Math_Arrays;
#include<type_traits>
#include "../../core/CalculusFieldElement.h"

/**
 * This class : an interpolator for Adams integrators using Nordsieck representation.
 *
 * <p>This interpolator computes dense output around the current point.
 * The interpolation equation is based on Taylor series formulas.
 *
 * @see Adams_Bashforth_Field_Integrator
 * @see Adams_moultonFieldIntegrator
 * @param <T> the type of the field elements
 */

template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Adams_Field_State_Interpolator : public AbstractFieldODE_StateInterpolator<T> 
{

    /** Step size used in the first scaled derivative and Nordsieck vector. */
    private T scaling_h;

    /** Reference state.
     * <p>Sometimes, the reference state is the same as global_previous_state, * sometimes it is the same as global_current_state, so we use a separate
     * field to avoid any confusion.
     * </p>
     */
    private const Field_ODE_State_And_Derivative<T> reference;

    /** First scaled derivative. */
    private const std::vector<T> scaled;

    /** Nordsieck vector. */
    private const Array2DRowField_Matrix<T> nordsieck;

    /** Simple constructor.
     * @param step_size step size used in the scaled and Nordsieck arrays
     * @param reference reference state from which Taylor expansion are estimated
     * @param scaled first scaled derivative
     * @param nordsieck Nordsieck vector
     * @param is_forward integration direction indicator
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param equations_mapper mapper for ODE equations primary and secondary components
     */
    Adams_Field_State_Interpolator(const T step_size, const Field_ODE_State_And_Derivative<T> reference, const std::vector<T> scaled, const Array2DRowField_Matrix<T> nordsieck, const bool is_forward, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const FieldEquations_mapper<T> equations_mapper) 
    {
        this(step_size, reference, scaled, nordsieck, is_forward, global_previous_state, global_current_state, global_previous_state, global_current_state, equations_mapper);
    }

    /** Simple constructor.
     * @param step_size step size used in the scaled and Nordsieck arrays
     * @param reference reference state from which Taylor expansion are estimated
     * @param scaled first scaled derivative
     * @param nordsieck Nordsieck vector
     * @param is_forward integration direction indicator
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param soft_previous_state start of the restricted step
     * @param soft_current_state end of the restricted step
     * @param equations_mapper mapper for ODE equations primary and secondary components
     */
    private Adams_Field_State_Interpolator(const T step_size, const Field_ODE_State_And_Derivative<T> reference, const std::vector<T> scaled, const Array2DRowField_Matrix<T> nordsieck, const bool is_forward, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const Field_ODE_State_And_Derivative<T> soft_previous_state, const Field_ODE_State_And_Derivative<T> soft_current_state, const FieldEquations_mapper<T> equations_mapper) 
    {
        super(is_forward, global_previous_state, global_current_state, soft_previous_state, soft_current_state, equations_mapper);
        this.scaling_h  = step_size;
        this.reference = reference;
        this.scaled    = scaled.clone();
        this.nordsieck = Array2DRowField_Matrix<>(nordsieck.get_data(), false);
    }

    /** Create a instance.
     * @param new_forward integration direction indicator
     * @param new_global_previous_state start of the global step
     * @param new_global_current_state end of the global step
     * @param new_soft_previous_state start of the restricted step
     * @param new_soft_current_state end of the restricted step
     * @param new_mapper equations mapper for the all equations
     * @return a instance
     */
    //override
    protected Adams_Field_State_Interpolator<T> create(bool new_forward, Field_ODE_State_And_Derivative<T> new_global_previous_state, Field_ODE_State_And_Derivative<T> new_global_current_state, Field_ODE_State_And_Derivative<T> new_soft_previous_state, Field_ODE_State_And_Derivative<T> new_soft_current_state, FieldEquations_mapper<T> new_mapper) 
    {
        return Adams_Field_State_Interpolator<T>(scaling_h, reference, scaled, nordsieck, new_forward, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);

    }

    /** Get the first scaled derivative.
     * @return first scaled derivative
     */
    public std::vector<T> get_scaled() 
    {
        return scaled.clone();
    }

    /** Get the Nordsieck vector.
     * @return Nordsieck vector
     */
    public Array2DRowField_Matrix<T> get_nordsieck() 
    {
        return nordsieck;
    }

    /** {@inherit_doc} */
    //override
    protected Field_ODE_State_And_Derivative<T> compute_interpolated_state_and_derivatives(const FieldEquations_mapper<T> equations_mapper, const T time, const T theta, const T theta_h, const T one_minus_theta_h) 
    {
        return taylor(equations_mapper, reference, time, scaling_h, scaled, nordsieck);
    }

    /** Estimate state by applying Taylor formula.
     * @param equations_mapper mapper for ODE equations primary and secondary components
     * @param reference reference state
     * @param time time at which state must be estimated
     * @param step_size step size used in the scaled and Nordsieck arrays
     * @param scaled first scaled derivative
     * @param nordsieck Nordsieck vector
     * @return estimated state
     * @param <S> the type of the field elements
     */
    public static <S extends Calculus_Field_Element<S>> Field_ODE_State_And_Derivative<S> taylor(const FieldEquations_mapper<S> equations_mapper, const Field_ODE_State_And_Derivative<S> reference, const S time, const S step_size, const S[] scaled, const Array2DRowField_Matrix<S> nordsieck) 
    {

        const S x = time.subtract(reference.get_time());
        const S normalized_abscissa = x.divide(step_size);

        S[] state_variation = Math_Arrays::build_array(time.get_field(), scaled.size());
        Arrays.fill(state_variation, time.get_field().get_zero());
        S[] estimated_derivatives = Math_Arrays::build_array(time.get_field(), scaled.size());
        Arrays.fill(estimated_derivatives, time.get_field().get_zero());

        // apply Taylor formula from high order to low order, // for the sake of numerical accuracy
        const S[][] n_data = nordsieck.get_data_ref();
        for (int i = n_data.size() - 1; i >= 0; --i) 
        {
            const int order = i + 2;
            const S[] n_data_i = n_data[i];
            const S power = normalized_abscissa.pow(order);
            for (int j{}; j < n_data_i.size(); ++j) 
            {
                const S d = n_data_i[j].multiply(power);
                state_variation[j]          = state_variation[j].add(d);
                estimated_derivatives[j] = estimated_derivatives[j].add(d.multiply(order));
            }
        }

        S[] estimated_state = reference.get_complete_state();
        for (int j{}; j < state_variation.size(); ++j) 
        {
            state_variation[j] = state_variation[j].add(scaled[j].multiply(normalized_abscissa));
            estimated_state[j] = estimated_state[j].add(state_variation[j]);
            estimated_derivatives[j] =
                estimated_derivatives[j].add(scaled[j].multiply(normalized_abscissa)).divide(x);
        }

        return equations_mapper.map_state_and_derivative(time, estimated_state, estimated_derivatives);

    }

};