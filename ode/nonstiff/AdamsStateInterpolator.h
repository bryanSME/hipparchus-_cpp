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

//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.ode.Equations_mapper;
//import org.hipparchus.ode.ODE_State_And_Derivative;
//import org.hipparchus.ode.sampling.AbstractODE_StateInterpolator;
//import org.hipparchus.util.FastMath;

/**
 * This class : an interpolator for integrators using Nordsieck representation.
 *
 * <p>This interpolator computes dense output around the current point.
 * The interpolation equation is based on Taylor series formulas.
 *
 * @see org.hipparchus.ode.nonstiff.AdamsBashforth_integrator
 * @see org.hipparchus.ode.nonstiff.Adams_moultonIntegrator
 */

class Adams_State_Interpolator extends AbstractODE_StateInterpolator 
{

    /** Serializable version identifier */
    private static const long serial_version_uid = 20160402L;

    /** Step size used in the first scaled derivative and Nordsieck vector. */
    private double scaling_h;

    /** Reference state.
     * <p>Sometimes, the reference state is the same as global_previous_state, * sometimes it is the same as global_current_state, so we use a separate
     * field to avoid any confusion.
     * </p>
     */
    private const ODE_State_And_Derivative reference;

    /** First scaled derivative. */
    private std::vector<double> scaled;

    /** Nordsieck vector. */
    private Array_2D_Row_Real_Matrix nordsieck;

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
    Adams_State_Interpolator(const double step_size, const ODE_State_And_Derivative reference, const std::vector<double> scaled, const Array_2D_Row_Real_Matrix nordsieck, const bool is_forward, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const Equations_mapper equations_mapper) 
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
    private Adams_State_Interpolator(const double step_size, const ODE_State_And_Derivative reference, const std::vector<double> scaled, const Array_2D_Row_Real_Matrix nordsieck, const bool is_forward, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const ODE_State_And_Derivative soft_previous_state, const ODE_State_And_Derivative soft_current_state, const Equations_mapper equations_mapper) 
    {
        super(is_forward, global_previous_state, global_current_state, soft_previous_state, soft_current_state, equations_mapper);
        this.scaling_h  = step_size;
        this.reference = reference;
        this.scaled    = scaled.clone();
        this.nordsieck = Array_2D_Row_Real_Matrix(nordsieck.get_data(), false);
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
    protected Adams_State_Interpolator create(bool new_forward, ODE_State_And_Derivative new_global_previous_state, ODE_State_And_Derivative new_global_current_state, ODE_State_And_Derivative new_soft_previous_state, ODE_State_And_Derivative new_soft_current_state, Equations_mapper new_mapper) 
    {
        return Adams_State_Interpolator(scaling_h, reference, scaled, nordsieck, new_forward, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);

    }

    /** Get the first scaled derivative.
     * @return first scaled derivative
     */
    public std::vector<double> get_scaled() 
    {
        return scaled.clone();
    }

    /** Get the Nordsieck vector.
     * @return Nordsieck vector
     */
    public Array_2D_Row_Real_Matrix get_nordsieck() 
    {
        return nordsieck;
    }

    /** {@inherit_doc} */
    //override
    protected ODE_State_And_Derivative compute_interpolated_state_and_derivatives(const Equations_mapper equations_mapper, const double time, const double& theta, const double theta_h, const double one_minus_theta_h) 
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
     */
    public static ODE_State_And_Derivative taylor(const Equations_mapper equations_mapper, const ODE_State_And_Derivative reference, const double time, const double step_size, const std::vector<double> scaled, const Array_2D_Row_Real_Matrix nordsieck) 
    {

        const double x = time - reference.get_time();
        const double normalized_abscissa = x / step_size;

        std::vector<double> state_variation       = std::vector<double>(scaled.size()];
        std::vector<double> estimated_derivatives = std::vector<double>(scaled.size()];

        // apply Taylor formula from high order to low order, // for the sake of numerical accuracy
        const std::vector<std::vector<double>> n_data = nordsieck.get_data_ref();
        for (int i = n_data.size() - 1; i >= 0; --i) 
        {
            const int order = i + 2;
            const std::vector<double> n_data_i = n_data[i];
            const double power = std::pow(normalized_abscissa, order);
            for (int j{}; j < n_data_i.size(); ++j) 
            {
                const double d = n_data_i[j] * power;
                state_variation[j]          = state_variation[j] + d;
                estimated_derivatives[j] = estimated_derivatives[j] + d * order;
            }
        }

        std::vector<double> estimated_state = reference.get_complete_state();
        for (int j{}; j < state_variation.size(); ++j) 
        {
            state_variation[j]       = state_variation[j] + scaled[j] * normalized_abscissa;
            estimated_state[j]       = estimated_state[j] + state_variation[j];
            estimated_derivatives[j] = (estimated_derivatives[j] + scaled[j] * normalized_abscissa) / x;
        }

        return equations_mapper.map_state_and_derivative(time, estimated_state, estimated_derivatives);

    }

}


