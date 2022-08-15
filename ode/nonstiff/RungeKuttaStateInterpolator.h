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
//import org.hipparchus.ode.sampling.AbstractODE_StateInterpolator;

/** This class represents an interpolator over the last step during an
 * ODE integration for Runge-Kutta and embedded Runge-Kutta integrators.
 *
 * @see Runge_Kutta_Integrator
 * @see EmbeddedRunge_Kutta_Integrator
 *
 */

virtual class Runge_Kutta_State_Interpolator extends AbstractODE_StateInterpolator 
{

    /** Serializable UID. */
    private static const long serial_version_uid = 20160328L;

    /** Slopes at the intermediate points */
    protected std::vector<std::vector<double>> y_dot_k;

    /** Simple constructor.
     * @param forward integration direction indicator
     * @param y_dot_k slopes at the intermediate points
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param soft_previous_state start of the restricted step
     * @param soft_current_state end of the restricted step
     * @param mapper equations mapper for the all equations
     */
    protected Runge_Kutta_State_Interpolator(const bool forward, const std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const ODE_State_And_Derivative soft_previous_state, const ODE_State_And_Derivative soft_current_state, const Equations_mapper mapper) 
    {
        super(forward, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
        this.y_dot_k = std::vector<double>(y_dot_k.size()][];
        for (int i{}; i < y_dot_k.size(); ++i) 
        {
            this.y_dot_k[i] = y_dot_k[i].clone();
        }
    }

    /** {@inherit_doc} */
    //override
    protected Runge_Kutta_State_Interpolator create(bool new_forward, ODE_State_And_Derivative new_global_previous_state, ODE_State_And_Derivative new_global_current_state, ODE_State_And_Derivative new_soft_previous_state, ODE_State_And_Derivative new_soft_current_state, Equations_mapper new_mapper) 
    {
        return create(new_forward, y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** Create a instance.
     * @param new_forward integration direction indicator
     * @param new_y_dot_k slopes at the intermediate points
     * @param new_global_previous_state start of the global step
     * @param new_global_current_state end of the global step
     * @param new_soft_previous_state start of the restricted step
     * @param new_soft_current_state end of the restricted step
     * @param new_mapper equations mapper for the all equations
     * @return a instance
     */
    protected virtual Runge_Kutta_State_Interpolator create(bool new_forward, std::vector<std::vector<double>> new_y_dot_k, ODE_State_And_Derivative new_global_previous_state, ODE_State_And_Derivative new_global_current_state, ODE_State_And_Derivative new_soft_previous_state, ODE_State_And_Derivative new_soft_current_state, Equations_mapper new_mapper);

    /** Compute a state by linear combination added to previous state.
     * @param coefficients coefficients to apply to the method staged derivatives
     * @return combined state
     */
    protected const std::vector<double> previous_state_linear_combination(const double ... coefficients) 
    {
        return combine(get_global_previous_state().get_complete_state(), coefficients);
    }

    /** Compute a state by linear combination added to current state.
     * @param coefficients coefficients to apply to the method staged derivatives
     * @return combined state
     */
    protected std::vector<double> current_state_linear_combination(const double ... coefficients) 
    {
        return combine(get_global_current_state().get_complete_state(), coefficients);
    }

    /** Compute a state derivative by linear combination.
     * @param coefficients coefficients to apply to the method staged derivatives
     * @return combined state
     */
    protected std::vector<double> derivative_linear_combination(const double ... coefficients) 
    {
        return combine(std::vector<double>(y_dot_k[0].size()], coefficients);
    }

    /** Linearly combine arrays.
     * @param a array to add to
     * @param coefficients coefficients to apply to the method staged derivatives
     * @return a itself, as a conveniency for fluent API
     */
    private std::vector<double> combine(const std::vector<double> a, const double ... coefficients) 
    {
        for (int i{}; i < a.size(); ++i) 
        {
            for (int k{}; k < coefficients.size(); ++k) 
            {
                a[i] += coefficients[k] * y_dot_k[k][i];
            }
        }
        return a;
    }

}


