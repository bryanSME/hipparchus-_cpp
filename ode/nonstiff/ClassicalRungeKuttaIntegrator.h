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
 * @see Euler_Integrator
 * @see Gill_Integrator
 * @see Midpoint_Integrator
 * @see Three_Eighthes_Integrator
 * @see Luther_Integrator
 */

class Classical_Runge_Kutta_Integrator extends Runge_Kutta_Integrator 
{

    /** Simple constructor.
     * Build a fourth-order Runge-Kutta integrator with the given
     * step.
     * @param step integration step
     */
    public Classical_Runge_Kutta_Integrator(const double step) 
    {
        super("classical Runge-Kutta", step);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_c() 
    {
        return std::vector<double> 
        {
            1.0 / 2.0, 1.0 / 2.0, 1.0
        };
    }

    /** {@inherit_doc} */
    //override
    public std::vector<std::vector<double>> get_a() 
    {
        return std::vector<std::vector<double>> 
        {
            { 1.0 / 2.0 }, { 0.0, 1.0 / 2.0 }, { 0.0, 0.0, 1.0 }
        };
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_b() 
    {
        return std::vector<double> 
        {
            1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0
        };
    }

    /** {@inherit_doc} */
    //override
    protected ClassicalRunge_Kutta_State_Interpolator
    create_interpolator(const bool forward, std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const Equations_mapper mapper) 
    {
        return ClassicalRunge_Kutta_State_Interpolator(forward, y_dot_k, global_previous_state, global_current_state, global_previous_state, global_current_state, mapper);
    }

}


