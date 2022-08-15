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
 * This class : the 3/8 fourth order Runge-Kutta
 * integrator for Ordinary Differential Equations.
 *
 * <p>This method is an explicit Runge-Kutta method, its Butcher-array
 * is the following one :
 * <pre>
 *    0  |  0    0    0    0
 *   1/3 | 1/3   0    0    0
 *   2/3 |-1/3   1    0    0
 *    1  |  1   -1    1    0
 *       |--------------------
 *       | 1/8  3/8  3/8  1/8
 * </pre>
 * </p>
 *
 * @see Euler_Integrator
 * @see Classical_Runge_Kutta_Integrator
 * @see Gill_Integrator
 * @see Midpoint_Integrator
 * @see Luther_Integrator
 */

class Three_Eighthes_Integrator extends Runge_Kutta_Integrator 
{

    /** Simple constructor.
     * Build a 3/8 integrator with the given step.
     * @param step integration step
     */
    public Three_Eighthes_Integrator(const double step) 
    {
        super("3/8", step);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_c() 
    {
        return std::vector<double> 
        {
            1.0 / 3.0, 2.0 / 3.0, 1.0
        };
    }

    /** {@inherit_doc} */
    //override
    public std::vector<std::vector<double>> get_a() 
    {
        return std::vector<std::vector<double>> 
        {
            {  1.0 / 3.0 }, { -1.0 / 3.0, 1.0 }, {  1.0, -1.0, 1.0 }
        };
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_b() 
    {
        return std::vector<double> 
        {
            1.0 / 8.0, 3.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0
        };
    }

    /** {@inherit_doc} */
    //override
    protected Three_Eighthes_State_Interpolator
    create_interpolator(const bool forward, std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const Equations_mapper mapper) 
    {
        return Three_Eighthes_State_Interpolator(forward, y_dot_k, global_previous_state, global_current_state, global_previous_state, global_current_state, mapper);
    }

}


