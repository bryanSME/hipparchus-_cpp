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
 * @see Midpoint_Integrator
 * @see Classical_Runge_Kutta_Integrator
 * @see Gill_Integrator
 * @see Three_Eighthes_Integrator
 * @see Luther_Integrator
 */

class Euler_Integrator extends Runge_Kutta_Integrator 
{

    /** Simple constructor.
     * Build an Euler integrator with the given step.
     * @param step integration step
     */
    public Euler_Integrator(const double step) 
    {
        super("Euler", step);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_c() 
    {
        return std::vector<double>(0];
    }

    /** {@inherit_doc} */
    //override
    public std::vector<std::vector<double>> get_a() 
    {
        return std::vector<double>(0][];
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_b() 
    {
        return std::vector<double> { 1 };
    }

    /** {@inherit_doc} */
    //override
    protected Euler_State_Interpolator
        create_interpolator(const bool forward, std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const Equations_mapper mapper) 
        {
        return Euler_State_Interpolator(forward, y_dot_k, global_previous_state, global_current_state, global_previous_state, global_current_state, mapper);
    }

}


