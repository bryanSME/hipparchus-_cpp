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

//package org.hipparchus.ode;

//import org.hipparchus.complex.std::complex<double>;

/** This class converts {@link std::complex<double>Ordinary_Differential_Equation complex Ordinary
 * Differential Equations} into {@link Ordinary_Differential_Equation real ones}.
 *
 * <p>This class is a wrapper around a {@link std::complex<double>Ordinary_Differential_Equation} which
 * allow to use a {@link ODE_Integrator} to integrate it.</p>
 *
 * <p>The transformation is done by changing the n dimension state
 * vector to a 2n dimension vector, where the even components are
 * real parts and odd components are imaginary parts.</p>
 *
 * <p>One should be aware that the data is duplicated during the
 * transformation process and that for each call to {@link
 * Ordinary_Differential_Equation#compute_derivatives(double, std::vector<double>)
 * compute_derivatives}, this wrapper does copy 4n scalars : 2n before
 * the call to {@link
 * Ordinary_Differential_Equation#compute_derivatives(double, std::vector<double>)
 * compute_derivatives} in order to dispatch the y state vector, * and 2n after the call to gather z_dot. sin_ce the underlying problem
 * by itself perhaps also needs to copy data and dispatch the arrays
 * into domain objects, this has an impact on both memory and CPU usage.
 * The only way to avoid this duplication is to perform the transformation
 * at the problem level, i.e. to implement the problem as a first order one
 * and then avoid using this class.</p>
 *
 * <p>
 * The proper way to use the converter is as follows:
 * </p>
 * <pre>
 *   ODE_Integrator                       integrator       = ...build some integrator...;
 *   std::complex<double>Ordinary_Differential_Equation complex_equations = ...set up the complex problem...;
 *   std::complex<double>ODE_State                     initial_state     = ...set up initial state...;
 *   std::complex<double>_ODE_Converter                 converter        = std::complex<double>_ODE_Converter();
 *   std::complex<double>ODE_State_And_Derivative        conststate       =
 *      converter.convert_state_and_derivative(integrator.integrate(converter.convert_equations(complex_equations), *                                                               converter.convert_state(initial_state), *                                                               t);
 * </pre>
 * <p>
 * If there are {@link std::complex<double>Secondary_ODE complex secondary equations}, they must be converted
 * too and both the converted primary equations and converted secondary equations must be
 * combined together using {@link Expandable_ODE Expandable_ODE} as usual for regular real equations.
 * </p>
 *
 * @see std::complex<double>Ordinary_Differential_Equation
 * @see Ordinary_Differential_Equation
 * @since 1.4
 */

class std::complex<double>_ODE_Converter 
{

    /** Convert an equations set.
     * @param equations equations to convert
     * @return converted equations
     */
    public Ordinary_Differential_Equation convert_equations(const std::complex<double>Ordinary_Differential_Equation equations) 
    {
        return Ordinary_Differential_Equation() 
        {

            /** {@inherit_doc}
             * <p>The dimension of the real problem is twice the
             * dimension of the underlying complex problem.</p>
             * @return dimension of the problem
             */
            //override
            public int get_dimension() 
            {
                return 2 * equations.get_dimension();
            }

            /** {@inherit_doc} */
            //override
            public void init(const double t0, const std::vector<double> y0, const double const_time) 
            {
                equations.init(t0, convert(y0), const_time);
            }

            /** {@inherit_doc} */
            //override
            public std::vector<double> compute_derivatives(const double t, const std::vector<double> y) 
            {
                return convert(equations.compute_derivatives(t, convert(y)));
            }

        };
    }

    /** Convert a secondary equations set.
     * @param equations equations to convert
     * @return converted equations
     */
    public Secondary_ODE convert_secondary_equations(const std::complex<double>Secondary_ODE equations) 
    {
        return Secondary_ODE() 
        {

            /** {@inherit_doc}
             * <p>The dimension of the real problem is twice the
             * dimension of the underlying complex problem.</p>
             * @return dimension of the problem
             */
            //override
            public int get_dimension() 
            {
                return 2 * equations.get_dimension();
            }

            /** {@inherit_doc} */
            //override
            public void init(const double t0, const std::vector<double> primary0, const std::vector<double> secondary0, const double const_time) 
            {
                equations.init(t0, convert(primary0), convert(secondary0), const_time);
            }

            /** {@inherit_doc} */
            //override
            public std::vector<double> compute_derivatives(const double t, const std::vector<double> primary, const std::vector<double> primary_dot, const std::vector<double> secondary) 
            {
                return convert(equations.compute_derivatives(t, convert(primary), convert(primary_dot), convert(secondary)));
            }

        };
    }

    /** Convert a complex state (typically the initial state).
     * @param state state to convert
     * @return converted state
     */
    public ODE_State convert_state(const std::complex<double>ODE_State state) 
    {
        const std::vector<std::vector<double>> secondary = std::vector<double>(state.get_number_of_secondary_states()][];
        for (int index = 0; index < secondary.size(); ++index) 
        {
            secondary[index] = convert(state.get_secondary_state(index + 1));
        }
        return ODE_State(state.get_time(), convert(state.get_primary_state()), secondary);
    }

    /** Convert a real state and derivatives (typically the const state or some intermediate state for
     * step handling or event handling).
     * @param state state to convert
     * @return converted state
     */
    public std::complex<double>ODE_State_And_Derivative convert_state(const ODE_State_And_Derivative state) 
    {
        const std::complex<double>[][] secondary           = std::complex<double>[state.get_number_of_secondary_states()][];
        const std::complex<double>[][] secondary_derivative = std::complex<double>[state.get_number_of_secondary_states()][];
        for (int index = 0; index < secondary.size(); ++index) 
        {
            secondary[index]           = convert(state.get_secondary_state(index + 1));
            secondary_derivative[index] = convert(state.get_secondary_derivative(index + 1));
        }
        return std::complex<double>ODE_State_And_Derivative(state.get_time(), convert(state.get_primary_state()), convert(state.get_primary_derivative()), secondary, secondary_derivative);
    }

    /** Convert a real array into a complex array.
     * @param a array to convert
     * @return converted array
     */
    private std::vector<std::complex<double>>convert(const std::vector<double> a) 
    {
        const std::vector<std::complex<double>>converted = std::complex<double>[a.size() / 2];
        for (int i{}; i < converted.size(); ++i) 
        {
            converted[i] = std::complex<double>(a[2 * i], a[2 * i + 1]);
        }
        return converted;
    }

    /** Convert a complex array into a real array.
     * @param a array to convert
     * @return converted array
     */
    private std::vector<double> convert(const std::vector<std::complex<double>>a) 
    {
        const std::vector<double> converted = std::vector<double>(a.size() * 2];
        for (int i{}; i < a.size(); ++i) 
        {
            converted[2 * i]     = a[i].get_real();
            converted[2 * i + 1] = a[i].get_imaginary();
        }
        return converted;
    }

}


