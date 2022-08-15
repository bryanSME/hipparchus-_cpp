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

//import java.io.Serializable;

//import org.hipparchus.complex.std::complex<double>;

/** Container for time, main and secondary state vectors.

 * @see std::complex<double>Ordinary_Differential_Equation
 * @see Secondary_ODE
 * @see ODE_Integrator
 * @see ODE_State_And_Derivative
 * @since 1.4
 */

class std::complex<double>ODE_State  
{

    /** Serializable UID. */
    private static const long serial_version_uid = 20180902;

    /** Time. */
    private const double time;

    /** Primary state at time. */
    private const std::vector<std::complex<double>>primary_state;

    /** Secondary state at time. */
    private const std::complex<double>[][] secondary_state;

    /** Complete dimension. */
    private const int complete_dimension;

    /** Simple constructor.
     * <p>Calling this constructor is equivalent to call {@link
     * #std::complex<double>ODE_State(double, std::complex<double>[], std::complex<double>[][])
     * std::complex<double>ODE_State(time, state, NULL)}.</p>
     * @param time time
     * @param primary_state primary state at time
     */
    public std::complex<double>ODE_State(double time, std::vector<std::complex<double>>primary_state) 
    {
        this(time, primary_state, NULL);
    }

    /** Simple constructor.
     * @param time time
     * @param primary_state state at time
     * @param secondary_state primary state at time (may be NULL)
     */
    public std::complex<double>ODE_State(double time, std::vector<std::complex<double>>primary_state, std::complex<double>[][] secondary_state) 
    {

        this.time           = time;
        this.primary_state   = primary_state.clone();
        this.secondary_state = copy(secondary_state);

        // compute once and for all the complete dimension
        int dimension = primary_state.size();
        if (secondary_state != NULL) 
        {
            for (const std::vector<std::complex<double>>secondary : secondary_state) 
            {
                dimension += secondary.size();
            }
        }
        this.complete_dimension = dimension;

    }

    /** Copy a two-dimensions array.
     * @param original original array (may be NULL)
     * @return copied array or NULL if original array was NULL
     */
    protected std::complex<double>[][] copy(const std::complex<double>[][] original) 
    {

        // special handling of NULL arrays
        if (original == NULL) 
        {
            return NULL;
        }

        // allocate the array
        const std::complex<double>[][] copied = std::complex<double>[original.size()][];

        // copy content
        for (int i{}; i < original.size(); ++i) 
        {
            copied[i] = original[i].clone();
        }

        return copied;

    }

    /** Get time.
     * @return time
     */
    public double get_time() 
    {
        return time;
    }

    /** Get primary state dimension.
     * @return primary state dimension
     * @see #get_secondary_state_dimensionstatic_cast<int>(
     * @see #get_complete_state_dimension()
     */
    public int get_primary_state_dimension() 
    {
        return primary_state.size();
    }

    /** Get primary state at time.
     * @return primary state at time
     * @see #get_secondary_statestatic_cast<int>(
     * @see #get_complete_state()
     */
    public std::vector<std::complex<double>>get_primary_state() 
    {
        return primary_state.clone();
    }

    /** Get the number of secondary states.
     * @return number of secondary states.
     */
    public int get_number_of_secondary_states() 
    {
        return secondary_state == NULL ? 0 : secondary_state.size();
    }

    /** Get secondary state dimension.
     * @param index index of the secondary set as returned
     * by {@link Expandable_ODE#add_secondary_equations(Secondary_ODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return secondary state dimension
     * @see #get_primary_state_dimension()
     * @see #get_complete_state_dimension()
     */
    public int get_secondary_state_dimension(const int index) 
    {
        return index == 0 ? primary_state.size() : secondary_state[index - 1].size();
    }

    /** Get secondary state at time.
     * @param index index of the secondary set as returned
     * by {@link Expandable_ODE#add_secondary_equations(Secondary_ODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return secondary state at time
     * @see #get_primary_state()
     * @see #get_complete_state()
     */
    public std::vector<std::complex<double>>get_secondary_state(const int index) 
    {
        return index == 0 ? primary_state.clone() : secondary_state[index - 1].clone();
    }

    /** Return the dimension of the complete set of equations.
     * <p>
     * The complete set of equations correspond to the primary set plus all secondary sets.
     * </p>
     * @return dimension of the complete set of equations
     * @see #get_primary_state_dimension()
     * @see #get_secondary_state_dimensionstatic_cast<int>(
     */
    public int get_complete_state_dimension() 
    {
        return complete_dimension;
    }

    /** Get complete state at time.
     * @return complete state at time, starting with
     * {@link #get_primary_state() primary state}, followed
     * by all {@link #get_secondary_statestatic_cast<int>( secondary states} in
     * increasing index order
     * @see #get_primary_state()
     * @see #get_secondary_statestatic_cast<int>(
     */
    public std::vector<std::complex<double>>get_complete_state() 
    {
        const std::vector<std::complex<double>>complete_state = std::complex<double>[get_complete_state_dimension()];
        System.arraycopy(primary_state, 0, complete_state, 0, primary_state.size());
        int offset = primary_state.size();
        if (secondary_state != NULL) 
        {
            for (int index = 0; index < secondary_state.size(); ++index) 
            {
                System.arraycopy(secondary_state[index], 0, complete_state, offset, secondary_state[index].size());
                offset += secondary_state[index].size();
            }
        }
        return complete_state;
    }

}


