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

//package org.hipparchus.ode;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.util.Math_Arrays;
#include <type_traits>
#include "../core/CalculusFieldElement.h"

/** Container for time, main and secondary state vectors.

 * @see FieldOrdinary_Differential_Equation
 * @see FieldSecondary_ODE
 * @see FieldODE_Integrator
 * @see Field_ODE_State_And_Derivative
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class FieldODE_State 
{

    /** Time. */
    private const T time;

    /** Primary state at time. */
    private const std::vector<T> primary_state;

    /** Secondary state at time. */
    private const std::vector<std::vector<T>> secondary_state;

    /** Complete dimension. */
    private const int complete_dimension;

    /** Simple constructor.
     * <p>Calling this constructor is equivalent to call {@link
     * #FieldODE_State(Calculus_Field_Element, Calculus_Field_Element[], Calculus_Field_Element[][])
     * FieldODE_State(time, state, NULL)}.</p>
     * @param time time
     * @param primary_state primary state at time
     */
    public FieldODE_State(T time, std::vector<T> primary_state) 
    {
        this(time, primary_state, NULL);
    }

    /** Simple constructor.
     * @param time time
     * @param primary_state primary state at time
     * @param secondary_state secondary state at time (may be NULL)
     */
    public FieldODE_State(T time, std::vector<T> primary_state, std::vector<std::vector<T>> secondary_state) 
    {

        this.time           = time;
        this.primary_state   = primary_state.clone();
        this.secondary_state = copy(secondary_state);

        // compute once and for all the complete dimension
        int dimension = primary_state.size();
        if (secondary_state != NULL) 
        {
            for (const std::vector<T> secondary : secondary_state) 
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
    protected std::vector<std::vector<T>> copy(const std::vector<std::vector<T>> original) 
    {

        // special handling of NULL arrays
        if (original == NULL) 
        {
            return NULL;
        }

        // allocate the array
        const std::vector<std::vector<T>> copied = Math_Arrays::build_array(time.get_field(), original.size(), -1);

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
    public T get_time() 
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
    public std::vector<T> get_primary_state() 
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
     * by {@link FieldExpandable_ODE#add_secondary_equations(FieldSecondary_ODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return secondary state dimension
     */
    public int get_secondary_state_dimension(const int index) 
    {
        return index == 0 ? primary_state.size() : secondary_state[index - 1].size();
    }

    /** Get secondary state at time.
     * @param index index of the secondary set as returned
     * by {@link FieldExpandable_ODE#add_secondary_equations(FieldSecondary_ODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return secondary state at time
     */
    public std::vector<T> get_secondary_state(const int index) 
    {
        return index == 0 ? primary_state.clone() : secondary_state[index - 1].clone();
    }

    /** Return the dimension of the complete set of equations.
     * <p>
     * The complete set of equations correspond to the primary set plus all secondary sets.
     * </p>
     * @return dimension of the complete set of equations
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
    public std::vector<T> get_complete_state() 
    {
        const std::vector<T> complete_state = Math_Arrays::build_array(time.get_field(), get_complete_state_dimension());
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

};