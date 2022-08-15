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
#include <vector>

/** Container for time, main and secondary state vectors.

 * @see Ordinary_Differential_Equation
 * @see Secondary_ODE
 * @see ODE_Integrator
 * @see ODE_State_And_Derivative
 */

class ODE_State
{
private:

    /** Time. */
    const double my_time;

    /** Primary state at time. */
    const std::vector<double> my_primary_state;

    /** Secondary state at time. */
    const std::vector<std::vector<double>> my_secondary_state;

    /** Complete dimension. */
    const int my_complete_dimension;

public:
    /** Simple constructor.
     * <p>Calling this constructor is equivalent to call {@link
     * #ODE_State(double, std::vector<double>, std::vector<std::vector<double>>)
     * ODE_State(time, state, NULL)}.</p>
     * @param time time
     * @param primary_state primary state at time
     */
    ODE_State(const double& time, const std::vector<double>& primary_state)
    {
        ODE_State(time, primary_state, std::vector<std::vector<double>>{});
    }

    /** Simple constructor.
     * @param time time
     * @param primary_state state at time
     * @param secondary_state primary state at time (may be NULL)
     */
    ODE_State(const double& time, const std::vector<double>& primary_state, const std::vector<std::vector<double>>& secondary_state)
        :
        my_time{ time },
        my_primary_state{ primary_state },
        my_secondary_state{ copy(secondary_state) }
    {
        // compute once and for all the complete dimension
        auto dimension = primary_state.size();

        for (const auto& secondary : secondary_state)
        {
            dimension += secondary.size();
        }
        my_complete_dimension = dimension;
    }

    /** Get time.
     * @return time
     */
    double get_time() const
    {
        return my_time;
    }

    /** Get primary state dimension.
     * @return primary state dimension
     * @see #get_secondary_state_dimensionstatic_cast<int>(
     * @see #get_complete_state_dimension()
     */
    int get_primary_state_dimension() const
    {
        return my_primary_state.size();
    }

    /** Get primary state at time.
     * @return primary state at time
     * @see #get_secondary_statestatic_cast<int>(
     * @see #get_complete_state()
     */
    std::vector<double> get_primary_state() const
    {
        return my_primary_state;
    }

    /** Get the number of secondary states.
     * @return number of secondary states.
     */
    int get_number_of_secondary_states() const
    {
        return my_secondary_state.size();
    }

    /** Get secondary state dimension.
     * @param index index of the secondary set as returned
     * by {@link Expandable_ODE#add_secondary_equations(Secondary_ODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return secondary state dimension
     * @see #get_primary_state_dimension()
     * @see #get_complete_state_dimension()
     */
    int get_secondary_state_dimension(const int& index) const
    {
        return index == 0
            ? my_primary_state.size()
            : my_secondary_state[index - 1].size();
    }

    /** Get secondary state at time.
     * @param index index of the secondary set as returned
     * by {@link Expandable_ODE#add_secondary_equations(Secondary_ODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return secondary state at time
     * @see #get_primary_state()
     * @see #get_complete_state()
     */
    std::vector<double> get_secondary_state(const int& index) const
    {
        return index == 0
            ? my_primary_state
            : my_secondary_state[index - 1];
    }

    /** Return the dimension of the complete set of equations.
     * <p>
     * The complete set of equations correspond to the primary set plus all secondary sets.
     * </p>
     * @return dimension of the complete set of equations
     * @see #get_primary_state_dimension()
     * @see #get_secondary_state_dimensionstatic_cast<int>(
     */
    int get_complete_state_dimension() const
    {
        return my_complete_dimension;
    }

    /** Get complete state at time.
     * @return complete state at time, starting with
     * {@link #get_primary_state() primary state}, followed
     * by all {@link #get_secondary_statestatic_cast<int>( secondary states} in
     * increasing index order
     * @see #get_primary_state()
     * @see #get_secondary_statestatic_cast<int>(
     */
    std::vector<double> get_complete_state()
    {
        const auto complete_state = std::vector<double>(get_complete_state_dimension()];
        System.arraycopy(my_primary_state, 0, complete_state, 0, my_primary_state.size());
        auto offset = my_primary_state.size();
        for (int index{}; index < my_secondary_state.size(); ++index)
        {
            System.arraycopy(my_secondary_state[index], 0, complete_state, offset, my_secondary_state[index].size());
            offset += my_secondary_state[index].size();
        }
        return complete_state;
    }

protected:
    /** Copy a two-dimensions array.
     * @param original original array (may be NULL)
     * @return copied array or NULL if original array was NULL
     */
    std::vector<std::vector<double>> copy(const std::vector<std::vector<double>> original)
    {
        // allocate the array
        auto copied = std::vector<std::vector<double>>(original.size());

        // copy content
        for (int i{}; i < original.size(); ++i)
        {
            copied[i] = original[i];
        }

        return copied;
    }
};