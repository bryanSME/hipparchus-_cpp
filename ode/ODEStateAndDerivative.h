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

/** Container for time, main and secondary state vectors as well as their derivatives.

 * @see Ordinary_Differential_Equation
 * @see Secondary_ODE
 * @see ODE_Integrator
 */

class ODE_State_And_Derivative extends ODE_State 
{

    /** Serializable UID. */
    private static const long serial_version_uid = 20160408L;

    /** Derivative of the primary state at time. */
    private const std::vector<double> primary_derivative;

    /** Derivative of the secondary state at time. */
    private const std::vector<std::vector<double>> secondary_derivative;

    /** Simple constructor.
     * <p>Calling this constructor is equivalent to call {@link
     * #ODE_State_And_Derivative(double, std::vector<double>, std::vector<double>, * std::vector<std::vector<double>>, std::vector<std::vector<double>>) ODE_State_And_Derivative(time, state, * derivative, NULL, NULL)}.</p>
     * @param time time
     * @param primary_state primary state at time
     * @param primary_derivative derivative of the primary state at time
     */
    public ODE_State_And_Derivative(double time, std::vector<double> primary_state, std::vector<double> primary_derivative) 
    {
        this(time, primary_state, primary_derivative, NULL, NULL);
    }

    /** Simple constructor.
     * @param time time
     * @param primary_state primary state at time
     * @param primary_derivative derivative of the primary state at time
     * @param secondary_state state at time (may be NULL)
     * @param secondary_derivative derivative of the state at time (may be NULL)
     */
    public ODE_State_And_Derivative(double time, std::vector<double> primary_state, std::vector<double> primary_derivative, std::vector<std::vector<double>> secondary_state, std::vector<std::vector<double>> secondary_derivative) 
    {
        super(time, primary_state, secondary_state);
        this.primary_derivative   = primary_derivative.clone();
        this.secondary_derivative = copy(secondary_derivative);
    }

    /** Get derivative of the primary state at time.
     * @return derivative of the primary state at time
     * @see #get_secondary_derivativestatic_cast<int>(
     * @see #get_complete_derivative()
     */
    public std::vector<double> get_primary_derivative() 
    {
        return primary_derivative.clone();
    }

    /** Get derivative of the secondary state at time.
     * @param index index of the secondary set as returned
     * by {@link Expandable_ODE#add_secondary_equations(Secondary_ODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return derivative of the secondary state at time
     * @see #get_primary_derivative()
     * @see #get_complete_derivative()
     */
    public std::vector<double> get_secondary_derivative(const int index) 
    {
        return index == 0 ? primary_derivative.clone() : secondary_derivative[index - 1].clone();
    }

    /** Get complete derivative at time.
     * @return complete derivative at time, starting with
     * {@link #get_primary_derivative() primary derivative}, followed
     * by all {@link #get_secondary_derivativestatic_cast<int>( secondary derivatives} in
     * increasing index order
     * @see #get_primary_derivative()
     * @see #get_secondary_derivativestatic_cast<int>(
     */
    public std::vector<double> get_complete_derivative() 
    {
        const std::vector<double> complete_derivative = std::vector<double>(get_complete_state_dimension()];
        System.arraycopy(primary_derivative, 0, complete_derivative, 0, primary_derivative.size());
        int offset = primary_derivative.size();
        if (secondary_derivative != NULL) 
        {
            for (int index = 0; index < secondary_derivative.size(); ++index) 
            {
                System.arraycopy(secondary_derivative[index], 0, complete_derivative, offset, secondary_derivative[index].size());
                offset += secondary_derivative[index].size();
            }
        }
        return complete_derivative;
    }

}


