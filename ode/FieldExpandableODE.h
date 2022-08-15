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

//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.Math_Arrays;
#include <type_traits>
#include "../core/CalculusFieldElement.h"


/**
 * This class represents a combined set of first order differential equations, * with at least a primary set of equations expandable by some sets of secondary
 * equations.
 * <p>
 * One typical use case is the computation of the Jacobian matrix for some ODE.
 * In this case, the primary set of equations corresponds to the raw ODE, and we
 * add to this set another bunch of secondary equations which represent the Jacobian
 * matrix of the primary set.
 * </p>
 * <p>
 * We want the integrator to use <em>only</em> the primary set to estimate the
 * errors and hence the step sizes. It should <em>not</em> use the secondary
 * equations in this computation. The {@link FieldODE_Integrator integrator} will
 * be able to know where the primary set ends and so where the secondary sets begin.
 * </p>
 *
 * @see FieldOrdinary_Differential_Equation
 * @see FieldSecondary_ODE
 *
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class FieldExpandable_ODE 
{

    /** Primary differential equation. */
    private const FieldOrdinary_Differential_Equation<T> primary;

    /** Components of the expandable ODE. */
    private List<FieldSecondary_ODE<T>> components;

    /** Mapper for all equations. */
    private FieldEquations_mapper<T> mapper;

    /** Build an expandable set from its primary ODE set.
     * @param primary the primary set of differential equations to be integrated.
     */
    public FieldExpandable_ODE(const FieldOrdinary_Differential_Equation<T> primary) 
    {
        this.primary    = primary;
        this.components = Array_list<>();
        this.mapper     = FieldEquations_mapper<>(null, primary.get_dimension());
    }

    /** Get the primary set of differential equations to be integrated.
     * @return primary set of differential equations to be integrated
     * @since 2.2
     */
    public FieldOrdinary_Differential_Equation<T> get_primary() 
    {
        return primary;
    }

    /** Get the mapper for the set of equations.
     * @return mapper for the set of equations
     */
    public FieldEquations_mapper<T> get_mapper() 
    {
        return mapper;
    }

    /** Add a set of secondary equations to be integrated along with the primary set.
     * @param secondary secondary equations set
     * @return index of the secondary equation in the expanded state, to be used
     * as the parameter to {@link FieldODE_State#get_secondary_statestatic_cast<int>(} and
     * {@link Field_ODE_State_And_Derivative#get_secondary_derivativestatic_cast<int>(} (beware index
     * 0 corresponds to primary state, secondary states start at 1)
     */
    public int add_secondary_equations(const FieldSecondary_ODE<T> secondary) 
    {

        components.add(secondary);
        mapper = FieldEquations_mapper<>(mapper, secondary.get_dimension());

        return components.size();

    }

    /** Initialize equations at the start of an ODE integration.
     * @param s0 state at integration start
     * @param const_time target time for the integration
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if arrays dimensions do not match equations settings
     */
    public void init(const FieldODE_State<T> s0, const T const_time) 
    {

        const T t0 = s0.get_time();

        // initialize primary equations
        const std::vector<T> primary0 = s0.get_primary_state();
        primary.init(t0, primary0, const_time);

        // initialize secondary equations
        for (const int& index = 1; index < mapper.get_number_of_equations(); ++index) 
        {
            const std::vector<T> secondary0 = s0.get_secondary_state(index);
            components.get(index - 1).init(t0, primary0, secondary0, const_time);
        }

    }

    /** Get the current time derivative of the complete state vector.
     * @param t current value of the independent <I>time</I> variable
     * @param y array containing the current value of the complete state vector
     * @return time derivative of the complete state vector
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if arrays dimensions do not match equations settings
     */
    public std::vector<T> compute_derivatives(const T t, const std::vector<T> y)
        , Math_Illegal_State_Exception 
        {

        const std::vector<T> y_dot = Math_Arrays::build_array(t.get_field(), mapper.get_total_dimension());

        // compute derivatives of the primary equations
        const std::vector<T> primary_state    = mapper.extract_equation_data(0, y);
        const std::vector<T> primary_state_dot = primary.compute_derivatives(t, primary_state);

        // Add contribution for secondary equations
        for (const int& index = 1; index < mapper.get_number_of_equations(); ++index) 
        {
            const std::vector<T> component_state    = mapper.extract_equation_data(index, y);
            const std::vector<T> component_state_dot = components.get(index - 1).compute_derivatives(t, primary_state, primary_state_dot, component_state);
            mapper.insert_equation_data(index, component_state_dot, y_dot);
        }

        // we retrieve the primary_state_dot array after the secondary equations have
        // been computed in case they change the main state derivatives; this happens
        // for example in optimal control when the secondary equations handle co-state, // which changes control, and the control changes the primary state
        mapper.insert_equation_data(0, primary_state_dot, y_dot);

        return y_dot;

    }

};