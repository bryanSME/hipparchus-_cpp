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

//package org.hipparchus.ode.sampling;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.ode.Field_ODE_State_And_Derivative;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"

/** This interface represents an interpolator over the last step
 * during an ODE integration.
 *
 * <p>The various ODE integrators provide objects implementing this
 * interface to the step handlers. These objects are often custom
 * objects tightly bound to the integrator internal algorithms. The
 * handlers can use these objects to retrieve the state vector at
 * intermediate times between the previous and the current grid points
 * (this feature is often called dense output).</p>
 *
 * @param <T> the type of the field elements
 * @see org.hipparchus.ode.FieldODE_Integrator
 * @see FieldODE_Step_Handler
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class FieldODE_StateInterpolator
{

    /**
     * Get the state at previous grid point time.
     * @return state at previous grid point time
     */
    Field_ODE_State_And_Derivative<T> get_previous_state();

    /**
     * Determines if the {@link #get_previous_state() previous state} is computed directly
     * by the integrator, or if it is calculated using {@link #get_interpolated_state(Calculus_Field_Element)
     * interpolation}.
     *
     * <p> Typically the previous state is directly computed by the integrator, but when
     * events are detected the steps are shortened so that events occur on step boundaries
     * which means the previous state may be computed by the interpolator.
     *
     * @return {@code true} if the previous state was calculated by the interpolator and
     * false if it was computed directly by the integrator.
     */
    bool is_previous_state_interpolated();

    /**
     * Get the state at current grid point time.
     * @return state at current grid point time
     */
    Field_ODE_State_And_Derivative<T> get_current_state();

    /**
     * Determines if the {@link #get_current_state() current state} is computed directly by
     * the integrator, or if it is calculated using {@link #get_interpolated_state(Calculus_Field_Element)
     * interpolation}.
     *
     * <p> Typically the current state is directly computed by the integrator, but when
     * events are detected the steps are shortened so that events occur on step boundaries
     * which means the current state may be computed by the interpolator.
     *
     * @return {@code true} if the current state was calculated by the interpolator and
     * false if it was computed directly by the integrator.
     */
    bool is_current_state_interpolated();

    /**
     * Get the state at interpolated time.
     * <p>Setting the time outside of the current step is allowed, but
     * should be used with care since the accuracy of the interpolator will
     * probably be very poor far from this step. This allowance has been
     * added to simplify implementation of search algorithms near the
     * step endpoints.</p>
     * @param time time of the interpolated point
     * @return state at interpolated time
     */
    Field_ODE_State_And_Derivative<T> get_interpolated_state(T time);

    /** Check if the natural integration direction is forward.
     * <p>This method provides the integration direction as specified by
     * the integrator itself, it avoid some nasty problems in
     * degenerated cases like NULL steps due to cancellation at step
     * initialization, step control or discrete events
     * triggering.</p>
     * @return true if the integration variable (time) increases during
     * integration
     */
    bool is_forward();

};