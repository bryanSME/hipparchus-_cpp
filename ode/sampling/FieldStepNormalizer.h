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

//package org.hipparchus.ode.sampling;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.ode.Field_ODE_State_And_Derivative;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Precision;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"

/**
 * This class wraps an object implementing {@link FieldODE_Fixed_Step_Handler}
 * into a {@link FieldODE_Step_Handler}.

 * <p>This wrapper allows to use fixed step handlers with general
 * integrators which cannot guaranty their integration steps will
 * remain constant and therefore only accept general step
 * handlers.</p>
 *
 * <p>The stepsize used is selected at construction time. The {@link
 * FieldODE_Fixed_Step_Handler#handle_step handle_step} method of the underlying
 * {@link FieldODE_Fixed_Step_Handler} object is called at normalized times. The
 * normalized times can be influenced by the {@link _step_normalizer_mode} and
 * {@link Step_Normalizer_Bounds}.</p>
 *
 * <p>There is no constraint on the integrator, it can use any time step
 * it needs (time steps longer or shorter than the fixed time step and
 * non-integer ratios are all allowed).</p>
 *
 * <p>
 * <table border="1" align="center">
 * <tr BGCOLOR="#CCCCFF"><td colspan=6><font size="+2">Examples (step size = 0.5)</font></td></tr>
 * <tr BGCOLOR="#EEEEFF"><font size="+1"><td>Start time</td><td>End time</td>
 *  <td>Direction</td><td>{@link _step_normalizer_mode Mode}</td>
 *  <td>{@link Step_Normalizer_Bounds Bounds}</td><td>Output</td></font></tr>
 * <tr><td>0.3</td><td>3.1</td><td>forward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#NEITHER NEITHER}</td><td>0.8, 1.3, 1.8, 2.3, 2.8</td></tr>
 * <tr><td>0.3</td><td>3.1</td><td>forward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#FIRST FIRST}</td><td>0.3, 0.8, 1.3, 1.8, 2.3, 2.8</td></tr>
 * <tr><td>0.3</td><td>3.1</td><td>forward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#LAST LAST}</td><td>0.8, 1.3, 1.8, 2.3, 2.8, 3.1</td></tr>
 * <tr><td>0.3</td><td>3.1</td><td>forward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#BOTH BOTH}</td><td>0.3, 0.8, 1.3, 1.8, 2.3, 2.8, 3.1</td></tr>
 * <tr><td>0.3</td><td>3.1</td><td>forward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#NEITHER NEITHER}</td><td>0.5, 1.0, 1.5, 2.0, 2.5, 3.0</td></tr>
 * <tr><td>0.3</td><td>3.1</td><td>forward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#FIRST FIRST}</td><td>0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0</td></tr>
 * <tr><td>0.3</td><td>3.1</td><td>forward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#LAST LAST}</td><td>0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.1</td></tr>
 * <tr><td>0.3</td><td>3.1</td><td>forward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#BOTH BOTH}</td><td>0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.1</td></tr>
 * <tr><td>0.0</td><td>3.0</td><td>forward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#NEITHER NEITHER}</td><td>0.5, 1.0, 1.5, 2.0, 2.5, 3.0</td></tr>
 * <tr><td>0.0</td><td>3.0</td><td>forward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#FIRST FIRST}</td><td>0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0</td></tr>
 * <tr><td>0.0</td><td>3.0</td><td>forward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#LAST LAST}</td><td>0.5, 1.0, 1.5, 2.0, 2.5, 3.0</td></tr>
 * <tr><td>0.0</td><td>3.0</td><td>forward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#BOTH BOTH}</td><td>0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0</td></tr>
 * <tr><td>0.0</td><td>3.0</td><td>forward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#NEITHER NEITHER}</td><td>0.5, 1.0, 1.5, 2.0, 2.5, 3.0</td></tr>
 * <tr><td>0.0</td><td>3.0</td><td>forward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#FIRST FIRST}</td><td>0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0</td></tr>
 * <tr><td>0.0</td><td>3.0</td><td>forward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#LAST LAST}</td><td>0.5, 1.0, 1.5, 2.0, 2.5, 3.0</td></tr>
 * <tr><td>0.0</td><td>3.0</td><td>forward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#BOTH BOTH}</td><td>0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0</td></tr>
 * <tr><td>3.1</td><td>0.3</td><td>backward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#NEITHER NEITHER}</td><td>2.6, 2.1, 1.6, 1.1, 0.6</td></tr>
 * <tr><td>3.1</td><td>0.3</td><td>backward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#FIRST FIRST}</td><td>3.1, 2.6, 2.1, 1.6, 1.1, 0.6</td></tr>
 * <tr><td>3.1</td><td>0.3</td><td>backward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#LAST LAST}</td><td>2.6, 2.1, 1.6, 1.1, 0.6, 0.3</td></tr>
 * <tr><td>3.1</td><td>0.3</td><td>backward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#BOTH BOTH}</td><td>3.1, 2.6, 2.1, 1.6, 1.1, 0.6, 0.3</td></tr>
 * <tr><td>3.1</td><td>0.3</td><td>backward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#NEITHER NEITHER}</td><td>3.0, 2.5, 2.0, 1.5, 1.0, 0.5</td></tr>
 * <tr><td>3.1</td><td>0.3</td><td>backward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#FIRST FIRST}</td><td>3.1, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5</td></tr>
 * <tr><td>3.1</td><td>0.3</td><td>backward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#LAST LAST}</td><td>3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.3</td></tr>
 * <tr><td>3.1</td><td>0.3</td><td>backward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#BOTH BOTH}</td><td>3.1, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.3</td></tr>
 * <tr><td>3.0</td><td>0.0</td><td>backward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#NEITHER NEITHER}</td><td>2.5, 2.0, 1.5, 1.0, 0.5, 0.0</td></tr>
 * <tr><td>3.0</td><td>0.0</td><td>backward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#FIRST FIRST}</td><td>3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0</td></tr>
 * <tr><td>3.0</td><td>0.0</td><td>backward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#LAST LAST}</td><td>2.5, 2.0, 1.5, 1.0, 0.5, 0.0</td></tr>
 * <tr><td>3.0</td><td>0.0</td><td>backward</td><td>{@link _step_normalizer_mode#INCREMENT INCREMENT}</td><td>{@link Step_Normalizer_Bounds#BOTH BOTH}</td><td>3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0</td></tr>
 * <tr><td>3.0</td><td>0.0</td><td>backward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#NEITHER NEITHER}</td><td>2.5, 2.0, 1.5, 1.0, 0.5, 0.0</td></tr>
 * <tr><td>3.0</td><td>0.0</td><td>backward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#FIRST FIRST}</td><td>3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0</td></tr>
 * <tr><td>3.0</td><td>0.0</td><td>backward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#LAST LAST}</td><td>2.5, 2.0, 1.5, 1.0, 0.5, 0.0</td></tr>
 * <tr><td>3.0</td><td>0.0</td><td>backward</td><td>{@link _step_normalizer_mode#MULTIPLES MULTIPLES}</td><td>{@link Step_Normalizer_Bounds#BOTH BOTH}</td><td>3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0</td></tr>
 * </table>
 * </p>
 *
 * @param <T> the type of the field elements
 * @see FieldODE_Step_Handler
 * @see FieldODE_Fixed_Step_Handler
 * @see _step_normalizer_mode
 * @see Step_Normalizer_Bounds
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class FieldStep_Normalizer : public FieldODE_Step_Handler<T> 
{
    /** Fixed time step. */
    private double h;

    /** Underlying step handler. */
    private const FieldODE_Fixed_Step_Handler<T> handler;

    /** First step state. */
    private Field_ODE_State_And_Derivative<T> first;

    /** Last step step. */
    private Field_ODE_State_And_Derivative<T> last;

    /** Integration direction indicator. */
    private bool forward;

    /** The step normalizer bounds settings to use. */
    private const Step_Normalizer_Bounds bounds;

    /** The step normalizer mode to use. */
    private const _step_normalizer_mode mode;

    /** Simple constructor. Uses {@link _step_normalizer_mode#INCREMENT INCREMENT}
     * mode, and {@link Step_Normalizer_Bounds#FIRST FIRST} bounds setting, for
     * backwards compatibility.
     * @param h fixed time step (sign is not used)
     * @param handler fixed time step handler to wrap
     */
    public FieldStep_Normalizer(const double h, const FieldODE_Fixed_Step_Handler<T> handler) 
    {
        this(h, handler, _step_normalizer_mode.INCREMENT, Step_Normalizer_Bounds.FIRST);
    }

    /** Simple constructor. Uses {@link Step_Normalizer_Bounds#FIRST FIRST}
     * bounds setting.
     * @param h fixed time step (sign is not used)
     * @param handler fixed time step handler to wrap
     * @param mode step normalizer mode to use
     */
    public FieldStep_Normalizer(const double h, const FieldODE_Fixed_Step_Handler<T> handler, const _step_normalizer_mode mode) 
    {
        this(h, handler, mode, Step_Normalizer_Bounds.FIRST);
    }

    /** Simple constructor. Uses {@link _step_normalizer_mode#INCREMENT INCREMENT}
     * mode.
     * @param h fixed time step (sign is not used)
     * @param handler fixed time step handler to wrap
     * @param bounds step normalizer bounds setting to use
     */
    public FieldStep_Normalizer(const double h, const FieldODE_Fixed_Step_Handler<T> handler, const Step_Normalizer_Bounds bounds) 
    {
        this(h, handler, _step_normalizer_mode.INCREMENT, bounds);
    }

    /** Simple constructor.
     * @param h fixed time step (sign is not used)
     * @param handler fixed time step handler to wrap
     * @param mode step normalizer mode to use
     * @param bounds step normalizer bounds setting to use
     */
    public FieldStep_Normalizer(const double h, const FieldODE_Fixed_Step_Handler<T> handler, const _step_normalizer_mode mode, const Step_Normalizer_Bounds bounds) 
    {
        this.h       = std::abs(h);
        this.handler = handler;
        this.mode    = mode;
        this.bounds  = bounds;
        first        = NULL;
        last         = NULL;
        forward      = true;
    }

    /** {@inherit_doc} */
    //override
    public void init(const Field_ODE_State_And_Derivative<T> initial_state, const T const_time) 
    {

        first   = NULL;
        last    = NULL;
        forward = true;

        // initialize the underlying handler
        handler.init(initial_state, const_time);

    }

    /** {@inherit_doc} */
    //override
    public void handle_step(const FieldODE_StateInterpolator<T> interpolator) 
    {
        // The first time, update the last state with the start information.
        if (last == NULL) 
        {

            first   = interpolator.get_previous_state();
            last    = first;

            // Take the integration direction into account.
            forward = interpolator.is_forward();
            if (!forward) 
            {
                h = -h;
            }
        }

        // Calculate next normalized step time.
        T next_time = (mode == _step_normalizer_mode.INCREMENT) ?
                     last.get_time().add(h) :
                     last.get_time().get_field().get_zero().add((std::floor(last.get_time().get_real() / h) + 1) * h);
        if (mode == _step_normalizer_mode.MULTIPLES &&
            Precision.equals(next_time.get_real(), last.get_time().get_real(), 1)) 
            {
            next_time = next_time.add(h);
        }

        // Process normalized steps as long as they are in the current step.
        bool next_in_step = is_next_in_step(next_time, interpolator);
        while (next_in_step) 
        {
            // Output the stored previous step.
            do_normalized_step(false);

            // Store the next step as last step.
            last = interpolator.get_interpolated_state(next_time);

            // Move on to the next step.
            next_time = next_time.add(h);
            next_in_step = is_next_in_step(next_time, interpolator);
        }
    }

    /** {@inherit_doc} */
    //override
    public void finish(const Field_ODE_State_And_Derivative<T> const_state) 
    {
        // There will be no more steps. The stored one should be given to
        // the handler. We may have to output one more step. Only the last
        // one of those should be flagged as being the last.
        const bool add_last = bounds.last_included() && last.get_time().get_real() != const_state.get_time().get_real();
        do_normalized_step(!add_last);
        if (add_last) 
        {
            last = const_state;
            do_normalized_step(true);
        }
    }

    /**
     * Returns a value indicating whether the next normalized time is in the
     * current step.
     * @param next_time the next normalized time
     * @param interpolator interpolator for the last accepted step, to use to
     * get the end time of the current step
     * @return value indicating whether the next normalized time is in the
     * current step
     */
    private bool is_next_in_step(const T next_time, const FieldODE_StateInterpolator<T> interpolator) 
    {
        return forward ?
               next_time.get_real() <= interpolator.get_current_state().get_time().get_real() :
               next_time.get_real() >= interpolator.get_current_state().get_time().get_real();
    }

    /**
     * Invokes the underlying step handler for the current normalized step.
     * @param is_last true if the step is the last one
     */
    private void do_normalized_step(const bool is_last) 
    {
        if (!bounds.first_included() && first.get_time().get_real() == last.get_time().get_real()) 
        {
            return;
        }
        handler.handle_step(last, is_last);
    }

};