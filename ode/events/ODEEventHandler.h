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

//package org.hipparchus.ode.events;

//import org.hipparchus.ode.ODE_State;
//import org.hipparchus.ode.ODE_State_And_Derivative;

/** This interface represents a handler for discrete events triggered
 * during ODE integration.
 *
 * <p>Some events can be triggered at discrete times as an ODE problem
 * is solved. This occurs for example when the integration process
 * should be stopped as some state is reached (G-stop facility) when the
 * precise date is unknown a priori, or when the derivatives have
 * discontinuities, or simply when the user wants to monitor some
 * states boundaries crossings.
 * </p>
 *
 * <p>These events are defined as occurring when a <code>g</code>
 * switching function sign changes.</p>
 *
 * <p>sin_ce events are only problem-dependent and are triggered by the
 * independent <i>time</i> variable and the state vector, they can
 * occur at virtually any time, unknown in advance. The integrators will
 * take care to avoid sign changes inside the steps, they will reduce
 * the step size when such an event is detected in order to put this
 * event exactly at the end of the current step. This guarantees that
 * step interpolation (which always has a one step scope) is relevant
 * even in presence of discontinuities. This is independent from the
 * stepsize control provided by integrators that monitor the local
 * error (this event handling feature is available for all integrators, * including fixed step ones).</p>
 *
 * @see org.hipparchus.ode.events
 */
class ODE_Event_Handler  
{

    /** Initialize event handler at the start of an ODE integration.
     * <p>
     * This method is called once at the start of the integration. It
     * may be used by the event handler to initialize some internal data
     * if needed.
     * </p>
     * <p>
     * The default implementation does nothing
     * </p>
     * @param initial_state initial time, state vector and derivative
     * @param const_time target time for the integration
     */
    default void init(ODE_State_And_Derivative initial_state, double const_time) 
    {
        // nothing by default
    }

    /** Compute the value of the switching function.

     * <p>The discrete events are generated when the sign of this
     * switching function changes. The integrator will take care to change
     * the stepsize in such a way these events occur exactly at step boundaries.
     * The switching function must be continuous in its roots neighborhood
     * (but not necessarily smooth), as the integrator will need to find its
     * roots to locate precisely the events.</p>
     *
     * <p>Also note that for the integrator to detect an event the sign of the switching
     * function must have opposite signs just before and after the event. If this
     * consistency is not preserved the integrator may not detect any events.
     *
     * <p>This need for consistency is sometimes tricky to achieve. A typical
     * example is using an event to model a ball bouncing on the floor. The first
     * idea to represent this would be to have {@code g(state) = h(state)} where h is the
     * height above the floor at time {@code state.get_time()}. When {@code g(state)}
     * reaches 0, the ball is on the floor, so it should bounce and the typical way to do this is
     * to reverse its vertical velocity. However, this would mean that before the
     * event {@code g(state)} was decreasing from positive values to 0, and after the
     * event {@code g(state)} would be increasing from 0 to positive values again.
     * Consistency is broken here! The solution here is to have {@code g(state) = sign
     * * h(state)}, where sign is a variable with initial value set to {@code +1}. Each
     * time {@link #event_occurred(ODE_State_And_Derivative, bool) event_occurred} is called, * {@code sign} is reset to {@code -sign}. This allows the {@code g(state)}
     * function to remain continuous (and even smooth) even across events, despite
     * {@code h(state)} is not. Basically, the event is used to <em>fold</em> {@code h(state)}
     * at bounce points, and {@code sign} is used to <em>unfold</em> it back, so the
     * solvers sees a {@code g(state)} function which behaves smoothly even across events.</p>
     *
     * <p>This method is idempotent, that is calling this multiple times with the same
     * state will result in the same value, with two exceptions. First, the definition of
     * the g function may change when an {@link #event_occurred(ODE_State_And_Derivative, * bool) event occurs} on this handler, as in the above example. Second, the
     * definition of the g function may change when the {@link
     * #event_occurred(ODE_State_And_Derivative, bool) event_occurred} method of any other
     * event handler in the same integrator returns {@link Action#RESET_EVENTS}, {@link
     * Action#RESET_DERIVATIVES}, or {@link Action#RESET_STATE}.
     *
     * @param state current value of the independent <i>time</i> variable, state vector
     * and derivative
     * @return value of the g switching function
     * @see org.hipparchus.ode.events
     */
    double g(ODE_State_And_Derivative state);

    /** Handle an event and choose what to do next.

     * <p>This method is called when the integrator has accepted a step
     * ending exactly on a sign change of the function, just <em>after</em>
     * the step handler itself is called (see below for scheduling). It
     * allows the user to update his internal data to acknowledge the fact
     * the event has been handled (for example setting a flag in the {@link
     * org.hipparchus.ode.Ordinary_Differential_Equation
     * differential equations} to switch the derivatives computation in
     * case of discontinuity), or to direct the integrator to either stop
     * or continue integration, possibly with a reset state or derivatives.</p>
     *
     * <ul>
     *   <li>if {@link Action#STOP} is returned, the integration will be stopped,</li>
     *   <li>if {@link Action#RESET_STATE} is returned, the {@link #reset_state
     *   reset_state} method will be called once the step handler has
     *   finished its task, and the integrator will also recompute the
     *   derivatives,</li>
     *   <li>if {@link Action#RESET_DERIVATIVES} is returned, the integrator
     *   will recompute the derivatives, *   <li>if {@link Action#RESET_EVENTS} is returned, the integrator
     *   will recheck all event handlers, *   <li>if {@link Action#CONTINUE} is returned, no specific action will
     *   be taken (apart from having called this method) and integration
     *   will continue.</li>
     * </ul>
     *
     * <p>The scheduling between this method and the {@link
     * org.hipparchus.ode.sampling.ODE_Step_Handler ODE_Step_Handler} method {@link
     * org.hipparchus.ode.sampling.ODE_Step_Handler#handle_step(org.hipparchus.ode.sampling.ODE_StateInterpolator)
     * handle_step(interpolator)} is to call {@code handle_step} first and this method afterwards
     * (this scheduling changed as of Hipparchus 2.0). This scheduling allows user code
     * called by this method and user code called by step handlers to get values
     * of the independent time variable consistent with integration direction.</p>
     *
     * @param state current value of the independent <i>time</i> variable, state vector
     * and derivative
     * @param increasing if true, the value of the switching function increases
     * when times increases around event (note that increase is measured with respect
     * to physical time, not with respect to integration which may go backward in time)
     * @return indication of what the integrator should do next, this
     * value must be one of {@link Action#STOP}, {@link Action#RESET_STATE}, * {@link Action#RESET_DERIVATIVES}, {@link Action#RESET_EVENTS}, or
     * {@link Action#CONTINUE}
     */
    Action event_occurred(ODE_State_And_Derivative state, bool increasing);

    /** Reset the state prior to continue the integration.
     *
     * <p>This method is called after the step handler has returned and
     * before the next step is started, but only when {@link
     * #event_occurred} has itself returned the {@link Action#RESET_STATE}
     * indicator. It allows the user to reset the state vector for the
     * next step, without perturbing the step handler of the finishing
     * step.</p>
     * <p>The default implementation returns its argument.</p>
     * @param state current value of the independent <i>time</i> variable, state vector
     * and derivative
     * @return reset state (note that it does not include the derivatives, they will
     * be added automatically by the integrator afterwards)
     */
    default ODE_State reset_state(ODE_State_And_Derivative state) 
    {
        return state;
    }

}


