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

//import java.util.Collection;

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.analysis.solvers.Bracketed_Univariate_Solver;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.ode.events.Event_Handler_Configuration;
//import org.hipparchus.ode.events.ODE_Event_Handler;
//import org.hipparchus.ode.sampling.ODE_Step_Handler;

/** This interface represents a first order integrator for
 * differential equations.

 * <p>The classes which are devoted to solve first order differential
 * equations should implement this interface. The problems which can
 * be handled should implement the {@link
 * Ordinary_Differential_Equation} interface.</p>
 *
 * @see Ordinary_Differential_Equation
 * @see org.hipparchus.ode.sampling.ODE_Step_Handler
 * @see org.hipparchus.ode.events.ODE_Event_Handler
 */
class ODE_Integrator  
{

    /** Get the name of the method.
     * @return name of the method
     */
    std::string get_name();

    /** Add a step handler to this integrator.
     * <p>The handler will be called by the integrator for each accepted
     * step.</p>
     * @param handler handler for the accepted steps
     * @see #get_step_handlers()
     * @see #clear_step_handlers()
     */
    void add_step_handler(ODE_Step_Handler handler);

    /** Get all the step handlers that have been added to the integrator.
     * @return an unmodifiable collection of the added events handlers
     * @see #add_step_handler(ODE_Step_Handler)
     * @see #clear_step_handlers()
     */
    Collection<ODE_Step_Handler> get_step_handlers();

    /** Remove all the step handlers that have been added to the integrator.
     * @see #add_step_handler(ODE_Step_Handler)
     * @see #get_step_handlers()
     */
    void clear_step_handlers();

    /**
     * Add an event handler to the integrator.
     *
     * <p> Uses a default {@link org.hipparchus.analysis.solvers.Univariate_Solver} with an absolute accuracy equal to the
     * given convergence threshold, as root-finding algorithm to detect the state events.
     *
     * @param handler           event handler
     * @param max_check_interval  maximal time interval between switching function checks
     *                          (this interval prevents missing sign changes in case the
     *                          integration steps becomes very large)
     * @param convergence       convergence threshold in the event time search. Must be
     *                          smaller than {@code max_check_interval} and should be small
     *                          compared to time scale of the ODE dynamics.
     * @param max_iteration_count upper limit of the iteration count in the event time
     *                          search
     * @see #get_event_handlers()
     * @see #get_event_handlers_configurations()
     * @see #clear_event_handlers()
     */
    void add_event_handler(ODE_Event_Handler handler, double max_check_interval, double convergence, int max_iteration_count);

    /**
     * Add an event handler to the integrator.
     *
     * @param handler           event handler
     * @param max_check_interval  maximal time interval between switching function checks
     *                          (this interval prevents missing sign changes in case the
     *                          integration steps becomes very large)
     * @param convergence       convergence threshold in the event time search. Must be
     *                          smaller than {@code max_check_interval} and should be small
     *                          compared to time scale of the ODE dynamics.
     * @param max_iteration_count upper limit of the iteration count in the event time
     *                          search
     * @param solver            The root-finding algorithm to use to detect the state
     *                          events.
     * @see #get_event_handlers()
     * @see #get_event_handlers_configurations()
     * @see #clear_event_handlers()
     */
    void add_event_handler(ODE_Event_Handler handler, double max_check_interval, double convergence, int max_iteration_count, Bracketed_Univariate_Solver<Univariate_Function> solver);

    /** Get all the event handlers that have been added to the integrator.
     * @return an unmodifiable collection of the added events handlers
     * @see #add_event_handler(ODE_Event_Handler, double, double, int)
     * @see #add_event_handler(ODE_Event_Handler, double, double, int, Bracketed_Univariate_Solver)
     * @see #get_event_handlers_configurations()
     * @see #clear_event_handlers()
     */
    Collection<ODE_Event_Handler> get_event_handlers();

    /** Get all the event handlers configurations that have been added to the integrator.
     * @return an unmodifiable collection of the added events handlers configurations
     * @see #add_event_handler(ODE_Event_Handler, double, double, int)
     * @see #add_event_handler(ODE_Event_Handler, double, double, int, Bracketed_Univariate_Solver)
     * @see #get_event_handlers()
     * @see #clear_event_handlers()
     * @since 2.0
     */
    Collection<Event_Handler_Configuration> get_event_handlers_configurations();

    /** Remove all the event handlers that have been added to the integrator.
     * @see #add_event_handler(ODE_Event_Handler, double, double, int)
     * @see #add_event_handler(ODE_Event_Handler, double, double, int, Bracketed_Univariate_Solver)
     * @see #get_event_handlers()
     * @see #get_event_handlers_configurations()
     */
    void clear_event_handlers();

    /** Get the state at step start time t<sub>i</sub>.
     * <p>This method can be called during integration (typically by
     * the object implementing the {@link Ordinary_Differential_Equation
     * differential equations} problem) if the value of the current step that
     * is attempted is needed.</p>
     * <p>The result is undefined if the method is called outside of
     * calls to <code>integrate</code>.</p>
     * @return state at step start time t<sub>i</sub>
     */
    ODE_State_And_Derivative get_step_start();

    /** Get the current signed value of the integration stepsize.
     * <p>This method can be called during integration (typically by
     * the object implementing the {@link Ordinary_Differential_Equation
     * differential equations} problem) if the signed value of the current stepsize
     * that is tried is needed.</p>
     * <p>The result is undefined if the method is called outside of
     * calls to <code>integrate</code>.</p>
     * @return current signed value of the stepsize
     */
    double get_current_signed_stepsize();

    /** Set the maximal number of differential equations function evaluations.
     * <p>The purpose of this method is to avoid infinite loops which can occur
     * for example when stringent error constraints are set or when lots of
     * discrete events are triggered, thus leading to many rejected steps.</p>
     * @param max_evaluations maximal number of function evaluations (negative
     * values are silently converted to maximal integer value, thus representing
     * almost unlimited evaluations)
     */
    void set_max_evaluations(const int& max_evaluations);

    /** Get the maximal number of functions evaluations.
     * @return maximal number of functions evaluations
     */
    int get_max_evaluations();

    /** Get the number of evaluations of the differential equations function.
     * <p>
     * The number of evaluations corresponds to the last call to the
     * <code>integrate</code> method. It is 0 if the method has not been called yet.
     * </p>
     * @return number of evaluations of the differential equations function
     */
    int get_evaluations();

    /** Integrate the differential equations up to the given time.
     * <p>This method solves an Initial Value Problem (IVP).</p>
     * <p>sin_ce this method stores some internal state variables made
     * available in its class during integration ({@link
     * #get_current_signed_stepsize()}), it is <em>not</em> thread-safe.</p>
     * @param equations differential equations to integrate
     * @param initial_state initial state (time, primary and secondary state vectors)
     * @param const_time target time for the integration
     * (can be set to a value smaller than {@code t0} for backward integration)
     * @return const state, its time will be the same as {@code const_time} if
     * integration reached its target, but may be different if some {@link
     * org.hipparchus.ode.events.ODE_Event_Handler} stops it at some point.
     * @exception  if integration step is too small
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if the location of an event cannot be bracketed
     */
    ODE_State_And_Derivative integrate(Expandable_ODE equations, ODE_State initial_state, double const_time)
        , Math_Illegal_State_Exception;

    /** Integrate the differential equations up to the given time.
     * <p>This method solves an Initial Value Problem (IVP).</p>
     * <p>sin_ce this method stores some internal state variables made
     * available in its class during integration ({@link
     * #get_current_signed_stepsize()}), it is <em>not</em> thread-safe.</p>
     * @param equations differential equations to integrate
     * @param initial_state initial state (time, primary and secondary state vectors)
     * @param const_time target time for the integration
     * (can be set to a value smaller than {@code t0} for backward integration)
     * @return const state, its time will be the same as {@code const_time} if
     * integration reached its target, but may be different if some {@link
     * org.hipparchus.ode.events.ODE_Event_Handler} stops it at some point.
     * @exception  if integration step is too small
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if the location of an event cannot be bracketed
     */
    default ODE_State_And_Derivative integrate(Ordinary_Differential_Equation equations, ODE_State initial_state, double const_time)
        , Math_Illegal_State_Exception 
        {
        return integrate(new Expandable_ODE(equations), initial_state, const_time);
    }

}


