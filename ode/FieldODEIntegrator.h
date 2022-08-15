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

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.analysis.solvers.Bracketed_Real_Field_Univariate_Solver;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.ode.events.Field_Event_Handler_Configuration;
//import org.hipparchus.ode.events.FieldODE_Event_Handler;
//import org.hipparchus.ode.sampling.FieldODE_Step_Handler;
#include <type_traits>
#include "../core/CalculusFieldElement.h"

/** This interface represents a first order integrator for
 * differential equations.

 * <p>The classes which are devoted to solve first order differential
 * equations should implement this interface. The problems which can
 * be handled should implement the {@link
 * FieldOrdinary_Differential_Equation} interface.</p>
 *
 * @see FieldOrdinary_Differential_Equation
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class FieldODE_Integrator 
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
    void add_step_handler(FieldODE_Step_Handler<T> handler);

    /** Get all the step handlers that have been added to the integrator.
     * @return an unmodifiable collection of the added events handlers
     * @see #add_step_handler(FieldODE_Step_Handler)
     * @see #clear_step_handlers()
     */
    Collection<FieldODE_Step_Handler<T>> get_step_handlers();

    /** Remove all the step handlers that have been added to the integrator.
     * @see #add_step_handler(FieldODE_Step_Handler)
     * @see #get_step_handlers()
     */
    void clear_step_handlers();

    /** Add an event handler to the integrator.
     * <p>
     * The default solver is a 5<sup>th</sup> order {@link
     * org.hipparchus.analysis.solvers.FieldBracketing_Nth_Order_Brent_Solver}.
     * </p>
     * @param handler event handler
     * @param max_check_interval maximal time interval between switching
     * function checks (this interval prevents missing sign changes in
     * case the integration steps becomes very large)
     * @param convergence convergence threshold in the event time search
     * @param max_iteration_count upper limit of the iteration count in
     * the event time search events.
     * @see #add_event_handler(FieldODE_Event_Handler, double, double, int, Bracketed_Real_Field_Univariate_Solver)
     * @see #get_event_handlers()
     * @see #get_event_handlers_configurations()
     * @see #clear_event_handlers()
     */
    void add_event_handler(FieldODE_Event_Handler<T>  handler, double max_check_interval, double convergence, int max_iteration_count);

    /** Add an event handler to the integrator.
     * @param handler event handler
     * @param max_check_interval maximal time interval between switching
     * function checks (this interval prevents missing sign changes in
     * case the integration steps becomes very large)
     * @param convergence convergence threshold in the event time search
     * @param max_iteration_count upper limit of the iteration count in
     * the event time search events.
     * @param solver solver to use to locate the event
     * @see #add_event_handler(FieldODE_Event_Handler, double, double, int)
     * @see #get_event_handlers()
     * @see #get_event_handlers_configurations()
     * @see #clear_event_handlers()
     */
    void add_event_handler(FieldODE_Event_Handler<T>  handler, double max_check_interval, double convergence, int max_iteration_count, Bracketed_Real_Field_Univariate_Solver<T> solver);

    /** Get all the event handlers that have been added to the integrator.
     * @return an unmodifiable collection of the added events handlers
     * @see #add_event_handler(FieldODE_Event_Handler, double, double, int)
     * @see #add_event_handler(FieldODE_Event_Handler, double, double, int, Bracketed_Real_Field_Univariate_Solver)
     * @see #get_event_handlers_configurations()
     * @see #clear_event_handlers()
     */
    Collection<FieldODE_Event_Handler<T> > get_event_handlers();

    /** Get all the event handlers configurations that have been added to the integrator.
     * @return an unmodifiable collection of the added events handlers configurations
     * @see #add_event_handler(FieldODE_Event_Handler, double, double, int)
     * @see #add_event_handler(FieldODE_Event_Handler, double, double, int, Bracketed_Real_Field_Univariate_Solver)
     * @see #get_event_handlers()
     * @see #clear_event_handlers()
     * @since 2.0
     */
    Collection<Field_Event_Handler_Configuration<T>> get_event_handlers_configurations();

    /** Remove all the event handlers that have been added to the integrator.
     * @see #add_event_handler(FieldODE_Event_Handler, double, double, int)
     * @see #add_event_handler(FieldODE_Event_Handler, double, double, int, Bracketed_Real_Field_Univariate_Solver)
     * @see #get_event_handlers()
     * @see #get_event_handlers_configurations()
     */
    void clear_event_handlers();

    /** Get the state at step start time t<sub>i</sub>.
     * <p>This method can be called during integration (typically by
     * the object implementing the {@link FieldOrdinary_Differential_Equation
     * differential equations} problem) if the value of the current step that
     * is attempted is needed.</p>
     * <p>The result is undefined if the method is called outside of
     * calls to <code>integrate</code>.</p>
     * @return state at step start time t<sub>i</sub>
     */
    Field_ODE_State_And_Derivative<T> get_step_start();

    /** Get the current signed value of the integration stepsize.
     * <p>This method can be called during integration (typically by
     * the object implementing the {@link FieldOrdinary_Differential_Equation
     * differential equations} problem) if the signed value of the current stepsize
     * that is tried is needed.</p>
     * <p>The result is undefined if the method is called outside of
     * calls to <code>integrate</code>.</p>
     * @return current signed value of the stepsize
     */
    T get_current_signed_stepsize();

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
     * org.hipparchus.ode.events.FieldODE_Event_Handler} stops it at some point.
     * @exception  if integration step is too small
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if the location of an event cannot be bracketed
     */
    Field_ODE_State_And_Derivative<T> integrate(FieldExpandable_ODE<T> equations, FieldODE_State<T> initial_state, T const_time)
        , Math_Illegal_State_Exception;

};