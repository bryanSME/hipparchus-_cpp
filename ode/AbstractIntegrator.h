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

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.analysis.solvers.Bracketed_Univariate_Solver;
//import org.hipparchus.analysis.solvers.Bracketing_Nth_Order_Brent_Solver;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.ode.events.Action;
//import org.hipparchus.ode.events.Event_Handler_Configuration;
//import org.hipparchus.ode.events.Event_State;
//import org.hipparchus.ode.events.Event_State.Event_Occurrence;
//import org.hipparchus.ode.events.ODE_Event_Handler;
//import org.hipparchus.ode.sampling.AbstractODE_StateInterpolator;
//import org.hipparchus.ode.sampling.ODE_Step_Handler;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Incrementor;

//import java.util.Array_list;
//import java.util.Collection;
//import java.util.Collections;
//import java.util.Comparator;
//import java.util.List;
//import java.util.Priority_Queue;
//import java.util.Queue;

/**
 * Base class managing common boilerplate for all integrators.
 */
class Abstract_Integrator : ODE_Integrator 
{

    /** Default relative accuracy. */
    private static const double DEFAULT_RELATIVE_ACCURACY = 0;

    /** Default function value accuracy. */
    private static const double DEFAULT_FUNCTION_VALUE_ACCURACY = 0;

    /** Step handler. */
    private Collection<ODE_Step_Handler> step_handlers;

    /** Current step start time. */
    private ODE_State_And_Derivative step_start;

    /** Current stepsize. */
    private double step_size;

    /** Indicator for last step. */
    private bool is_last_step;

    /** Indicator that a state or derivative reset was triggered by some event. */
    private bool reset_occurred;

    /** Events states. */
    private Collection<Event_State> events_states;

    /** Initialization indicator of events states. */
    private bool states_initialized;

    /** Name of the method. */
    private const std::string name;

    /** Counter for number of evaluations. */
    private Incrementor evaluations;

    /** Differential equations to integrate. */
    private transient Expandable_ODE equations;

    /** Build an instance.
     * @param name name of the method
     */
    protected Abstract_Integrator(const std::string name) 
    {
        this.name         = name;
        step_handlers      = Array_list<>();
        step_start         = NULL;
        step_size          = std::numeric_limits<double>::quiet_NaN();
        events_states      = Array_list<>();
        states_initialized = false;
        evaluations       = Incrementor();
    }

    /** {@inherit_doc} */
    //override
    public std::string get_name() 
    {
        return name;
    }

    /** {@inherit_doc} */
    //override
    public void add_step_handler(const ODE_Step_Handler handler) 
    {
        step_handlers.add(handler);
    }

    /** {@inherit_doc} */
    //override
    public Collection<ODE_Step_Handler> get_step_handlers() 
    {
        return Collections.unmodifiable_collection(step_handlers);
    }

    /** {@inherit_doc} */
    //override
    public void clear_step_handlers() 
    {
        step_handlers.clear();
    }

    /** {@inherit_doc} */
    //override
    public void add_event_handler(const ODE_Event_Handler handler, const double max_check_interval, const double convergence, const int max_iteration_count) 
    {
        add_event_handler(handler, max_check_interval, convergence, max_iteration_count, Bracketing_Nth_Order_Brent_Solver(DEFAULT_RELATIVE_ACCURACY, convergence, DEFAULT_FUNCTION_VALUE_ACCURACY, 5));
    }

    /** {@inherit_doc} */
    //override
    public void add_event_handler(const ODE_Event_Handler handler, const double max_check_interval, const double convergence, const int max_iteration_count, const Bracketed_Univariate_Solver<Univariate_Function> solver) 
    {
        events_states.add(new Event_State(handler, max_check_interval, convergence, max_iteration_count, solver));
    }

    /** {@inherit_doc} */
    //override
    public Collection<ODE_Event_Handler> get_event_handlers() 
    {
        const List<ODE_Event_Handler> list = Array_list<>(events_states.size());
        for (Event_State state : events_states) 
        {
            list.add(state.get_event_handler());
        }
        return Collections.unmodifiable_collection(list);
    }

    /** {@inherit_doc} */
    //override
    public Collection<Event_Handler_Configuration> get_event_handlers_configurations() 
    {
        return Collections.unmodifiable_collection(events_states);
    }

    /** {@inherit_doc} */
    //override
    public void clear_event_handlers() 
    {
        events_states.clear();
    }

    /** {@inherit_doc} */
    //override
    public double get_current_signed_stepsize() 
    {
        return step_size;
    }

    /** {@inherit_doc} */
    //override
    public void set_max_evaluations(const int& max_evaluations) 
    {
        evaluations = evaluations.with_maximal_count((max_evaluations < 0) ? Integer.MAX_VALUE : max_evaluations);
    }

    /** {@inherit_doc} */
    //override
    public int get_max_evaluations() 
    {
        return evaluations.get_maximal_count();
    }

    /** {@inherit_doc} */
    //override
    public int get_evaluations() 
    {
        return evaluations.get_count();
    }

    /**
     * Prepare the start of an integration.
     *
     * @param eqn equations to integrate
     * @param s0  initial state vector
     * @param t   target time for the integration
     * @return Initial state with computed derivatives.
     */
    protected ODE_State_And_Derivative init_integration(const Expandable_ODE eqn, const ODE_State s0, const double t) 
    {

        this.equations = eqn;
        evaluations    = evaluations.with_count(0);

        // initialize ODE
        eqn.init(s0, t);

        // set up derivatives of initial state (including primary and secondary components)
        const double   t0    = s0.get_time();
        const std::vector<double> y0    = s0.get_complete_state();
        const std::vector<double> y0_dot = compute_derivatives(t0, y0);

        // built the state
        const ODE_State_And_Derivative s0_with_derivatives =
                        eqn.get_mapper().map_state_and_derivative(t0, y0, y0_dot);

        // initialize event handlers
        for (const Event_State state : events_states) 
        {
            state.get_event_handler().init(s0_with_derivatives, t);
        }

        // initialize step handlers
        for (ODE_Step_Handler handler : step_handlers) 
        {
            handler.init(s0_with_derivatives, t);
        }

        set_state_initialized(false);

        return s0_with_derivatives;

    }

    /** Get the differential equations to integrate.
     * @return differential equations to integrate
     */
    protected Expandable_ODE get_equations() 
    {
        return equations;
    }

    /** Get the evaluations counter.
     * @return evaluations counter
     */
    protected Incrementor get_evaluations_counter() 
    {
        return evaluations;
    }

    /** Compute the derivatives and check the number of evaluations.
     * @param t current value of the independent <I>time</I> variable
     * @param y array containing the current value of the state vector
     * @return state completed with derivatives
     * @exception  if arrays dimensions do not match equations settings
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception Null_Pointer_Exception if the ODE equations have not been set (i.e. if this method
     * is called outside of a call to {@link #integrate(Expandable_ODE, ODE_State, double) integrate}
     */
    public std::vector<double> compute_derivatives(const double t, const std::vector<double> y)
        , Math_Illegal_State_Exception, Null_Pointer_Exception 
        {
        evaluations.increment();
        return equations.compute_derivatives(t, y);
    }

    /** Set the state_initialized flag.
     * <p>This method must be called by integrators with the value
     * {@code false} before they start integration, so a proper lazy
     * initialization is done automatically on the first step.</p>
     * @param state_initialized value for the flag
     */
    protected void set_state_initialized(const bool state_initialized) 
    {
        this.states_initialized = state_initialized;
    }

    /** Accept a step, triggering events and step handlers.
     * @param interpolator step interpolator
     * @param t_end const integration time
     * @return state at end of step
     * @exception Math_Illegal_State_Exception if the interpolator one because
     * the number of functions evaluations is exceeded
     * @exception  if the location of an event cannot be bracketed
     * @exception  if arrays dimensions do not match equations settings
     */
    protected ODE_State_And_Derivative accept_step(const AbstractODE_StateInterpolator interpolator, const double t_end)
            , Math_Illegal_State_Exception 
            {

        ODE_State_And_Derivative previous_state = interpolator.get_global_previous_state();
        const ODE_State_And_Derivative current_state = interpolator.get_global_current_state();
        AbstractODE_StateInterpolator restricted = interpolator;


        // initialize the events states if needed
        if (!states_initialized) 
        {
            for (Event_State state : events_states) 
            {
                state.reinitialize_begin(interpolator);
            }
            states_initialized = true;
        }

        // search for next events that may occur during the step
        const int ordering_sign = interpolator.is_forward() ? +1 : -1;
        const Queue<Event_State> occurring_events = Priority_Queue<>(new Comparator<Event_State>() 
        {
            /** {@inherit_doc} */
            //override
            public int compare(const Event_State es0, const Event_State es1) 
            {
                return ordering_sign * Double.compare(es0.get_event_time(), es1.get_event_time());
            }
        });

        reset_occurred = false;
        bool done_with_step = false;
        reset_events:
        do 
        {

            // Evaluate all event detectors for events
            occurring_events.clear();
            for (const Event_State state : events_states) 
            {
                if (state.evaluate_step(restricted)) 
                {
                    // the event occurs during the current step
                    occurring_events.add(state);
                }
            }

            do 
            {

                event_loop:
                while (!occurring_events.is_empty()) 
                {

                    // handle the chronologically first event
                    const Event_State current_event = occurring_events.poll();

                    // get state at event time
                    ODE_State_And_Derivative event_state = restricted.get_interpolated_state(current_event.get_event_time());

                    // restrict the interpolator to the first part of the step, up to the event
                    restricted = restricted.restrict_step(previous_state, event_state);

                    // try to advance all event states to current time
                    for (const Event_State state : events_states) 
                    {
                        if (state != current_event && state.try_advance(event_state, interpolator)) 
                        {
                            // we need to handle another event first
                            // remove event we just updated to prevent heap corruption
                            occurring_events.remove(state);
                            // add it back to update its position in the heap
                            occurring_events.add(state);
                            // re-queue the event we were processing
                            occurring_events.add(current_event);
                            continue event_loop;
                        }
                    }
                    // all event detectors agree we can advance to the current event time

                    // handle the first part of the step, up to the event
                    for (const ODE_Step_Handler handler : step_handlers) 
                    {
                        handler.handle_step(restricted);
                    }

                    // acknowledge event occurrence
                    const Event_Occurrence occurrence = current_event.do_event(event_state);
                    const Action action = occurrence.get_action();
                    is_last_step = action == Action.STOP;

                    if (is_last_step) 
                    {

                        // ensure the event is after the root if it is returned STOP
                        // this lets the user integrate to a STOP event and then restart
                        // integration from the same time.
                        const ODE_State_And_Derivative saved_state = event_state;
                        event_state = interpolator.get_interpolated_state(occurrence.get_stop_time());
                        restricted = interpolator.restrict_step(saved_state, event_state);

                        // handle the almost zero size last part of the const step, at event time
                        for (const ODE_Step_Handler handler : step_handlers) 
                        {
                            handler.handle_step(restricted);
                            handler.finish(restricted.get_current_state());
                        }

                    }

                    if (is_last_step) 
                    {
                        // the event asked to stop integration
                        return event_state;
                    }

                    if (action == Action.RESET_DERIVATIVES || action == Action.RESET_STATE) 
                    {
                        // some event handler has triggered changes that
                        // invalidate the derivatives, we need to recompute them
                        const ODE_State new_state = occurrence.get_new_state();
                        const std::vector<double> y = new_state.get_complete_state();
                        const std::vector<double> y_dot = compute_derivatives(new_state.get_time(), y);
                        reset_occurred = true;
                        return equations.get_mapper().map_state_and_derivative(new_state.get_time(), y, y_dot);
                    }
                    // at this point action == Action.CONTINUE or Action.RESET_EVENTS

                    // prepare handling of the remaining part of the step
                    previous_state = event_state;
                    restricted = restricted.restrict_step(event_state, current_state);

                    if (action == Action.RESET_EVENTS) 
                    {
                        continue reset_events;
                    }

                    // at this point action == Action.CONTINUE
                    // check if the same event occurs again in the remaining part of the step
                    if (current_event.evaluate_step(restricted)) 
                    {
                        // the event occurs during the current step
                        occurring_events.add(current_event);
                    }

                }

                // last part of the step, after the last event. Advance all event
                // detectors to the end of the step. Should never find events unless
                // a previous event modified the g function of another event detector when
                // it returned Action.CONTINUE. Detecting such events here is unreliable
                // and RESET_EVENTS should be used instead. Other option is to replace
                // try_advance(...) with a do_advance(...) that an exception when
                // the g function sign is not as expected.
                for (const Event_State state : events_states) 
                {
                    if (state.try_advance(current_state, interpolator)) 
                    {
                        occurring_events.add(state);
                    }
                }

            } while (!occurring_events.is_empty());

            done_with_step = true;
        } while (!done_with_step);

        is_last_step = is_last_step || std::abs(current_state.get_time() - t_end) <= FastMath.ulp(t_end);

        // handle the remaining part of the step, after all events if any
        for (ODE_Step_Handler handler : step_handlers) 
        {
            handler.handle_step(restricted);
            if (is_last_step) 
            {
                handler.finish(restricted.get_current_state());
            }
        }

        return current_state;

    }

    /** Check the integration span.
     * @param initial_state initial state
     * @param t target time for the integration
     * @exception  if integration span is too small
     * @exception  if adaptive step size integrators
     * tolerance arrays dimensions are not compatible with equations settings
     */
    protected void sanity_checks(const ODE_State initial_state, const double t)
         
        {

        const double threshold = 1000 * FastMath.ulp(std::max(std::abs(initial_state.get_time()), std::abs(t)));
        const double dt = std::abs(initial_state.get_time() - t);
        if (dt <= threshold) 
        {
            throw (Localized_ODE_Formats.TOO_SMALL_INTEGRATION_INTERVAL, dt, threshold, false);
        }

    }

    /** Check if a reset occurred while last step was accepted.
     * @return true if a reset occurred while last step was accepted
     */
    protected bool reset_occurred() 
    {
        return reset_occurred;
    }

    /** Set the current step size.
     * @param step_size step size to set
     */
    protected void set_step_size(const double step_size) 
    {
        this.step_size = step_size;
    }

    /** Get the current step size.
     * @return current step size
     */
    protected double get_step_size() 
    {
        return step_size;
    }
    /** Set current step start.
     * @param step_start step start
     */
    protected void set_step_start(const ODE_State_And_Derivative step_start) 
    {
        this.step_start = step_start;
    }

    /**  {@inherit_doc} */
    //override
    public ODE_State_And_Derivative get_step_start() 
    {
        return step_start;
    }

    /** Set the last state flag.
     * @param is_last_step if true, this step is the last one
     */
    protected void set_is_last_step(const bool is_last_step) 
    {
        this.is_last_step = is_last_step;
    }

    /** Check if this step is the last one.
     * @return true if this step is the last one
     */
    protected bool is_last_step() 
    {
        return is_last_step;
    }

}


