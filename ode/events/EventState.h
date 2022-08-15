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

//package org.hipparchus.ode.events;

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.analysis.solvers.Bracketed_Univariate_Solver;
//import org.hipparchus.analysis.solvers.Bracketed_Univariate_Solver.Interval;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.ode.ODE_State;
//import org.hipparchus.ode.ODE_State_And_Derivative;
//import org.hipparchus.ode.sampling.ODE_StateInterpolator;
//import org.hipparchus.util.FastMath;

/** This class handles the state for one {@link ODE_Event_Handler
 * event handler} during integration steps.
 *
 * <p>Each time the integrator proposes a step, the event handler
 * switching function should be checked. This class handles the state
 * of one handler during one integration step, with references to the
 * state at the end of the preceding step. This information is used to
 * decide if the handler should trigger an event or not during the
 * proposed step.</p>
 *
 */
class Event_State : Event_Handler_Configuration 
{

    /** Event handler. */
    private const ODE_Event_Handler handler;

    /** Maximal time interval between events handler checks. */
    private const double max_check_interval;

    /** Convergence threshold for event localization. */
    private const double convergence;

    /** Upper limit in the iteration count for event localization. */
    private const int max_iteration_count;

    /** Time at the beginning of the step. */
    private double t0;

    /** Value of the events handler at the beginning of the step. */
    private double g0;

    /** Sign of g0. */
    private bool g0_positive;

    /** Indicator of event expected during the step. */
    private bool pending_event;

    /** Occurrence time of the pending event. */
    private double pending_event_time;

    /**
     * Time to stop propagation if the event is a stop event. Used to enable stopping at
     * an event and then restarting after that event.
     */
    private double stop_time;

    /** Time after the current event. */
    private double after_event;

    /** Value of the g function after the current event. */
    private double after_g;

    /** The earliest time considered for events. */
    private double earliest_time_considered;

    /** Integration direction. */
    private bool forward;

    /**
     * Direction of g(t) in the propagation direction for the pending event, or if there
     * is no pending event the direction of the previous event.
     */
    private bool increasing;

    /** Root-finding algorithm to use to detect state events. */
    private const Bracketed_Univariate_Solver<Univariate_Function> solver;

    /** Simple constructor.
     * @param handler event handler
     * @param max_check_interval maximal time interval between switching
     * function checks (this interval prevents missing sign changes in
     * case the integration steps becomes very large)
     * @param convergence convergence threshold in the event time search
     * @param max_iteration_count upper limit of the iteration count in
     * the event time search
     * @param solver Root-finding algorithm to use to detect state events
     */
    public Event_State(const ODE_Event_Handler handler, const double max_check_interval, const double convergence, const int max_iteration_count, const Bracketed_Univariate_Solver<Univariate_Function> solver) 
    {
        this.handler           = handler;
        this.max_check_interval  = max_check_interval;
        this.convergence       = std::abs(convergence);
        this.max_iteration_count = max_iteration_count;
        this.solver            = solver;

        // some dummy values ...
        t0                = std::numeric_limits<double>::quiet_NaN();
        g0                = std::numeric_limits<double>::quiet_NaN();
        g0_positive        = true;
        pending_event      = false;
        pending_event_time  = std::numeric_limits<double>::quiet_NaN();
        increasing        = true;
        earliest_time_considered = std::numeric_limits<double>::quiet_NaN();
        after_event = std::numeric_limits<double>::quiet_NaN();
        after_g = std::numeric_limits<double>::quiet_NaN();
    }

    /** {@inherit_doc} */
    //override
    public ODE_Event_Handler get_event_handler() 
    {
        return handler;
    }

    /** {@inherit_doc} */
    //override
    public double get_max_check_interval() 
    {
        return max_check_interval;
    }

    /** {@inherit_doc} */
    //override
    public double get_convergence() 
    {
        return convergence;
    }

    /** {@inherit_doc} */
    //override
    public int get_max_iteration_count() 
    {
        return max_iteration_count;
    }

    /** {@inherit_doc} */
    //override
    public Bracketed_Univariate_Solver<Univariate_Function> get_solver() 
    {
        return solver;
    }

    /** Reinitialize the beginning of the step.
     * @param interpolator valid for the current step
     * @exception Math_Illegal_State_Exception if the interpolator one because
     * the number of functions evaluations is exceeded
     */
    public void reinitialize_begin(const ODE_StateInterpolator interpolator)
        Math_Illegal_State_Exception 
        {

        forward = interpolator.is_forward();
        const ODE_State_And_Derivative s0 = interpolator.get_previous_state();
        t0 = s0.get_time();
        g0 = handler.g(s0);
        while (g0 == 0) 
        {
            // excerpt from MATH-421 issue:
            // If an ODE solver is setup with an ODE_Event_Handler that return STOP
            // when the even is triggered, the integrator stops (which is exactly
            // the expected behavior). If however the user wants to restart the
            // solver from the const state reached at the event with the same
            // configuration (expecting the event to be triggered again at a
            // later time), then the integrator may fail to start. It can get stuck
            // at the previous event. The use case for the bug MATH-421 is fairly
            // general, so events occurring exactly at start in the first step should
            // be ignored. Some g functions may be zero for multiple adjacent values of t
            // so keep skipping roots while g(t) is zero.

            // extremely rare case: there is a zero EXACTLY at interval start
            // we will use the sign slightly after step beginning to force ignoring this zero
            const double epsilon = std::max(solver.get_absolute_accuracy(), std::abs(solver.get_relative_accuracy() * t0));
            double t_start = t0 + (forward ? 0.5 : -0.5) * epsilon;
            // check for case where tolerance is too small to make a difference
            if (t_start == t0) 
            {
                t_start = next_after(t0);
            }
            t0 = t_start;
            g0 = handler.g(interpolator.get_interpolated_state(t_start));
        }
        g0_positive = g0 > 0;
        // "last" event was increasing
        increasing = g0_positive;

    }

    /**
     * Evaluate the impact of the proposed step on the event handler.
     *
     * @param interpolator step interpolator for the proposed step
     * @return true if the event handler triggers an event before the end of the proposed
     * step
     * @Math_Illegal_State_Exception    if the interpolator one because the
     *                                      number of functions evaluations is exceeded
     * @ if the event cannot be bracketed
     */
    public bool evaluate_step(const ODE_StateInterpolator interpolator)
            , Math_Illegal_State_Exception 
            {

        forward = interpolator.is_forward();
        const ODE_State_And_Derivative s1 = interpolator.get_current_state();
        const double t1 = s1.get_time();
        const double dt = t1 - t0;
        if (std::abs(dt) < convergence) 
        {
            // we cannot do anything on such a small step, don't trigger any events
            return false;
        }
        // number of points to check in the current step
        const int n = std::max(1, static_cast<int>( FastMath.ceil(std::abs(dt) / max_check_interval));
        const double h = dt / n;


        double ta = t0;
        double ga = g0;
        for (int i{}; i < n; ++i) 
        {

            // evaluate handler value at the end of the substep
            const double tb = (i == n - 1) ? t1 : t0 + (i + 1) * h;
            const double gb = handler.g(interpolator.get_interpolated_state(tb));

            // check events occurrence
            if (gb == 0.0 || (g0_positive ^ (gb > 0))) 
            {
                // there is a sign change: an event is expected during this step
                if (find_root(interpolator, ta, ga, tb, gb)) 
                {
                    return true;
                }
            }
else 
            {
                // no sign change: there is no event for now
                ta = tb;
                ga = gb;
            }

        }

        // no event during the whole step
        pending_event     = false;
        pending_event_time = std::numeric_limits<double>::quiet_NaN();
        return false;

    }

    /**
     * Find a root in a bracketing interval.
     *
     * <p> When calling this method one of the following must be true. Either ga == 0, gb
     * == 0, (ga < 0  and gb > 0), or (ga > 0 and gb < 0).
     *
     * @param interpolator that covers the interval.
     * @param ta           earliest possible time for root.
     * @param ga           g(ta).
     * @param tb           latest possible time for root.
     * @param gb           g(tb).
     * @return if a zero crossing was found.
     */
    private bool find_root(const ODE_StateInterpolator interpolator, const double ta, const double ga, const double tb, const double gb) 
    {
        // check there appears to be a root in [ta, tb]
        check(ga == 0.0 || gb == 0.0 || (ga > 0.0 && gb < 0.0) || (ga < 0.0 && gb > 0.0));

        const Univariate_Function f = t -> handler.g(interpolator.get_interpolated_state(t));

        // prepare loop below
        double loop_t = ta;
        double loop_g = ga;

        // event time, just at or before the actual root.
        double before_root_t = std::numeric_limits<double>::quiet_NaN();
        double before_root_g = std::numeric_limits<double>::quiet_NaN();
        // time on the other side of the root.
        // Initialized the the loop below executes once.
        double after_root_t = ta;
        double after_root_g = 0.0;

        // check for some conditions that the root finders don't like
        // these conditions cannot not happen in the loop below
        // the ga == 0.0 case is handled by the loop below
        if (ta == tb) 
        {
            // both non-zero but times are the same. Probably due to reset state
            before_root_t = ta;
            before_root_g = ga;
            after_root_t = shifted_by(before_root_t, convergence);
            after_root_g = f.value(after_root_t);
        }
else if (ga != 0.0 && gb == 0.0) 
        {
            // hard: ga != 0.0 and gb == 0.0
            // look past gb by up to convergence to find next sign
            // throw an exception if g(t) = 0.0 in [tb, tb + convergence]
            before_root_t = tb;
            before_root_g = gb;
            after_root_t = shifted_by(before_root_t, convergence);
            after_root_g = f.value(after_root_t);
        }
else if (ga != 0.0) 
        {
            const double new_ga = f.value(ta);
            if (ga > 0 != new_ga > 0) 
            {
                // both non-zero, step sign change at ta, possibly due to reset state
                const double next_t = min_time(shifted_by(ta, convergence), tb);
                const double next_g = f.value(next_t);
                if (next_g > 0.0 == g0_positive) 
                {
                    // the sign change between ga and Ga just moved the root less than one convergence
                    // threshold later, we are still in a regular search for another root before tb, // we just need to fix the bracketing interval
                    // (see issue https://github.com/Hipparchus-Math/hipparchus/issues/184)
                    loop_t = next_t;
                    loop_g = next_g;
                }
else 
                {
                    before_root_t = ta;
                    before_root_g = new_ga;
                    after_root_t = next_t;
                    after_root_g = next_g;
                }

            }
        }

        // loop to skip through "fake" roots, i.e. where g(t) = g'(t) = 0.0
        // executed once if we didn't hit a special case above
        while ((after_root_g == 0.0 || after_root_g > 0.0 == g0_positive) &&
               strictly_after(after_root_t, tb)) 
               {
            if (loop_g == 0.0) 
            {
                // ga == 0.0 and gb may or may not be 0.0
                // handle the root at ta first
                before_root_t = loop_t;
                before_root_g = loop_g;
                after_root_t = min_time(shifted_by(before_root_t, convergence), tb);
                after_root_g = f.value(after_root_t);
            }
else 
            {
                // both non-zero, the usual case, use a root finder.
                if (forward) 
                {
                    const Interval interval =
                            solver.solve_interval(max_iteration_count, f, loop_t, tb);
                    before_root_t = interval.get_left_abscissa();
                    before_root_g = interval.get_left_value();
                    after_root_t = interval.get_right_abscissa();
                    after_root_g = interval.get_right_value();
                }
else 
                {
                    const Interval interval =
                            solver.solve_interval(max_iteration_count, f, tb, loop_t);
                    before_root_t = interval.get_right_abscissa();
                    before_root_g = interval.get_right_value();
                    after_root_t = interval.get_left_abscissa();
                    after_root_g = interval.get_left_value();
                }
            }
            // tolerance is set to less than 1 ulp
            // assume tolerance is 1 ulp
            if (before_root_t == after_root_t) 
            {
                after_root_t = next_after(after_root_t);
                after_root_g = f.value(after_root_t);
            }
            // check loop is making some progress
            check((forward && after_root_t > before_root_t) || (!forward && after_root_t < before_root_t));
            // setup next iteration
            loop_t = after_root_t;
            loop_g = after_root_g;
        }

        // figure out the result of root finding, and return accordingly
        if (after_root_g == 0.0 || after_root_g > 0.0 == g0_positive) 
        {
            // loop gave up and didn't find any crossing within this step
            return false;
        }
else 
        {
            // real crossing
            check(!std::isnan(before_root_t) && !std::isnan(before_root_g));
            // variation direction, with respect to the integration direction
            increasing = !g0_positive;
            pending_event_time = before_root_t;
            stop_time = before_root_g == 0.0 ? before_root_t : after_root_t;
            pending_event = true;
            after_event = after_root_t;
            after_g = after_root_g;

            // check increasing set correctly
            check(after_g > 0 == increasing);
            check(increasing == gb >= ga);

            return true;
        }

    }

    /**
     * Get the next number after the given number in the current propagation direction.
     *
     * @param t input time
     * @return t +/- 1 ulp depending on the direction.
     */
    private double next_after(const double t) 
    {
        // direction
        const double dir = forward ? INFINITY : -INFINITY;
        return FastMath.next_after(t, dir);
    }

    /**
     * Get the occurrence time of the event triggered in the current step.
     *
     * @return occurrence time of the event triggered in the current step or infinity if
     * no events are triggered
     */
    public double get_event_time() 
    {
        return pending_event ?
               pending_event_time :
               (forward ? INFINITY : -INFINITY);
    }

    /**
     * Try to accept the current history up to the given time.
     *
     * <p> It is not necessary to call this method before calling {@link
     * #do_event(ODE_State_And_Derivative)} with the same state. It is necessary to call this
     * method before you call {@link #do_event(ODE_State_And_Derivative)} on some other event
     * detector.
     *
     * @param state        to try to accept.
     * @param interpolator to use to find the root, if any.
     * @return if the event detector has an event it has not detected before that is on or
     * before the same time as {@code state}. In other words {@code false} means continue
     * on while {@code true} means stop and handle my event first.
     */
    public bool try_advance(const ODE_State_And_Derivative state, const ODE_StateInterpolator interpolator) 
    {
        const double t = state.get_time();
        // check this is only called before a pending event.
        check(!pending_event || !strictly_after(pending_event_time, t));

        const bool me_first;

        if (strictly_after(t, earliest_time_considered)) 
        {
            // just found an event and we know the next time we want to search again
            me_first = false;
        }
else 
        {
            // check g function to see if there is a event
            const double g = handler.g(state);
            const bool positive = g > 0;

            if (positive == g0_positive) 
            {
                // g function has expected sign
                g0 = g; // g0_positive is the same
                me_first = false;
            }
else 
            {
                // found a root we didn't expect -> find precise location
                const double old_pending_event_time = pending_event_time;
                const bool found_root = find_root(interpolator, t0, g0, t, g);
                // make sure the root is not the same as the old root, if one exists
                me_first = found_root &&
                          (std::isnan(old_pending_event_time) || old_pending_event_time != pending_event_time);
            }
        }

        if (!me_first) 
        {
            // advance t0 to the current time so we can't find events that occur before t
            t0 = t;
        }

        return me_first;
    }

    /**
     * Notify the user's listener of the event. The event occurs wholly within this method
     * call including a call to {@link ODE_Event_Handler#reset_state(ODE_State_And_Derivative)}
     * if necessary.
     *
     * @param state the state at the time of the event. This must be at the same time as
     *              the current value of {@link #get_event_time()}.
     * @return the user's requested action and the state if the action is {@link
     * Action#RESET_STATE}. Otherwise the state is {@code state}. The stop time
     * indicates what time propagation should stop if the action is {@link Action#STOP}.
     * This guarantees the integration will stop on or after the root, so that integration
     * may be restarted safely.
     */
    public Event_Occurrence do_event(const ODE_State_And_Derivative state) 
    {
        // check event is pending and is at the same time
        check(pending_event);
        check(state.get_time() == this.pending_event_time);

        const Action action = handler.event_occurred(state, increasing == forward);
        const ODE_State new_state;
        if (action == Action.RESET_STATE) 
        {
            new_state = handler.reset_state(state);
        }
else 
        {
            new_state = state;
        }
        // clear pending event
        pending_event = false;
        pending_event_time = std::numeric_limits<double>::quiet_NaN();
        // setup for next search
        earliest_time_considered = after_event;
        t0 = after_event;
        g0 = after_g;
        g0_positive = increasing;
        // check g0_positive set correctly
        check(g0 == 0.0 || g0_positive == (g0 > 0));
        return Event_Occurrence(action, new_state, stop_time);
    }

    /**
     * Shift a time value along the current integration direction: {@link #forward}.
     *
     * @param t     the time to shift.
     * @param delta the amount to shift.
     * @return t + delta if forward, else t - delta. If the result has to be rounded it
     * will be rounded to be before the true value of t + delta.
     */
    private double shifted_by(const double t, const double delta) 
    {
        if (forward) 
        {
            const double ret = t + delta;
            if (ret - t > delta) 
            {
                return FastMath.next_down(ret);
            }
else 
            {
                return ret;
            }
        }
else 
        {
            const double ret = t - delta;
            if (t - ret > delta) 
            {
                return FastMath.next_up(ret);
            }
else 
            {
                return ret;
            }
        }
    }

    /**
     * Get the time that happens first along the current propagation direction: {@link
     * #forward}.
     *
     * @param a first time
     * @param b second time
     * @return min(a, b) if forward, else max (a, b)
     */
    private double min_time(const double& a, const double& b) 
    {
        return forward ? std::min(a, b) : std::max(a, b);
    }

    /**
     * Check the ordering of two times.
     *
     * @param t1 the first time.
     * @param t2 the second time.
     * @return true if {@code t2} is strictly after {@code t1} in the propagation
     * direction.
     */
    private bool strictly_after(const double t1, const double t2) 
    {
        return forward ? t1 < t2 : t2 < t1;
    }

    /**
     * Same as keyword assert, but throw a {@link Math_Runtime_Exception}.
     *
     * @param condition to check
     * @Math_Runtime_Exception if {@code condition} is false.
     */
    private void check(const bool condition) Math_Runtime_Exception 
    {
        if (!condition) 
        {
            throw Math_Runtime_Exception.create_internal_error();
        }
    }

    /**
     * Class to hold the data related to an event occurrence that is needed to decide how
     * to modify integration.
     */
    public static class Event_Occurrence 
    {

        /** User requested action. */
        private const Action action;
        /** New state for a reset action. */
        private const ODE_State new_state;
        /** The time to stop propagation if the action is a stop event. */
        private const double stop_time;

        /**
         * Create a occurrence of an event.
         *
         * @param action   the user requested action.
         * @param new_state for a reset event. Should be the current state unless the
         *                 action is {@link Action#RESET_STATE}.
         * @param stop_time to stop propagation if the action is {@link Action#STOP}. Used
         *                 to move the stop time to just after the root.
         */
        Event_Occurrence(const Action action, const ODE_State new_state, const double stop_time) 
        {
            this.action = action;
            this.new_state = new_state;
            this.stop_time = stop_time;
        }

        /**
         * Get the user requested action.
         *
         * @return the action.
         */
        public Action get_action() 
        {
            return action;
        }

        /**
         * Get the state for a reset action.
         *
         * @return the state.
         */
        public ODE_State get_new_state() 
        {
            return new_state;
        }

        /**
         * Get the time for a stop action.
         *
         * @return when to stop propagation.
         */
        public double get_stop_time() 
        {
            return stop_time;
        }

    }

}


