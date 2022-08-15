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

  //import java.util.Array_list;
  //import java.util.Collection;
  //import java.util.Collections;
  //import java.util.Comparator;
  //import java.util.List;
  //import java.util.Priority_Queue;
  //import java.util.Queue;

  //import org.hipparchus.Calculus_Field_Element;
  //import org.hipparchus.Field;
  //import org.hipparchus.analysis.solvers.Bracketed_Real_Field_Univariate_Solver;
  //import org.hipparchus.analysis.solvers.FieldBracketing_Nth_Order_Brent_Solver;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.ode.events.Action;
  //import org.hipparchus.ode.events.Field_Event_Handler_Configuration;
  //import org.hipparchus.ode.events.Field_Event_State;
  //import org.hipparchus.ode.events.Field_Event_State.Event_Occurrence;
  //import org.hipparchus.ode.events.FieldODE_Event_Handler;
  //import org.hipparchus.ode.sampling.AbstractFieldODE_StateInterpolator;
  //import org.hipparchus.ode.sampling.FieldODE_Step_Handler;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Incrementor;
#include <type_traits>

  /**
   * Base class managing common boilerplate for all integrators.
   * @param <T> the type of the field elements
   */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Abstract_Field_Integrator : FieldODE_Integrator<T>
{
	/** Default relative accuracy. */
	private static const double DEFAULT_RELATIVE_ACCURACY = 0;

	/** Default function value accuracy. */
	private static const double DEFAULT_FUNCTION_VALUE_ACCURACY = 0;

	/** Step handler. */
	private Collection<FieldODE_Step_Handler<T>> step_handlers;

	/** Current step start. */
	private Field_ODE_State_And_Derivative<T> step_start;

	/** Current stepsize. */
	private T step_size;

	/** Indicator for last step. */
	private bool is_last_step;

	/** Indicator that a state or derivative reset was triggered by some event. */
	private bool reset_occurred;

	/** Field to which the time and state vector elements belong. */
	private const Field<T> field;

	/** Events states. */
	private Collection<Field_Event_State<T>> events_states;

	/** Initialization indicator of events states. */
	private bool states_initialized;

	/** Name of the method. */
	private const std::string name;

	/** Counter for number of evaluations. */
	private Incrementor evaluations;

	/** Differential equations to integrate. */
	private transient FieldExpandable_ODE<T> equations;

	/** Build an instance.
	 * @param field field to which the time and state vector elements belong
	 * @param name name of the method
	 */
	protected Abstract_Field_Integrator(const Field<T> field, const std::string name)
	{
		this.field = field;
		this.name = name;
		step_handlers = Array_list<>();
		step_start = NULL;
		step_size = NULL;
		events_states = Array_list<>();
		states_initialized = false;
		evaluations = Incrementor();
	}

	/** Get the field to which state vector elements belong.
	 * @return field to which state vector elements belong
	 */
	public Field<T> get_field()
	{
		return field;
	}

	/** {@inherit_doc} */
	//override
	public std::string get_name()
	{
		return name;
	}

	/** {@inherit_doc} */
	//override
	public void add_step_handler(const FieldODE_Step_Handler<T> handler)
	{
		step_handlers.add(handler);
	}

	/** {@inherit_doc} */
	//override
	public Collection<FieldODE_Step_Handler<T>> get_step_handlers()
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
	public void add_event_handler(const FieldODE_Event_Handler<T> handler, const double max_check_interval, const double convergence, const int max_iteration_count)
	{
		add_event_handler(handler, max_check_interval, convergence, max_iteration_count, FieldBracketing_Nth_Order_Brent_Solver<T>(field.get_zero().add(DEFAULT_RELATIVE_ACCURACY), field.get_zero().add(convergence), field.get_zero().add(DEFAULT_FUNCTION_VALUE_ACCURACY), 5));
	}

	/** {@inherit_doc} */
	//override
	public void add_event_handler(const FieldODE_Event_Handler<T> handler, const double max_check_interval, const double convergence, const int max_iteration_count, const Bracketed_Real_Field_Univariate_Solver<T> solver)
	{
		events_states.add(new Field_Event_State<T>(handler, max_check_interval, field.get_zero().new_instance(convergence), max_iteration_count, solver));
	}

	/** {@inherit_doc} */
	//override
	public Collection<FieldODE_Event_Handler<T>> get_event_handlers()
	{
		const List<FieldODE_Event_Handler<T>> list = Array_list<>(events_states.size());
		for (Field_Event_State<T> state : events_states)
		{
			list.add(state.get_event_handler());
		}
		return Collections.unmodifiable_collection(list);
	}

	/** {@inherit_doc} */
	//override
	public Collection<Field_Event_Handler_Configuration<T>> get_event_handlers_configurations()
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
	public T get_current_signed_stepsize()
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

	/** Prepare the start of an integration.
	 * @param eqn equations to integrate
	 * @param s0 initial state vector
	 * @param t target time for the integration
	 * @return initial state with derivatives added
	 */
	protected Field_ODE_State_And_Derivative<T> init_integration(const FieldExpandable_ODE<T> eqn, const FieldODE_State<T> s0, const T t)
	{
		this.equations = eqn;
		evaluations = evaluations.with_count(0);

		// initialize ODE
		eqn.init(s0, t);

		// set up derivatives of initial state (including primary and secondary components)
		const T   t0 = s0.get_time();
		const std::vector<T> y0 = s0.get_complete_state();
		const std::vector<T> y0_dot = compute_derivatives(t0, y0);

		// built the state
		const Field_ODE_State_And_Derivative<T> s0_with_derivatives =
			eqn.get_mapper().map_state_and_derivative(t0, y0, y0_dot);

		// initialize event handlers
		for (const Field_Event_State<T> state : events_states)
		{
			state.get_event_handler().init(s0_with_derivatives, t);
		}

		// initialize step handlers
		for (FieldODE_Step_Handler<T> handler : step_handlers)
		{
			handler.init(s0_with_derivatives, t);
		}

		set_state_initialized(false);

		return s0_with_derivatives;
	}

	/** Get the differential equations to integrate.
	 * @return differential equations to integrate
	 */
	protected FieldExpandable_ODE<T> get_equations()
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
	 * is called outside of a call to {@link #integrate(FieldExpandable_ODE, FieldODE_State, * Calculus_Field_Element) integrate}
	 */
	public std::vector<T> compute_derivatives(const T t, const std::vector<T> y)
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

	/**
	 * Accept a step, triggering events and step handlers.
	 *
	 * @param interpolator step interpolator
	 * @param t_end         const integration time
	 * @return state at end of step
	 * @Math_Illegal_State_Exception    if the interpolator one because the
	 *                                      number of functions evaluations is exceeded
	 * @ if the location of an event cannot be
	 *                                      bracketed
	 * @ if arrays dimensions do not match equations
	 *                                      settings
	 */
	protected Field_ODE_State_And_Derivative<T> accept_step(const AbstractFieldODE_StateInterpolator<T> interpolator, const T t_end)
		, Math_Illegal_State_Exception
	{
Field_ODE_State_And_Derivative<T> previous_state = interpolator.get_global_previous_state();
const Field_ODE_State_And_Derivative<T> current_state = interpolator.get_global_current_state();
AbstractFieldODE_StateInterpolator<T> restricted = interpolator;

// initialize the events states if needed
if (!states_initialized)
{
	for (Field_Event_State<T> state : events_states)
	{
		state.reinitialize_begin(interpolator);
	}
	states_initialized = true;
}

// search for next events that may occur during the step
const int ordering_sign = interpolator.is_forward() ? +1 : -1;
const Queue<Field_Event_State<T>> occurring_events = Priority_Queue<>(new Comparator<Field_Event_State<T>>()
{
	/** {@inherit_doc} */
	//override
	public int compare(Field_Event_State<T> es0, Field_Event_State<T> es1)
	{
		return ordering_sign * Double.compare(es0.get_event_time().get_real(), es1.get_event_time().get_real());
	}
});

reset_occurred = false;
bool done_with_step = false;
reset_events:
do
{
	// Evaluate all event detectors for events
	occurring_events.clear();
	for (const Field_Event_State<T> state : events_states)
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
			const Field_Event_State<T> current_event = occurring_events.poll();

			// get state at event time
			Field_ODE_State_And_Derivative<T> event_state =
					restricted.get_interpolated_state(current_event.get_event_time());

			// restrict the interpolator to the first part of the step, up to the event
			restricted = restricted.restrict_step(previous_state, event_state);

			// try to advance all event states to current time
			for (const Field_Event_State<T> state : events_states)
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
			for (const FieldODE_Step_Handler<T> handler : step_handlers)
			{
				handler.handle_step(restricted);
			}

			// acknowledge event occurrence
			const Event_Occurrence<T> occurrence = current_event.do_event(event_state);
			const Action action = occurrence.get_action();
			is_last_step = action == Action.STOP;

			if (is_last_step)
			{
				// ensure the event is after the root if it is returned STOP
				// this lets the user integrate to a STOP event and then restart
				// integration from the same time.
				const Field_ODE_State_And_Derivative<T> saved_state = event_state;
				event_state = interpolator.get_interpolated_state(occurrence.get_stop_time());
				restricted = interpolator.restrict_step(saved_state, event_state);

				// handle the almost zero size last part of the const step, at event time
				for (const FieldODE_Step_Handler<T> handler : step_handlers)
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
				const FieldODE_State<T> new_state = occurrence.get_new_state();
				const std::vector<T> y = new_state.get_complete_state();
				const std::vector<T> y_dot = compute_derivatives(new_state.get_time(), y);
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

		// last part of the step, after the last event
		// may be a event here if the last event modified the g function of
		// another event detector.
		for (const Field_Event_State<T> state : events_states)
		{
			if (state.try_advance(current_state, interpolator))
			{
				occurring_events.add(state);
			}
		}
	} while (!occurring_events.is_empty());

	done_with_step = true;
} while (!done_with_step);

is_last_step = is_last_step || current_state.get_time().subtract(t_end).norm() <= FastMath.ulp(t_end.get_real());

// handle the remaining part of the step, after all events if any
for (FieldODE_Step_Handler<T> handler : step_handlers)
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
		protected void sanity_checks(const FieldODE_State<T> initial_state, const T t)

	{
		const double threshold = 1000 * FastMath.ulp(std::max(std::abs(initial_state.get_time().get_real()), std::abs(t.get_real())));
		const double dt = initial_state.get_time().subtract(t).norm();
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
	protected void set_step_size(const T step_size)
	{
		this.step_size = step_size;
	}

	/** Get the current step size.
	 * @return current step size
	 */
	protected T get_step_size()
	{
		return step_size;
	}

	/** Set current step start.
	 * @param step_start step start
	 */
	protected void set_step_start(const Field_ODE_State_And_Derivative<T> step_start)
	{
		this.step_start = step_start;
	}

	/**  {@inherit_doc} */
	//override
	public Field_ODE_State_And_Derivative<T> get_step_start()
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