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

//import java.util.Arrays;

//import org.hipparchus.ode.ODE_State;
//import org.hipparchus.ode.ODE_State_And_Derivative;

/** Wrapper used to detect only increasing or decreasing events.
 *
 * <p>General {@link ODE_Event_Handler events} are defined implicitly
 * by a {@link ODE_Event_Handler#g(ODE_State_And_Derivative) g function} crossing
 * zero. This function needs to be continuous in the event neighborhood, * and its sign must remain consistent between events. This implies that
 * during an ODE integration, events triggered are alternately events
 * for which the function increases from negative to positive values, * and events for which the function decreases from positive to
 * negative values.
 * </p>
 *
 * <p>Sometimes, users are only interested in one type of event (say
 * increasing events for example) and not in the other type. In these
 * cases, looking precisely for all events location and triggering
 * events that will later be ignored is a waste of computing time.</p>
 *
 * <p>Users can wrap a regular {@link ODE_Event_Handler event handler} in
 * an instance of this class and provide this wrapping instance to
 * the {@link org.hipparchus.ode.ODE_Integrator ODE solver}
 * in order to avoid wasting time looking for uninteresting events.
 * The wrapper will intercept the calls to the {@link
 * ODE_Event_Handler#g(ODE_State_And_Derivative) g function} and to the {@link
 * ODE_Event_Handler#event_occurred(ODE_State_And_Derivative, bool)
 * event_occurred} method in order to ignore uninteresting events. The
 * wrapped regular {@link ODE_Event_Handler event handler} will the see only
 * the interesting events, i.e. either only {@code increasing} events or
 * {@code decreasing} events. the number of calls to the {@link
 * ODE_Event_Handler#g(ODE_State_And_Derivative) g function} will also be reduced.</p>
 *
 */

class Event_Filter : ODE_Event_Handler 
{

    /** Number of past transformers updates stored. */
    private static const int HISTORY_SIZE = 100;

    /** Wrapped event handler. */
    private const ODE_Event_Handler raw_handler;

    /** Filter to use. */
    private const Filter_Type filter;

    /** Transformers of the g function. */
    private const Transformer[] transformers;

    /** Update time of the transformers. */
    private const std::vector<double> updates;

    /** Indicator for forward integration. */
    private bool forward;

    /** Extreme time encountered so far. */
    private double extreme_t;

    /** Wrap an {@link ODE_Event_Handler event handler}.
     * @param raw_handler event handler to wrap
     * @param filter filter to use
     */
    public Event_Filter(const ODE_Event_Handler raw_handler, const Filter_Type filter) 
    {
        this.raw_handler   = raw_handler;
        this.filter       = filter;
        this.transformers = Transformer[HISTORY_SIZE];
        this.updates      = std::vector<double>(HISTORY_SIZE];
    }

    /**  {@inherit_doc} */
    //override
    public void init(const ODE_State_And_Derivative initial_state, double const_time) 
    {

        // delegate to raw handler
        raw_handler.init(initial_state, const_time);

        // initialize events triggering logic
        forward  = const_time >= initial_state.get_time();
        extreme_t = forward ? -INFINITY : INFINITY;
        Arrays.fill(transformers, Transformer.UNINITIALIZED);
        Arrays.fill(updates, extreme_t);

    }

    /**  {@inherit_doc} */
    //override
    public double g(const ODE_State_And_Derivative state) 
    {

        const double raw_g = raw_handler.g(state);

        // search which transformer should be applied to g
        if (forward) 
        {
            const int last = transformers.size() - 1;
            if (extreme_t < state.get_time()) 
            {
                // we are at the forward end of the history

                // check if a rough root has been crossed
                const Transformer previous = transformers[last];
                const Transformer next     = filter.select_transformer(previous, raw_g, forward);
                if (next != previous) 
                {
                    // there is a root somewhere between extreme_t and t.
                    // the transformer is valid for t (this is how we have just computed
                    // it above), but it is in fact valid on both sides of the root, so
                    // it was already valid before t and even up to previous time. We store
                    // the switch at extreme_t for safety, to ensure the previous transformer
                    // is not applied too close of the root
                    System.arraycopy(updates,      1, updates,      0, last);
                    System.arraycopy(transformers, 1, transformers, 0, last);
                    updates[last]      = extreme_t;
                    transformers[last] = next;
                }

                extreme_t = state.get_time();

                // apply the transform
                return next.transformed(raw_g);

            }
else 
            {
                // we are in the middle of the history

                // select the transformer
                for (int i = last; i > 0; --i) 
                {
                    if (updates[i] <= state.get_time()) 
                    {
                        // apply the transform
                        return transformers[i].transformed(raw_g);
                    }
                }

                return transformers[0].transformed(raw_g);

            }
        }
else 
        {
            if (state.get_time() < extreme_t) 
            {
                // we are at the backward end of the history

                // check if a rough root has been crossed
                const Transformer previous = transformers[0];
                const Transformer next     = filter.select_transformer(previous, raw_g, forward);
                if (next != previous) 
                {
                    // there is a root somewhere between extreme_t and t.
                    // the transformer is valid for t (this is how we have just computed
                    // it above), but it is in fact valid on both sides of the root, so
                    // it was already valid before t and even up to previous time. We store
                    // the switch at extreme_t for safety, to ensure the previous transformer
                    // is not applied too close of the root
                    System.arraycopy(updates,      0, updates,      1, updates.size() - 1);
                    System.arraycopy(transformers, 0, transformers, 1, transformers.size() - 1);
                    updates[0]      = extreme_t;
                    transformers[0] = next;
                }

                extreme_t = state.get_time();

                // apply the transform
                return next.transformed(raw_g);

            }
else 
            {
                // we are in the middle of the history

                // select the transformer
                for (int i{}; i < updates.size() - 1; ++i) 
                {
                    if (state.get_time() <= updates[i]) 
                    {
                        // apply the transform
                        return transformers[i].transformed(raw_g);
                    }
                }

                return transformers[updates.size() - 1].transformed(raw_g);

            }
       }

    }

    /**  {@inherit_doc} */
    //override
    public Action event_occurred(const ODE_State_And_Derivative state, const bool increasing) 
    {
        // delegate to raw handler, fixing increasing status on the fly
        return raw_handler.event_occurred(state, filter.is_triggered_on_increasing());
    }

    /**  {@inherit_doc} */
    //override
    public ODE_State reset_state(const ODE_State_And_Derivative state) 
    {
        // delegate to raw handler
        return raw_handler.reset_state(state);
    }

}


