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
//import java.util.List;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.ode.sampling.FieldODE_Step_Handler;
//import org.hipparchus.ode.sampling.FieldODE_StateInterpolator;
//import org.hipparchus.util.FastMath;
#include <type_traits>
#include "../core/CalculusFieldElement.h"

/**
 * This class stores all information provided by an ODE integrator
 * during the integration process and build a continuous model of the
 * solution from this.
 *
 * <p>This class act as a step handler from the integrator point of
 * view. It is called iteratively during the integration process and
 * stores a copy of all steps information in a sorted collection for
 * later use. Once the integration process is over, the user can use
 * the {@link #get_interpolated_state(Calculus_Field_Element) get_interpolated_state}
 * method to retrieve this information at any time. It is important to wait
 * for the integration to be over before attempting to call {@link
 * #get_interpolated_state(Calculus_Field_Element)} because some internal
 * variables are set only once the last step has been handled.</p>
 *
 * <p>This is useful for example if the main loop of the user
 * application should remain independent from the integration process
 * or if one needs to mimic the behaviour of an analytical model
 * despite a numerical model is used (i.e. one needs the ability to
 * get the model value at any time or to navigate through the
 * data).</p>
 *
 * <p>If problem modeling is done with several separate
 * integration phases for contiguous intervals, the same
 * FieldDense_Output_Model can be used as step handler for all
 * integration phases as long as they are performed in order and in
 * the same direction. As an example, one can extrapolate the
 * trajectory of a satellite with one model (i.e. one set of
 * differential equations) up to the beginning of a maneuver, use
 * another more complex model including thrusters modeling and
 * accurate attitude control during the maneuver, and revert to the
 * first model after the end of the maneuver. If the same continuous
 * output model handles the steps of all integration phases, the user
 * do not need to bother when the maneuver begins or ends, he has all
 * the data available in a transparent manner.</p>
 *
 * <p>One should be aware that the amount of data stored in a
 * FieldDense_Output_Model instance can be important if the state vector
 * is large, if the integration interval is long or if the steps are
 * small (which can result from small tolerance settings in {@link
 * org.hipparchus.ode.nonstiff.Adaptive_Stepsize_Field_Integrator adaptive
 * step size integrators}).</p>
 *
 * @see FieldODE_Step_Handler
 * @see FieldODE_StateInterpolator
 * @param <T> the type of the field elements
 */

template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class FieldDense_Output_Model : public FieldODE_Step_Handler<T> 
{

    /** Initial integration time. */
    private T initial_time;

    /** Final integration time. */
    private T const_time;

    /** Integration direction indicator. */
    private bool forward;

    /** Current interpolator index. */
    private int index;

    /** Steps table. */
    private List<FieldODE_StateInterpolator<T>> steps;

    /** Simple constructor.
     * Build an empty continuous output model.
     */
    public FieldDense_Output_Model() 
    {
        steps       = Array_list<>();
        initial_time = NULL;
        const_time   = NULL;
        forward     = true;
        index       = 0;
    }

    /** Append another model at the end of the instance.
     * @param model model to add at the end of the instance
     * @exception  if the model to append is not
     * compatible with the instance (dimension of the state vector, * propagation direction, hole between the dates)
     * @exception  if the dimensions of the states or
     * the number of secondary states do not match
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * during step constization
     */
    public void append(const FieldDense_Output_Model<T> model)
        , Math_Illegal_State_Exception 
        {

        if (model.steps.is_empty()) 
        {
            return;
        }

        if (steps.is_empty()) 
        {
            initial_time = model.initial_time;
            forward     = model.forward;
        }
else 
        {

            // safety checks
            const Field_ODE_State_And_Derivative<T> s1 = steps.get(0).get_previous_state();
            const Field_ODE_State_And_Derivative<T> s2 = model.steps.get(0).get_previous_state();
            check_dimensions_equality(s1.get_primary_state_dimension(), s2.get_primary_state_dimension());
            check_dimensions_equality(s1.get_number_of_secondary_states(), s2.get_number_of_secondary_states());
            for (int i{}; i < s1.get_number_of_secondary_states(); ++i) 
            {
                check_dimensions_equality(s1.get_secondary_state_dimension(i), s2.get_secondary_state_dimension(i));
            }

            if (forward ^ model.forward) 
            {
                throw (Localized_ODE_Formats.PROPAGATION_DIRECTION_MISMATCH);
            }

            const FieldODE_StateInterpolator<T> last_interpolator = steps.get(index);
            const T current  = last_interpolator.get_current_state().get_time();
            const T previous = last_interpolator.get_previous_state().get_time();
            const T step = current.subtract(previous);
            const T gap = model.get_initial_time().subtract(current);
            if (gap.abs().subtract(step.abs().multiply(1.0e-3)).get_real() > 0) 
            {
                throw (Localized_ODE_Formats.HOLE_BETWEEN_MODELS_TIME_RANGES, gap.norm());
            }

        }

        for (FieldODE_StateInterpolator<T> interpolator : model.steps) 
        {
            steps.add(interpolator);
        }

        index = steps.size() - 1;
        const_time = (steps.get(index)).get_current_state().get_time();

    }

    /** Check dimensions equality.
     * @param d1 first dimension
     * @param d2 second dimansion
     * @exception  if dimensions do not match
     */
    private void check_dimensions_equality(const int d1, const int d2)
         
        {
        if (d1 != d2) 
        {
            throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, d2, d1);
        }
    }

    /** {@inherit_doc} */
    //override
    public void init(const Field_ODE_State_And_Derivative<T> initial_state, const T t) 
    {
        initial_time = initial_state.get_time();
        const_time   = t;
        forward     = true;
        index       = 0;
        steps.clear();
    }

    /** {@inherit_doc} */
    //override
    public void handle_step(const FieldODE_StateInterpolator<T> interpolator) 
    {

        if (steps.is_empty()) 
        {
            initial_time = interpolator.get_previous_state().get_time();
            forward     = interpolator.is_forward();
        }

        steps.add(interpolator);
    }

    /** {@inherit_doc} */
    //override
    public void finish(Field_ODE_State_And_Derivative<T> const_state) 
    {
        const_time = const_state.get_time();
        index     = steps.size() - 1;
    }

    /**
     * Get the initial integration time.
     * @return initial integration time
     */
    public T get_initial_time() 
    {
        return initial_time;
    }

    /**
     * Get the const integration time.
     * @return const integration time
     */
    public T get_final_time() 
    {
        return const_time;
    }

    /**
     * Get the state at interpolated time.
     * @param time time of the interpolated point
     * @return state at interpolated time
     */
    public Field_ODE_State_And_Derivative<T> get_interpolated_state(const T time) 
    {

        // initialize the search with the complete steps table
        int i_min = 0;
        const FieldODE_StateInterpolator<T> s_min = steps.get(i_min);
        T t_min = s_min.get_previous_state().get_time().add(s_min.get_current_state().get_time()).multiply(0.5);

        int i_max = steps.size() - 1;
        const FieldODE_StateInterpolator<T> s_max = steps.get(i_max);
        T t_max = s_max.get_previous_state().get_time().add(s_max.get_current_state().get_time()).multiply(0.5);

        // handle points outside of the integration interval
        // or in the first and last step
        if (locate_point(time, s_min) <= 0) 
        {
            index = i_min;
            return s_min.get_interpolated_state(time);
        }
        if (locate_point(time, s_max) >= 0) 
        {
            index = i_max;
            return s_max.get_interpolated_state(time);
        }

        // reduction of the table slice size
        while (i_max - i_min > 5) 
        {

            // use the last estimated index as the splitting index
            const FieldODE_StateInterpolator<T> si = steps.get(index);
            const int location = locate_point(time, si);
            if (location < 0) 
            {
                i_max = index;
                t_max = si.get_previous_state().get_time().add(si.get_current_state().get_time()).multiply(0.5);
            }
else if (location > 0) 
            {
                i_min = index;
                t_min = si.get_previous_state().get_time().add(si.get_current_state().get_time()).multiply(0.5);
            }
else 
            {
                // we have found the target step, no need to continue searching
                return si.get_interpolated_state(time);
            }

            // compute a estimate of the index in the reduced table slice
            const int i_med = (i_min + i_max) / 2;
            const FieldODE_StateInterpolator<T> s_med = steps.get(i_med);
            const T t_med = s_med.get_previous_state().get_time().add(s_med.get_current_state().get_time()).multiply(0.5);

            if (t_med.subtract(t_min).abs().subtract(1.0e-6).get_real() < 0 ||
                t_max.subtract(t_med).abs().subtract(1.0e-6).get_real() < 0) 
                {
                // too close to the bounds, we estimate using a simple dichotomy
                index = i_med;
            }
else 
            {
                // estimate the index using a reverse quadratic polynomial
                // (reverse means we have i = P(t), thus allowing to simply
                // compute index = P(time) rather than solving a quadratic equation)
                const T d12 = t_max.subtract(t_med);
                const T d23 = t_med.subtract(t_min);
                const T d13 = t_max.subtract(t_min);
                const T dt1 = time.subtract(t_max);
                const T dt2 = time.subtract(t_med);
                const T dt3 = time.subtract(t_min);
                const T i_lagrange =           dt2.multiply(dt3).multiply(d23).multiply(i_max).
                                     subtract(dt1.multiply(dt3).multiply(d13).multiply(i_med)).
                                     add(     dt1.multiply(dt2).multiply(d12).multiply(i_min)).
                                     divide(d12.multiply(d23).multiply(d13));
                index = static_cast<int>( FastMath.rint(i_lagrange.get_real());
            }

            // force the next size reduction to be at least one tenth
            const int low  = std::max(i_min + 1, (9 * i_min + i_max) / 10);
            const int high = std::min(i_max - 1, (i_min + 9 * i_max) / 10);
            if (index < low) 
            {
                index = low;
            }
else if (index > high) 
            {
                index = high;
            }

        }

        // now the table slice is very small, we perform an iterative search
        index = i_min;
        while (index <= i_max && locate_point(time, steps.get(index)) > 0) 
        {
            ++index;
        }

        return steps.get(index).get_interpolated_state(time);

    }

    /** Compare a step interval and a double.
     * @param time point to locate
     * @param interval step interval
     * @return -1 if the double is before the interval, 0 if it is in
     * the interval, and +1 if it is after the interval, according to
     * the interval direction
     */
    private int locate_point(const T time, const FieldODE_StateInterpolator<T> interval) 
    {
        if (forward) 
        {
            if (time.subtract(interval.get_previous_state().get_time()).get_real() < 0) 
            {
                return -1;
            }
else if (time.subtract(interval.get_current_state().get_time()).get_real() > 0) 
            {
                return +1;
            }
else 
            {
                return 0;
            }
        }
        if (time.subtract(interval.get_previous_state().get_time()).get_real() > 0) 
        {
            return -1;
        }
else if (time.subtract(interval.get_current_state().get_time()).get_real() < 0) 
        {
            return +1;
        }
else 
        {
            return 0;
        }
    }

};