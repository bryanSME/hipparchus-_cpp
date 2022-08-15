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

//import java.io.Serializable;
//import java.util.Array_list;
//import java.util.List;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.ode.sampling.ODE_StateInterpolator;
//import org.hipparchus.ode.sampling.ODE_Step_Handler;
//import org.hipparchus.util.FastMath;

/**
 * This class stores all information provided by an ODE integrator
 * during the integration process and build a continuous model of the
 * solution from this.
 *
 * <p>This class act as a step handler from the integrator point of
 * view. It is called iteratively during the integration process and
 * stores a copy of all steps information in a sorted collection for
 * later use. Once the integration process is over, the user can use
 * the {@link #get_interpolated_statestatic_cast<double>( get_interpolated_state}
 * method to retrieve this information at any time. It is important
 * to wait for the integration to be over before attempting to call
 * {@link #get_interpolated_statestatic_cast<double>( get_interpolated_state} because
 * some internal variables are set only once the last step has been
 * handled.</p>
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
 * Dense_Output_Model can be used as step handler for all
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
 * <p>An important feature of this class is that it : the
 * <code>Serializable</code> interface. This means that the result of
 * an integration can be serialized and reused later (if stored into a
 * persistent medium like a filesystem or a database) or elsewhere (if
 * sent to another application). Only the result of the integration is
 * stored, there is no reference to the integrated problem by
 * itself.</p>
 *
 * <p>One should be aware that the amount of data stored in a
 * Dense_Output_Model instance can be important if the state vector
 * is large, if the integration interval is long or if the steps are
 * small (which can result from small tolerance settings in {@link
 * org.hipparchus.ode.nonstiff.Adaptive_Stepsize_Integrator adaptive
 * step size integrators}).</p>
 *
 * @see ODE_Step_Handler
 * @see ODE_StateInterpolator
 */

class Dense_Output_Model : ODE_Step_Handler
{

    /** Serializable version identifier */
    private static const long serial_version_uid = 20160328L;

    /** Initial integration time. */
    private double initial_time;

    /** Final integration time. */
    private double const_time;

    /** Integration direction indicator. */
    private bool forward;

    /** Current interpolator index. */
    private int index;

    /** Steps table. */
    private List<ODE_StateInterpolator> steps;

    /** Simple constructor.
     * Build an empty continuous output model.
     */
    public Dense_Output_Model() 
    {
        steps       = Array_list<>();
        initial_time = std::numeric_limits<double>::quiet_NaN();
        const_time   = std::numeric_limits<double>::quiet_NaN();
        forward     = true;
        index       = 0;
    }

    /** Append another model at the end of the instance.
     * @param model model to add at the end of the instance
     * @exception  if the model to append is not
     * compatible with the instance (dimension of the state vector, * propagation direction, hole between the dates)
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * during step constization
     */
    public void append(const Dense_Output_Model model)
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

            const ODE_State_And_Derivative s1 = steps.get(0).get_previous_state();
            const ODE_State_And_Derivative s2 = model.steps.get(0).get_previous_state();
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

            const ODE_StateInterpolator last_interpolator = steps.get(index);
            const double current  = last_interpolator.get_current_state().get_time();
            const double previous = last_interpolator.get_previous_state().get_time();
            const double step = current - previous;
            const double gap = model.get_initial_time() - current;
            if (std::abs(gap) > 1.0e-3 * std::abs(step)) 
            {
                throw (Localized_ODE_Formats.HOLE_BETWEEN_MODELS_TIME_RANGES, std::abs(gap));
            }

        }

        for (ODE_StateInterpolator interpolator : model.steps) 
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
    public void init(const ODE_State_And_Derivative initial_state, const double target_time) 
    {
        initial_time    = initial_state.get_time();
        this.const_time = target_time;
        forward        = true;
        index          = 0;
        steps.clear();
    }

    /** {@inherit_doc} */
    //override
    public void handle_step(const ODE_StateInterpolator interpolator) 
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
    public void finish(const ODE_State_And_Derivative const_state) 
    {
        const_time = const_state.get_time();
        index     = steps.size() - 1;
    }

    /**
     * Get the initial integration time.
     * @return initial integration time
     */
    public double get_initial_time() 
    {
        return initial_time;
    }

    /**
     * Get the const integration time.
     * @return const integration time
     */
    public double get_final_time() 
    {
        return const_time;
    }

    /**
     * Get the state at interpolated time.
     * @param time time of the interpolated point
     * @return state at interpolated time
     */
    public ODE_State_And_Derivative get_interpolated_state(const double time) 
    {

        // initialize the search with the complete steps table
        int i_min = 0;
        const ODE_StateInterpolator s_min = steps.get(i_min);
        double t_min = 0.5 * (s_min.get_previous_state().get_time() + s_min.get_current_state().get_time());

        int i_max = steps.size() - 1;
        const ODE_StateInterpolator s_max = steps.get(i_max);
        double t_max = 0.5 * (s_max.get_previous_state().get_time() + s_max.get_current_state().get_time());

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
            const ODE_StateInterpolator si = steps.get(index);
            const int location = locate_point(time, si);
            if (location < 0) 
            {
                i_max = index;
                t_max = 0.5 * (si.get_previous_state().get_time() + si.get_current_state().get_time());
            }
else if (location > 0) 
            {
                i_min = index;
                t_min = 0.5 * (si.get_previous_state().get_time() + si.get_current_state().get_time());
            }
else 
            {
                // we have found the target step, no need to continue searching
                return si.get_interpolated_state(time);
            }

            // compute a estimate of the index in the reduced table slice
            const int i_med = (i_min + i_max) / 2;
            const ODE_StateInterpolator s_med = steps.get(i_med);
            const double t_med = 0.5 * (s_med.get_previous_state().get_time() + s_med.get_current_state().get_time());

            if ((std::abs(t_med - t_min) < 1e-6) || (std::abs(t_max - t_med) < 1e-6)) 
            {
                // too close to the bounds, we estimate using a simple dichotomy
                index = i_med;
            }
else 
            {
                // estimate the index using a reverse quadratic polynom
                // (reverse means we have i = P(t), thus allowing to simply
                // compute index = P(time) rather than solving a quadratic equation)
                const double d12 = t_max - t_med;
                const double d23 = t_med - t_min;
                const double d13 = t_max - t_min;
                const double dt1 = time - t_max;
                const double dt2 = time - t_med;
                const double dt3 = time - t_min;
                const double i_lagrange = ((dt2 * dt3 * d23) * i_max -
                                (dt1 * dt3 * d13) * i_med +
                                (dt1 * dt2 * d12) * i_min) /
                                (d12 * d23 * d13);
                index = static_cast<int>( FastMath.rint(i_lagrange);
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
        while ((index <= i_max) && (locate_point(time, steps.get(index)) > 0)) 
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
    private int locate_point(const double time, const ODE_StateInterpolator interval) 
    {
        if (forward) 
        {
            if (time < interval.get_previous_state().get_time()) 
            {
                return -1;
            }
else if (time > interval.get_current_state().get_time()) 
            {
                return +1;
            }
else 
            {
                return 0;
            }
        }
        if (time > interval.get_previous_state().get_time()) 
        {
            return -1;
        }
else if (time < interval.get_current_state().get_time()) 
        {
            return +1;
        }
else 
        {
            return 0;
        }
    }

}


