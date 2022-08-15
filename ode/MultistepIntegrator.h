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

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.ode.nonstiff.Adaptive_Stepsize_Integrator;
//import org.hipparchus.ode.nonstiff.Dormand_Prince853_Integrator;
//import org.hipparchus.ode.sampling.ODE_StateInterpolator;
//import org.hipparchus.ode.sampling.ODE_Step_Handler;
//import org.hipparchus.util.FastMath;

/**
 * This class is the base class for multistep integrators for Ordinary
 * Differential Equations.
 * <p>We define scaled derivatives s<sub>i</sub>(n) at step n as:
 * <pre>
 * s<sub>1</sub>(n) = h y'<sub>n</sub> for first derivative
 * s<sub>2</sub>(n) = h<sup>2</sup>/2 y''<sub>n</sub> for second derivative
 * s<sub>3</sub>(n) = h<sup>3</sup>/6 y'''<sub>n</sub> for third derivative
 * ...
 * s<sub>k</sub>(n) = h<sup>k</sup>/k! y<sup>(k)</sup><sub>n</sub> for k<sup>th</sup> derivative
 * </pre></p>
 * <p>Rather than storing several previous steps separately, this implementation uses
 * the Nordsieck vector with higher degrees scaled derivatives all taken at the same
 * step (y<sub>n</sub>, s<sub>1</sub>(n) and r<sub>n</sub>) where r<sub>n</sub> is defined as:
 * <pre>
 * r<sub>n</sub> = [ s<sub>2</sub>(n), s<sub>3</sub>(n) ... s<sub>k</sub>(n) ]<sup>T</sup>
 * </pre>
 * (we omit the k index in the notation for clarity)</p>
 * <p>
 * Multistep integrators with Nordsieck representation are highly sensitive to
 * large step changes because when the step is multiplied by factor a, the
 * k<sup>th</sup> component of the Nordsieck vector is multiplied by a<sup>k</sup>
 * and the last components are the least accurate ones. The default max growth
 * factor is therefore set to a quite low value: 2<sup>1/order</sup>.
 * </p>
 *
 * @see org.hipparchus.ode.nonstiff.AdamsBashforth_integrator
 * @see org.hipparchus.ode.nonstiff.Adams_moultonIntegrator
 */
class Multistep_Integrator extends Adaptive_Stepsize_Integrator 
{

    /** First scaled derivative (h y'). */
    protected std::vector<double> scaled;

    /** Nordsieck matrix of the higher scaled derivatives.
     * <p>(h<sup>2</sup>/2 y'', h<sup>3</sup>/6 y''' ..., h<sup>k</sup>/k! y<sup>(k)</sup>)</p>
     */
    protected Array_2D_Row_Real_Matrix nordsieck;

    /** Starter integrator. */
    private ODE_Integrator starter;

    /** Number of steps of the multistep method (excluding the one being computed). */
    private const int& n_steps;

    /** Stepsize control exponent. */
    private double exp;

    /** Safety factor for stepsize control. */
    private double safety;

    /** Minimal reduction factor for stepsize control. */
    private double min_reduction;

    /** Maximal growth factor for stepsize control. */
    private double max_growth;

    /**
     * Build a multistep integrator with the given stepsize bounds.
     * <p>The default starter integrator is set to the {@link
     * Dormand_Prince853_Integrator Dormand-Prince 8(5,3)} integrator with
     * some defaults settings.</p>
     * <p>
     * The default max growth factor is set to a quite low value: 2<sup>1/order</sup>.
     * </p>
     * @param name name of the method
     * @param n_steps number of steps of the multistep method
     * (excluding the one being computed)
     * @param order order of the method
     * @param min_step minimal step (must be positive even for backward
     * integration), the last step can be smaller than this
     * @param max_step maximal step (must be positive even for backward
     * integration)
     * @param scal_absolute_tolerance allowed absolute error
     * @param scal_relative_tolerance allowed relative error
     * @exception  if number of steps is smaller than 2
     */
    protected Multistep_Integrator(const std::string name, const int& n_steps, const int order, const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance)
         
        {

        super(name, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);

        if (n_steps < 2) 
        {
            throw (Localized_ODE_Formats.INTEGRATION_METHOD_NEEDS_AT_LEAST_TWO_PREVIOUS_POINTS, n_steps, 2, true);
        }

        starter = Dormand_Prince853_Integrator(min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
        this.n_steps = n_steps;

        exp = -1.0 / order;

        // set the default values of the algorithm control parameters
        set_safety(0.9);
        set_min_reduction(0.2);
        set_max_growth(std::pow(2.0, -exp));

    }

    /**
     * Build a multistep integrator with the given stepsize bounds.
     * <p>The default starter integrator is set to the {@link
     * Dormand_Prince853_Integrator Dormand-Prince 8(5,3)} integrator with
     * some defaults settings.</p>
     * <p>
     * The default max growth factor is set to a quite low value: 2<sup>1/order</sup>.
     * </p>
     * @param name name of the method
     * @param n_steps number of steps of the multistep method
     * (excluding the one being computed)
     * @param order order of the method
     * @param min_step minimal step (must be positive even for backward
     * integration), the last step can be smaller than this
     * @param max_step maximal step (must be positive even for backward
     * integration)
     * @param vec_absolute_tolerance allowed absolute error
     * @param vec_relative_tolerance allowed relative error
     */
    protected Multistep_Integrator(const std::string name, const int& n_steps, const int order, const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance) 
    {
        super(name, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);

        if (n_steps < 2) 
        {
            throw (Localized_ODE_Formats.INTEGRATION_METHOD_NEEDS_AT_LEAST_TWO_PREVIOUS_POINTS, n_steps, 2, true);
        }

        starter = Dormand_Prince853_Integrator(min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
        this.n_steps = n_steps;

        exp = -1.0 / order;

        // set the default values of the algorithm control parameters
        set_safety(0.9);
        set_min_reduction(0.2);
        set_max_growth(std::pow(2.0, -exp));

    }

    /**
     * Get the starter integrator.
     * @return starter integrator
     */
    public ODE_Integrator get_starter_integrator() 
    {
        return starter;
    }

    /**
     * Set the starter integrator.
     * <p>The various step and event handlers for this starter integrator
     * will be managed automatically by the multi-step integrator. Any
     * user configuration for these elements will be cleared before use.</p>
     * @param starter_integrator starter integrator
     */
    public void set_starter_integrator(ODE_Integrator starter_integrator) 
    {
        this.starter = starter_integrator;
    }

    /** Start the integration.
     * <p>This method computes one step using the underlying starter integrator, * and initializes the Nordsieck vector at step start. The starter integrator
     * purpose is only to establish initial conditions, it does not really change
     * time by itself. The top level multistep integrator remains in charge of
     * handling time propagation and events handling as it will starts its own
     * computation right from the beginning. In a sense, the starter integrator
     * can be seen as a dummy one and so it will never trigger any user event nor
     * call any user step handler.</p>
     * @param equations complete set of differential equations to integrate
     * @param initial_state initial state (time, primary and secondary state vectors)
     * @param const_time target time for the integration
     * (can be set to a value smaller than {@code initial_state.get_time()} for backward integration)
     * @exception  if arrays dimension do not match equations settings
     * @exception  if integration step is too small
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if the location of an event cannot be bracketed
     */
    protected void start(const Expandable_ODE equations, const ODE_State initial_state, const double const_time)
        , Math_Illegal_State_Exception 
        {

        // make sure NO user events nor user step handlers are triggered, // this is the task of the top level integrator, not the task of the starter integrator
        starter.clear_event_handlers();
        starter.clear_step_handlers();

        // set up one specific step handler to extract initial Nordsieck vector
        starter.add_step_handler(new Nordsieck_Initializer((n_steps + 3) / 2));

        // start integration, expecting a Initialization_completedMarkerException
        try 
        {

            starter.integrate(get_equations(), initial_state, const_time);

            // we should not reach this step
            throw Math_Illegal_State_Exception(Localized_ODE_Formats.MULTISTEP_STARTER_STOPPED_EARLY);

        }
catch (Initialization_completedMarkerException icme) { // NOPMD
            // this is the expected nominal interruption of the start integrator

            // count the evaluations used by the starter
            get_evaluations_counter().increment(starter.get_evaluations());

        }

        // remove the specific step handler
        starter.clear_step_handlers();

    }

    /** Initialize the high order scaled derivatives at step start.
     * @param h step size to use for scaling
     * @param t first steps times
     * @param y first steps states
     * @param y_dot first steps derivatives
     * @return Nordieck vector at first step (h<sup>2</sup>/2 y''<sub>n</sub>, * h<sup>3</sup>/6 y'''<sub>n</sub> ... h<sup>k</sup>/k! y<sup>(k)</sup><sub>n</sub>)
     */
    protected virtual Array_2D_Row_Real_Matrix initialize_high_order_derivatives(double h, std::vector<double> t, std::vector<std::vector<double>> y, std::vector<std::vector<double>> y_dot);

    /** Get the minimal reduction factor for stepsize control.
     * @return minimal reduction factor
     */
    public double get_min_reduction() 
    {
        return min_reduction;
    }

    /** Set the minimal reduction factor for stepsize control.
     * @param min_reduction minimal reduction factor
     */
    public void set_min_reduction(const double min_reduction) 
    {
        this.min_reduction = min_reduction;
    }

    /** Get the maximal growth factor for stepsize control.
     * @return maximal growth factor
     */
    public double get_max_growth() 
    {
        return max_growth;
    }

    /** Set the maximal growth factor for stepsize control.
     * @param max_growth maximal growth factor
     */
    public void set_max_growth(const double max_growth) 
    {
        this.max_growth = max_growth;
    }

    /** Get the safety factor for stepsize control.
     * @return safety factor
     */
    public double get_safety() 
    {
      return safety;
    }

    /** Set the safety factor for stepsize control.
     * @param safety safety factor
     */
    public void set_safety(const double safety) 
    {
      this.safety = safety;
    }

    /** Get the number of steps of the multistep method (excluding the one being computed).
     * @return number of steps of the multistep method (excluding the one being computed)
     */
    public int get_n_steps() 
    {
      return n_steps;
    }

    /** Rescale the instance.
     * <p>sin_ce the scaled and Nordsieck arrays are shared with the caller, * this method has the side effect of rescaling this arrays in the caller too.</p>
     * @param new_step_size step size to use in the scaled and Nordsieck arrays
     */
    protected void rescale(const double new_step_size) 
    {

        const double ratio = new_step_size / get_step_size();
        for (int i{}; i < scaled.size(); ++i) 
        {
            scaled[i] = scaled[i] * ratio;
        }

        const std::vector<std::vector<double>> n_data = nordsieck.get_data_ref();
        double power = ratio;
        for (int i{}; i < n_data.size(); ++i) 
        {
            power = power * ratio;
            const std::vector<double> n_data_i = n_data[i];
            for (int j{}; j < n_data_i.size(); ++j) 
            {
                n_data_i[j] = n_data_i[j] * power;
            }
        }

        set_step_size(new_step_size);

    }

    /** Compute step grow/shrink factor according to normalized error.
     * @param error normalized error of the current step
     * @return grow/shrink factor for next step
     */
    protected double compute_step_grow_shrink_factor(const double error) 
    {
        return std::min(max_growth, std::max(min_reduction, safety * std::pow(error, exp)));
    }

    /** Specialized step handler storing the first step. */
    private class Nordsieck_Initializer : ODE_Step_Handler 
    {

        /** Steps counter. */
        private int count;

        /** Start of the integration. */
        private ODE_State_And_Derivative saved_start;

        /** First steps times. */
        private const std::vector<double> t;

        /** First steps states. */
        private const std::vector<std::vector<double>> y;

        /** First steps derivatives. */
        private const std::vector<std::vector<double>> y_dot;

        /** Simple constructor.
         * @param nb_start_points number of start points (including the initial point)
         */
        Nordsieck_Initializer(const int& nb_start_points) 
        {
            this.count  = 0;
            this.t      = std::vector<double>(nb_start_points];
            this.y      = std::vector<double>(nb_start_points][];
            this.y_dot   = std::vector<double>(nb_start_points][];
        }

        /** {@inherit_doc} */
        //override
        public void handle_step(ODE_StateInterpolator interpolator) 
        {

            if (count == 0) 
            {
                // first step, we need to store also the point at the beginning of the step
                saved_start   = interpolator.get_previous_state();
                t[0]    = saved_start.get_time();
                y[0]    = saved_start.get_complete_state();
                y_dot[0] = saved_start.get_complete_derivative();
            }

            // store the point at the end of the step
            ++count;
            const ODE_State_And_Derivative curr = interpolator.get_current_state();
            t[count]    = curr.get_time();
            y[count]    = curr.get_complete_state();
            y_dot[count] = curr.get_complete_derivative();

            if (count == t.size() - 1) 
            {

                // this was the last point we needed, we can compute the derivatives
                set_step_start(saved_start);
                const double raw_step = (t[t.size() - 1] - t[0]) / (t.size() - 1);
                set_step_size(get_step_size_helper().filter_step(raw_step, raw_step >= 0, true));

                // first scaled derivative
                scaled = y_dot[0].clone();
                for (int j{}; j < scaled.size(); ++j) 
                {
                    scaled[j] *= get_step_size();
                }

                // higher order derivatives
                nordsieck = initialize_high_order_derivatives(get_step_size(), t, y, y_dot);

                // stop the integrator now that all needed steps have been handled
                throw Initialization_completedMarkerException();

            }

        }

    }

    /** Marker exception used ONLY to stop the starter integrator after first step. */
    private static class Initialization_completedMarkerException
        extends Runtime_Exception 
        {

        /** Serializable version identifier. */
        private static const long serial_version_uid = -1914085471038046418L;

        /** Simple constructor. */
        Initialization_completedMarkerException() 
        {
            super((Throwable) NULL);
        }

    }

}



