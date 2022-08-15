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

//package org.hipparchus.ode.nonstiff;

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.ode.Equations_mapper;
//import org.hipparchus.ode.Expandable_ODE;
//import org.hipparchus.ode.Localized_ODE_Formats;
//import org.hipparchus.ode.Multistep_Integrator;
//import org.hipparchus.ode.ODE_State;
//import org.hipparchus.ode.ODE_State_And_Derivative;


/** Base class for {@link AdamsBashforth_integrator Adams-Bashforth} and
 * {@link Adams_moultonIntegrator Adams-Moulton} integrators.
 */
class Adams_Integrator : public Multistep_Integrator 
{

    /** Transformer. */
    private const Adams_Nordsieck_Transformer transformer;

    /**
     * Build an Adams integrator with the given order and step control parameters.
     * @param name name of the method
     * @param n_steps number of steps of the method excluding the one being computed
     * @param order order of the method
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scal_absolute_tolerance allowed absolute error
     * @param scal_relative_tolerance allowed relative error
     * @exception  if order is 1 or less
     */
    public Adams_Integrator(const std::string name, const int& n_steps, const int order, const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance)
         
        {
        super(name, n_steps, order, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
        transformer = Adams_Nordsieck_Transformer.get_instance(n_steps);
    }

    /**
     * Build an Adams integrator with the given order and step control parameters.
     * @param name name of the method
     * @param n_steps number of steps of the method excluding the one being computed
     * @param order order of the method
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vec_absolute_tolerance allowed absolute error
     * @param vec_relative_tolerance allowed relative error
     * @exception Illegal_Argument_Exception if order is 1 or less
     */
    public Adams_Integrator(const std::string name, const int& n_steps, const int order, const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance)
        Illegal_Argument_Exception 
        {
        super(name, n_steps, order, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
        transformer = Adams_Nordsieck_Transformer.get_instance(n_steps);
    }

    /** {@inherit_doc} */
    //override
    public ODE_State_And_Derivative integrate(const Expandable_ODE equations, const ODE_State initial_state, const double const_time)
        , Math_Illegal_State_Exception 
        {

        sanity_checks(initial_state, const_time);
        set_step_start(init_integration(equations, initial_state, const_time));
        const bool forward = const_time > initial_state.get_time();

        // compute the initial Nordsieck vector using the configured starter integrator
        start(equations, get_step_start(), const_time);

        // reuse the step that was chosen by the starter integrator
        ODE_State_And_Derivative step_end   =
                        Adams_State_Interpolator.taylor(equations.get_mapper(), get_step_start(), get_step_start().get_time() + get_step_size(), get_step_size(), scaled, nordsieck);

        // main integration loop
        set_is_last_step(false);
        const std::vector<double> y  = get_step_start().get_complete_state();
        do 
        {

            std::vector<double> predicted_y  = NULL;
            const std::vector<double> predicted_scaled = std::vector<double>(y.size()];
            Array_2D_Row_Real_Matrix predicted_nordsieck = NULL;
            double error = 10;
            while (error >= 1.0) 
            {

                // predict a first estimate of the state at step end
                predicted_y = step_end.get_complete_state();

                // evaluate the derivative
                const std::vector<double> y_dot = compute_derivatives(step_end.get_time(), predicted_y);

                // predict Nordsieck vector at step end
                for (int j{}; j < predicted_scaled.size(); ++j) 
                {
                    predicted_scaled[j] = get_step_size() * y_dot[j];
                }
                predicted_nordsieck = update_high_order_derivatives_phase_1(nordsieck);
                update_high_order_derivatives_phase_2(scaled, predicted_scaled, predicted_nordsieck);

                // evaluate error
                error = error_estimation(y, step_end.get_time(), predicted_y, predicted_scaled, predicted_nordsieck);
                if (std::isnan(error)) 
                {
                    throw Math_Illegal_State_Exception(Localized_ODE_Formats.NAN_APPEARING_DURING_INTEGRATION, step_end.get_time());
                }

                if (error >= 1.0) 
                {
                    // reject the step and attempt to reduce error by stepsize control
                    const double factor = compute_step_grow_shrink_factor(error);
                    rescale(get_step_size_helper().filter_step(get_step_size() * factor, forward, false));
                    step_end = Adams_State_Interpolator.taylor(equations.get_mapper(), get_step_start(), get_step_start().get_time() + get_step_size(), get_step_size(), scaled, nordsieck);

                }
            }

            const Adams_State_Interpolator interpolator =
                            constize_step(get_step_size(), predicted_y, predicted_scaled, predicted_nordsieck, forward, get_step_start(), step_end, equations.get_mapper());

            // discrete events handling
            set_step_start(accept_step(interpolator, const_time));
            scaled    = interpolator.get_scaled();
            nordsieck = interpolator.get_nordsieck();

            if (!is_last_step()) 
            {

                if (reset_occurred()) 
                {

                    // some events handler has triggered changes that
                    // invalidate the derivatives, we need to restart from scratch
                    start(equations, get_step_start(), const_time);

                    const double  next_t      = get_step_start().get_time() + get_step_size();
                    const bool next_is_last = forward ?
                                               (next_t >= const_time) :
                                               (next_t <= const_time);
                    const double h_new = next_is_last ? const_time - get_step_start().get_time() : get_step_size();

                    rescale(h_new);
                    System.arraycopy(get_step_start().get_complete_state(), 0, y, 0, y.size());

                }
else 
                {

                    // stepsize control for next step
                    const double  factor     = compute_step_grow_shrink_factor(error);
                    const double  scaled_h    = get_step_size() * factor;
                    const double  next_t      = get_step_start().get_time() + scaled_h;
                    const bool next_is_last = forward ?
                                               (next_t >= const_time) :
                                               (next_t <= const_time);
                    double h_new = get_step_size_helper().filter_step(scaled_h, forward, next_is_last);

                    const double  filtered_next_t      = get_step_start().get_time() + h_new;
                    const bool filtered_next_is_last = forward ? (filtered_next_t >= const_time) : (filtered_next_t <= const_time);
                    if (filtered_next_is_last) 
                    {
                        h_new = const_time - get_step_start().get_time();
                    }

                    rescale(h_new);
                    System.arraycopy(predicted_y, 0, y, 0, y.size());

                }

                step_end = Adams_State_Interpolator.taylor(equations.get_mapper(), get_step_start(), get_step_start().get_time() + get_step_size(), get_step_size(), scaled, nordsieck);

            }

        } while (!is_last_step());

        const ODE_State_And_Derivative const_state = get_step_start();
        set_step_start(null);
        set_step_size(Double.NaN);
        return const_state;

    }

    /** {@inherit_doc} */
    //override
    protected Array_2D_Row_Real_Matrix initialize_high_order_derivatives(const double h, const std::vector<double> t, const std::vector<std::vector<double>> y, const std::vector<std::vector<double>> y_dot) 
    {
        return transformer.initialize_high_order_derivatives(h, t, y, y_dot);
    }

    /** Update the high order scaled derivatives for Adams integrators (phase 1).
     * <p>The complete update of high order derivatives has a form similar to:
     * <pre>
     * r<sub>n+1</sub> = (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u + P<sup>-1</sup> A P r<sub>n</sub>
     * </pre>
     * this method computes the P<sup>-1</sup> A P r<sub>n</sub> part.</p>
     * @param high_order high order scaled derivatives
     * (h<sup>2</sup>/2 y'', ... h<sup>k</sup>/k! y(k))
     * @return updated high order derivatives
     * @see #update_high_order_derivatives_phase_2(std::vector<double>, std::vector<double>, Array_2D_Row_Real_Matrix)
     */
    public Array_2D_Row_Real_Matrix update_high_order_derivatives_phase_1(const Array_2D_Row_Real_Matrix high_order) 
    {
        return transformer.update_high_order_derivatives_phase_1(high_order);
    }

    /** Update the high order scaled derivatives Adams integrators (phase 2).
     * <p>The complete update of high order derivatives has a form similar to:
     * <pre>
     * r<sub>n+1</sub> = (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u + P<sup>-1</sup> A P r<sub>n</sub>
     * </pre>
     * this method computes the (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u part.</p>
     * <p>Phase 1 of the update must already have been performed.</p>
     * @param start first order scaled derivatives at step start
     * @param end first order scaled derivatives at step end
     * @param high_order high order scaled derivatives, will be modified
     * (h<sup>2</sup>/2 y'', ... h<sup>k</sup>/k! y(k))
     * @see #update_high_order_derivatives_phase_1(Array_2D_Row_Real_Matrix)
     */
    public void update_high_order_derivatives_phase_2(const std::vector<double> start, const std::vector<double> end, const Array_2D_Row_Real_Matrix high_order) 
    {
        transformer.update_high_order_derivatives_phase_2(start, end, high_order);
    }

    /** Estimate error.
     * @param previous_state state vector at step start
     * @param predicted_time time at step end
     * @param predicted_state predicted state vector at step end
     * @param predicted_scaled predicted value of the scaled derivatives at step end
     * @param predicted_nordsieck predicted value of the Nordsieck vector at step end
     * @return estimated normalized local discretization error
     * @since 2.0
     */
    protected virtual double error_estimation(std::vector<double> previous_state, double predicted_time, std::vector<double> predicted_state, std::vector<double> predicted_scaled, Real_Matrix predicted_nordsieck);

    /** Finalize the step.
     * @param step_size step size used in the scaled and Nordsieck arrays
     * @param predicted_state predicted state at end of step
     * @param predicted_scaled predicted first scaled derivative
     * @param predicted_nordsieck predicted Nordsieck vector
     * @param is_forward integration direction indicator
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param equations_mapper mapper for ODE equations primary and secondary components
     * @return step interpolator
     * @since 2.0
     */
    protected virtual Adams_State_Interpolator constize_step(double step_size, std::vector<double> predicted_state, std::vector<double> predicted_scaled, Array_2D_Row_Real_Matrix predicted_nordsieck, bool is_forward, ODE_State_And_Derivative global_previous_state, ODE_State_And_Derivative global_current_state, Equations_mapper equations_mapper);

}


