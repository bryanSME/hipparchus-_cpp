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

//package org.hipparchus.ode.nonstiff;

//import org.hipparchus.Field;
//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.linear.Array2DRowField_Matrix;
//import org.hipparchus.linear.Field_Matrix;
//import org.hipparchus.ode.FieldEquations_mapper;
//import org.hipparchus.ode.FieldExpandable_ODE;
//import org.hipparchus.ode.FieldODE_State;
//import org.hipparchus.ode.Field_ODE_State_And_Derivative;
//import org.hipparchus.ode.Localized_ODE_Formats;
//import org.hipparchus.ode.Multistep_Field_Integrator;
//import org.hipparchus.util.Math_Arrays;
#include <type_traits>
#include <vector>
#include <cmath>
#include <string>
#include "../../core/CalculusFieldElement.h"
#include "../../ode/MultistepFieldIntegrator.hpp"

/** Base class for {@link Adams_Bashforth_Field_Integrator Adams-Bashforth} and
 * {@link Adams_moultonFieldIntegrator Adams-Moulton} integrators.
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Adams_Field_Integrator : public Multistep_Field_Integrator<T> 
{
private:
    /** Transformer. */
    const Adams_Nordsieck_Field_Transformer<T> transformer;

public:
    /**
     * Build an Adams integrator with the given order and step control parameters.
     * @param field field to which the time and state vector elements belong
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
    Adams_Field_Integrator(const Field<T> field, const std::string name, const int& n_steps, const int order, const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance)
         
        {
        super(field, name, n_steps, order, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
        transformer = Adams_Nordsieck_Field_Transformer.get_instance(field, n_steps);
    }

    /**
     * Build an Adams integrator with the given order and step control parameters.
     * @param field field to which the time and state vector elements belong
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
    public Adams_Field_Integrator(const Field<T> field, const std::string name, const int& n_steps, const int order, const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance)
        Illegal_Argument_Exception 
        {
        super(field, name, n_steps, order, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
        transformer = Adams_Nordsieck_Field_Transformer.get_instance(field, n_steps);
    }

    /** {@inherit_doc} */
    //override
    public Field_ODE_State_And_Derivative<T> integrate(const FieldExpandable_ODE<T> equations, const FieldODE_State<T> initial_state, const T const_time)
        , Math_Illegal_State_Exception 
        {

        sanity_checks(initial_state, const_time);
        set_step_start(init_integration(equations, initial_state, const_time));
        const bool forward = const_time.subtract(initial_state.get_time()).get_real() > 0;

        // compute the initial Nordsieck vector using the configured starter integrator
        start(equations, get_step_start(), const_time);

        // reuse the step that was chosen by the starter integrator
        Field_ODE_State_And_Derivative<T> step_start = get_step_start();
        Field_ODE_State_And_Derivative<T> step_end   =
                        Adams_Field_State_Interpolator.taylor(equations.get_mapper(), step_start, step_start.get_time().add(get_step_size()), get_step_size(), scaled, nordsieck);

        // main integration loop
        set_is_last_step(false);
        const std::vector<T> y = step_start.get_complete_state();
        do 
        {

            std::vector<T> predicted_y = NULL;
            const std::vector<T> predicted_scaled = Math_Arrays::build_array(get_field(), y.size());
            Array2DRowField_Matrix<T> predicted_nordsieck = NULL;
            double error = 10;
            while (error >= 1.0) 
            {

                // predict a first estimate of the state at step end
                predicted_y = step_end.get_complete_state();

                // evaluate the derivative
                const std::vector<T> y_dot = compute_derivatives(step_end.get_time(), predicted_y);

                // predict Nordsieck vector at step end
                for (int j{}; j < predicted_scaled.size(); ++j) 
                {
                    predicted_scaled[j] = get_step_size().multiply(y_dot[j]);
                }
                predicted_nordsieck = update_high_order_derivatives_phase_1(nordsieck);
                update_high_order_derivatives_phase_2(scaled, predicted_scaled, predicted_nordsieck);

                // evaluate error
                error = error_estimation(y, step_end.get_time(), predicted_y, predicted_scaled, predicted_nordsieck);
                if (std::isnan(error)) 
                {
                    throw Math_Illegal_State_Exception(Localized_ODE_Formats.NAN_APPEARING_DURING_INTEGRATION, step_end.get_time().get_real());
                }

                if (error >= 1.0) 
                {
                    // reject the step and attempt to reduce error by stepsize control
                    const double factor = compute_step_grow_shrink_factor(error);
                    rescale(get_step_size_helper().filter_step(get_step_size().multiply(factor), forward, false));
                    step_end = Adams_Field_State_Interpolator.taylor(equations.get_mapper(), get_step_start(), get_step_start().get_time().add(get_step_size()), get_step_size(), scaled, nordsieck);

                }
            }

            const Adams_Field_State_Interpolator<T> interpolator =
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

                    const T  next_t      = get_step_start().get_time().add(get_step_size());
                    const bool next_is_last = forward ?
                                               next_t.subtract(const_time).get_real() >= 0 :
                                               next_t.subtract(const_time).get_real() <= 0;
                    const T h_new = next_is_last ? const_time.subtract(get_step_start().get_time()) : get_step_size();

                    rescale(h_new);
                    System.arraycopy(get_step_start().get_complete_state(), 0, y, 0, y.size());

                }
                else 
                {

                    // stepsize control for next step
                    const double  factor     = compute_step_grow_shrink_factor(error);
                    const T       scaled_h    = get_step_size().multiply(factor);
                    const T       next_t      = get_step_start().get_time().add(scaled_h);
                    const bool next_is_last = forward ?
                                               next_t.subtract(const_time).get_real() >= 0 :
                                               next_t.subtract(const_time).get_real() <= 0;
                    T h_new = get_step_size_helper().filter_step(scaled_h, forward, next_is_last);

                    const T       filtered_next_t      = get_step_start().get_time().add(h_new);
                    const bool filtered_next_is_last = forward ?
                                                       filtered_next_t.subtract(const_time).get_real() >= 0 :
                                                       filtered_next_t.subtract(const_time).get_real() <= 0;
                    if (filtered_next_is_last) 
                    {
                        h_new = const_time.subtract(get_step_start().get_time());
                    }

                    rescale(h_new);
                    System.arraycopy(predicted_y, 0, y, 0, y.size());

                }

                step_end = Adams_Field_State_Interpolator.taylor(equations.get_mapper(), get_step_start(), get_step_start().get_time().add(get_step_size()), get_step_size(), scaled, nordsieck);

            }

        } while (!is_last_step());

        const Field_ODE_State_And_Derivative<T> const_state = get_step_start();
        set_step_start(null);
        set_step_size(null);
        return const_state;

    }

    /** {@inherit_doc} */
    //override
    protected Array2DRowField_Matrix<T> initialize_high_order_derivatives(const T h, const std::vector<T> t, const std::vector<std::vector<T>> y, const std::vector<std::vector<T>> y_dot) 
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
     * @see #update_high_order_derivatives_phase_2(Calculus_Field_Element[], Calculus_Field_Element[], Array2DRowField_Matrix)
     */
    public Array2DRowField_Matrix<T> update_high_order_derivatives_phase_1(const Array2DRowField_Matrix<T> high_order) 
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
     * @see #update_high_order_derivatives_phase_1(Array2DRowField_Matrix)
     */
    public void update_high_order_derivatives_phase_2(const std::vector<T> start, const std::vector<T> end, const Array2DRowField_Matrix<T> high_order) 
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
    protected virtual double error_estimation(std::vector<T> previous_state, T predicted_time, std::vector<T> predicted_state, std::vector<T> predicted_scaled, Field_Matrix<T> predicted_nordsieck);

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
    protected virtual Adams_Field_State_Interpolator<T> constize_step(T step_size, std::vector<T> predicted_state, std::vector<T> predicted_scaled, Array2DRowField_Matrix<T> predicted_nordsieck, bool is_forward, Field_ODE_State_And_Derivative<T> global_previous_state, Field_ODE_State_And_Derivative<T> global_current_state, FieldEquations_mapper<T> equations_mapper);

};