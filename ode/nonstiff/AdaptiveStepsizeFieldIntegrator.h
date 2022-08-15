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

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.ode.Abstract_Field_Integrator;
//import org.hipparchus.ode.FieldEquations_mapper;
//import org.hipparchus.ode.FieldODE_State;
//import org.hipparchus.ode.Field_ODE_State_And_Derivative;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
#include <vector>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include "../../core/CalculusFieldElement.h"
#include "StepsizeHelper.h"

/**
 * This virtual class holds the common part of all adaptive
 * stepsize integrators for Ordinary Differential Equations.
 *
 * <p>These algorithms perform integration with stepsize control, which
 * means the user does not specify the integration step but rather a
 * tolerance on error. The error threshold is computed as
 * <pre>
 * threshold_i = abs_tol_i + rel_tol_i * max (abs (ym), abs (ym+1))
 * </pre>
 * where abs_tol_i is the absolute tolerance for component i of the
 * state vector and rel_tol_i is the relative tolerance for the same
 * component. The user can also use only two scalar values abs_tol and
 * rel_tol which will be used for all components.
 * </p>
 * <p>
 * Note that <em>only</em> the {@link FieldODE_State#get_primary_state() main part}
 * of the state vector is used for stepsize control. The {@link
 * FieldODE_State#get_secondary_statestatic_cast<int>( secondary parts} of the state
 * vector are explicitly ignored for stepsize control.
 * </p>
 *
 * <p>If the estimated error for ym+1 is such that
 * <pre>
 * sqrt((sum (err_est_i / threshold_i)^2 ) / n) &lt; 1
 * </pre>
 *
 * (where n is the main set dimension) then the step is accepted, * otherwise the step is rejected and a attempt is made with a new
 * stepsize.</p>
 *
 * @param <T> the type of the field elements
 *
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Adaptive_Stepsize_Field_Integrator : public Abstract_Field_Integrator<T> 
{
private:
    /** Helper for step size control. */
    Stepsize_Helper my_stepsize_helper;

protected:
    /** Get the stepsize helper.
     * @return stepsize helper
     * @since 2.0
     */
    Stepsize_Helper get_step_size_helper()
    {
        return my_stepsize_helper;
    }

    /** {@inherit_doc} */
    //override
    void sanity_checks(const FieldODE_State<T> initial_state, const T t)

    {
        super.sanity_checks(initial_state, t);
        my_stepsize_helper.set_main_set_dimension(initial_state.get_primary_state_dimension());
    }

    /** Reset internal state to dummy values. */
    void reset_internal_state()
    {
        set_step_start(null);
        set_step_size(get_field().get_zero().add(stepsize_helper.get_dummy_stepsize()));
    }

public:

    /** Build an integrator with the given stepsize bounds.
     * The default step handler does nothing.
     * @param field field to which the time and state vector elements belong
     * @param name name of the method
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scal_absolute_tolerance allowed absolute error
     * @param scal_relative_tolerance allowed relative error
     */
    Adaptive_Stepsize_Field_Integrator(const Field<T> field, const std::string name, const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance) 
    {
        super(field, name);
        my_stepsize_helper = Stepsize_Helper(min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
        reset_internal_state();
    }

    /** Build an integrator with the given stepsize bounds.
     * The default step handler does nothing.
     * @param field field to which the time and state vector elements belong
     * @param name name of the method
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vec_absolute_tolerance allowed absolute error
     * @param vec_relative_tolerance allowed relative error
     */
    Adaptive_Stepsize_Field_Integrator(const Field<T> field, const std::string name, const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance) 
    {
        super(field, name);
        my_stepsize_helper = Stepsize_Helper(min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
        reset_internal_state();
    }

    /** Set the adaptive step size control parameters.
     * <p>
     * A side effect of this method is to also reset the initial
     * step so it will be automatically computed by the integrator
     * if {@link #set_initial_step_sizestatic_cast<double>( set_initial_step_size}
     * is not called by the user.
     * </p>
     * @param minimal_step minimal step (must be positive even for backward
     * integration), the last step can be smaller than this
     * @param maximal_step maximal step (must be positive even for backward
     * integration)
     * @param absolute_tolerance allowed absolute error
     * @param relative_tolerance allowed relative error
     */
    void set_step_size_control(const double minimal_step, const double maximal_step, const double& absolute_tolerance, const double relative_tolerance) 
    {
        my_stepsize_helper = Stepsize_Helper(minimal_step, maximal_step, absolute_tolerance, relative_tolerance);
    }

    /** Set the adaptive step size control parameters.
     * <p>
     * A side effect of this method is to also reset the initial
     * step so it will be automatically computed by the integrator
     * if {@link #set_initial_step_sizestatic_cast<double>( set_initial_step_size}
     * is not called by the user.
     * </p>
     * @param minimal_step minimal step (must be positive even for backward
     * integration), the last step can be smaller than this
     * @param maximal_step maximal step (must be positive even for backward
     * integration)
     * @param absolute_tolerance allowed absolute error
     * @param relative_tolerance allowed relative error
     */
    void set_step_size_control(const double minimal_step, const double maximal_step, const std::vector<double> absolute_tolerance, const std::vector<double> relative_tolerance) 
    {
        my_stepsize_helper = Stepsize_Helper(minimal_step, maximal_step, absolute_tolerance, relative_tolerance);
    }

    /** Set the initial step size.
     * <p>This method allows the user to specify an initial positive
     * step size instead of letting the integrator guess it by
     * itself. If this method is not called before integration is
     * started, the initial step size will be estimated by the
     * integrator.</p>
     * @param initial_step_size initial step size to use (must be positive even
     * for backward integration ; providing a negative value or a value
     * outside of the min/max step interval will lead the integrator to
     * ignore the value and compute the initial step size by itself)
     */
    void set_initial_step_size(const double initial_step_size) 
    {
        my_stepsize_helper.set_initial_step_size(initial_step_size);
    }

    /** Initialize the integration step.
     * @param forward forward integration indicator
     * @param order order of the method
     * @param scale scaling vector for the state vector (can be shorter than state vector)
     * @param state0 state at integration start time
     * @param mapper mapper for all the equations
     * @return first integration step
     * @exception Math_Illegal_State_Exception if the number of functions evaluations is exceeded
     * @exception  if arrays dimensions do not match equations settings
     */
    double initialize_step(const bool forward, const int order, const std::vector<T> scale, const Field_ODE_State_And_Derivative<T> state0, const FieldEquations_mapper<T> mapper)
    {
        if (my_stepsize_helper.get_initial_step() > 0) 
        {
            // use the user provided value
            return forward
                ? my_stepsize_helper.get_initial_step() 
                : -stepsize_helper.get_initial_step();
        }

        // very rough first guess : h = 0.01 * ||y/scale|| / ||y'/scale||
        // this guess will be used to perform an Euler step
        const std::vector<T> y0    = state0.get_complete_state();
        const std::vector<T> y_dot_0 = state0.get_complete_derivative();
        double y_on_scale_2     = 0;
        double y_dot_on_scale_2  = 0;
        for (int j{}; j < scale.size(); ++j) 
        {
            const double ratio    = y0[j].get_real() / scale[j].get_real();
            y_on_scale_2            += ratio * ratio;
            const double ratio_dot = y_dot_0[j].get_real() / scale[j].get_real();
            y_dot_on_scale_2         += ratio_dot * ratio_dot;
        }

        double h = ((y_on_scale_2 < 1.0e-10) || (y_dot_on_scale_2 < 1.0e-10)) ?
                   1.0e-6 : (0.01 * std::sqrt(y_on_scale_2 / y_dot_on_scale_2));
        if (! forward) 
        {
            h = -h;
        }

        // perform an Euler step using the preceding rough guess
        const std::vector<T> y1 = Math_Arrays::build_array(get_field(), y0.size());
        for (int j{}; j < y0.size(); ++j) 
        {
            y1[j] = y0[j].add(y_dot_0[j].multiply(h));
        }
        const std::vector<T> y_dot1 = compute_derivatives(state0.get_time().add(h), y1);

        // estimate the second derivative of the solution
        double y_d_dot_on_scale{};
        for (int j{}; j < scale.size(); ++j) 
        {
            const double ratio_dot_dot = (y_dot1[j].get_real() - y_dot_0[j].get_real()) / scale[j].get_real();
            y_d_dot_on_scale += ratio_dot_dot * ratio_dot_dot;
        }
        y_d_dot_on_scale = std::sqrt(y_d_dot_on_scale) / h;

        // step size is computed such that
        // h^order * max (||y'/tol||, ||y''/tol||) = 0.01
        const double max_inv2 = std::max(std::sqrt(y_dot_on_scale_2), y_d_dot_on_scale);
        const double h1 = (max_inv2 < 1.0e-15)
            ? std::max(1.0e-6, 0.001 * std::abs(h))
            : std::pow(0.01 / max_inv2, 1.0 / order);

        h = std::min(100.0 * std::abs(h), h1);
        h = std::max(h, 1.0e-12 * std::abs(state0.get_time().get_real()));  // avoids cancellation when computing t1 - t0
        if (h < get_min_step()) 
        {
            h = get_min_step();
        }
        if (h > get_max_step()) 
        {
            h = get_max_step();
        }

        if (! forward) 
        {
            h = -h;
        }

        return h;
    }

    /** Get the minimal step.
     * @return minimal step
     */
    double get_min_step() 
    {
        returnmy_stepsize_helper.get_min_step();
    }

    /** Get the maximal step.
     * @return maximal step
     */
    double get_max_step() 
    {
        returnmy_stepsize_helper.get_max_step();
    }

};