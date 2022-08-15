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
//import org.hipparchus.ode.Abstract_Integrator;
//import org.hipparchus.ode.Equations_mapper;
//import org.hipparchus.ode.ODE_State;
//import org.hipparchus.ode.ODE_State_And_Derivative;
//import org.hipparchus.util.FastMath;
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
 * If the Ordinary Differential Equations is an {@link org.hipparchus.ode.Expandable_ODE
 * extended ODE} rather than a {@link
 * org.hipparchus.ode.Ordinary_Differential_Equation basic ODE}, then
 * <em>only</em> the {@link org.hipparchus.ode.Expandable_ODE#get_primary() primary part}
 * of the state vector is used for stepsize control, not the complete state vector.
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
 *
 */

class Adaptive_Stepsize_Integrator
    extends Abstract_Integrator 
    {

    /** Helper for step size control. */
    private Stepsize_Helper stepsize_helper;

    /** Build an integrator with the given stepsize bounds.
     * The default step handler does nothing.
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
    public Adaptive_Stepsize_Integrator(const std::string name, const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance) 
    {
        super(name);
        stepsize_helper = Stepsize_Helper(min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
        reset_internal_state();
    }

    /** Build an integrator with the given stepsize bounds.
     * The default step handler does nothing.
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
    public Adaptive_Stepsize_Integrator(const std::string name, const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance) 
    {
        super(name);
        stepsize_helper = Stepsize_Helper(min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
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
    public void set_step_size_control(const double minimal_step, const double maximal_step, const double& absolute_tolerance, const double relative_tolerance) 
    {
        stepsize_helper = Stepsize_Helper(minimal_step, maximal_step, absolute_tolerance, relative_tolerance);
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
    public void set_step_size_control(const double minimal_step, const double maximal_step, const std::vector<double> absolute_tolerance, const std::vector<double> relative_tolerance) 
    {
        stepsize_helper = Stepsize_Helper(minimal_step, maximal_step, absolute_tolerance, relative_tolerance);
    }

    /** Get the stepsize helper.
     * @return stepsize helper
     * @since 2.0
     */
    protected Stepsize_Helper get_step_size_helper() 
    {
        return stepsize_helper;
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
    public void set_initial_step_size(const double initial_step_size) 
    {
        stepsize_helper.set_initial_step_size(initial_step_size);
    }

    /** {@inherit_doc} */
    //override
    protected void sanity_checks(const ODE_State initial_state, const double t)
                     
                    {
        super.sanity_checks(initial_state, t);
        stepsize_helper.set_main_set_dimension(initial_state.get_primary_state_dimension());
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
    public double initialize_step(const bool forward, const int order, const std::vector<double> scale, const ODE_State_And_Derivative state0, const Equations_mapper mapper)
        , Math_Illegal_State_Exception 
        {

        if (stepsize_helper.get_initial_step() > 0) 
        {
            // use the user provided value
            return forward ? stepsize_helper.get_initial_step() : -stepsize_helper.get_initial_step();
        }

        // very rough first guess : h = 0.01 * ||y/scale|| / ||y'/scale||
        // this guess will be used to perform an Euler step
        const std::vector<double> y0    = state0.get_complete_state();
        const std::vector<double> y_dot_0 = state0.get_complete_derivative();
        double y_on_scale_2 = 0;
        double y_dot_on_scale_2 = 0;
        for (int j{}; j < scale.size(); ++j) 
        {
            const double ratio    = y0[j] / scale[j];
            y_on_scale_2            += ratio * ratio;
            const double ratio_dot = y_dot_0[j] / scale[j];
            y_dot_on_scale_2         += ratio_dot * ratio_dot;
        }

        double h = ((y_on_scale_2 < 1.0e-10) || (y_dot_on_scale_2 < 1.0e-10)) ?
                   1.0e-6 : (0.01 * std::sqrt(y_on_scale_2 / y_dot_on_scale_2));
        if (! forward) 
        {
            h = -h;
        }

        // perform an Euler step using the preceding rough guess
        const std::vector<double> y1 = std::vector<double>(y0.size()];
        for (int j{}; j < y0.size(); ++j) 
        {
            y1[j] = y0[j] + h * y_dot_0[j];
        }
        const std::vector<double> y_dot1 = compute_derivatives(state0.get_time() + h, y1);

        // estimate the second derivative of the solution
        double y_d_dot_on_scale = 0;
        for (int j{}; j < scale.size(); ++j) 
        {
            const double ratio_dot_dot = (y_dot1[j] - y_dot_0[j]) / scale[j];
            y_d_dot_on_scale += ratio_dot_dot * ratio_dot_dot;
        }
        y_d_dot_on_scale = std::sqrt(y_d_dot_on_scale) / h;

        // step size is computed such that
        // h^order * max (||y'/tol||, ||y''/tol||) = 0.01
        const double max_inv2 = std::max(std::sqrt(y_dot_on_scale_2), y_d_dot_on_scale);
        const double h1 = (max_inv2 < 1.0e-15) ?
                           std::max(1.0e-6, 0.001 * std::abs(h)) :
                           std::pow(0.01 / max_inv2, 1.0 / order);
        h = std::min(100.0 * std::abs(h), h1);
        h = std::max(h, 1.0e-12 * std::abs(state0.get_time()));  // avoids cancellation when computing t1 - t0
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

    /** Reset internal state to dummy values. */
    protected void reset_internal_state() 
    {
        set_step_start(null);
        set_step_size(stepsize_helper.get_dummy_stepsize());
    }

    /** Get the minimal step.
     * @return minimal step
     */
    public double get_min_step() 
    {
        return stepsize_helper.get_min_step();
    }

    /** Get the maximal step.
     * @return maximal step
     */
    public double get_max_step() 
    {
        return stepsize_helper.get_max_step();
    }

}


