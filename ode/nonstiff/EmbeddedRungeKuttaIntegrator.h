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
//import org.hipparchus.ode.Equations_mapper;
//import org.hipparchus.ode.Expandable_ODE;
//import org.hipparchus.ode.Localized_ODE_Formats;
//import org.hipparchus.ode.ODE_State;
//import org.hipparchus.ode.ODE_State_And_Derivative;
//import org.hipparchus.util.FastMath;
#include "StepsizeHelper.h"
#include <vector>
#include <string>

/**
 * This class : the common part of all embedded Runge-Kutta
 * integrators for Ordinary Differential Equations.
 *
 * <p>These methods are embedded explicit Runge-Kutta methods with two
 * sets of coefficients allowing to estimate the error, their Butcher
 * arrays are as follows :
 * <pre>
 *    0  |
 *   c2  | a21
 *   c3  | a31  a32
 *   ... |        ...
 *   cs  | as1  as2  ...  ass-1
 *       |--------------------------
 *       |  b1   b2  ...   bs-1  bs
 *       |  b'1  b'2 ...   b's-1 b's
 * </pre>
 * </p>
 *
 * <p>In fact, we rather use the array defined by ej = bj - b'j to
 * compute directly the error rather than computing two estimates and
 * then comparing them.</p>
 *
 * <p>Some methods are qualified as <i>fsal</i> (first same as last)
 * methods. This means the last evaluation of the derivatives in one
 * step is the same as the first in the next step. Then, this
 * evaluation can be reused from one step to the next one and the cost
 * of such a method is really s-1 evaluations despite the method still
 * has s stages. This behaviour is true only for successful steps, if
 * the step is rejected after the error estimation phase, no
 * evaluation is saved. For an <i>fsal</i> method, we have cs = 1 and
 * asi = bi for all i.</p>
 *
 */

class EmbeddedRunge_Kutta_Integrator
    extends Adaptive_Stepsize_Integrator
    : Butcher_Array_Provider 
    {

    /** Index of the pre-computed derivative for <i>fsal</i> methods. */
    private const int fsal;

    /** Time steps from Butcher array (without the first zero). */
    private const std::vector<double> c;

    /** Internal weights from Butcher array (without the first empty row). */
    private const std::vector<std::vector<double>> a;

    /** External weights for the high order method from Butcher array. */
    private const std::vector<double> b;

    /** Stepsize control exponent. */
    private const double exp;

    /** Safety factor for stepsize control. */
    private double safety;

    /** Minimal reduction factor for stepsize control. */
    private double min_reduction;

    /** Maximal growth factor for stepsize control. */
    private double max_growth;

    /** Build a Runge-Kutta integrator with the given Butcher array.
     * @param name name of the method
     * @param fsal index of the pre-computed derivative for <i>fsal</i> methods
     * or -1 if method is not <i>fsal</i>
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scal_absolute_tolerance allowed absolute error
     * @param scal_relative_tolerance allowed relative error
     */
    protected EmbeddedRunge_Kutta_Integrator(const std::string name, const int fsal, const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance) 
    {

        super(name, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);

        this.fsal = fsal;
        this.c    = get_c();
        this.a    = get_a();
        this.b    = get_b();

        exp = -1.0 / get_order();

        // set the default values of the algorithm control parameters
        set_safety(0.9);
        set_min_reduction(0.2);
        set_max_growth(10.0);

    }

    /** Build a Runge-Kutta integrator with the given Butcher array.
     * @param name name of the method
     * @param fsal index of the pre-computed derivative for <i>fsal</i> methods
     * or -1 if method is not <i>fsal</i>
     * @param min_step minimal step (must be positive even for backward
     * integration), the last step can be smaller than this
     * @param max_step maximal step (must be positive even for backward
     * integration)
     * @param vec_absolute_tolerance allowed absolute error
     * @param vec_relative_tolerance allowed relative error
     */
    protected EmbeddedRunge_Kutta_Integrator(const std::string name, const int fsal, const double   min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance) 
    {

        super(name, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);

        this.fsal = fsal;
        this.c    = get_c();
        this.a    = get_a();
        this.b    = get_b();

        exp = -1.0 / get_order();

        // set the default values of the algorithm control parameters
        set_safety(0.9);
        set_min_reduction(0.2);
        set_max_growth(10.0);

    }

    /** Create an interpolator.
     * @param forward integration direction indicator
     * @param y_dot_k slopes at the intermediate points
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param mapper equations mapper for the all equations
     * @return external weights for the high order method from Butcher array
     */
    protected virtual Runge_Kutta_State_Interpolator create_interpolator(bool forward, std::vector<std::vector<double>> y_dot_k, ODE_State_And_Derivative global_previous_state, ODE_State_And_Derivative global_current_state, Equations_mapper mapper);
    /** Get the order of the method.
     * @return order of the method
     */
    public virtual int get_order();

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

    /** {@inherit_doc} */
    //override
    public ODE_State_And_Derivative integrate(const Expandable_ODE equations, const ODE_State initial_state, const double const_time)
        , Math_Illegal_State_Exception 
        {

        sanity_checks(initial_state, const_time);
        set_step_start(init_integration(equations, initial_state, const_time));
        const bool forward = const_time > initial_state.get_time();

        // create some internal working arrays
        const int        stages  = c.size() + 1;
        const std::vector<std::vector<double>> y_dot_k   = std::vector<double>(stages][];
        const std::vector<double>   y_tmp    = std::vector<double>(equations.get_mapper().get_total_dimension()];

        // set up integration control objects
        double  h_new      = 0;
        bool first_time = true;

        // main integration loop
        set_is_last_step(false);
        do 
        {

            // iterate over step size, ensuring local normalized error is smaller than 1
            double error = 10;
            while (error >= 1.0) 
            {

                // first stage
                const std::vector<double> y = get_step_start().get_complete_state();
                y_dot_k[0] = get_step_start().get_complete_derivative();

                if (first_time) 
                {
                    const Stepsize_Helper helper = get_step_size_helper();
                    const std::vector<double> scale = std::vector<double>(helper.get_main_set_dimension()];
                    for (int i{}; i < scale.size(); ++i) 
                    {
                        scale[i] = helper.get_tolerance(i, std::abs(y[i]));
                    }
                    h_new = initialize_step(forward, get_order(), scale, get_step_start(), equations.get_mapper());
                    first_time = false;
                }

                set_step_size(h_new);
                if (forward) 
                {
                    if (get_step_start().get_time() + get_step_size() >= const_time) 
                    {
                        set_step_size(const_time - get_step_start().get_time());
                    }
                }
else 
                {
                    if (get_step_start().get_time() + get_step_size() <= const_time) 
                    {
                        set_step_size(const_time - get_step_start().get_time());
                    }
                }

                // next stages
                for (int k{ 1 }; k < stages; ++k) 
                {

                    for (int j{}; j < y.size(); ++j) 
                    {
                        double sum = a[k-1][0] * y_dot_k[0][j];
                        for (const int& l = 1; l < k; ++l) 
                        {
                            sum += a[k-1][l] * y_dot_k[l][j];
                        }
                        y_tmp[j] = y[j] + get_step_size() * sum;
                    }

                    y_dot_k[k] = compute_derivatives(get_step_start().get_time() + c[k-1] * get_step_size(), y_tmp);

                }

                // estimate the state at the end of the step
                for (int j{}; j < y.size(); ++j) 
                {
                    double sum    = b[0] * y_dot_k[0][j];
                    for (const int& l = 1; l < stages; ++l) 
                    {
                        sum    += b[l] * y_dot_k[l][j];
                    }
                    y_tmp[j] = y[j] + get_step_size() * sum;
                }

                // estimate the error at the end of the step
                error = estimate_error(y_dot_k, y, y_tmp, get_step_size());
                if (std::isnan(error)) 
                {
                    throw Math_Illegal_State_Exception(Localized_ODE_Formats.NAN_APPEARING_DURING_INTEGRATION, get_step_start().get_time() + get_step_size());
                }
                if (error >= 1.0) 
                {
                    // reject the step and attempt to reduce error by stepsize control
                    const double factor =
                                    std::min(max_growth, std::max(min_reduction, safety * std::pow(error, exp)));
                    h_new = get_step_size_helper().filter_step(get_step_size() * factor, forward, false);
                }

            }
            const double   step_end = get_step_start().get_time() + get_step_size();
            const std::vector<double> y_dot_tmp = (fsal >= 0) ? y_dot_k[fsal] : compute_derivatives(step_end, y_tmp);
            const ODE_State_And_Derivative state_tmp = equations.get_mapper().map_state_and_derivative(step_end, y_tmp, y_dot_tmp);

            // local error is small enough: accept the step, trigger events and step handlers
            set_step_start(accept_step(create_interpolator(forward, y_dot_k, get_step_start(), state_tmp, equations.get_mapper()), const_time));

            if (!is_last_step()) 
            {

                // stepsize control for next step
                const double factor =
                                std::min(max_growth, std::max(min_reduction, safety * std::pow(error, exp)));
                const double  scaled_h    = get_step_size() * factor;
                const double  next_t      = get_step_start().get_time() + scaled_h;
                const bool next_is_last = forward ? (next_t >= const_time) : (next_t <= const_time);
                h_new = get_step_size_helper().filter_step(scaled_h, forward, next_is_last);

                const double  filtered_next_t      = get_step_start().get_time() + h_new;
                const bool filtered_next_is_last = forward ? (filtered_next_t >= const_time) : (filtered_next_t <= const_time);
                if (filtered_next_is_last) 
                {
                    h_new = const_time - get_step_start().get_time();
                }

            }

        } while (!is_last_step());

        const ODE_State_And_Derivative const_state = get_step_start();
        reset_internal_state();
        return const_state;

    }

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

    /** Compute the error ratio.
     * @param y_dot_k derivatives computed during the first stages
     * @param y0 estimate of the step at the start of the step
     * @param y1 estimate of the step at the end of the step
     * @param h  current step
     * @return error ratio, greater than 1 if step should be rejected
     */
    protected virtual double estimate_error(std::vector<std::vector<double>> y_dot_k, std::vector<double> y0, std::vector<double> y1, double h);

}


