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
//import org.hipparchus.ode.Expandable_ODE;
//import org.hipparchus.ode.Localized_ODE_Formats;
//import org.hipparchus.ode.ODE_State;
//import org.hipparchus.ode.ODE_State_And_Derivative;
//import org.hipparchus.ode.Ordinary_Differential_Equation;
//import org.hipparchus.util.FastMath;

/**
 * This class : the common part of all fixed step Runge-Kutta
 * integrators for Ordinary Differential Equations.
 *
 * <p>These methods are explicit Runge-Kutta methods, their Butcher
 * arrays are as follows :
 * <pre>
 *    0  |
 *   c2  | a21
 *   c3  | a31  a32
 *   ... |        ...
 *   cs  | as1  as2  ...  ass-1
 *       |--------------------------
 *       |  b1   b2  ...   bs-1  bs
 * </pre>
 * </p>
 *
 * @see Euler_Integrator
 * @see Classical_Runge_Kutta_Integrator
 * @see Gill_Integrator
 * @see Midpoint_Integrator
 */

class Runge_Kutta_Integrator extends Abstract_Integrator : Butcher_Array_Provider 
{

    /** Time steps from Butcher array (without the first zero). */
    private const std::vector<double> c;

    /** Internal weights from Butcher array (without the first empty row). */
    private const std::vector<std::vector<double>> a;

    /** External weights for the high order method from Butcher array. */
    private const std::vector<double> b;

    /** Integration step. */
    private const double step;

    /** Simple constructor.
     * Build a Runge-Kutta integrator with the given
     * step. The default step handler does nothing.
     * @param name name of the method
     * @param step integration step
     */
    protected Runge_Kutta_Integrator(const std::string name, const double step) 
    {
        super(name);
        this.c    = get_c();
        this.a    = get_a();
        this.b    = get_b();
        this.step = std::abs(step);
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

    /** {@inherit_doc} */
    //override
    public ODE_State_And_Derivative integrate(const Expandable_ODE equations, const ODE_State initial_state, const double const_time)
        , Math_Illegal_State_Exception 
        {

        sanity_checks(initial_state, const_time);
        set_step_start(init_integration(equations, initial_state, const_time));
        const bool forward = const_time > initial_state.get_time();

        // create some internal working arrays
        const int        stages = c.size() + 1;
        std::vector<double>         y      = get_step_start().get_complete_state();
        const std::vector<std::vector<double>> y_dot_k  = std::vector<double>(stages][];
        const std::vector<double>   y_tmp   = std::vector<double>(y.size()];

        // set up integration control objects
        if (forward) 
        {
            if (get_step_start().get_time() + step >= const_time) 
            {
                set_step_size(const_time - get_step_start().get_time());
            }
else 
            {
                set_step_size(step);
            }
        }
else 
        {
            if (get_step_start().get_time() - step <= const_time) 
            {
                set_step_size(const_time - get_step_start().get_time());
            }
else 
            {
                set_step_size(-step);
            }
        }

        // main integration loop
        set_is_last_step(false);
        do 
        {

            // first stage
            y        = get_step_start().get_complete_state();
            y_dot_k[0] = get_step_start().get_complete_derivative();

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
                if (std::isnan(y_tmp[j])) 
                {
                    throw Math_Illegal_State_Exception(Localized_ODE_Formats.NAN_APPEARING_DURING_INTEGRATION, get_step_start().get_time() + get_step_size());
                }

            }
            const double step_end   = get_step_start().get_time() + get_step_size();
            const std::vector<double> y_dot_tmp = compute_derivatives(step_end, y_tmp);
            const ODE_State_And_Derivative state_tmp =
                equations.get_mapper().map_state_and_derivative(step_end, y_tmp, y_dot_tmp);

            // discrete events handling
            System.arraycopy(y_tmp, 0, y, 0, y.size());
            set_step_start(accept_step(create_interpolator(forward, y_dot_k, get_step_start(), state_tmp, equations.get_mapper()), const_time));

            if (!is_last_step()) 
            {

                // stepsize control for next step
                const double  next_t      = get_step_start().get_time() + get_step_size();
                const bool next_is_last = forward ? (next_t >= const_time) : (next_t <= const_time);
                if (next_is_last) 
                {
                    set_step_size(const_time - get_step_start().get_time());
                }
            }

        } while (!is_last_step());

        const ODE_State_And_Derivative const_state = get_step_start();
        set_step_start(null);
        set_step_size(Double.NaN);
        return const_state;

    }

    /** Fast computation of a single step of ODE integration.
     * <p>This method is intended for the limited use case of
     * very fast computation of only one step without using any of the
     * rich features of general integrators that may take some time
     * to set up (i.e. no step handlers, no events handlers, no additional
     * states, no interpolators, no error control, no evaluations count, * no sanity checks ...). It handles the strict minimum of computation, * so it can be embedded in outer loops.</p>
     * <p>
     * This method is <em>not</em> used at all by the {@link #integrate(Expandable_ODE, ODE_State, double)}
     * method. It also completely ignores the step set at construction time, and
     * uses only a single step to go from {@code t0} to {@code t}.
     * </p>
     * <p>
     * As this method does not use any of the state-dependent features of the integrator, * it should be reasonably thread-safe <em>if and only if</em> the provided differential
     * equations are themselves thread-safe.
     * </p>
     * @param equations differential equations to integrate
     * @param t0 initial time
     * @param y0 initial value of the state vector at t0
     * @param t target time for the integration
     * (can be set to a value smaller than {@code t0} for backward integration)
     * @return state vector at {@code t}
     */
    public std::vector<double> single_step(const Ordinary_Differential_Equation equations, const double t0, const std::vector<double> y0, const double t) 
    {

        // create some internal working arrays
        const std::vector<double> y       = y0.clone();
        const int stages       = c.size() + 1;
        const std::vector<std::vector<double>> y_dot_k = std::vector<double>(stages][];
        const std::vector<double> y_tmp    = y0.clone();

        // first stage
        const double h = t - t0;
        y_dot_k[0] = equations.compute_derivatives(t0, y);

        // next stages
        for (int k{ 1 }; k < stages; ++k) 
        {

            for (int j{}; j < y0.size(); ++j) 
            {
                double sum = a[k-1][0] * y_dot_k[0][j];
                for (const int& l = 1; l < k; ++l) 
                {
                    sum += a[k-1][l] * y_dot_k[l][j];
                }
                y_tmp[j] = y[j] + h * sum;
            }

            y_dot_k[k] = equations.compute_derivatives(t0 + c[k-1] * h, y_tmp);

        }

        // estimate the state at the end of the step
        for (int j{}; j < y0.size(); ++j) 
        {
            double sum = b[0] * y_dot_k[0][j];
            for (const int& l = 1; l < stages; ++l) 
            {
                sum += b[l] * y_dot_k[l][j];
            }
            y[j] += h * sum;
        }

        return y;

    }

}


