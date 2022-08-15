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
#include <type_traits>
#include "../../core/CalculusFieldElement.h"

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.ode.Abstract_Field_Integrator;
//import org.hipparchus.ode.FieldEquations_mapper;
//import org.hipparchus.ode.FieldExpandable_ODE;
//import org.hipparchus.ode.FieldODE_State;
//import org.hipparchus.ode.Field_ODE_State_And_Derivative;
//import org.hipparchus.ode.FieldOrdinary_Differential_Equation;
//import org.hipparchus.util.Math_Arrays;

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
 * @see Euler_fieldIntegrator
 * @see ClassicalRunge_Kutta_Field_Integrator
 * @see Gill_Field_Integrator
 * @see Midpoint_Field_Integrator
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Runge_Kutta_Field_Integrator : public Abstract_Field_Integrator<T>, FieldButcher_Array_Provider<T> 
{

    /** Time steps from Butcher array (without the first zero). */
    private const std::vector<T> c;

    /** Internal weights from Butcher array (without the first empty row). */
    private const std::vector<std::vector<T>> a;

    /** External weights for the high order method from Butcher array. */
    private const std::vector<T> b;

    /** Integration step. */
    private const T step;

    /** Simple constructor.
     * Build a Runge-Kutta integrator with the given
     * step. The default step handler does nothing.
     * @param field field to which the time and state vector elements belong
     * @param name name of the method
     * @param step integration step
     */
    protected Runge_Kutta_Field_Integrator(const Field<T> field, const std::string name, const T step) 
    {
        super(field, name);
        this.c    = get_c();
        this.a    = get_a();
        this.b    = get_b();
        this.step = step.abs();
    }

    /** Create a fraction.
     * @param p numerator
     * @param q denominator
     * @return p/q computed in the instance field
     */
    protected T fraction(const int p, const int q) 
    {
        return get_field().get_zero().add(p).divide(q);
    }

    /** Create an interpolator.
     * @param forward integration direction indicator
     * @param y_dot_k slopes at the intermediate points
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param mapper equations mapper for the all equations
     * @return external weights for the high order method from Butcher array
     */
    protected virtual Runge_Kutta_Field_State_Interpolator<T> create_interpolator(bool forward, std::vector<std::vector<T>> y_dot_k, Field_ODE_State_And_Derivative<T> global_previous_state, Field_ODE_State_And_Derivative<T> global_current_state, FieldEquations_mapper<T> mapper);

    /** {@inherit_doc} */
    //override
    public Field_ODE_State_And_Derivative<T> integrate(const FieldExpandable_ODE<T> equations, const FieldODE_State<T> initial_state, const T const_time)
        , Math_Illegal_State_Exception 
        {

        sanity_checks(initial_state, const_time);
        set_step_start(init_integration(equations, initial_state, const_time));
        const bool forward = const_time.subtract(initial_state.get_time()).get_real() > 0;

        // create some internal working arrays
        const int   stages = c.size() + 1;
        const std::vector<std::vector<T>> y_dot_k  = Math_Arrays::build_array(get_field(), stages, -1);
        const std::vector<T>   y_tmp   = Math_Arrays::build_array(get_field(), equations.get_mapper().get_total_dimension());

        // set up integration control objects
        if (forward) 
        {
            if (get_step_start().get_time().add(step).subtract(const_time).get_real() >= 0) 
            {
                set_step_size(const_time.subtract(get_step_start().get_time()));
            }
else 
            {
                set_step_size(step);
            }
        }
else 
        {
            if (get_step_start().get_time().subtract(step).subtract(const_time).get_real() <= 0) 
            {
                set_step_size(const_time.subtract(get_step_start().get_time()));
            }
else 
            {
                set_step_size(step.negate());
            }
        }

        // main integration loop
        set_is_last_step(false);
        do 
        {

            // first stage
            const std::vector<T> y = get_step_start().get_complete_state();
            y_dot_k[0]    = get_step_start().get_complete_derivative();

            // next stages
            for (int k{ 1 }; k < stages; ++k) 
            {

                for (int j{}; j < y.size(); ++j) 
                {
                    T sum = y_dot_k[0][j].multiply(a[k-1][0]);
                    for (const int& l = 1; l < k; ++l) 
                    {
                        sum = sum.add(y_dot_k[l][j].multiply(a[k-1][l]));
                    }
                    y_tmp[j] = y[j].add(get_step_size().multiply(sum));
                }

                y_dot_k[k] = compute_derivatives(get_step_start().get_time().add(get_step_size().multiply(c[k-1])), y_tmp);

            }

            // estimate the state at the end of the step
            for (int j{}; j < y.size(); ++j) 
            {
                T sum = y_dot_k[0][j].multiply(b[0]);
                for (const int& l = 1; l < stages; ++l) 
                {
                    sum = sum.add(y_dot_k[l][j].multiply(b[l]));
                }
                y_tmp[j] = y[j].add(get_step_size().multiply(sum));
            }
            const T step_end   = get_step_start().get_time().add(get_step_size());
            const std::vector<T> y_dot_tmp = compute_derivatives(step_end, y_tmp);
            const Field_ODE_State_And_Derivative<T> state_tmp = equations.get_mapper().map_state_and_derivative(step_end, y_tmp, y_dot_tmp);

            // discrete events handling
            set_step_start(accept_step(create_interpolator(forward, y_dot_k, get_step_start(), state_tmp, equations.get_mapper()), const_time));

            if (!is_last_step()) 
            {

                // stepsize control for next step
                const T  next_t      = get_step_start().get_time().add(get_step_size());
                const bool next_is_last = forward ?
                                           (next_t.subtract(const_time).get_real() >= 0) :
                                           (next_t.subtract(const_time).get_real() <= 0);
                if (next_is_last) 
                {
                    set_step_size(const_time.subtract(get_step_start().get_time()));
                }
            }

        } while (!is_last_step());

        const Field_ODE_State_And_Derivative<T> const_state = get_step_start();
        set_step_start(null);
        set_step_size(null);
        return const_state;

    }

    /** Fast computation of a single step of ODE integration.
     * <p>This method is intended for the limited use case of
     * very fast computation of only one step without using any of the
     * rich features of general integrators that may take some time
     * to set up (i.e. no step handlers, no events handlers, no additional
     * states, no interpolators, no error control, no evaluations count, * no sanity checks ...). It handles the strict minimum of computation, * so it can be embedded in outer loops.</p>
     * <p>
     * This method is <em>not</em> used at all by the {@link #integrate(FieldExpandable_ODE, * FieldODE_State, Calculus_Field_Element)} method. It also completely ignores the step set at
     * construction time, and uses only a single step to go from {@code t0} to {@code t}.
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
    public std::vector<T> single_step(const FieldOrdinary_Differential_Equation<T> equations, const T t0, const std::vector<T> y0, const T t) 
    {

        // create some internal working arrays
        const std::vector<T> y       = y0.clone();
        const int stages  = c.size() + 1;
        const std::vector<std::vector<T>> y_dot_k = Math_Arrays::build_array(get_field(), stages, -1);
        const std::vector<T> y_tmp    = y0.clone();

        // first stage
        const T h = t.subtract(t0);
        y_dot_k[0] = equations.compute_derivatives(t0, y);

        // next stages
        for (int k{ 1 }; k < stages; ++k) 
        {

            for (int j{}; j < y0.size(); ++j) 
            {
                T sum = y_dot_k[0][j].multiply(a[k-1][0]);
                for (const int& l = 1; l < k; ++l) 
                {
                    sum = sum.add(y_dot_k[l][j].multiply(a[k-1][l]));
                }
                y_tmp[j] = y[j].add(h.multiply(sum));
            }

            y_dot_k[k] = equations.compute_derivatives(t0.add(h.multiply(c[k-1])), y_tmp);

        }

        // estimate the state at the end of the step
        for (int j{}; j < y0.size(); ++j) 
        {
            T sum = y_dot_k[0][j].multiply(b[0]);
            for (const int& l = 1; l < stages; ++l) 
            {
                sum = sum.add(y_dot_k[l][j].multiply(b[l]));
            }
            y[j] = y[j].add(h.multiply(sum));
        }

        return y;

    }

};