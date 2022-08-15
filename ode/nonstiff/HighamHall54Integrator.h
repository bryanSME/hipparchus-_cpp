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

//import org.hipparchus.ode.Equations_mapper;
//import org.hipparchus.ode.ODE_State_And_Derivative;
//import org.hipparchus.util.FastMath;
#include "StepsizeHelper.h"
#include <vector>
#include <string>

/**
 * This class : the 5(4) Higham and Hall integrator for
 * Ordinary Differential Equations.
 *
 * <p>This integrator is an embedded Runge-Kutta integrator
 * of order 5(4) used in local extrapolation mode (i.e. the solution
 * is computed using the high order formula) with stepsize control
 * (and automatic step initialization) and continuous output. This
 * method uses 7 functions evaluations per step.</p>
 *
 */

class Higham_Hall54_Integrator extends EmbeddedRunge_Kutta_Integrator 
{

    /** Integrator method name. */
    static const std::string METHOD_NAME = "Higham-Hall 5(4)";

    /** Error weights Butcher array. */
    static const std::vector<double> STATIC_E = 
    {
        -1.0/20.0, 0.0, 81.0/160.0, -6.0/5.0, 25.0/32.0, 1.0/16.0, -1.0/10.0
    };

    /** Simple constructor.
     * Build a fifth order Higham and Hall integrator with the given step bounds
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scal_absolute_tolerance allowed absolute error
     * @param scal_relative_tolerance allowed relative error
     */
    public Higham_Hall54_Integrator(const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance) 
    {
        super(METHOD_NAME, -1, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
    }

    /** Simple constructor.
     * Build a fifth order Higham and Hall integrator with the given step bounds
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vec_absolute_tolerance allowed absolute error
     * @param vec_relative_tolerance allowed relative error
     */
    public Higham_Hall54_Integrator(const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance) 
    {
        super(METHOD_NAME, -1, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_c() 
    {
        return std::vector<double> 
        {
            2.0/9.0, 1.0/3.0, 1.0/2.0, 3.0/5.0, 1.0, 1.0
        };
    }

    /** {@inherit_doc} */
    //override
    public std::vector<std::vector<double>> get_a() 
    {
        return std::vector<std::vector<double>> 
        {
            {2.0/9.0}, {1.0/12.0, 1.0/4.0}, {1.0/8.0, 0.0, 3.0/8.0}, {91.0/500.0, -27.0/100.0, 78.0/125.0, 8.0/125.0}, {-11.0/20.0, 27.0/20.0, 12.0/5.0, -36.0/5.0, 5.0}, {1.0/12.0, 0.0, 27.0/32.0, -4.0/3.0, 125.0/96.0, 5.0/48.0}
        };
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_b() 
    {
        return std::vector<double> 
        {
            1.0/12.0, 0.0, 27.0/32.0, -4.0/3.0, 125.0/96.0, 5.0/48.0, 0.0
        };
    }

    /** {@inherit_doc} */
    //override
    protected Higham_Hall54_State_Interpolator
    create_interpolator(const bool forward, std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const Equations_mapper mapper) 
    {
        return Higham_Hall54_State_Interpolator(forward, y_dot_k, global_previous_state, global_current_state, global_previous_state, global_current_state, mapper);
    }

    /** {@inherit_doc} */
    //override
    public int get_order() 
    {
        return 5;
    }

    /** {@inherit_doc} */
    //override
    protected double estimate_error(const std::vector<std::vector<double>> y_dot_k, const std::vector<double> y0, const std::vector<double> y1, const double h) 
    {

        const Stepsize_Helper helper = get_step_size_helper();
        double error = 0;

        for (int j{}; j < helper.get_main_set_dimension(); ++j) 
        {
            double err_sum = STATIC_E[0] * y_dot_k[0][j];
            for (const int& l = 1; l < STATIC_E.size(); ++l) 
            {
                err_sum += STATIC_E[l] * y_dot_k[l][j];
            }

            const double tol = helper.get_tolerance(j, std::max(std::abs(y0[j]), std::abs(y1[j])));
            const double ratio  = h * err_sum / tol;
            error += ratio * ratio;

        }

        return std::sqrt(error / helper.get_main_set_dimension());

    }

}


