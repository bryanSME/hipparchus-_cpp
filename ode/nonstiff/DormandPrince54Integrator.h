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
#include <string>
#include <vector>

/**
 * This class : the 5(4) Dormand-Prince integrator for Ordinary
 * Differential Equations.

 * <p>This integrator is an embedded Runge-Kutta integrator
 * of order 5(4) used in local extrapolation mode (i.e. the solution
 * is computed using the high order formula) with stepsize control
 * (and automatic step initialization) and continuous output. This
 * method uses 7 functions evaluations per step. However, since this
 * is an <i>fsal</i>, the last evaluation of one step is the same as
 * the first evaluation of the next step and hence can be avoided. So
 * the cost is really 6 functions evaluations per step.</p>
 *
 * <p>This method has been published (whithout the continuous output
 * that was added by Shampine in 1986) in the following article :
 * <pre>
 *  A family of embedded Runge-Kutta formulae
 *  J. R. Dormand and P. J. Prince
 *  Journal of Computational and Applied Mathematics
 *  volume 6, no 1, 1980, pp. 19-26
 * </pre></p>
 *
 */

class Dormand_Prince54_Integrator : public EmbeddedRunge_Kutta_Integrator 
{

    /** Integrator method name. */
    static const std::string METHOD_NAME = "Dormand-Prince 5(4)";

    /** Error array, element 1. */
    static const double E1 =     71.0 / 57600.0;

    // element 2 is zero, so it is neither stored nor used

    /** Error array, element 3. */
    static const double E3 =    -71.0 / 16695.0;

    /** Error array, element 4. */
    static const double E4 =     71.0 / 1920.0;

    /** Error array, element 5. */
    static const double E5 = -17253.0 / 339200.0;

    /** Error array, element 6. */
    static const double E6 =     22.0 / 525.0;

    /** Error array, element 7. */
    static const double E7 =     -1.0 / 40.0;

    /** Simple constructor.
     * Build a fifth order Dormand-Prince integrator with the given step bounds
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scal_absolute_tolerance allowed absolute error
     * @param scal_relative_tolerance allowed relative error
     */
    public Dormand_Prince54_Integrator(const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance) 
    {
        super(METHOD_NAME, 6, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
    }

    /** Simple constructor.
     * Build a fifth order Dormand-Prince integrator with the given step bounds
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vec_absolute_tolerance allowed absolute error
     * @param vec_relative_tolerance allowed relative error
     */
    public Dormand_Prince54_Integrator(const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance) 
    {
        super(METHOD_NAME, 6, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_c() 
    {
        return std::vector<double> 
        {
            1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0
        };
    }

    /** {@inherit_doc} */
    //override
    public std::vector<std::vector<double>> get_a() 
    {
        return std::vector<std::vector<double>> 
        {
            { 1.0 / 5.0 }, { 3.0 / 40.0, 9.0 / 40.0 }, { 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0 }, { 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0,  -212.0 / 729.0 }, { 9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0 }, { 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0 }
        };
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> get_b() 
    {
        return std::vector<double> 
        {
            35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0
        };
    }

    /** {@inherit_doc} */
    //override
    protected Dormand_Prince54_State_Interpolator
    create_interpolator(const bool forward, std::vector<std::vector<double>> y_dot_k, const ODE_State_And_Derivative global_previous_state, const ODE_State_And_Derivative global_current_state, const Equations_mapper mapper) 
    {
        return Dormand_Prince54_State_Interpolator(forward, y_dot_k, global_previous_state, global_current_state, global_previous_state, global_current_state, mapper);
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
            const double err_sum = E1 * y_dot_k[0][j] +  E3 * y_dot_k[2][j] +
                                  E4 * y_dot_k[3][j] +  E5 * y_dot_k[4][j] +
                                  E6 * y_dot_k[5][j] +  E7 * y_dot_k[6][j];

            const double tol = helper.get_tolerance(j, std::max(std::abs(y0[j]), std::abs(y1[j])));
            const double ratio  = h * err_sum / tol;
            error += ratio * ratio;

        }

        return std::sqrt(error / helper.get_main_set_dimension());

    }

}


