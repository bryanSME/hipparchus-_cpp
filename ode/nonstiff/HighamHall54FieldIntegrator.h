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
//import org.hipparchus.ode.FieldEquations_mapper;
//import org.hipparchus.ode.Field_ODE_State_And_Derivative;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"
#include "StepsizeHelper.h"
#include <vector>


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
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Higham_Hall54_Field_Integrator : public EmbeddedRunge_Kutta_Field_Integrator<T> 
{

    /** Simple constructor.
     * Build a fifth order Higham and Hall integrator with the given step bounds
     * @param field field to which the time and state vector elements belong
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scal_absolute_tolerance allowed absolute error
     * @param scal_relative_tolerance allowed relative error
     */
    public Higham_Hall54_Field_Integrator(const Field<T> field, const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance) 
    {
        super(field, Higham_Hall54_Integrator.METHOD_NAME, -1, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
    }

    /** Simple constructor.
     * Build a fifth order Higham and Hall integrator with the given step bounds
     * @param field field to which the time and state vector elements belong
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vec_absolute_tolerance allowed absolute error
     * @param vec_relative_tolerance allowed relative error
     */
    public Higham_Hall54_Field_Integrator(const Field<T> field, const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance) 
    {
        super(field, Higham_Hall54_Integrator.METHOD_NAME, -1, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> get_c() 
    {
        const std::vector<T> c = Math_Arrays::build_array(get_field(), 6);
        c[0] = fraction(2, 9);
        c[1] = fraction(1, 3);
        c[2] = fraction(1, 2);
        c[3] = fraction(3, 5);
        c[4] = get_field().get_one();
        c[5] = get_field().get_one();
        return c;
    }

    /** {@inherit_doc} */
    //override
    public std::vector<std::vector<T>> get_a() 
    {
        const std::vector<std::vector<T>> a = Math_Arrays::build_array(get_field(), 6, -1);
        for (int i{}; i < a.size(); ++i) 
        {
            a[i] = Math_Arrays::build_array(get_field(), i + 1);
        }
        a[0][0] = fraction(     2,     9);
        a[1][0] = fraction(     1,    12);
        a[1][1] = fraction(     1,     4);
        a[2][0] = fraction(     1,     8);
        a[2][1] = get_field().get_zero();
        a[2][2] = fraction(     3,     8);
        a[3][0] = fraction(    91,   500);
        a[3][1] = fraction(   -27,   100);
        a[3][2] = fraction(    78,   125);
        a[3][3] = fraction(     8,   125);
        a[4][0] = fraction(   -11,    20);
        a[4][1] = fraction(    27,    20);
        a[4][2] = fraction(    12,     5);
        a[4][3] = fraction(   -36,     5);
        a[4][4] = fraction(     5,     1);
        a[5][0] = fraction(     1,    12);
        a[5][1] = get_field().get_zero();
        a[5][2] = fraction(    27,    32);
        a[5][3] = fraction(    -4,     3);
        a[5][4] = fraction(   125,    96);
        a[5][5] = fraction(     5,    48);
        return a;
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> get_b() 
    {
        const std::vector<T> b = Math_Arrays::build_array(get_field(), 7);
        b[0] = fraction(  1, 12);
        b[1] = get_field().get_zero();
        b[2] = fraction( 27, 32);
        b[3] = fraction( -4,  3);
        b[4] = fraction(125, 96);
        b[5] = fraction(  5, 48);
        b[6] = get_field().get_zero();
        return b;
    }

    /** {@inherit_doc} */
    //override
    protected Higham_Hall54_Field_State_Interpolator<T>
        create_interpolator(const bool forward, std::vector<std::vector<T>> y_dot_k, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const FieldEquations_mapper<T> mapper) 
        {
        return Higham_Hall54_Field_State_Interpolator<T>(get_field(), forward, y_dot_k, global_previous_state, global_current_state, global_previous_state, global_current_state, mapper);
    }

    /** {@inherit_doc} */
    //override
    public int get_order() 
    {
        return 5;
    }

    /** {@inherit_doc} */
    //override
    protected double estimate_error(const std::vector<std::vector<T>> y_dot_k, const std::vector<T> y0, const std::vector<T> y1, const T h) 
    {

        const Stepsize_Helper helper = get_step_size_helper();
        double error = 0;

        for (int j{}; j < helper.get_main_set_dimension(); ++j) 
        {
            double err_sum = Higham_Hall54_Integrator.STATIC_E[0] * y_dot_k[0][j].get_real();
            for (const int& l = 1; l < Higham_Hall54_Integrator.STATIC_E.size(); ++l) 
            {
                err_sum += Higham_Hall54_Integrator.STATIC_E[l] * y_dot_k[l][j].get_real();
            }
            const double tol   = helper.get_tolerance(j, std::max(std::abs(y0[j].get_real()), std::abs(y1[j].get_real())));
            const double ratio = h.get_real() * err_sum / tol;
            error += ratio * ratio;

        }

        return std::sqrt(error / helper.get_main_set_dimension());

    }

};