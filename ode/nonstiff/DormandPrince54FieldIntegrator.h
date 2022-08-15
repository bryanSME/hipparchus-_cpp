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
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Dormand_Prince54_Field_Integrator : public EmbeddedRunge_Kutta_Field_Integrator<T> 
{

    /** Simple constructor.
     * Build a fifth order Dormand-Prince integrator with the given step bounds
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
    public Dormand_Prince54_Field_Integrator(const Field<T> field, const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance) 
    {
        super(field, Dormand_Prince54_Integrator.METHOD_NAME, 6, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
    }

    /** Simple constructor.
     * Build a fifth order Dormand-Prince integrator with the given step bounds
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
    public Dormand_Prince54_Field_Integrator(const Field<T> field, const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance) 
    {
        super(field, Dormand_Prince54_Integrator.METHOD_NAME, 6, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> get_c() 
    {
        const std::vector<T> c = Math_Arrays::build_array(get_field(), 6);
        c[0] = fraction(1,  5);
        c[1] = fraction(3, 10);
        c[2] = fraction(4,  5);
        c[3] = fraction(8,  9);
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
        a[0][0] = fraction(     1,     5);
        a[1][0] = fraction(     3,    40);
        a[1][1] = fraction(     9,    40);
        a[2][0] = fraction(    44,    45);
        a[2][1] = fraction(   -56,    15);
        a[2][2] = fraction(    32,     9);
        a[3][0] = fraction( 19372,  6561);
        a[3][1] = fraction(-25360,  2187);
        a[3][2] = fraction( 64448,  6561);
        a[3][3] = fraction(  -212,   729);
        a[4][0] = fraction(  9017,  3168);
        a[4][1] = fraction(  -355,    33);
        a[4][2] = fraction( 46732,  5247);
        a[4][3] = fraction(    49,   176);
        a[4][4] = fraction( -5103, 18656);
        a[5][0] = fraction(    35,   384);
        a[5][1] = get_field().get_zero();
        a[5][2] = fraction(   500,  1113);
        a[5][3] = fraction(   125,   192);
        a[5][4] = fraction( -2187,  6784);
        a[5][5] = fraction(    11,    84);
        return a;
    }

    /** {@inherit_doc} */
    //override
    public std::vector<T> get_b() 
    {
        const std::vector<T> b = Math_Arrays::build_array(get_field(), 7);
        b[0] = fraction(   35,   384);
        b[1] = get_field().get_zero();
        b[2] = fraction(  500, 1113);
        b[3] = fraction(  125,  192);
        b[4] = fraction(-2187, 6784);
        b[5] = fraction(   11,   84);
        b[6] = get_field().get_zero();
        return b;
    }

    /** {@inherit_doc} */
    //override
    protected Dormand_Prince54_Field_State_Interpolator<T>
        create_interpolator(const bool forward, std::vector<std::vector<T>> y_dot_k, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const FieldEquations_mapper<T> mapper) 
        {
        return Dormand_Prince54_Field_State_Interpolator<T>(get_field(), forward, y_dot_k, global_previous_state, global_current_state, global_previous_state, global_current_state, mapper);
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
            const double err_sum = Dormand_Prince54_Integrator.E1 * y_dot_k[0][j].get_real() +  Dormand_Prince54_Integrator.E3 * y_dot_k[2][j].get_real() +
                                  Dormand_Prince54_Integrator.E4 * y_dot_k[3][j].get_real() +  Dormand_Prince54_Integrator.E5 * y_dot_k[4][j].get_real() +
                                  Dormand_Prince54_Integrator.E6 * y_dot_k[5][j].get_real() +  Dormand_Prince54_Integrator.E7 * y_dot_k[6][j].get_real();
            const double tol = helper.get_tolerance(j, std::max(std::abs(y0[j].get_real()), std::abs(y1[j].get_real())));
            const double ratio  = h.get_real() * err_sum / tol;
            error += ratio * ratio;
        }

        return std::sqrt(error / helper.get_main_set_dimension());

    }

};