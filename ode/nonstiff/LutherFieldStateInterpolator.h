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
#include <type_traits>
#include "../../core/CalculusFieldElement.h"

/**
 * This class represents an interpolator over the last step during an
 * ODE integration for the 6th order Luther integrator.
 *
 * <p>This interpolator computes dense output inside the last
 * step computed. The interpolation equation is consistent with the
 * integration scheme.</p>
 *
 * @see Luther_fieldIntegrator
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Luther_fieldStateInterpolator : public Runge_Kutta_Field_State_Interpolator<T> 
{

    /** -49 - 49 q. */
    private const T c5a;

    /** 392 + 287 q. */
    private const T c5b;

    /** -637 - 357 q. */
    private const T c5c;

    /** 833 + 343 q. */
    private const T c5d;

    /** -49 + 49 q. */
    private const T c6a;

    /** -392 - 287 q. */
    private const T c6b;

    /** -637 + 357 q. */
    private const T c6c;

    /** 833 - 343 q. */
    private const T c6d;

    /** 49 + 49 q. */
    private const T d5a;

    /** -1372 - 847 q. */
    private const T d5b;

    /** 2254 + 1029 q */
    private const T d5c;

    /** 49 - 49 q. */
    private const T d6a;

    /** -1372 + 847 q. */
    private const T d6b;

    /** 2254 - 1029 q */
    private const T d6c;

    /** Simple constructor.
     * @param field field to which the time and state vector elements belong
     * @param forward integration direction indicator
     * @param y_dot_k slopes at the intermediate points
     * @param global_previous_state start of the global step
     * @param global_current_state end of the global step
     * @param soft_previous_state start of the restricted step
     * @param soft_current_state end of the restricted step
     * @param mapper equations mapper for the all equations
     */
    Luther_fieldStateInterpolator(const Field<T> field, const bool forward, const std::vector<std::vector<T>> y_dot_k, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const Field_ODE_State_And_Derivative<T> soft_previous_state, const Field_ODE_State_And_Derivative<T> soft_current_state, const FieldEquations_mapper<T> mapper) 
    {
        super(field, forward, y_dot_k, global_previous_state, global_current_state, soft_previous_state, soft_current_state, mapper);
        const T q = field.get_zero().add(21).sqrt();
        c5a = q.multiply(  -49).add(  -49);
        c5b = q.multiply(  287).add(  392);
        c5c = q.multiply( -357).add( -637);
        c5d = q.multiply(  343).add(  833);
        c6a = q.multiply(   49).add(  -49);
        c6b = q.multiply( -287).add(  392);
        c6c = q.multiply(  357).add( -637);
        c6d = q.multiply( -343).add(  833);
        d5a = q.multiply(   49).add(   49);
        d5b = q.multiply( -847).add(-1372);
        d5c = q.multiply( 1029).add( 2254);
        d6a = q.multiply(  -49).add(   49);
        d6b = q.multiply(  847).add(-1372);
        d6c = q.multiply(-1029).add( 2254);
    }

    /** {@inherit_doc} */
    //override
    protected Luther_fieldStateInterpolator<T> create(const Field<T> new_field, const bool new_forward, const std::vector<std::vector<T>> new_y_dot_k, const Field_ODE_State_And_Derivative<T> new_global_previous_state, const Field_ODE_State_And_Derivative<T> new_global_current_state, const Field_ODE_State_And_Derivative<T> new_soft_previous_state, const Field_ODE_State_And_Derivative<T> new_soft_current_state, const FieldEquations_mapper<T> new_mapper) 
    {
        return Luther_fieldStateInterpolator<T>(new_field, new_forward, new_y_dot_k, new_global_previous_state, new_global_current_state, new_soft_previous_state, new_soft_current_state, new_mapper);
    }

    /** {@inherit_doc} */
    //@Suppress_Warnings("unchecked")
    //override
    protected Field_ODE_State_And_Derivative<T> compute_interpolated_state_and_derivatives(const FieldEquations_mapper<T> mapper, const T time, const T theta, const T theta_h, const T one_minus_theta_h) 
    {

        // the coefficients below have been computed by solving the
        // order conditions from a theorem from Butcher (1963), using
        // the method explained in Folkmar Bornemann paper "Runge-Kutta
        // Methods, Trees, and Maple", Center of Mathematical Sciences, Munich
        // University of Technology, February 9, 2001
        //<http://wwwzenger.informatik.tu-muenchen.de/selcuk/sjam012101.html>

        // the method is implemented in the rkcheck tool
        // <https://www.spaceroots.org/software/rkcheck/index.html>.
        // Running it for order 5 gives the following order conditions
        // for an interpolator:
        // order 1 conditions
        // \\sum_{i=1}^{i=s}\\left(b_{i} 
ight) =1
        // order 2 conditions
        // \\sum_{i=1}^{i=s}\\left(b_{i} c_{i}
ight) = \x0crac{	heta}{2}
        // order 3 conditions
        // \\sum_{i=2}^{i=s}\\left(b_{i} \\sum_{j=1}^{j=i-1}{\\left(a_{i,j} c_{j} 
ight)}
ight) = \x0crac{	heta^{2}}{6}
        // \\sum_{i=1}^{i=s}\\left(b_{i} c_{i}^{2}
ight) = \x0crac{	heta^{2}}{3}
        // order 4 conditions
        // \\sum_{i=3}^{i=s}\\left(b_{i} \\sum_{j=2}^{j=i-1}{\\left(a_{i,j} \\sum_{k=1}^{k=j-1}{\\left(a_{j,k} c_{k} 
ight)} 
ight)}
ight) = \x0crac{	heta^{3}}{24}
        // \\sum_{i=2}^{i=s}\\left(b_{i} \\sum_{j=1}^{j=i-1}{\\left(a_{i,j} c_{j}^{2} 
ight)}
ight) = \x0crac{	heta^{3}}{12}
        // \\sum_{i=2}^{i=s}\\left(b_{i} c_{i}\\sum_{j=1}^{j=i-1}{\\left(a_{i,j} c_{j} 
ight)}
ight) = \x0crac{	heta^{3}}{8}
        // \\sum_{i=1}^{i=s}\\left(b_{i} c_{i}^{3}
ight) = \x0crac{	heta^{3}}{4}
        // order 5 conditions
        // \\sum_{i=4}^{i=s}\\left(b_{i} \\sum_{j=3}^{j=i-1}{\\left(a_{i,j} \\sum_{k=2}^{k=j-1}{\\left(a_{j,k} \\sum_{l=1}^{l=k-1}{\\left(a_{k,l} c_{l} 
ight)} 
ight)} 
ight)}
ight) = \x0crac{	heta^{4}}{120}
        // \\sum_{i=3}^{i=s}\\left(b_{i} \\sum_{j=2}^{j=i-1}{\\left(a_{i,j} \\sum_{k=1}^{k=j-1}{\\left(a_{j,k} c_{k}^{2} 
ight)} 
ight)}
ight) = \x0crac{	heta^{4}}{60}
        // \\sum_{i=3}^{i=s}\\left(b_{i} \\sum_{j=2}^{j=i-1}{\\left(a_{i,j} c_{j}\\sum_{k=1}^{k=j-1}{\\left(a_{j,k} c_{k} 
ight)} 
ight)}
ight) = \x0crac{	heta^{4}}{40}
        // \\sum_{i=2}^{i=s}\\left(b_{i} \\sum_{j=1}^{j=i-1}{\\left(a_{i,j} c_{j}^{3} 
ight)}
ight) = \x0crac{	heta^{4}}{20}
        // \\sum_{i=3}^{i=s}\\left(b_{i} c_{i}\\sum_{j=2}^{j=i-1}{\\left(a_{i,j} \\sum_{k=1}^{k=j-1}{\\left(a_{j,k} c_{k} 
ight)} 
ight)}
ight) = \x0crac{	heta^{4}}{30}
        // \\sum_{i=2}^{i=s}\\left(b_{i} c_{i}\\sum_{j=1}^{j=i-1}{\\left(a_{i,j} c_{j}^{2} 
ight)}
ight) = \x0crac{	heta^{4}}{15}
        // \\sum_{i=2}^{i=s}\\left(b_{i} \\left(\\sum_{j=1}^{j=i-1}{\\left(a_{i,j} c_{j} 
ight)} 
ight)^{2}
ight) = \x0crac{	heta^{4}}{20}
        // \\sum_{i=2}^{i=s}\\left(b_{i} c_{i}^{2}\\sum_{j=1}^{j=i-1}{\\left(a_{i,j} c_{j} 
ight)}
ight) = \x0crac{	heta^{4}}{10}
        // \\sum_{i=1}^{i=s}\\left(b_{i} c_{i}^{4}
ight) = \x0crac{	heta^{4}}{5}

        // The a_{j,k} and c_{k} are given by the integrator Butcher arrays. What remains to solve
        // are the b_i for the interpolator. They are found by solving the above equations.
        // For a given interpolator, some equations are redundant, so in our case when we select
        // all equations from order 1 to 4, we still don't have enough independent equations
        // to solve from b_1 to b_7. We need to also select one equation from order 5. Here, // we selected the last equation. It appears this choice implied at least the last 3 equations
        // are fulfilled, but some of the former ones are not, so the resulting interpolator is order 5.
        // At the end, we get the b_i as polynomials in theta.

        const T coeff_dot_1 =  theta.multiply(theta.multiply(theta.multiply(theta.multiply(   21        ).add( -47          )).add(   36         )).add( -54     /   5.0)).add(1);
        const T coeff_dot_2 =  time.get_field().get_zero();
        const T coeff_dot_3 =  theta.multiply(theta.multiply(theta.multiply(theta.multiply(  112        ).add(-608    /  3.0)).add(  320   / 3.0 )).add(-208    /  15.0));
        const T coeff_dot_4 =  theta.multiply(theta.multiply(theta.multiply(theta.multiply( -567  /  5.0).add( 972    /  5.0)).add( -486   / 5.0 )).add( 324    /  25.0));
        const T coeff_dot5 =  theta.multiply(theta.multiply(theta.multiply(theta.multiply(c5a.divide(5)).add(c5b.divide(15))).add(c5c.divide(30))).add(c5d.divide(150)));
        const T coeff_dot6 =  theta.multiply(theta.multiply(theta.multiply(theta.multiply(c6a.divide(5)).add(c6b.divide(15))).add(c6c.divide(30))).add(c6d.divide(150)));
        const T coeff_dot7 =  theta.multiply(theta.multiply(theta.multiply(                                             3.0 ).add(   -3         )).add(   3   /   5.0));
        const std::vector<T> interpolated_state;
        const std::vector<T> interpolated_derivatives;

        if (get_global_previous_state() != NULL && theta.get_real() <= 0.5) 
        {

            const T s         = theta_h;
            const T coeff1    = s.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply(  21    /  5.0).add( -47    /  4.0)).add(   12         )).add( -27    /   5.0)).add(1));
            const T coeff2    = time.get_field().get_zero();
            const T coeff3    = s.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply( 112    /  5.0).add(-152    /  3.0)).add(  320   / 9.0 )).add(-104    /  15.0)));
            const T coeff4    = s.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply(-567    / 25.0).add( 243    /  5.0)).add( -162   / 5.0 )).add( 162    /  25.0)));
            const T coeff5    = s.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply(c5a.divide(25)).add(c5b.divide(60))).add(c5c.divide(90))).add(c5d.divide(300))));
            const T coeff6    = s.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply(c6a.divide(25)).add(c6b.divide(60))).add(c6c.divide(90))).add(c6d.divide(300))));
            const T coeff7    = s.multiply(theta.multiply(theta.multiply(theta.multiply(                                      3    /  4.0 ).add(   -1         )).add(   3    /  10.0)));
            interpolated_state       = previous_state_linear_combination(coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7);
            interpolated_derivatives = derivative_linear_combination(coeff_dot_1, coeff_dot_2, coeff_dot_3, coeff_dot_4, coeff_dot5, coeff_dot6, coeff_dot7);
        }
else 
        {

            const T s         = one_minus_theta_h;
            const T coeff1    = s.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply( -21   /   5.0).add(   151  /  20.0)).add( -89   /   20.0)).add(  19 /  20.0)).add(- 1 / 20.0));
            const T coeff2    = time.get_field().get_zero();
            const T coeff3    = s.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply(-112   /   5.0).add(   424  /  15.0)).add( -328  /   45.0)).add( -16 /  45.0)).add(-16 /  45.0));
            const T coeff4    = s.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply( 567   /  25.0).add(  -648  /  25.0)).add(  162  /   25.0))));
            const T coeff5    = s.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply(d5a.divide(25)).add(d5b.divide(300))).add(d5c.divide(900))).add( -49 / 180.0)).add(-49 / 180.0));
            const T coeff6    = s.multiply(theta.multiply(theta.multiply(theta.multiply(theta.multiply(d6a.divide(25)).add(d6b.divide(300))).add(d6c.divide(900))).add( -49 / 180.0)).add(-49 / 180.0));
            const T coeff7    = s.multiply(               theta.multiply(theta.multiply(theta.multiply(                        -3  /   4.0 ).add(   1   /    4.0)).add(  -1 /  20.0)).add( -1 /  20.0));
            interpolated_state       = current_state_linear_combination(coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7);
            interpolated_derivatives = derivative_linear_combination(coeff_dot_1, coeff_dot_2, coeff_dot_3, coeff_dot_4, coeff_dot5, coeff_dot6, coeff_dot7);
        }

        return mapper.map_state_and_derivative(time, interpolated_state, interpolated_derivatives);

    }

};