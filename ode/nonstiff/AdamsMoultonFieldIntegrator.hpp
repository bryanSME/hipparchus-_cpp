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

//import java.util.Arrays;

//import org.hipparchus.Field;
//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.linear.Array2DRowField_Matrix;
//import org.hipparchus.linear.Field_Matrix;
//import org.hipparchus.linear.Field_Matrix_Preserving_Visitor;
//import org.hipparchus.ode.FieldEquations_mapper;
//import org.hipparchus.ode.Field_ODE_State_And_Derivative;
//import org.hipparchus.ode.Localized_ODE_Formats;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"
#include "StepsizeHelper.h"

/**
 * This class : implicit Adams-Moulton integrators for Ordinary
 * Differential Equations.
 *
 * <p>Adams-Moulton methods (in fact due to Adams alone) are implicit
 * multistep ODE solvers. This implementation is a variation of the classical
 * one: it uses adaptive stepsize to implement error control, whereas
 * classical implementations are fixed step size. The value of state vector
 * at step n+1 is a simple combination of the value at step n and of the
 * derivatives at steps n+1, n, n-1 ... sin_ce y'<sub>n+1</sub> is needed to
 * compute y<sub>n+1</sub>, another method must be used to compute a first
 * estimate of y<sub>n+1</sub>, then compute y'<sub>n+1</sub>, then compute
 * a const estimate of y<sub>n+1</sub> using the following formulas. Depending
 * on the number k of previous steps one wants to use for computing the next
 * value, different formulas are available for the const estimate:</p>
 * <ul>
 *   <li>k = 1: y<sub>n+1</sub> = y<sub>n</sub> + h y'<sub>n+1</sub></li>
 *   <li>k = 2: y<sub>n+1</sub> = y<sub>n</sub> + h (y'<sub>n+1</sub>+y'<sub>n</sub>)/2</li>
 *   <li>k = 3: y<sub>n+1</sub> = y<sub>n</sub> + h (5y'<sub>n+1</sub>+8y'<sub>n</sub>-y'<sub>n-1</sub>)/12</li>
 *   <li>k = 4: y<sub>n+1</sub> = y<sub>n</sub> + h (9y'<sub>n+1</sub>+19y'<sub>n</sub>-5y'<sub>n-1</sub>+y'<sub>n-2</sub>)/24</li>
 *   <li>...</li>
 * </ul>
 *
 * <p>A k-steps Adams-Moulton method is of order k+1.</p>
 *
 * <p> There must be sufficient time for the {@link #set_starter_integrator(org.hipparchus.ode.FieldODE_Integrator)
 * starter integrator} to take several steps between the the last reset event, and the end
 * of integration, otherwise an exception may be thrown during integration. The user can
 * adjust the end date of integration, or the step size of the starter integrator to
 * ensure a sufficient number of steps can be completed before the end of integration.
 * </p>
 *
 * <h3>Implementation details</h3>
 *
 * <p>We define scaled derivatives s<sub>i</sub>(n) at step n as:
 * <pre>
 * s<sub>1</sub>(n) = h y'<sub>n</sub> for first derivative
 * s<sub>2</sub>(n) = h<sup>2</sup>/2 y''<sub>n</sub> for second derivative
 * s<sub>3</sub>(n) = h<sup>3</sup>/6 y'''<sub>n</sub> for third derivative
 * ...
 * s<sub>k</sub>(n) = h<sup>k</sup>/k! y<sup>(k)</sup><sub>n</sub> for k<sup>th</sup> derivative
 * </pre></p>
 *
 * <p>The definitions above use the classical representation with several previous first
 * derivatives. Lets define
 * <pre>
 *   q<sub>n</sub> = [ s<sub>1</sub>(n-1) s<sub>1</sub>(n-2) ... s<sub>1</sub>(n-(k-1)) ]<sup>T</sup>
 * </pre>
 * (we omit the k index in the notation for clarity). With these definitions, * Adams-Moulton methods can be written:
 * <ul>
 *   <li>k = 1: y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n+1)</li>
 *   <li>k = 2: y<sub>n+1</sub> = y<sub>n</sub> + 1/2 s<sub>1</sub>(n+1) + [ 1/2 ] q<sub>n+1</sub></li>
 *   <li>k = 3: y<sub>n+1</sub> = y<sub>n</sub> + 5/12 s<sub>1</sub>(n+1) + [ 8/12 -1/12 ] q<sub>n+1</sub></li>
 *   <li>k = 4: y<sub>n+1</sub> = y<sub>n</sub> + 9/24 s<sub>1</sub>(n+1) + [ 19/24 -5/24 1/24 ] q<sub>n+1</sub></li>
 *   <li>...</li>
 * </ul></p>
 *
 * <p>Instead of using the classical representation with first derivatives only (y<sub>n</sub>, * s<sub>1</sub>(n+1) and q<sub>n+1</sub>), our implementation uses the Nordsieck vector with
 * higher degrees scaled derivatives all taken at the same step (y<sub>n</sub>, s<sub>1</sub>(n)
 * and r<sub>n</sub>) where r<sub>n</sub> is defined as:
 * <pre>
 * r<sub>n</sub> = [ s<sub>2</sub>(n), s<sub>3</sub>(n) ... s<sub>k</sub>(n) ]<sup>T</sup>
 * </pre>
 * (here again we omit the k index in the notation for clarity)
 * </p>
 *
 * <p>Taylor series formulas show that for any index offset i, s<sub>1</sub>(n-i) can be
 * computed from s<sub>1</sub>(n), s<sub>2</sub>(n) ... s<sub>k</sub>(n), the formula being exact
 * for degree k polynomials.
 * <pre>
 * s<sub>1</sub>(n-i) = s<sub>1</sub>(n) + &sum;<sub>j&gt;0</sub> (j+1) (-i)<sup>j</sup> s<sub>j+1</sub>(n)
 * </pre>
 * The previous formula can be used with several values for i to compute the transform between
 * classical representation and Nordsieck vector. The transform between r<sub>n</sub>
 * and q<sub>n</sub> resulting from the Taylor series formulas above is:
 * <pre>
 * q<sub>n</sub> = s<sub>1</sub>(n) u + P r<sub>n</sub>
 * </pre>
 * where u is the [ 1 1 ... 1 ]<sup>T</sup> vector and P is the (k-1)&times;(k-1) matrix built
 * with the (j+1) (-i)<sup>j</sup> terms with i being the row number starting from 1 and j being
 * the column number starting from 1:
 * <pre>
 *        [  -2   3   -4    5  ... ]
 *        [  -4  12  -32   80  ... ]
 *   P =  [  -6  27 -108  405  ... ]
 *        [  -8  48 -256 1280  ... ]
 *        [          ...           ]
 * </pre></p>
 *
 * <p>Using the Nordsieck vector has several advantages:
 * <ul>
 *   <li>it greatly simplifies step interpolation as the interpolator mainly applies
 *   Taylor series formulas,</li>
 *   <li>it simplifies step changes that occur when discrete events that truncate
 *   the step are triggered,</li>
 *   <li>it allows to extend the methods in order to support adaptive stepsize.</li>
 * </ul></p>
 *
 * <p>The predicted Nordsieck vector at step n+1 is computed from the Nordsieck vector at step
 * n as follows:
 * <ul>
 *   <li>Y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n) + u<sup>T</sup> r<sub>n</sub></li>
 *   <li>S<sub>1</sub>(n+1) = h f(t<sub>n+1</sub>, Y<sub>n+1</sub>)</li>
 *   <li>R<sub>n+1</sub> = (s<sub>1</sub>(n) - S<sub>1</sub>(n+1)) P<sup>-1</sup> u + P<sup>-1</sup> A P r<sub>n</sub></li>
 * </ul>
 * where A is a rows shifting matrix (the lower left part is an identity matrix):
 * <pre>
 *        [ 0 0   ...  0 0 | 0 ]
 *        [ ---------------+---]
 *        [ 1 0   ...  0 0 | 0 ]
 *    A = [ 0 1   ...  0 0 | 0 ]
 *        [       ...      | 0 ]
 *        [ 0 0   ...  1 0 | 0 ]
 *        [ 0 0   ...  0 1 | 0 ]
 * </pre>
 * From this predicted vector, the corrected vector is computed as follows:
 * <ul>
 *   <li>y<sub>n+1</sub> = y<sub>n</sub> + S<sub>1</sub>(n+1) + [ -1 +1 -1 +1 ... &plusmn;1 ] r<sub>n+1</sub></li>
 *   <li>s<sub>1</sub>(n+1) = h f(t<sub>n+1</sub>, y<sub>n+1</sub>)</li>
 *   <li>r<sub>n+1</sub> = R<sub>n+1</sub> + (s<sub>1</sub>(n+1) - S<sub>1</sub>(n+1)) P<sup>-1</sup> u</li>
 * </ul>
 * where the upper case Y<sub>n+1</sub>, S<sub>1</sub>(n+1) and R<sub>n+1</sub> represent the
 * predicted states whereas the lower case y<sub>n+1</sub>, s<sub>n+1</sub> and r<sub>n+1</sub>
 * represent the corrected states.</p>
 *
 * <p>The P<sup>-1</sup>u vector and the P<sup>-1</sup> A P matrix do not depend on the state, * they only depend on k and therefore are precomputed once for all.</p>
 *
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class Adams_moultonFieldIntegrator extends Adams_Field_Integrator<T> 
{

    /** Integrator method name. */
    private static const std::string METHOD_NAME = "Adams-Moulton";

    /**
     * Build an Adams-Moulton integrator with the given order and error control parameters.
     * @param field field to which the time and state vector elements belong
     * @param n_steps number of steps of the method excluding the one being computed
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scal_absolute_tolerance allowed absolute error
     * @param scal_relative_tolerance allowed relative error
     * @exception  if order is 1 or less
     */
    public Adams_moultonFieldIntegrator(const Field<T> field, const int& n_steps, const double min_step, const double max_step, const double scal_absolute_tolerance, const double scal_relative_tolerance)
         
        {
        super(field, METHOD_NAME, n_steps, n_steps + 1, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
    }

    /**
     * Build an Adams-Moulton integrator with the given order and error control parameters.
     * @param field field to which the time and state vector elements belong
     * @param n_steps number of steps of the method excluding the one being computed
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vec_absolute_tolerance allowed absolute error
     * @param vec_relative_tolerance allowed relative error
     * @exception Illegal_Argument_Exception if order is 1 or less
     */
    public Adams_moultonFieldIntegrator(const Field<T> field, const int& n_steps, const double min_step, const double max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance)
        Illegal_Argument_Exception 
        {
        super(field, METHOD_NAME, n_steps, n_steps + 1, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
    }

    /** {@inherit_doc} */
    //override
    protected double error_estimation(const std::vector<T> previous_state, const T predicted_time, const std::vector<T> predicted_state, const std::vector<T> predicted_scaled, const Field_Matrix<T> predicted_nordsieck) 
    {
        const double error = predicted_nordsieck.walk_in_optimized_order(new Corrector(previous_state, predicted_scaled, predicted_state)).get_real();
        if (std::isnan(error)) 
        {
            throw Math_Illegal_State_Exception(Localized_ODE_Formats.NAN_APPEARING_DURING_INTEGRATION, predicted_time.get_real());
        }
        return error;
    }

    /** {@inherit_doc} */
    //override
    protected Adams_Field_State_Interpolator<T> constize_step(const T step_size, const std::vector<T> predicted_y, const std::vector<T> predicted_scaled, const Array2DRowField_Matrix<T> predicted_nordsieck, const bool is_forward, const Field_ODE_State_And_Derivative<T> global_previous_state, const Field_ODE_State_And_Derivative<T> global_current_state, const FieldEquations_mapper<T> equations_mapper) 
    {

        const std::vector<T> corrected_y_dot = compute_derivatives(global_current_state.get_time(), predicted_y);

        // update Nordsieck vector
        const std::vector<T> corrected_scaled = Math_Arrays::build_array(get_field(), predicted_y.size());
        for (int j{}; j < corrected_scaled.size(); ++j) 
        {
            corrected_scaled[j] = get_step_size().multiply(corrected_y_dot[j]);
        }
        update_high_order_derivatives_phase_2(predicted_scaled, corrected_scaled, predicted_nordsieck);

        const Field_ODE_State_And_Derivative<T> updated_step_end =
                        equations_mapper.map_state_and_derivative(global_current_state.get_time(), predicted_y, corrected_y_dot);
        return Adams_Field_State_Interpolator<>(get_step_size(), updated_step_end, corrected_scaled, predicted_nordsieck, is_forward, get_step_start(), updated_step_end, equations_mapper);

    }

    /** Corrector for current state in Adams-Moulton method.
     * <p>
     * This visitor : the Taylor series formula:
     * <pre>
     * Y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n+1) + [ -1 +1 -1 +1 ... &plusmn;1 ] r<sub>n+1</sub>
     * </pre>
     * </p>
     */
    private class Corrector : Field_Matrix_Preserving_Visitor<T> 
    {

        /** Previous state. */
        private const std::vector<T> previous;

        /** Current scaled first derivative. */
        private const std::vector<T> scaled;

        /** Current state before correction. */
        private const std::vector<T> before;

        /** Current state after correction. */
        private const std::vector<T> after;

        /** Simple constructor.
         * <p>
         * All arrays will be stored by reference to caller arrays.
         * </p>
         * @param previous previous state
         * @param scaled current scaled first derivative
         * @param state state to correct (will be overwritten after visit)
         */
        Corrector(const std::vector<T> previous, const std::vector<T> scaled, const std::vector<T> state) { // NOPMD - array reference storage is intentional and documented here
            this.previous = previous;
            this.scaled   = scaled;
            this.after    = state;
            this.before   = state.clone();
        }

        /** {@inherit_doc} */
        //override
        public void start(const int& rows, int columns, int start_row, int end_row, int start_column, int end_column) 
        {
            Arrays.fill(after, get_field().get_zero());
        }

        /** {@inherit_doc} */
        //override
        public void visit(const int& row, const int& column, T value) 
        {
            if ((row & 0x1) == 0) 
            {
                after[column] = after[column].subtract(value);
            }
else 
            {
                after[column] = after[column].add(value);
            }
        }

        /**
         * End visiting the Nordsieck vector.
         * <p>The correction is used to control stepsize. So its amplitude is
         * considered to be an error, which must be normalized according to
         * error control settings. If the normalized value is greater than 1, * the correction was too large and the step must be rejected.</p>
         * @return the normalized correction, if greater than 1, the step
         * must be rejected
         */
        //override
        public T end() 
        {

            const Stepsize_Helper helper = get_step_size_helper();
            T error = get_field().get_zero();
            for (int i{}; i < after.size(); ++i) 
            {
                after[i] = after[i].add(previous[i].add(scaled[i]));
                if (i < helper.get_main_set_dimension()) 
                {
                    const T tol   = helper.get_tolerance(i, Math_Utils::max(previous[i].abs(), after[i].abs()));
                    const T ratio = after[i].subtract(before[i]).divide(tol); // (corrected-predicted)/tol
                    error = error.add(ratio.multiply(ratio));
                }
            }

            return error.divide(helper.get_main_set_dimension()).sqrt();

        }
    }
};