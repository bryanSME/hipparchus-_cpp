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

//import java.util.Arrays;

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.linear.Array_2D_Row_Real_Matrix;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.linear.Real_Matrix_Preserving_Visitor;
//import org.hipparchus.ode.Equations_mapper;
//import org.hipparchus.ode.Localized_ODE_Formats;
//import org.hipparchus.ode.ODE_State_And_Derivative;
//import org.hipparchus.util.FastMath;
#include "AdamsIntegrator.h"
#include <vector>
#include <string>
#include <cmath>
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
 * <p> There must be sufficient time for the {@link #set_starter_integrator(org.hipparchus.ode.ODE_Integrator)
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
 */
class Adams_moultonIntegrator : public Adams_Integrator
{
private:
    /** Integrator method name. */
    static const std::string METHOD_NAME{ "Adams-Moulton" };

    /**
     * Build an Adams-Moulton integrator with the given order and error control parameters.
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
    public Adams_moultonIntegrator(const int& n_steps, const double& min_step, const double& max_step, const double& scal_absolute_tolerance, const double& scal_relative_tolerance)
    {
        Adams_Integrator(METHOD_NAME, n_steps, n_steps + 1, min_step, max_step, scal_absolute_tolerance, scal_relative_tolerance);
    }

    /**
     * Build an Adams-Moulton integrator with the given order and error control parameters.
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
    public Adams_moultonIntegrator(const int& n_steps, const double& min_step, const double& max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance)
    {
        Adams_Integrator(METHOD_NAME, n_steps, n_steps + 1, min_step, max_step, vec_absolute_tolerance, vec_relative_tolerance);
    }

    /** {@inherit_doc} */
    //override
    protected double error_estimation(const std::vector<double>& previous_state, const double& predicted_time, const std::vector<double>& predicted_state, const std::vector<double>& predicted_scaled, const Real_Matrix& predicted_nordsieck)
    {
        const double error = predicted_nordsieck.walk_in_optimized_order(new Corrector(previous_state, predicted_scaled, predicted_state));
        if (std::isnan(error))
        {
            throw std::exception("not fully implemented -- Adams_moultonIntegrator");
            //throw Math_Illegal_State_Exception(Localized_ODE_Formats.NAN_APPEARING_DURING_INTEGRATION, predicted_time);
        }
        return error;
    }

    /** {@inherit_doc} */
    //override
    protected Adams_State_Interpolator constize_step(const double& step_size, const std::vector<double>& predicted_state, const std::vector<double>& predicted_scaled, const Array_2D_Row_Real_Matrix& predicted_nordsieck, const bool is_forward, const ODE_State_And_Derivative& global_previous_state, const ODE_State_And_Derivative& global_current_state, const Equations_mapper& equations_mapper)
    {
        const std::vector<double> corrected_y_dot = compute_derivatives(global_current_state.get_time(), predicted_state);

        // update Nordsieck vector
        const std::vector<double> corrected_scaled = std::vector<double>(predicted_state.size()];
        for (int j{}; j < corrected_scaled.size(); ++j)
        {
            corrected_scaled[j] = get_step_size() * corrected_y_dot[j];
        }
        update_high_order_derivatives_phase_2(predicted_scaled, corrected_scaled, predicted_nordsieck);

        const ODE_State_And_Derivative updated_step_end =
            equations_mapper.map_state_and_derivative(global_current_state.get_time(), predicted_state, corrected_y_dot);
        return Adams_State_Interpolator(get_step_size(), updated_step_end, corrected_scaled, predicted_nordsieck, is_forward, get_step_start(), updated_step_end, equations_mapper);

    }

    /** Corrector for current state in Adams-Moulton method.
     * <p>
     * This visitor : the Taylor series formula:
     * <pre>
     * Y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n+1) + [ -1 +1 -1 +1 ... &plusmn;1 ] r<sub>n+1</sub>
     * </pre>
     * </p>
     */
    class Corrector : public Real_Matrix_Preserving_Visitor
    {
    private:
        /** Previous state. */
        const std::vector<double> my_previous;

        /** Current scaled first derivative. */
        const std::vector<double> my_scaled;

        /** Current state before correction. */
        const std::vector<double> my_before;

        /** Current state after correction. */
        std::vector<double> my_after;

    public:
        /** Simple constructor.
         * <p>
         * All arrays will be stored by reference to caller arrays.
         * </p>
         * @param previous previous state
         * @param scaled current scaled first derivative
         * @param state state to correct (will be overwritten after visit)
         */
        Corrector(const std::vector<double>& previous, const std::vector<double>& scaled, const std::vector<double>& state)
            :
            my_previous{previous },
            my_scaled{ scaled},
            my_after{state},
            my_before{state}
        { // NOPMD - array reference storage is intentional and documented here
        }

        /** {@inherit_doc} */
        //override
        void start(const int& rows, [[maybe_unused]]const int& columns, [[maybe_unused]]const int& start_row, [[maybe_unused]]const int& end_row, [[maybe_unused]]const int& start_column, [[maybe_unused]]const int& end_column)
        {
            my_after = std::vector<double>(my_after.size(), 0.0);
        }

        /** {@inherit_doc} */
        //override
        void visit(const int& row, const int& column, const double& value)
        {
            if ((row & 0x1) == 0)
            {
                my_after[column] -= value;
            }
            else
            {
                my_after[column] += value;
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
        double end()
        {
            const Stepsize_Helper helper = get_step_size_helper();
            double error = 0;
            for (int i{}; i < my_after.size(); ++i)
            {
                my_after[i] += my_previous[i] + my_scaled[i];
                if (i < helper.get_main_set_dimension())
                {
                    const double tol = helper.get_tolerance(i, std::max(std::abs(my_previous[i]), std::abs(my_after[i])));
                    const double ratio = (my_after[i] - my_before[i]) / tol; // (corrected-predicted)/tol
                    error += ratio * ratio;
                }
            }

            return std::sqrt(error / helper.get_main_set_dimension());

        }
    }
};