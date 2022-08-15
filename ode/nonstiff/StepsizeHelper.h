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

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.ode.Localized_ODE_Formats;
//import org.hipparchus.util.FastMath;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"
#include  <cmath>
#include <vector>


/** Helper for adaptive stepsize control.
 * @since 2.0
 */

class Stepsize_Helper 
{
private:
    /** Allowed absolute scalar error. */
    double my_scal_absolute_tolerance;

    /** Allowed relative scalar error. */
    double my_scal_relative_tolerance;

    /** Allowed absolute vectorial error. */
    std::vector<double> my_vec_absolute_tolerance;

    /** Allowed relative vectorial error. */
    std::vector<double> my_vec_relative_tolerance;

    /** Main set dimension. */
    int my_main_set_dimension;

    /** User supplied initial step. */
    double my_initial_step;

    /** Minimal step. */
    double my_min_step;

    /** Maximal step. */
    double my_max_step;

protected:
    /** Set main set dimension.
     * @param main_set_dimension dimension of the main set
     * @exception  if adaptive step size integrators
     * tolerance arrays dimensions are not compatible with equations settings
     */
    void set_main_set_dimension(const int& main_set_dimension)
    {
        my_main_set_dimension = main_set_dimension;
        throw std::exception("error not implemented");
        //if (my_vec_absolute_tolerance != NULL && my_vec_absolute_tolerance.size() != main_set_dimension)
        //{
        //    throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, main_set_dimension, my_vec_absolute_tolerance.size());
        //}

        //if (my_vec_absolute_tolerance != NULL && my_vec_absolute_tolerance.size() != main_set_dimension)
        //{
        //    throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, main_set_dimension, my_vec_absolute_tolerance.size());
        //}
    }

public:
    /** Simple constructor.
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scal_absolute_tolerance allowed absolute error
     * @param scal_relative_tolerance allowed relative error
     */
    Stepsize_Helper(const double& min_step, const double& max_step, const double& scal_absolute_tolerance, const double& scal_relative_tolerance) 
        :
        my_min_step{ std::abs(min_step) },
        my_max_step{ std::abs(max_step) },
        my_initial_step{ -1 },
        my_scal_absolute_tolerance{ scal_absolute_tolerance },
        my_scal_relative_tolerance{ scal_relative_tolerance },
        my_vec_absolute_tolerance{ NULL },
        my_vec_relative_tolerance{ NULL }
    {
    }

    /** Simple constructor..
     * @param min_step minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param max_step maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vec_absolute_tolerance allowed absolute error
     * @param vec_relative_tolerance allowed relative error
     */
    Stepsize_Helper(const double& min_step, const double& max_step, const std::vector<double>& vec_absolute_tolerance, const std::vector<double>& vec_relative_tolerance) 
        :
        my_min_step{ std::abs(min_step) },
        my_max_step{ std::abs(max_step) },
        my_initial_step{ -1 },
        my_scal_absolute_tolerance{},
        my_scal_relative_tolerance{},
        my_vec_absolute_tolerance{ vec_absolute_tolerance },
        my_vec_relative_tolerance{ vec_relative_tolerance }
    {}

    /** Get the main set dimension.
     * @return main set dimension
     */
    int get_main_set_dimension() const
    {
        return my_main_set_dimension;
    }

    /** Get the relative tolerance for one component.
     * @param i component to select
     * @return relative tolerance for selected component
     */
    double get_relative_tolerance(const int& i) const
    {
        return my_vec_absolute_tolerance == NULL
            ? my_scal_relative_tolerance
            : my_vec_relative_tolerance[i];
    }

    /** Get the tolerance for one component.
     * @param i component to select
     * @param scale scale factor for relative tolerance (i.e. y[i])
     * @return tolerance for selected component
     */
    double get_tolerance(const int& i, const double& scale) 
    {
        return my_vec_absolute_tolerance == NULL
            ? my_scal_absolute_tolerance + my_scal_relative_tolerance * scale
            : my_vec_absolute_tolerance[i] + my_vec_relative_tolerance[i] * scale;
    }

    /** Get the tolerance for one component.
     * @param i component to select
     * @param scale scale factor for relative tolerance (i.e. y[i])
     * @param <T> type of the field elements
     * @return tolerance for selected component
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
     T get_tolerance(const int& i, const T& scale) 
    {
        return my_vec_absolute_tolerance == NULL
            ? scale.multiply(my_scal_relative_tolerance).add(my_scal_absolute_tolerance)
            : scale.multiply(my_vec_relative_tolerance[i]).add(my_vec_absolute_tolerance[i]);
    }

    /** Filter the integration step.
     * @param h signed step
     * @param forward forward integration indicator
     * @param accept_small if true, steps smaller than the minimal value
     * are silently increased up to this value, if false such small
     * steps generate an exception
     * @return a bounded integration step (h if no bound is reach, or a bounded value)
     * @exception  if the step is too small and accept_small is false
     */
    double filter_step(const double& h, const bool forward, const bool accept_small)
    {

        double filtered_h = h;
        if (std::abs(h) < my_min_step) 
        {
            if (accept_small) 
            {
                filtered_h = forward ? my_min_step : -my_min_step;
            }
            else 
            {
                throw std::exception("error not implemented");
                //throw (Localized_ODE_Formats.MINIMAL_STEPSIZE_REACHED_DURING_INTEGRATION, std::abs(h), min_step, true);
            }
        }

        if (filtered_h > my_max_step) 
        {
            filtered_h = my_max_step;
        }
        else if (filtered_h < -my_max_step) 
        {
            filtered_h = -my_max_step;
        }

        return filtered_h;

    }

    /** Filter the integration step.
     * @param h signed step
     * @param forward forward integration indicator
     * @param accept_small if true, steps smaller than the minimal value
     * are silently increased up to this value, if false such small
     * steps generate an exception
     * @param <T> type of the field elements
     * @return a bounded integration step (h if no bound is reach, or a bounded value)
     * @exception  if the step is too small and accept_small is false
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
    T filter_step(const T h, const bool forward, const bool accept_small)
    {

        T filtered_h = h;
        if (h.abs().subtract(min_step).get_real() < 0) 
        {
            if (accept_small) 
            {
                filtered_h = h.get_field().get_zero().add(forward ? min_step : -min_step);
            }
            else 
            {
                throw std::exception("error not implemented");
                //throw (Localized_ODE_Formats.MINIMAL_STEPSIZE_REACHED_DURING_INTEGRATION, std::abs(h.get_real()), min_step, true);
            }
        }

        if (filtered_h.subtract(max_step).get_real() > 0) 
        {
            filtered_h = h.get_field().get_zero().add(max_step);
        }
        else if (filtered_h.add(max_step).get_real() < 0) 
        {
            filtered_h = h.get_field().get_zero().add(-max_step);
        }

        return filtered_h;

    }

    /** Set the initial step size.
     * <p>This method allows the user to specify an initial positive
     * step size instead of letting the integrator guess it by
     * itself. If this method is not called before integration is
     * started, the initial step size will be estimated by the
     * integrator.</p>
     * @param initial_step_size initial step size to use (must be positive even
     * for backward integration ; providing a negative value or a value
     * outside of the min/max step interval will lead the integrator to
     * ignore the value and compute the initial step size by itself)
     */
    void set_initial_step_size(const double& initial_step_size) 
    {
        if ((initial_step_size < my_min_step) || (initial_step_size > my_max_step)) 
        {
            my_initial_step = -1.0;
        }
        else 
        {
            my_initial_step = initial_step_size;
        }
    }

    /** Get the initial step.
     * @return initial step
     */
    double get_initial_step() const
    {
        return my_initial_step;
    }

    /** Get the minimal step.
     * @return minimal step
     */
    double get_min_step() const
    {
        return my_min_step;
    }

    /** Get the maximal step.
     * @return maximal step
     */
    double get_max_step() const 
    {
        return my_max_step;
    }

    /** Get a dummy step size.
     * @return geometric mean of {@link #get_min_step()} and {@link #get_max_step()}
     */
    double get_dummy_stepsize() const
    {
        return std::sqrt(my_min_step * my_max_step);
    }

};