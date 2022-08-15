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

//package org.hipparchus.ode;

//import java.io.Serializable;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <type_traits>
#include "../core/CalculusFieldElement.h"

/**
 * Class mapping the part of a complete state or derivative that pertains
 * to a set of differential equations.
 * <p>
 * Instances of this class are guaranteed to be immutable.
 * </p>
 * @see FieldExpandable_ODE
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class FieldEquations_mapper  
{

    /** Serializable UID. */
    private static const long serial_version_uid = 20151114L;

    /** Start indices of the components. */
    private const std::vector<int> start;

    /** Create a mapper by adding a equation to another mapper.
     * <p>
     * The equation will have index {@code mapper.}{@link #get_number_of_equations()}, * or 0 if {@code mapper} is NULL.
     * </p>
     * @param mapper former mapper, with one equation less (null for first equation)
     * @param dimension dimension of the equation state vector
     */
    FieldEquations_mapper(const FieldEquations_mapper<T> mapper, const int& dimension) 
    {
        const int index = (mapper == NULL) ? 0 : mapper.get_number_of_equations();
        this.start = int[index + 2];
        if (mapper == NULL) 
        {
            start[0] = 0;
        }
else 
        {
            System.arraycopy(mapper.start, 0, start, 0, index + 1);
        }
        start[index + 1] = start[index] + dimension;
    }

    /** Get the number of equations mapped.
     * @return number of equations mapped
     */
    public int get_number_of_equations() 
    {
        return start.size() - 1;
    }

    /** Return the dimension of the complete set of equations.
     * <p>
     * The complete set of equations correspond to the primary set plus all secondary sets.
     * </p>
     * @return dimension of the complete set of equations
     */
    public int get_total_dimension() 
    {
        return start[start.size() - 1];
    }

    /** Map flat arrays to a state and derivative.
     * @param t time
     * @param y state array to map, including primary and secondary components
     * @param y_dot state derivative array to map, including primary and secondary components
     * @return mapped state
     * @exception  if an array does not match total dimension
     */
    public Field_ODE_State_And_Derivative<T> map_state_and_derivative(const T t, const std::vector<T> y, const std::vector<T> y_dot)
         
        {

        if (y.size() != get_total_dimension()) 
        {
            throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, y.size(), get_total_dimension());
        }

        if (y_dot.size() != get_total_dimension()) 
        {
            throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, y_dot.size(), get_total_dimension());
        }

        const int n = get_number_of_equations();
        const std::vector<T> state      = extract_equation_data(0, y);
        const std::vector<T> derivative = extract_equation_data(0, y_dot);
        if (n < 2) 
        {
            return Field_ODE_State_And_Derivative<T>(t, state, derivative);
        }
else 
        {
            const std::vector<std::vector<T>> secondary_state      = Math_Arrays::build_array(t.get_field(), n - 1, -1);
            const std::vector<std::vector<T>> secondary_derivative = Math_Arrays::build_array(t.get_field(), n - 1, -1);
            for (const int& index = 1; index < get_number_of_equations(); ++index) 
            {
                secondary_state[index - 1]      = extract_equation_data(index, y);
                secondary_derivative[index - 1] = extract_equation_data(index, y_dot);
            }
            return Field_ODE_State_And_Derivative<T>(t, state, derivative, secondary_state, secondary_derivative);
        }
    }

    /** Extract equation data from a complete state or derivative array.
     * @param index index of the equation, must be between 0 included and
     * {@link #get_number_of_equations()} (excluded)
     * @param complete complete state or derivative array from which
     * equation data should be retrieved
     * @return equation data
     * @exception  if index is out of range
     * @exception  if complete state has not enough elements
     */
    public std::vector<T> extract_equation_data(const int index, const std::vector<T> complete)
         
        {
        check_index(index);
        const int begin     = start[index];
        const int end       = start[index + 1];
        if (complete.size() < end) 
        {
            throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, complete.size(), end);
        }
        const int dimension = end - begin;
        const std::vector<T> equation_data = Math_Arrays::build_array(complete[0].get_field(), dimension);
        System.arraycopy(complete, begin, equation_data, 0, dimension);
        return equation_data;
    }

    /** Insert equation data into a complete state or derivative array.
     * @param index index of the equation, must be between 0 included and
     * {@link #get_number_of_equations()} (excluded)
     * @param equation_data equation data to be inserted into the complete array
     * @param complete placeholder where to put equation data (only the
     * part corresponding to the equation will be overwritten)
     * @exception  if either array has not enough elements
     */
    public void insert_equation_data(const int index, std::vector<T> equation_data, std::vector<T> complete)
         
        {
        check_index(index);
        const int begin     = start[index];
        const int end       = start[index + 1];
        const int dimension = end - begin;
        if (complete.size() < end) 
        {
            throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, complete.size(), end);
        }
        if (equation_data.size() != dimension) 
        {
            throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, equation_data.size(), dimension);
        }
        System.arraycopy(equation_data, 0, complete, begin, dimension);
    }

    /** Check equation index.
     * @param index index of the equation, must be between 0 included and
     * {@link #get_number_of_equations()} (excluded)
     * @exception  if index is out of range
     */
    private void check_index(const int index)  
    {
        Math_Utils::check_range_inclusive(index, 0, start.size() - 2);
    }

};