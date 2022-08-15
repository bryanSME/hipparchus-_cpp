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

//package org.hipparchus.ode;

//import java.io.Serializable;
#include <vector>

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;

/**
 * Class mapping the part of a complete state or derivative that pertains
 * to a specific differential equation.
 * <p>
 * Instances of this class are guaranteed to be immutable.
 * </p>
 * @see Secondary_ODE
 */
class Equations_mapper  
{
private:

    /** Start indices of the components. */
    std::vector<int> my_start;

    /** Create a mapper by adding a equation to another mapper.
     * <p>
     * The equation will have index {@code mapper.}{@link #get_number_of_equations()}, * or 0 if {@code mapper} is NULL.
     * </p>
     * @param mapper former mapper, with one equation less (null for first equation)
     * @param dimension dimension of the equation state vector
     */
    Equations_mapper(const Equations_mapper& mapper, const int& dimension) 
    {
        const int index = (mapper == NULL) ? 0 : mapper.get_number_of_equations();
        my_start = std::vector<int>(index + 2);
        if (mapper == NULL) 
        {
            my_start[0] = 0;
        }
        else 
        {
            System.arraycopy(mapper.start, 0, my_start, 0, index + 1);
        }
        my_start[index + 1] = my_start[index] + dimension;
    }

    /** Check equation index.
     * @param index index of the equation, must be between 0 included and
     * {@link #get_number_of_equations()} (excluded)
     * @exception  if index is out of range
     */
    void check_index(const int& index)
    {
        if (index < 0 || index > my_start.size() - 2)
        {
            throw std::exception("not implemented");
            //throw (Localized_Core_Formats.OUT_OF_RANGE_SIMPLE, index, 0, my_start.size() - 2);
        }
    }

public:
    /** Get the number of equations mapped.
     * @return number of equations mapped
     */
    int get_number_of_equations() 
    {
        return my_start.size() - 1;
    }

    /** Return the dimension of the complete set of equations.
     * <p>
     * The complete set of equations correspond to the primary set plus all secondary sets.
     * </p>
     * @return dimension of the complete set of equations
     */
    int get_total_dimension() const
    {
        return my_start[my_start.size() - 1];
    }

    /** Map flat arrays to a state and derivative.
     * @param t time
     * @param y state array to map, including primary and secondary components
     * @param y_dot state derivative array to map, including primary and secondary components
     * @return mapped state
     * @exception  if an array does not match total dimension
     */
    ODE_State_And_Derivative map_state_and_derivative(const double t, const std::vector<double> y, const std::vector<double> y_dot)
    {
        if (y.size() != get_total_dimension()) 
        {
            throw std::exception("not implemented");
            //throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, y.size(), get_total_dimension());
        }

        if (y_dot.size() != get_total_dimension()) 
        {
            throw std::exception("not implemented");
            //throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, y_dot.size(), get_total_dimension());
        }

        const int n = get_number_of_equations();
        const std::vector<double> state      = extract_equation_data(0, y);
        const std::vector<double> derivative = extract_equation_data(0, y_dot);
        if (n < 2) 
        {
            return ODE_State_And_Derivative(t, state, derivative);
        }
        else 
        {
            auto secondary_state = std::vector<std::vector<double>>(n - 1);
            auto secondary_derivative = std::vector<std::vector<double>>(n - 1);
            for (int index{ 1 }; index < get_number_of_equations(); ++index)
            {
                secondary_state[index - 1]      = extract_equation_data(index, y);
                secondary_derivative[index - 1] = extract_equation_data(index, y_dot);
            }
            return ODE_State_And_Derivative(t, state, derivative, secondary_state, secondary_derivative);
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
    std::vector<double> extract_equation_data(const int index, const std::vector<double> complete)
         
        {
        check_index(index);
        const int begin     = my_start[index];
        const int end       = my_start[index + 1];
        if (complete.size() < end) 
        {
            throw std::exception("not implemented");
            //throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, complete.size(), end);
        }
        const auto dimension = end - begin;
        const auto equation_data = std::vector<double>(dimension);
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
    void insert_equation_data(const int index, std::vector<double> equation_data, std::vector<double> complete)
         
        {
        check_index(index);
        const int begin     = my_start[index];
        const int end       = my_start[index + 1];
        const int dimension = end - begin;
        if (complete.size() < end) 
        {
            throw std::exception("not implemented");
            //throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, complete.size(), end);
        }
        if (equation_data.size() != dimension) 
        {
            throw std::exception("not implemented");
            //throw (Localized_Core_Formats.DIMENSIONS_MISMATCH, equation_data.size(), dimension);
        }
        System.arraycopy(equation_data, 0, complete, begin, dimension);
    }
};