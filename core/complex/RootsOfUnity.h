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

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Sin_Cos;

#include <vector>
#include <numbers>
#include "../util/SinCos.h"

/**
 * A helper class for the computation and caching of the {@code n}-th roots
 * of unity.
 */
class Roots_Of_Unity 
{
private:
    /** Number of roots of unity. */
    int my_omega_count;

    /** Real part of the roots. */
    std::vector<double> my_omega_real;

    /**
     * Imaginary part of the {@code n}-th roots of unity, for positive values
     * of {@code n}. In this array, the roots are stored in counter-clockwise
     * order.
     */
    std::vector<double> my_omega_imaginary_counter_clockwise;

    /**
     * Imaginary part of the {@code n}-th roots of unity, for negative values
     * of {@code n}. In this array, the roots are stored in clockwise order.
     */
    std::vector<double> my_omega_imaginary_clockwise;

    /**
     * {@code true} if {@link #compute_rootsstatic_cast<int>(} was called with a positive
     * value of its argument {@code n}. In this case, counter-clockwise ordering
     * of the roots of unity should be used.
     */
    bool my_is_counter_clock_wise;

public:
    /**
     * Build an engine for computing the {@code n}-th roots of unity.
     */
    Roots_Of_Unity() 
    {
        omega_count = 0;
        my_omega_real = NULL;
        my_omega_imaginary_counter_clockwise = NULL;
        my_omega_imaginary_clockwise = NULL;
        my_is_counter_clock_wise = true;
    }

    /**
     * Returns {@code true} if {@link #compute_rootsstatic_cast<int>(} was called with a
     * positive value of its argument {@code n}. If {@code true}, then
     * counter-clockwise ordering of the roots of unity should be used.
     *
     * @return {@code true} if the roots of unity are stored in counter-clockwise order
     * @Math_Illegal_State_Exception if no roots of unity have been computed yet
     */
    bool my_is_counter_clock_wise() 
    {
        if (omega_count == 0) 
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::ROOTS_OF_UNITY_NOT_COMPUTED_YET);
        }
        return my_is_counter_clock_wise;
        
    }

    /**
     * Computes the {@code n}-th roots of unity.
     * <p>
     * The roots are stored in {@code omega[]}, such that {@code omega[k] = w ^ k}, * where {@code k = 0, ..., n - 1}, {@code w = exp(2 * pi * i / n)} and
     * {@code i = sqrt(-1)}.
     * <p>
     * Note that {@code n} can be positive of negative
     * <ul>
     * <li>{@code abs(n)} is always the number of roots of unity.</li>
     * <li>If {@code n > 0}, then the roots are stored in counter-clockwise order.</li>
     * <li>If {@code n < 0}, then the roots are stored in clockwise order.</li>
     * </ul>
     *
     * @param n the (signed) number of roots of unity to be computed
     * @ if {@code n = 0}
     */
    void compute_roots(const int& n)  
    {

        if (n == 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::CANNOT_COMPUTE_0TH_ROOT_OF_UNITY);
        }


        my_is_counter_clock_wise = n > 0;

        // avoid repetitive calculations
        const int abs_n = std::abs(n);

        if (abs_n == omega_count) 
        {
            return;
        }

        // calculate everything from scratch
        const double t  = 2.0 * std::numbers::pi / abs_n;
        const Sin_Cos sc = Sin_Cos(t);
        my_omega_real = std::vector<double>(abs_n);
        my_omega_imaginary_counter_clockwise = std::vector<double>(abs_n);
        my_omega_imaginary_clockwise = std::vector<double>(abs_n);
        my_omega_real[0] = 1.0;
        my_omega_imaginary_counter_clockwise[0] = 0.0;
        my_omega_imaginary_clockwise[0] = 0.0;
        for (int i{ 1 }; i < abs_n; i++) 
        {
            my_omega_real[i] = my_omega_real[i - 1] * sc.cos() - my_omega_imaginary_counter_clockwise[i - 1] * sc.sin();
            my_omega_imaginary_counter_clockwise[i] = my_omega_real[i - 1] * sc.sin() + my_omega_imaginary_counter_clockwise[i - 1] * sc.cos();
            my_omega_imaginary_clockwise[i] = -my_omega_imaginary_counter_clockwise[i];
        }
        my_omega_count = abs_n;
        
    }

    /**
     * Get the real part of the {@code k}-th {@code n}-th root of unity.
     *
     * @param k index of the {@code n}-th root of unity
     * @return real part of the {@code k}-th {@code n}-th root of unity
     * @Math_Illegal_State_Exception if no roots of unity have been computed yet
     * @ if {@code k} is out of range
     */
    double get_real(const int& k) const
    {

        if (my_omega_count == 0) 
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::ROOTS_OF_UNITY_NOT_COMPUTED_YET);
        }
        if ((k < 0) || (k >= my_omega_count))
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_ROOT_OF_UNITY_INDEX, static_cast<int>(k), static_cast<int>(0), static_cast<int>(my_omega_count - 1));
        }

        return my_omega_real[k];
        
    }

    /**
     * Get the imaginary part of the {@code k}-th {@code n}-th root of unity.
     *
     * @param k index of the {@code n}-th root of unity
     * @return imaginary part of the {@code k}-th {@code n}-th root of unity
     * @Math_Illegal_State_Exception if no roots of unity have been computed yet
     * @ if {@code k} is out of range
     */
    double get_imaginary(const int& k)
    {
        if (my_omega_count == 0)
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::ROOTS_OF_UNITY_NOT_COMPUTED_YET);
        }
        if ((k < 0) || (k >= my_omega_count))
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_ROOT_OF_UNITY_INDEX, static_cast<int>(k), static_cast<int>(0), static_cast<int>(my_omega_count - 1));
        }

        return my_is_counter_clock_wise 
            ? my_omega_imaginary_counter_clockwise[k]
            : my_omega_imaginary_clockwise[k];
    }

    /**
     * Returns the number of roots of unity currently stored.
     * <p>
     * If {@link #compute_rootsstatic_cast<int>(} was called with {@code n}, then this method
     * returns {@code abs(n)}. If no roots of unity have been computed yet, this
     * method returns 0.
     *
     * @return the number of roots of unity currently stored
     */
    int get_number_of_roots() const
    {
        return my_omega_count;
    }
};