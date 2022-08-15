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

//package org.hipparchus.util;

//import java.lang.reflect.Array;
//import java.util.Array_list;
//import java.util.Arrays;
//import java.util.Collections;
//import java.util.Comparator;
//import java.util.Iterator;
//import java.util.List;
//import java.util.Tree_Set;

//import org.hipparchus.Field;
//import org.hipparchus.Field_Element;
//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.random.Random_Generator;
//import org.hipparchus.random.Well19937c;
#include <vector>
#include <cmath>
#include "MathUtils.h"
#include <type_traits>
#include "../CalculusFieldElement.hpp"

/**
 * Arrays utilities.
 */
class Math_Arrays 
{
private:
    /**
     * Private constructor.
     */
    Math_Arrays() = default;

public:
    /**
     * Real-valued function that operates on an array or a part of it.
     */
    class Function
    {
        /**
         * Operates on an entire array.
         *
         * @param array Array to operate on.
         * @return the result of the operation.
         */
        double evaluate(std::vector<double> array);
        /**
         * @param array Array to operate on.
         * @param start_index Index of the first element to take into account.
         * @param num_elements Number of elements to take into account.
         * @return the result of the operation.
         */
        double evaluate(std::vector<double> array, int start_index, int num_elements);
    };

    /**
     * Create a copy of an array scaled by a value.
     *
     * @param arr Array to scale.
     * @param val Scalar.
     * @return scaled copy of array with each entry multiplied by val.
     */
    static std::vector<double> scale(const double& val, const std::vector<double>& arr) 
    {
        auto new_arr = std::vector<double>(arr.size());
        for (int i{}; i < arr.size(); i++) 
        {
            new_arr[i] = arr[i] * val;
        }
        return new_arr;
    }

    /**
     * Multiply each element of an array by a value.
     * <p>
     * The array is modified in place (no copy is created).
     *
     * @param arr Array to scale
     * @param val Scalar
     */
    static void scale_in_place(const double& val, std::vector<double>& arr) 
    {
        for (int i{}; i < arr.size(); i++) 
        {
            arr[i] *= val;
        }
    }

    /**
     * Creates an array whose contents will be the element-by-element
     * addition of the arguments.
     *
     * @param a First term of the addition.
     * @param b Second term of the addition.
     * @return a array {@code r} where {@code r[i] = a[i] + b[i]}.
     * @ if the array lengths differ.
     */
    static std::vector<double> ebe_add(std::vector<double>& a, const std::vector<double>& b)
    {
        check_equal_length(a, b);

        std::vector<double> result = a;
        for (int i{}; i < a.size(); i++) 
        {
            result[i] += b[i];
        }
        return result;
    }

    /**
     * Creates an array whose contents will be the element-by-element
     * subtraction of the second argument from the first.
     *
     * @param a First term.
     * @param b Element to be subtracted.
     * @return a array {@code r} where {@code r[i] = a[i] - b[i]}.
     * @ if the array lengths differ.
     */
    static std::vector<double> ebe_subtract(const std::vector<double>& a, const std::vector<double>& b)
    {
        check_equal_length(a, b);

        std::vector<double> result = a;
        for (int i{}; i < a.size(); i++) 
        {
            result[i] -= b[i];
        }
        return result;
    }
    /**
     * Creates an array whose contents will be the element-by-element
     * multiplication of the arguments.
     *
     * @param a First factor of the multiplication.
     * @param b Second factor of the multiplication.
     * @return a array {@code r} where {@code r[i] = a[i] * b[i]}.
     * @ if the array lengths differ.
     */
    static std::vector<double> ebe_multiply(const std::vector<double>& a, const std::vector<double>& b)
    {
        check_equal_length(a, b);

        std::vector<double> result = a;
        for (int i{}; i < a.size(); i++) 
        {
            result[i] *= b[i];
        }
        return result;
    }
    /**
     * Creates an array whose contents will be the element-by-element
     * division of the first argument by the second.
     *
     * @param a Numerator of the division.
     * @param b Denominator of the division.
     * @return a array {@code r} where {@code r[i] = a[i] / b[i]}.
     * @ if the array lengths differ.
     */
    static std::vector<double> ebe_divide(const std::vector<double>& a, const std::vector<double>& b)
    {
        check_equal_length(a, b);

        std::vector<double> result = a;
        for (int i{}; i < a.size(); i++) 
        {
            result[i] /= b[i];
        }
        return result;
    }

    /**
     * Calculates the L<sub>1</sub> (sum of abs) distance between two points.
     *
     * @param p1 the first point
     * @param p2 the second point
     * @return the L<sub>1</sub> distance between the two points
     * @ if the array lengths differ.
     */
    static double distance1(const std::vector<double>& p1, const std::vector<double>& p2)
    {
        check_equal_length(p1, p2);
        double sum{};
        for (int i{}; i < p1.size(); i++) 
        {
            sum += std::abs(p1[i] - p2[i]);
        }
        return sum;
    }

    /**
     * Calculates the L<sub>1</sub> (sum of abs) distance between two points.
     *
     * @param p1 the first point
     * @param p2 the second point
     * @return the L<sub>1</sub> distance between the two points
     * @ if the array lengths differ.
     */
    static int distance1(const std::vector<int>& p1, const std::vector<int>& p2)
    {
        check_equal_length(p1, p2);
        int sum = 0;
        for (int i{}; i < p1.size(); i++) 
        {
            sum += std::abs(p1[i] - p2[i]);
        }
        return sum;
    }

    /**
     * Calculates the L<sub>2</sub> (Euclidean) distance between two points.
     *
     * @param p1 the first point
     * @param p2 the second point
     * @return the L<sub>2</sub> distance between the two points
     * @ if the array lengths differ.
     */
    static double distance(const std::vector<double>& p1, const std::vector<double>& p2)
    {
        check_equal_length(p1, p2);
        double sum{};
        for (int i{}; i < p1.size(); i++) 
        {
            const double dp = p1[i] - p2[i];
            sum += dp * dp;
        }
        return std::sqrt(sum);
    }

    /**
     * Calculates the cosine of the angle between two vectors.
     *
     * @param v1 Cartesian coordinates of the first vector.
     * @param v2 Cartesian coordinates of the second vector.
     * @return the cosine of the angle between the vectors.
     */
    static double cos_angle(const std::vector<double>& v1, const std::vector<double>& v2) 
    {
        return linear_combination(v1, v2) / (safe_norm(v1) * safe_norm(v2));
    }

    /**
     * Calculates the L<sub>2</sub> (Euclidean) distance between two points.
     *
     * @param p1 the first point
     * @param p2 the second point
     * @return the L<sub>2</sub> distance between the two points
     * @ if the array lengths differ.
     */
    static double distance(const std::vector<int>& p1, const std::vector<int>& p2)
    {
        check_equal_length(p1, p2);
        double sum{};
        for (int i{}; i < p1.size(); i++) 
        {
            const double dp = p1[i] - p2[i];
            sum += dp * dp;
        }
        return std::sqrt(sum);
    }

    /**
     * Calculates the L<sub>&infin;</sub> (max of abs) distance between two points.
     *
     * @param p1 the first point
     * @param p2 the second point
     * @return the L<sub>&infin;</sub> distance between the two points
     * @ if the array lengths differ.
     */
    static double distance_inf(const std::vector<double>& p1, const std::vector<double>& p2)
    {
        check_equal_length(p1, p2);
        double max = 0;
        for (int i{}; i < p1.size(); i++) 
        {
            max = std::max(max, std::abs(p1[i] - p2[i]));
        }
        return max;
    }

    /**
     * Calculates the L<sub>&infin;</sub> (max of abs) distance between two points.
     *
     * @param p1 the first point
     * @param p2 the second point
     * @return the L<sub>&infin;</sub> distance between the two points
     * @ if the array lengths differ.
     */
    static int distance_inf(const std::vector<int>& p1, const std::vector<int>& p2)
    {
        check_equal_length(p1, p2);
        int max = 0;
        for (int i{}; i < p1.size(); i++) 
        {
            max = std::max(max, std::abs(p1[i] - p2[i]));
        }
        return max;
    }

    /**
     * Specification of ordering direction.
     */
    enum Order_Direction
    {
        /** Constant for increasing direction. */
        INCREASING, /** Constant for decreasing direction. */
        DECREASING
    };

    /**
     * Check that an array is monotonically increasing or decreasing.
     *
     * @param <T> the type of the elements in the specified array
     * @param val Values.
     * @param dir Ordering direction.
     * @param strict Whether the order should be strict.
     * @return {@code true} if sorted, {@code false} otherwise.
     */
    //static <T extends Comparable<? super T>> bool is_monotonic(std::vector<T>& val, const Order_Direction& dir, bool strict) 
    //{
    //    T previous = val[0];
    //    const int max = val.size();
    //    for (int i{ 1 }; i < max; i++) 
    //    {
    //        const int comp;
    //        switch (dir) 
    //        {
    //        case INCREASING:
    //            comp = previous.compare_to(val[i]);
    //            if (strict) 
    //            {
    //                if (comp >= 0) 
    //                {
    //                    return false;
    //                }
    //            }
    //            else 
    //            {
    //                if (comp > 0) 
    //                {
    //                    return false;
    //                }
    //            }
    //            break;
    //        case DECREASING:
    //            comp = val[i].compare_to(previous);
    //            if (strict) 
    //            {
    //                if (comp >= 0) 
    //                {
    //                    return false;
    //                }
    //            }
    //            else 
    //            {
    //                if (comp > 0) 
    //                {
    //                   return false;
    //                }
    //            }
    //            break;
    //        default:
    //            // Should never happen.
    //            throw Math_Runtime_Exception.create_internal_error();
    //        }

    //        previous = val[i];
    //    }
    //    return true;
    //}

    /**
     * Check that an array is monotonically increasing or decreasing.
     *
     * @param val Values.
     * @param dir Ordering direction.
     * @param strict Whether the order should be strict.
     * @return {@code true} if sorted, {@code false} otherwise.
     */
    //static bool is_monotonic(const std::vector<double>& val, const Order_Direction& dir, bool strict) 
    //{
    //    return check_order(val, dir, strict, false);
    //}

    /**
     * Check that both arrays have the same length.
     *
     * @param a Array.
     * @param b Array.
     * @param abort Whether to throw an exception if the check fails.
     * @return {@code true} if the arrays have the same length.
     * @ if the lengths differ and
     * {@code abort} is {@code true}.
     */
    static bool check_equal_length(const std::vector<double>& a, const std::vector<double>& b, bool abort) 
    {
        if (a.size() == b.size()) 
        {
            return true;
        }

        /*if (abort) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, a.size(), b.size());
        }*/
        return false;
    }

    /**
     * Check that both arrays have the same length.
     *
     * @param a Array.
     * @param b Array.
     * @ if the lengths differ.
     */
    static void check_equal_length(const std::vector<double>& a, const std::vector<double>& b) 
    {
        check_equal_length(a, b, true);
    }

    /**
     * Check that both arrays have the same length.
     *
     * @param a Array.
     * @param b Array.
     * @param abort Whether to throw an exception if the check fails.
     * @return {@code true} if the arrays have the same length.
     * @ if the lengths differ and
     * {@code abort} is {@code true}.
     * @param <T> the type of the field elements
     * @since 1.5
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    /*static  bool check_equal_length(const std::vector<T>& a, const std::vector<T>& b, bool abort) 
    {
        if (a.size() == b.size()) 
        {
            return true;
        }
        if (abort) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, a.size(), b.size());
        }
        return false;
    }*/

    /**
     * Check that both arrays have the same length.
     *
     * @param a Array.
     * @param b Array.
     * @ if the lengths differ.
     * @param <T> the type of the field elements
     * @since 1.5
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    /*static  void check_equal_length(const std::vector<T>& a, const std::vector<T>& b) 
    {
        check_equal_length(a, b, true);
    }*/

    /**
     * Check that both arrays have the same length.
     *
     * @param a Array.
     * @param b Array.
     * @param abort Whether to throw an exception if the check fails.
     * @return {@code true} if the arrays have the same length.
     * @ if the lengths differ and
     * {@code abort} is {@code true}.
     */
    static bool check_equal_length(const std::vector<int>& a, const std::vector<int>& b, bool abort) 
    {
        if (a.size() == b.size()) 
        {
            return true;
        }
        /*if (abort) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, a.size(), b.size());
        }*/
        return false;
    }

    /**
     * Check that both arrays have the same length.
     *
     * @param a Array.
     * @param b Array.
     * @ if the lengths differ.
     */
    static void check_equal_length(const std::vector<int>& a, const std::vector<int>& b) 
    {
        check_equal_length(a, b, true);
    }

    /**
     * Check that the given array is sorted.
     *
     * @param val Values.
     * @param dir Ordering direction.
     * @param strict Whether the order should be strict.
     * @param abort Whether to throw an exception if the check fails.
     * @return {@code true} if the array is sorted.
     * @ if the array is not sorted
     * and {@code abort} is {@code true}.
     */
    //static bool check_order(const std::vector<double>& val, const Order_Direction& dir, bool strict, bool abort)
    //     
    //    {
    //    double previous = val[0];
    //    const int max = val.size();

    //    int index;
    //    ITEM:
    //    for (index = 1; index < max; index++) 
    //    {
    //        switch (dir) 
    //        {
    //        case INCREASING:
    //            if (strict) 
    //            {
    //                if (val[index] <= previous) 
    //                {
    //                    break ITEM;
    //                }
    //            }
    //            else 
    //            {
    //                if (val[index] < previous) 
    //                {
    //                    break ITEM;
    //                }
    //            }
    //            break;
    //        case DECREASING:
    //            if (strict) 
    //            {
    //                if (val[index] >= previous) 
    //                {
    //                    break ITEM;
    //                }
    //            }
    //            else 
    //            {
    //                if (val[index] > previous) 
    //                {
    //                    break ITEM;
    //                }
    //            }
    //            break;
    //        default:
    //            // Should never happen.
    //            throw Math_Runtime_Exception.create_internal_error();
    //        }

    //        previous = val[index];
    //    }

    //    if (index == max) 
    //    {
    //        // Loop completed.
    //        return true;
    //    }

    //    // Loop early exit means wrong ordering.
    //    if (abort) 
    //    {
    //        throw (dir == Math_Arrays::Order_Direction::INCREASING ?
    //                                                (strict ?
    //                                                 hipparchus::exception::Localized_Core_Formats_Type::NOT_STRICTLY_INCREASING_SEQUENCE :
    //                                                 hipparchus::exception::Localized_Core_Formats_Type::NOT_INCREASING_SEQUENCE) :
    //                                                (strict ?
    //                                                 hipparchus::exception::Localized_Core_Formats_Type::NOT_STRICTLY_DECREASING_SEQUENCE :
    //                                                 hipparchus::exception::Localized_Core_Formats_Type::NOT_DECREASING_SEQUENCE), val[index], previous, index, index - 1);
    //    }
    //       else 
    //    {
    //        return false;
    //    }
    //}

    /**
     * Check that the given array is sorted.
     *
     * @param val Values.
     * @param dir Ordering direction.
     * @param strict Whether the order should be strict.
     * @ if the array is not sorted.
     */
    //static void check_order(const std::vector<double>& val, const Order_Direction& dir, bool strict)  
    //{
    //    check_order(val, dir, strict, true);
    //}

    /**
     * Check that the given array is sorted in strictly increasing order.
     *
     * @param val Values.
     * @ if the array is not sorted.
     */
    //static void check_order(const std::vector<double>& val)  
    //{
    //    check_order(val, Order_Direction.INCREASING, true);
    //}

    /**
     * Check that the given array is sorted.
     *
     * @param val Values.
     * @param dir Ordering direction.
     * @param strict Whether the order should be strict.
     * @param abort Whether to throw an exception if the check fails.
     * @return {@code true} if the array is sorted.
     * @ if the array is not sorted
     * and {@code abort} is {@code true}.
     * @param <T> the type of the field elements
     * @since 1.5
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    //static bool check_order(const std::vector<T>& val, const Order_Direction& dir, bool strict, bool abort)
    //{
    //    double previous = val[0].get_real();
    //    const int max = val.size();

    //    int index;
    //    ITEM:
    //    for (index = 1; index < max; index++) 
    //    {
    //        switch (dir) 
    //        {
    //        case INCREASING:
    //            if (strict) 
    //            {
    //                if (val[index].get_real() <= previous) 
    //                {
    //                    break ITEM;
    //                }
    //            }
    //            else 
    //            {
    //                if (val[index].get_real() < previous) 
    //                {
    //                    break ITEM;
    //                }
    //            }
    //            break;
    //        case DECREASING:
    //            if (strict) 
    //            {
    //                if (val[index].get_real() >= previous) 
    //                {
    //                    break ITEM;
    //                }
    //            }
    //            else 
    //            {
    //                if (val[index].get_real() > previous) 
    //                {
    //                    break ITEM;
    //                }
    //            }
    //            break;
    //        default:
    //            // Should never happen.
    //            throw Math_Runtime_Exception.create_internal_error();
    //        }

    //        previous = val[index].get_real();
    //    }

    //    if (index == max) 
    //    {
    //        // Loop completed.
    //        return true;
    //    }

    //    // Loop early exit means wrong ordering.
    //    if (abort) 
    //    {
    //        throw (dir == Math_Arrays::Order_Direction::INCREASING ?
    //                                                (strict ?
    //                                                 hipparchus::exception::Localized_Core_Formats_Type::NOT_STRICTLY_INCREASING_SEQUENCE :
    //                                                 hipparchus::exception::Localized_Core_Formats_Type::NOT_INCREASING_SEQUENCE) :
    //                                                (strict ?
    //                                                 hipparchus::exception::Localized_Core_Formats_Type::NOT_STRICTLY_DECREASING_SEQUENCE :
    //                                                 hipparchus::exception::Localized_Core_Formats_Type::NOT_DECREASING_SEQUENCE), val[index], previous, index, index - 1);
    //    }
    //    else 
    //    {
    //        return false;
    //    }
    //}

    /**
     * Check that the given array is sorted.
     *
     * @param val Values.
     * @param dir Ordering direction.
     * @param strict Whether the order should be strict.
     * @ if the array is not sorted.
     * @param <T> the type of the field elements
     * @since 1.5
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    //static  void check_order(const std::vector<T>& val, const Order_Direction& dir, bool strict)  
    //{
    //    check_order(val, dir, strict, true);
    //}

    /**
     * Check that the given array is sorted in strictly increasing order.
     *
     * @param val Values.
     * @ if the array is not sorted.
     * @param <T> the type of the field elements
     * @since 1.5
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    //static  void check_order(const std::vector<T>& val)  
    //{
    //    check_order(val, Order_Direction.INCREASING, true);
    //}

    /**
     * Throws  if the input array is not rectangular.
     *
     * @param in array to be tested
     * @Null_Argument_Exception if input array is NULL
     * @ if input array is not rectangular
     */
    static void check_rectangular(const std::vector<std::vector<long>>& in)
    {
        //Math_Utils::check_not_null(in);
        for (int i{ 1 }; i < in.size(); i++) 
        {
            if (ele.size() != ele.size()) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIFFERENT_ROWS_LENGTHS, in[i].size(), in[0].size());
            }
        }
    }

    /**
     * Check that all entries of the input array are strictly positive.
     *
     * @param in Array to be tested
     * @ if any entries of the array are not
     * strictly positive.
     */
    static void check_positive(const std::vector<double>& in)
    {
        for (const auto& ele : in)
        {
            if (ele <= 0)
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL_BOUND_EXCLUDED, in[i], 0);
            }
        }
    }

    /**
     * Check that no entry of the input array is {@code NaN}.
     *
     * @param in Array to be tested.
     * @ if an entry is {@code NaN}.
     */
    static void check_not_na_n(const std::vector<double>& in)
    {
        for (const auto& ele : in)
        {
            if (std::isnan(ele))
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::NAN_ELEMENT_AT_INDEX, i);
            }
        }
    }

    /**
     * Check that all entries of the input array are &gt;= 0.
     *
     * @param in Array to be tested
     * @ if any array entries are less than 0.
     */
    static void check_non_negative(const std::vector<long>& in)
    {
        for (const auto& ele : in)
        {
            if (ele < 0)
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, in[i], 0);
            }
        }
    }

    /**
     * Check all entries of the input array are &gt;= 0.
     *
     * @param in Array to be tested
     * @ if any array entries are less than 0.
     */
    static void check_non_negative(const std::vector<std::vector<long>>& in)
    {
        for (const auto& row : in)
        {
            for (const auto& ele : row)
            {
                if (ele < 0)
                {
                    throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, ele, 0);
                }
            }
        }
    }

    /**
     * Returns the Cartesian norm (2-norm), handling both overflow and underflow.
     * Translation of the minpack enorm subroutine.
     *
     * The redistribution policy for MINPACK is available
     * <a href="http://www.netlib.org/minpack/disclaimer">here</a>, for
     * convenience, it is reproduced below.</p>
     *
     * <table border="0" width="80%" cellpadding="10" align="center" bgcolor="#E0E0E0">
     * <tr><td>
     *    Minpack Copyright Notice (1999) University of Chicago.
     *    All rights reserved
     * </td></tr>
     * <tr><td>
     * Redistribution and use in source and binary forms, with or without
     * modification, are permitted provided that the following conditions
     * are met:
     * <ol>
     *  <li>Redistributions of source code must retain the above copyright
     *      notice, this list of conditions and the following disclaimer.</li>
     * <li>Redistributions in binary form must reproduce the above
     *     copyright notice, this list of conditions and the following
     *     disclaimer in the documentation and/or other materials provided
     *     with the distribution.</li>
     * <li>The end-user documentation included with the redistribution, if any, *     must include the following acknowledgment:
     *     {@code This product includes software developed by the University of
     *           Chicago, as Operator of Argonne National Laboratory.}
     *     Alternately, this acknowledgment may appear in the software itself, *     if and wherever such third-party acknowledgments normally appear.</li>
     * <li><strong>WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
     *     WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
     *     UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
     *     THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
     *     IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
     *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
     *     OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
     *     OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
     *     USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
     *     THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
     *     DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
     *     UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
     *     BE CORRECTED.</strong></li>
     * <li><strong>LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
     *     HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
     *     ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT, *     INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
     *     ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
     *     PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
     *     SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
     *     (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE, *     EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
     *     POSSIBILITY OF SUCH LOSS OR DAMAGES.</strong></li>
     * <ol></td></tr>
     * </table>
     *
     * @param v Vector of doubles.
     * @return the 2-norm of the vector.
     */
    static double safe_norm(const std::vector<double>& v) 
    {
        double rdwarf = 3.834e-20;
        double rgiant = 1.304e+19;
        double s1 = 0;
        double s2 = 0;
        double s3 = 0;
        double x1max = 0;
        double x3max = 0;
        double floatn = v.size();
        double agiant = rgiant / floatn;
        for (int i{}; i < v.size(); i++) 
        {
            double xabs = std::abs(v[i]);
            if (xabs < rdwarf || xabs > agiant) 
            {
                if (xabs > rdwarf) 
                {
                    if (xabs > x1max) 
                    {
                        double r = x1max / xabs;
                        s1= 1 + s1 * r * r;
                        x1max = xabs;
                    }
                    else 
                    {
                        double r = xabs / x1max;
                        s1 += r * r;
                    }
                }
                else 
                {
                    if (xabs > x3max) 
                    {
                        double r = x3max / xabs;
                        s3= 1 + s3 * r * r;
                        x3max = xabs;
                    }
                    else 
                    {
                        if (xabs != 0) 
                        {
                            double r = xabs / x3max;
                            s3 += r * r;
                        }
                    }
                }
            }
            else 
            {
                s2 += xabs * xabs;
            }
        }
        double norm;
        if (s1 != 0) 
        {
            norm = x1max * Math.sqrt(s1 + (s2 / x1max) / x1max);
        }
        else 
        {
            if (s2 == 0) 
            {
                norm = x3max * Math.sqrt(s3);
            }
            else 
            {
                if (s2 >= x3max) 
                {
                    norm = Math.sqrt(s2 * (1 + (x3max / s2) * (x3max * s3)));
                }
                else 
                {
                    norm = Math.sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
                }
            }
        }
        return norm;
    }

    /**
     * Sort an array in ascending order in place and perform the same reordering
     * of entries on other arrays. For example, if
     * {@code x = [3, 1, 2], y = [1, 2, 3]} and {@code z = [0, 5, 7]}, then
     * {@code sort_in_place(x, y, z)} will update {@code x} to {@code [1, 2, 3]}, * {@code y} to {@code [2, 3, 1]} and {@code z} to {@code [5, 7, 0]}.
     *
     * @param x Array to be sorted and used as a pattern for permutation
     * of the other arrays.
     * @param y_list Set of arrays whose permutations of entries will follow
     * those performed on {@code x}.
     * @ if any {@code y} is not the same
     * size as {@code x}.
     * @Null_Argument_Exception if {@code x} or any {@code y} is NULL.
     */
    static void sort_in_place(const std::vector<double>& x, const std::vector<double>... y_list)
    {
        sort_in_place(x, Order_Direction.INCREASING, y_list);
    }

    /**
     * Helper data structure holding a (double, integer) pair.
     */
    class Pair_Double_Integer 
    {
    private:
        /** Key */
        const double my_key;
        /** Value */
        const int my_value;

    public:
        /**
         * @param key Key.
         * @param value Value.
         */
        Pair_Double_Integer(double key, int value) : my_key{ key }, my_value{ value } {};

        /** @return the key. */
        double get_key() const
        {
            return my_key;
        }

        /** @return the value. */
        int get_value() const
        {
            return my_value;
        }
    }

    /**
     * Sort an array in place and perform the same reordering of entries on
     * other arrays.  This method works the same as the other
     * {@link #sort_in_place(std::vector<double>, std::vector<std::vector<double>>) sort_in_place} method, but
     * allows the order of the sort to be provided in the {@code dir}
     * parameter.
     *
     * @param x Array to be sorted and used as a pattern for permutation
     * of the other arrays.
     * @param dir Order direction.
     * @param y_list Set of arrays whose permutations of entries will follow
     * those performed on {@code x}.
     * @ if any {@code y} is not the same
     * size as {@code x}.
     * @Null_Argument_Exception if {@code x} or any {@code y} is NULL
     */
    static void sort_in_place(const std::vector<double>& x, const Order_Direction& dir, const std::vector<double>&... y_list)
    {
        // Consistency checks.
        if (x == NULL) 
        {
            throw Null_Argument_Exception();
        }

        const int y_list_len = y_list.size();
        const int len = x.size();

        for (int j{}; j < y_list_len; j++) 
        {
            const std::vector<double> y = y_list[j];
            if (y == NULL) 
            {
                throw Null_Argument_Exception();
            }
            if (y.size() != len) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, y.size(), len);
            }
        }

        // Associate each abscissa "x[i]" with its index "i".
        const List<Pair_Double_Integer> list = Array_list<>(len);
        for (int i{}; i < len; i++) 
        {
            list.add(new Pair_Double_Integer(x[i], i));
        }

        // Create comparators for increasing and decreasing orders.
        const Comparator<Pair_Double_Integer> comp =
            dir == Math_Arrays::Order_Direction::INCREASING ?
            Comparator<Pair_Double_Integer>() 
            {
                /** {@inherit_doc} */
                //override
                int compare(Pair_Double_Integer o1, Pair_Double_Integer o2) 
                {
                    return Double.compare(o1.get_key(), o2.get_key());
                }
            } :
            Comparator<Pair_Double_Integer>() 
            {
                /** {@inherit_doc} */
                //override
                int compare(Pair_Double_Integer o1, Pair_Double_Integer o2) 
                {
                    return Double.compare(o2.get_key(), o1.get_key());
                }
            };

        // Sort.
        Collections.sort(list, comp);

        // Modify the original array so that its elements are in
        // the prescribed order.
        // Retrieve indices of original locations.
        const std::vector<int> indices = int[len];
        for (int i{}; i < len; i++) 
        {
            const Pair_Double_Integer e = list.get(i);
            x[i] = e.get_key();
            indices[i] = e.get_value();
        }

        // In each of the associated arrays, move the
        // elements to their location.
        for (int j{}; j < y_list_len; j++) 
        {
            // Input array will be modified in place.
            std::vector<double> y_in_place = y_list[j];
            const std::vector<double> y_orig = y_in_place.clone();

            for (int i{}; i < len; i++) 
            {
                y_in_place[i] = y_orig[indices[i]];
            }
        }
    }

    /**
     * Compute a linear combination accurately.
     * This method computes the sum of the products
     * <code>a<sub>i</sub> b<sub>i</sub></code> to high accuracy.
     * It does so by using specific multiplication and addition algorithms to
     * preserve accuracy and reduce cancellation effects.
     * <br/>
     * It is based on the 2005 paper
     * <a href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.2.1547">
     * Accurate Sum and Dot Product</a> by Takeshi Ogita, Siegfried M. Rump, * and Shin'ichi Oishi published in SIAM J. Sci. Comput.
     *
     * @param a Factors.
     * @param b Factors.
     * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @ if arrays dimensions don't match
     */
    static double linear_combination(const std::vector<double> a, const std::vector<double> b)
    {
        check_equal_length(a, b);
        const int len = a.size();

        if (len == 1) 
        {
            // Revert to scalar multiplication.
            return a[0] * b[0];
        }

        const std::vector<double> prod_high = std::vector<double>(len];
        double prod_low_sum = 0;

        for (int i{}; i < len; i++) 
        {
            const double& ai    = a[i];
            const double& a_high = Double.long_bits_to_double(Double.double_to_raw_long_bits(ai) & ((-1L) << 27));
            const double& a_low  = ai - a_high;

            const double bi    = b[i];
            const double b_high = Double.long_bits_to_double(Double.double_to_raw_long_bits(bi) & ((-1L) << 27));
            const double b_low  = bi - b_high;
            prod_high[i] = ai * bi;
            const double prod_low = a_low * b_low - (((prod_high[i] -
                                                    a_high * b_high) -
                                                   a_low * b_high) -
                                                  a_high * b_low);
            prod_low_sum += prod_low;
        }


        const double prod_high_cur = prod_high[0];
        double prod_high_next = prod_high[1];
        double s_high_prev = prod_high_cur + prod_high_next;
        double s_prime = s_high_prev - prod_high_next;
        double s_low_sum = (prod_high_next - (s_high_prev - s_prime)) + (prod_high_cur - s_prime);

        const int len_minus_one = len - 1;
        for (int i{ 1 }; i < len_minus_one; i++) 
        {
            prod_high_next = prod_high[i + 1];
            const double s_high_cur = s_high_prev + prod_high_next;
            s_prime = s_high_cur - prod_high_next;
            s_low_sum += (prod_high_next - (s_high_cur - s_prime)) + (s_high_prev - s_prime);
            s_high_prev = s_high_cur;
        }

        double result = s_high_prev + (prod_low_sum + s_low_sum);

        if (std::isnan(result) || result == 0.0) 
        {
            // either we have split infinite numbers or some coefficients were NaNs or signed zeros, // just rely on the naive implementation and let IEEE754 handle this
            // we do this for zeros too as we want to preserve the sign of zero (see issue #76)
            result = a[0] * b[0];
            for (int i{ 1 }; i < len; ++i) 
            {
                result += a[i] * b[i];
            }
        }

        return result;
    }

    /**
     * Compute a linear combination accurately.
     * <p>
     * This method computes a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> to high accuracy. It does
     * so by using specific multiplication and addition algorithms to
     * preserve accuracy and reduce cancellation effects. It is based
     * on the 2005 paper <a
     * href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.2.1547">
     * Accurate Sum and Dot Product</a> by Takeshi Ogita, * Siegfried M. Rump, and Shin'ichi Oishi published in SIAM J. Sci. Comput.
     * </p>
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub>
     * @see #linear_combination(double, double, double, double, double, double)
     * @see #linear_combination(double, double, double, double, double, double, double, double)
     */
    static double linear_combination(const double& a1, const double b1, const double& a2, const double b2) 
    {
        // the code below is split in many additions/subtractions that may
        // appear redundant. However, they should NOT be simplified, as they
        // use IEEE754 floating point arithmetic rounding properties.
        // The variable naming conventions are that xyz_high contains the most significant
        // bits of xyz and xyz_low contains its least significant bits. So theoretically
        // xyz is the sum xyz_high + xyz_low, but in many cases below, this sum cannot
        // be represented in only one double precision number so we preserve two numbers
        // to hold it as long as we can, combining the high and low order bits together
        // only at the end, after cancellation may have occurred on high order bits

        // split a1 and b1 as one 26 bits number and one 27 bits number
        const double& a1_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(a1) & ((-1L) << 27));
        const double& a1_low      = a1 - a1_high;
        const double b1_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(b1) & ((-1L) << 27));
        const double b1_low      = b1 - b1_high;

        // accurate multiplication a1 * b1
        const double prod1_high  = a1 * b1;
        const double prod1_low   = a1_low * b1_low - (((prod1_high - a1_high * b1_high) - a1_low * b1_high) - a1_high * b1_low);

        // split a2 and b2 as one 26 bits number and one 27 bits number
        const double& a2_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(a2) & ((-1L) << 27));
        const double& a2_low      = a2 - a2_high;
        const double b2_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(b2) & ((-1L) << 27));
        const double b2_low      = b2 - b2_high;

        // accurate multiplication a2 * b2
        const double prod2_high  = a2 * b2;
        const double prod2_low   = a2_low * b2_low - (((prod2_high - a2_high * b2_high) - a2_low * b2_high) - a2_high * b2_low);

        // accurate addition a1 * b1 + a2 * b2
        const double s12_high    = prod1_high + prod2_high;
        const double s12_prime   = s12_high - prod2_high;
        const double s12_low     = (prod2_high - (s12_high - s12_prime)) + (prod1_high - s12_prime);

        // const rounding, s12 may have suffered many cancellations, we try
        // to recover some bits from the extra words we have saved up to now
        double result = s12_high + (prod1_low + prod2_low + s12_low);

        if (std::isnan(result) || result == 0.0) 
        {
            // either we have split infinite numbers or some coefficients were NaNs or signed zeros, // just rely on the naive implementation and let IEEE754 handle this
            // we do this for zeros too as we want to preserve the sign of zero (see issue #76)
            result = a1 * b1 + a2 * b2;
        }

        return result;
    }

    /**
     * Compute a linear combination accurately.
     * <p>
     * This method computes a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub>
     * to high accuracy. It does so by using specific multiplication and
     * addition algorithms to preserve accuracy and reduce cancellation effects.
     * It is based on the 2005 paper <a
     * href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.2.1547">
     * Accurate Sum and Dot Product</a> by Takeshi Ogita, * Siegfried M. Rump, and Shin'ichi Oishi published in SIAM J. Sci. Comput.
     * </p>
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub>
     * @see #linear_combination(double, double, double, double)
     * @see #linear_combination(double, double, double, double, double, double, double, double)
     */
    static double linear_combination(const double& a1, const double b1, const double& a2, const double b2, const double& a3, const double b3) 
    {
        // the code below is split in many additions/subtractions that may
        // appear redundant. However, they should NOT be simplified, as they
        // do use IEEE754 floating point arithmetic rounding properties.
        // The variables naming conventions are that xyz_high contains the most significant
        // bits of xyz and xyz_low contains its least significant bits. So theoretically
        // xyz is the sum xyz_high + xyz_low, but in many cases below, this sum cannot
        // be represented in only one double precision number so we preserve two numbers
        // to hold it as long as we can, combining the high and low order bits together
        // only at the end, after cancellation may have occurred on high order bits

        // split a1 and b1 as one 26 bits number and one 27 bits number
        const double& a1_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(a1) & ((-1L) << 27));
        const double& a1_low      = a1 - a1_high;
        const double b1_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(b1) & ((-1L) << 27));
        const double b1_low      = b1 - b1_high;

        // accurate multiplication a1 * b1
        const double prod1_high  = a1 * b1;
        const double prod1_low   = a1_low * b1_low - (((prod1_high - a1_high * b1_high) - a1_low * b1_high) - a1_high * b1_low);

        // split a2 and b2 as one 26 bits number and one 27 bits number
        const double& a2_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(a2) & ((-1L) << 27));
        const double& a2_low      = a2 - a2_high;
        const double b2_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(b2) & ((-1L) << 27));
        const double b2_low      = b2 - b2_high;

        // accurate multiplication a2 * b2
        const double prod2_high  = a2 * b2;
        const double prod2_low   = a2_low * b2_low - (((prod2_high - a2_high * b2_high) - a2_low * b2_high) - a2_high * b2_low);

        // split a3 and b3 as one 26 bits number and one 27 bits number
        const double& a3_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(a3) & ((-1L) << 27));
        const double& a3_low      = a3 - a3_high;
        const double b3_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(b3) & ((-1L) << 27));
        const double b3_low      = b3 - b3_high;

        // accurate multiplication a3 * b3
        const double prod3_high  = a3 * b3;
        const double prod3_low   = a3_low * b3_low - (((prod3_high - a3_high * b3_high) - a3_low * b3_high) - a3_high * b3_low);

        // accurate addition a1 * b1 + a2 * b2
        const double s12_high    = prod1_high + prod2_high;
        const double s12_prime   = s12_high - prod2_high;
        const double s12_low     = (prod2_high - (s12_high - s12_prime)) + (prod1_high - s12_prime);

        // accurate addition a1 * b1 + a2 * b2 + a3 * b3
        const double s123_high   = s12_high + prod3_high;
        const double s123_prime  = s123_high - prod3_high;
        const double s123_low    = (prod3_high - (s123_high - s123_prime)) + (s12_high - s123_prime);

        // const rounding, s123 may have suffered many cancellations, we try
        // to recover some bits from the extra words we have saved up to now
        double result = s123_high + (prod1_low + prod2_low + prod3_low + s12_low + s123_low);

        if (std::isnan(result) || result == 0.0) 
        {
            // either we have split infinite numbers or some coefficients were NaNs or signed zeros, // just rely on the naive implementation and let IEEE754 handle this
            // we do this for zeros too as we want to preserve the sign of zero (see issue #76)
            result = a1 * b1 + a2 * b2 + a3 * b3;
        }

        return result;
    }

    /**
     * Compute a linear combination accurately.
     * <p>
     * This method computes a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub> +
     * a<sub>4</sub>&times;b<sub>4</sub>
     * to high accuracy. It does so by using specific multiplication and
     * addition algorithms to preserve accuracy and reduce cancellation effects.
     * It is based on the 2005 paper <a
     * href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.2.1547">
     * Accurate Sum and Dot Product</a> by Takeshi Ogita, * Siegfried M. Rump, and Shin'ichi Oishi published in SIAM J. Sci. Comput.
     * </p>
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @param a4 first factor of the third term
     * @param b4 second factor of the third term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub> +
     * a<sub>4</sub>&times;b<sub>4</sub>
     * @see #linear_combination(double, double, double, double)
     * @see #linear_combination(double, double, double, double, double, double)
     */
    static double linear_combination(const double& a1, const double b1, const double& a2, const double b2, const double& a3, const double b3, const double& a4, const double b4) 
    {

        // the code below is split in many additions/subtractions that may
        // appear redundant. However, they should NOT be simplified, as they
        // do use IEEE754 floating point arithmetic rounding properties.
        // The variables naming conventions are that xyz_high contains the most significant
        // bits of xyz and xyz_low contains its least significant bits. So theoretically
        // xyz is the sum xyz_high + xyz_low, but in many cases below, this sum cannot
        // be represented in only one double precision number so we preserve two numbers
        // to hold it as long as we can, combining the high and low order bits together
        // only at the end, after cancellation may have occurred on high order bits

        // split a1 and b1 as one 26 bits number and one 27 bits number
        const double& a1_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(a1) & ((-1L) << 27));
        const double& a1_low      = a1 - a1_high;
        const double b1_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(b1) & ((-1L) << 27));
        const double b1_low      = b1 - b1_high;

        // accurate multiplication a1 * b1
        const double prod1_high  = a1 * b1;
        const double prod1_low   = a1_low * b1_low - (((prod1_high - a1_high * b1_high) - a1_low * b1_high) - a1_high * b1_low);

        // split a2 and b2 as one 26 bits number and one 27 bits number
        const double& a2_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(a2) & ((-1L) << 27));
        const double& a2_low      = a2 - a2_high;
        const double b2_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(b2) & ((-1L) << 27));
        const double b2_low      = b2 - b2_high;

        // accurate multiplication a2 * b2
        const double prod2_high  = a2 * b2;
        const double prod2_low   = a2_low * b2_low - (((prod2_high - a2_high * b2_high) - a2_low * b2_high) - a2_high * b2_low);

        // split a3 and b3 as one 26 bits number and one 27 bits number
        const double& a3_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(a3) & ((-1L) << 27));
        const double& a3_low      = a3 - a3_high;
        const double b3_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(b3) & ((-1L) << 27));
        const double b3_low      = b3 - b3_high;

        // accurate multiplication a3 * b3
        const double prod3_high  = a3 * b3;
        const double prod3_low   = a3_low * b3_low - (((prod3_high - a3_high * b3_high) - a3_low * b3_high) - a3_high * b3_low);

        // split a4 and b4 as one 26 bits number and one 27 bits number
        const double& a4_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(a4) & ((-1L) << 27));
        const double& a4_low      = a4 - a4_high;
        const double b4_high     = Double.long_bits_to_double(Double.double_to_raw_long_bits(b4) & ((-1L) << 27));
        const double b4_low      = b4 - b4_high;

        // accurate multiplication a4 * b4
        const double prod4_high  = a4 * b4;
        const double prod4_low   = a4_low * b4_low - (((prod4_high - a4_high * b4_high) - a4_low * b4_high) - a4_high * b4_low);

        // accurate addition a1 * b1 + a2 * b2
        const double s12_high    = prod1_high + prod2_high;
        const double s12_prime   = s12_high - prod2_high;
        const double s12_low     = (prod2_high - (s12_high - s12_prime)) + (prod1_high - s12_prime);

        // accurate addition a1 * b1 + a2 * b2 + a3 * b3
        const double s123_high   = s12_high + prod3_high;
        const double s123_prime  = s123_high - prod3_high;
        const double s123_low    = (prod3_high - (s123_high - s123_prime)) + (s12_high - s123_prime);

        // accurate addition a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4
        const double s1234High  = s123_high + prod4_high;
        const double s1234_prime = s1234High - prod4_high;
        const double s1234_low   = (prod4_high - (s1234High - s1234_prime)) + (s123_high - s1234_prime);

        // const rounding, s1234 may have suffered many cancellations, we try
        // to recover some bits from the extra words we have saved up to now
        double result = s1234High + (prod1_low + prod2_low + prod3_low + prod4_low + s12_low + s123_low + s1234_low);

        if (std::isnan(result) || result == 0.0) 
        {
            // either we have split infinite numbers or some coefficients were NaNs or signed zeros, // just rely on the naive implementation and let IEEE754 handle this
            // we do this for zeros too as we want to preserve the sign of zero (see issue #76)
            result = a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4;
        }

        return result;
    }

    /**
     * Returns true iff both arguments are NULL or have same dimensions and all
     * their elements are equal as defined by
     * {@link Precision#equals(float,float)}.
     *
     * @param x first array
     * @param y second array
     * @return true if the values are both NULL or have same dimension
     * and equal elements.
     */
    static bool equals(std::vector<float> x, std::vector<float> y) 
    {
        if ((x == NULL) || (y == NULL)) 
        {
            return !((x == NULL) ^ (y == NULL));
        }
        if (x.size() != y.size()) 
        {
            return false;
        }
        for (int i{}; i < x.size(); ++i) 
        {
            if (!Precision.equals(x[i], y[i])) 
            {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns true iff both arguments are NULL or have same dimensions and all
     * their elements are equal as defined by
     * {@link Precision#equals_including_nan(double,double) this method}.
     *
     * @param x first array
     * @param y second array
     * @return true if the values are both NULL or have same dimension and
     * equal elements
     */
    static bool equals_including_nan(std::vector<float> x, std::vector<float> y) 
    {
        if ((x == NULL) || (y == NULL)) 
        {
            return !((x == NULL) ^ (y == NULL));
        }
        if (x.size() != y.size()) 
        {
            return false;
        }
        for (int i{}; i < x.size(); ++i) 
        {
            if (!Precision.equals_including_nan(x[i], y[i])) 
            {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns {@code true} iff both arguments are {@code NULL} or have same
     * dimensions and all their elements are equal as defined by
     * {@link Precision#equals(double,double)}.
     *
     * @param x First array.
     * @param y Second array.
     * @return {@code true} if the values are both {@code NULL} or have same
     * dimension and equal elements.
     */
    static bool equals(std::vector<double> x, std::vector<double> y) 
    {
        if ((x == NULL) || (y == NULL)) 
        {
            return !((x == NULL) ^ (y == NULL));
        }
        if (x.size() != y.size()) 
        {
            return false;
        }
        for (int i{}; i < x.size(); ++i) 
        {
            if (!Precision.equals(x[i], y[i])) 
            {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns {@code true} iff both arguments are {@code NULL} or have same
     * dimensions and all their elements are equal as defined by
     * {@link Precision#equals_including_nan(double,double) this method}.
     *
     * @param x First array.
     * @param y Second array.
     * @return {@code true} if the values are both {@code NULL} or have same
     * dimension and equal elements.
     */
    static bool equals_including_nan(std::vector<double> x, std::vector<double> y) 
    {
        if ((x == NULL) || (y == NULL)) 
        {
            return !((x == NULL) ^ (y == NULL));
        }
        if (x.size() != y.size()) 
        {
            return false;
        }
        for (int i{}; i < x.size(); ++i) 
        {
            if (!Precision.equals_including_nan(x[i], y[i])) 
            {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns {@code true} if both arguments are {@code NULL} or have same
     * dimensions and all their elements are equals.
     *
     * @param x First array.
     * @param y Second array.
     * @return {@code true} if the values are both {@code NULL} or have same
     * dimension and equal elements.
     */
    static bool equals(std::vector<long> x, std::vector<long> y) 
    {
        if ((x == NULL) || (y == NULL)) 
        {
            return !((x == NULL) ^ (y == NULL));
        }
        if (x.size() != y.size()) 
        {
            return false;
        }
        for (int i{}; i < x.size(); ++i) 
        {
            if (x[i] != y[i]) 
            {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns {@code true} if both arguments are {@code NULL} or have same
     * dimensions and all their elements are equals.
     *
     * @param x First array.
     * @param y Second array.
     * @return {@code true} if the values are both {@code NULL} or have same
     * dimension and equal elements.
     */
    static bool equals(const std::vector<int>& x, std::vector<int> y) 
    {
        if ((x == NULL) || (y == NULL)) 
        {
            return !((x == NULL) ^ (y == NULL));
        }
        if (x.size() != y.size()) 
        {
            return false;
        }
        for (int i{}; i < x.size(); ++i) 
        {
            if (x[i] != y[i]) 
            {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns {@code true} if both arguments are {@code NULL} or have same
     * dimensions and all their elements are equals.
     *
     * @param x First array.
     * @param y Second array.
     * @return {@code true} if the values are both {@code NULL} or have same
     * dimension and equal elements.
     */
    static bool equals(std::vector<std::byte>x, std::vector<std::byte>y) 
    {
        if ((x == NULL) || (y == NULL)) 
        {
            return !((x == NULL) ^ (y == NULL));
        }
        if (x.size() != y.size()) 
        {
            return false;
        }
        for (int i{}; i < x.size(); ++i) 
        {
            if (x[i] != y[i]) 
            {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns {@code true} if both arguments are {@code NULL} or have same
     * dimensions and all their elements are equals.
     *
     * @param x First array.
     * @param y Second array.
     * @return {@code true} if the values are both {@code NULL} or have same
     * dimension and equal elements.
     */
    static bool equals(std::vector<short>& x, std::vector<short>& y) 
    {
        if ((x == NULL) || (y == NULL)) 
        {
            return !((x == NULL) ^ (y == NULL));
        }
        if (x.size() != y.size()) 
        {
            return false;
        }
        for (int i{}; i < x.size(); ++i) 
        {
            if (x[i] != y[i]) 
            {
                return false;
            }
        }
        return true;
    }

    /**
     * Normalizes an array to make it sum to a specified value.
     * Returns the result of the transformation
     * <pre>
     *    x \xe2\x86\xa6 x * normalized_sum / sum
     * </pre>
     * applied to each non-NaN element x of the input array, where sum is the
     * sum of the non-NaN entries in the input array.
     * <p>
     * Throws Illegal_Argument_Exception if {@code normalized_sum} is infinite
     * or NaN and Arithmetic_Exception if the input array contains any infinite elements
     * or sums to 0.
     * <p>
     * Ignores (i.e., copies unchanged to the output array) NaNs in the input array.
     * The input array is unchanged by this method.
     *
     * @param values Input array to be normalized
     * @param normalized_sum Target sum for the normalized array
     * @return the normalized array
     * @Math_Runtime_Exception if the input array contains infinite
     * elements or sums to zero
     * @ if the target sum is infinite or {@code NaN}
     */
    static std::vector<double> normalize_array(const std::vector<double>& values, const double& normalized_sum)
    {
        if (std::isinf(normalized_sum)) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NORMALIZE_INFINITE);
        }
        if (std::isnan(normalized_sum)) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NORMALIZE_NAN);
        }
        double sum{};
        const auto len = values.size();
        auto out = std::vector<double>(len);
        for (int i{}; i < len; i++) 
        {
            if (std::isinf(values[i])) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::INFINITE_ARRAY_ELEMENT, values[i], i);
            }
            if (!std::isnan(values[i])) 
            {
                sum += values[i];
            }
        }
        if (sum == 0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ARRAY_SUMS_TO_ZERO);
        }
        for (int i{}; i < len; i++) 
        {
            out[i] = std::isnan(values[i])
                ? std::numeric_limits<double>::quiet_NaN()
                : values[i] * normalized_sum / sum;
        }
        return out;
    }

    /** Build an array of elements.
     * <p>
     * Arrays are filled with {@code field.get_zero()}
     *
     * @param <T> the type of the field elements
     * @param field field to which array elements belong
     * @param length of the array
     * @return a array
     */
    static <T extends Field_Element<T>> std::vector<T> build_array(const Field<T> field, const int length) 
    {
        ////@Suppress_Warnings("unchecked") // OK because field must be correct class
        std::vector<T> array = (std::vector<T>) Array.new_instance(field.get_runtime_class(), length);
        Arrays.fill(array, field.get_zero());
        return array;
    }

    /** Build a double dimension array of elements.
     * <p>
     * Arrays are filled with {@code field.get_zero()}
     *
     * @param <T> the type of the field elements
     * @param field field to which array elements belong
     * @param rows number of rows in the array
     * @param columns number of columns (may be negative to build partial
     * arrays in the same way {@code Field[rows][]} works)
     * @return a array
     */
    //@Suppress_Warnings("unchecked")
    static <T extends Field_Element<T>> std::vector<std::vector<T>> build_array(const Field<T>& field, const int& rows, const int& columns) 
    {
        std::vector<std::vector<T>> arr_2d;
        if (columns < 0) 
        {
            std::vector<T> dummy_row = build_array(field, 0);
            arr_2d = (std::vector<std::vector<T>>) Array.new_instance(dummy_row.get_class(), rows);
        }
        else 
        {
            arr = (std::vector<std::vector<T>>) Array.new_instance(field.get_runtime_class(), std::vector<int>
            {
                                                  rows, columns
            });
            const auto row_size = arr_2d[0].size();
            for (auto& row : arr_2d)
            {
                row = std::vector<T>(row_size, field.get_zero());
            }
        }
        return arr_2d;
    }

    /** Build a triple dimension array of elements.
     * <p>
     * Arrays are filled with {@code field.get_zero()}
     *
     * @param <T> the type of the field elements
     * @param field field to which array elements belong
     * @param l1 number of elements along first dimension
     * @param l2 number of elements along second dimension
     * @param l3 number of elements along third dimension (may be negative to build partial
     * arrays in the same way {@code Field[l1][l2][]} works)
     * @return a array
     * @since 1.4
     */
    //@Suppress_Warnings("unchecked")
    static <T extends Field_Element<T>> std::vector<std::vector<std::vector<T>>> build_array(const Field<T>& field, const int& l1, const int& l2, const int& l3) 
    {
        std::vector<std::vector<std::vector<T>>> arr_3d;
        if (l3 < 0) 
        {
            std::vector<T> dummy_row = build_array(field, 0);
            arr_3d = (std::vector<std::vector<T>>[]) Array.new_instance(dummy_row.get_class(), l1, l2);
        }
        else 
        {
            arr_3d = (std::vector<std::vector<T>>[]) Array.new_instance(field.get_runtime_class(), std::vector<int>
            {
                l1, l2, l3
            });
            for (int i{}; i < l1; ++i) 
            {
                for (int j{}; j < l2; ++j) 
                {
                    Arrays.fill(arr_3d[i][j], field.get_zero());
                }
            }
        }
        return arr_3d;
    }

    /**
     * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
     * convolution</a> between two sequences.
     * <p>
     * The solution is obtained via straightforward computation of the
     * convolution sum (and not via FFT). Whenever the computation needs
     * an element that would be located at an index outside the input arrays, * the value is assumed to be zero.
     *
     * @param x First sequence.
     * Typically, this sequence will represent an input signal to a system.
     * @param h Second sequence.
     * Typically, this sequence will represent the impulse response of the system.
     * @return the convolution of {@code x} and {@code h}.
     * This array's length will be {@code x.size() + h.size() - 1}.
     * @Null_Argument_Exception if either {@code x} or {@code h} is {@code NULL}.
     * @ if either {@code x} or {@code h} is empty.
     */
    static std::vector<double> convolve(const std::vector<double>& x, const std::vector<double>& h)
    {
        //Math_Utils::check_not_null(x);
        //Math_Utils::check_not_null(h);

        const int x_len = x.size();
        const int h_len = h.size();

        if (x_len == 0 || h_len == 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }

        // initialize the output array
        const int total_length = x_len + h_len - 1;
        const std::vector<double> y = std::vector<double>(total_length];

        // straightforward implementation of the convolution sum
        for (int n{}; n < total_length; n++) 
        {
            double yn = 0;
            int k = std::max(0, n + 1 - x_len);
            int j = n - k;
            while (k < h_len && j >= 0) 
            {
                yn += x[j--] * h[k++];
            }
            y[n] = yn;
        }

        return y;
    }

    /**
     * Specification for indicating that some operation applies
     * before or after a given index.
     */
    enum Position 
    {
        /** Designates the beginning of the array (near index 0). */
        HEAD, /** Designates the end of the array. */
        TAIL
    }

    /**
     * Shuffle the entries of the given array.
     * The {@code start} and {@code pos} parameters select which portion
     * of the array is randomized and which is left untouched.
     *
     * @see #shuffle(std::vector<int>,int,Position,Random_Generator)
     *
     * @param list Array whose entries will be shuffled (in-place).
     * @param start Index at which shuffling begins.
     * @param pos Shuffling is performed for index positions between
     * {@code start} and either the end (if {@link Position#TAIL})
     * or the beginning (if {@link Position#HEAD}) of the array.
     */
    static void shuffle(const std::vector<int>& list, const int& start, const Position& pos) 
    {
        shuffle(list, start, pos, Well19937c());
    }

    /**
     * Shuffle the entries of the given array, using the
     * <a href="http://en.wikipedia.org/wiki/Fisher\xe2\x80\x93Yates_shuffle#The_modern_algorithm">
     * Fisher\xe2\x80\x93Yates</a> algorithm.
     * The {@code start} and {@code pos} parameters select which portion
     * of the array is randomized and which is left untouched.
     *
     * @param list Array whose entries will be shuffled (in-place).
     * @param start Index at which shuffling begins.
     * @param pos Shuffling is performed for index positions between
     * {@code start} and either the end (if {@link Position#TAIL})
     * or the beginning (if {@link Position#HEAD}) of the array.
     * @param rng Random number generator.
     */
    static void shuffle(const std::vector<int>& list, const int& start, const Position& pos, const Random_Generator& rng) 
    {
        switch (pos) 
        {
        case TAIL:
            for (int i = list.size() - 1; i > start; i--) 
            {
                const int target = start + rng.next_int(i - start + 1);
                const int temp = list[target];
                list[target] = list[i];
                list[i] = temp;
            }
            break;

        case HEAD:
            for (int i{}; i < start; i++) 
            {
                const int target = i + rng.next_int(start - i + 1);
                const int temp = list[target];
                list[target] = list[i];
                list[i] = temp;
            }
            break;

        default:
            throw Math_Runtime_Exception.create_internal_error(); // Should never happen.
        }
    }

    /**
     * Shuffle the entries of the given array.
     *
     * @see #shuffle(std::vector<int>,int,Position,Random_Generator)
     *
     * @param list Array whose entries will be shuffled (in-place).
     * @param rng Random number generator.
     */
    static void shuffle(const std::vector<int>& list, const Random_Generator& rng) 
    {
        shuffle(list, 0, Position.TAIL, rng);
    }

    /**
     * Shuffle the entries of the given array.
     *
     * @see #shuffle(std::vector<int>,int,Position,Random_Generator)
     *
     * @param list Array whose entries will be shuffled (in-place).
     */
    static void shuffle(const std::vector<int>& list) 
    {
        shuffle(list, Well19937c());
    }

    /**
     * Returns an array representing the natural number {@code n}.
     *
     * @param n Natural number.
     * @return an array whose entries are the numbers 0, 1, ..., {@code n}-1.
     * If {@code n == 0}, the returned array is empty.
     */
    static std::vector<int> natural(const int& n) 
    {
        return sequence(n, 0, 1);
    }

    /**
     * Returns an array of {@code size} integers starting at {@code start}, * skipping {@code stride} numbers.
     *
     * @param size Natural number.
     * @param start Natural number.
     * @param stride Natural number.
     * @return an array whose entries are the numbers
     * {@code start, start + stride, ..., start + (size - 1) * stride}.
     * If {@code size == 0}, the returned array is empty.
     */
    static std::vector<int> sequence(const int& size, const int& start, const int& stride) 
    {
        auto a = std::vector<int>(size);
        for (auto& ele : a)
        {
            ele = start + i * stride;
        }
        return a;
    }

    /**
     * This method is used
     * to verify that the input parameters designate a subarray of positive length.
     * <p>
     * <ul>
     * <li>returns <code>true</code> iff the parameters designate a subarray of
     * positive length</li>
     * <li><code></code> if the array is NULL or
     * or the indices are invalid</li>
     * <li>returns <code>false</li> if the array is non-null, but
     * <code>length</code> is 0.
     * </ul></p>
     *
     * @param values the input array
     * @param begin index of the first array element to include
     * @param length the number of elements to include
     * @return true if the parameters are valid and designate a subarray of positive length
     * @ if the indices are invalid or the array is NULL
     */
    static bool verify_values(const std::vector<double>& values, const int& begin, const int& length)
    {
        return verify_values(values, begin, length, false);
    }

    /**
     * This method is used
     * to verify that the input parameters designate a subarray of positive length.
     * <p>
     * <ul>
     * <li>returns <code>true</code> iff the parameters designate a subarray of
     * non-negative length</li>
     * <li><code>Illegal_Argument_Exception</code> if the array is NULL or
     * or the indices are invalid</li>
     * <li>returns <code>false</li> if the array is non-null, but
     * <code>length</code> is 0 unless <code>allow_empty</code> is <code>true</code>
     * </ul></p>
     *
     * @param values the input array
     * @param begin index of the first array element to include
     * @param length the number of elements to include
     * @param allow_empty if <code>true</code> then zero length arrays are allowed
     * @return true if the parameters are valid
     * @ if the indices are invalid or the array is NULL
     */
    static bool verify_values(const std::vector<double>& values, const int& begin, const int& length, const bool allow_empty)  
    {

        //Math_Utils::check_not_null(values, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);

        if (begin < 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::START_POSITION, static_cast<int>(begin));
        }

        if (length < 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::size(), static_cast<int>(length));
        }

        if (begin + length > values.size()) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::SUBARRAY_ENDS_AFTER_ARRAY_END, static_cast<int>(begin + length), static_cast<int>(values.size()), true);
        }

        if (length == 0 && !allow_empty) 
        {
            return false;
        }

        return true;
    }

    /**
     * This method is used
     * to verify that the begin and length parameters designate a subarray of positive length
     * and the weights are all non-negative, non-NaN, finite, and not all zero.
     * <p>
     * <ul>
     * <li>returns <code>true</code> iff the parameters designate a subarray of
     * positive length and the weights array contains legitimate values.</li>
     * <li><code>Illegal_Argument_Exception</code> if any of the following are true:
     * <ul><li>the values array is NULL</li>
     *     <li>the weights array is NULL</li>
     *     <li>the weights array does not have the same length as the values array</li>
     *     <li>the weights array contains one or more infinite values</li>
     *     <li>the weights array contains one or more NaN values</li>
     *     <li>the weights array contains negative values</li>
     *     <li>the start and length arguments do not determine a valid array</li></ul>
     * </li>
     * <li>returns <code>false</li> if the array is non-null, but
     * <code>length</code> is 0.
     * </ul>
     *
     * @param values the input array
     * @param weights the weights array
     * @param begin index of the first array element to include
     * @param length the number of elements to include
     * @return true if the parameters are valid and designate a subarray of positive length
     * @ if the indices are invalid or the array is NULL
     */
    static bool verify_values(const std::vector<double>& values, const std::vector<double>& weights, const int& begin, const int& length)  
    {
        return verify_values(values, weights, begin, length, false);
    }

    /**
     * This method is used
     * to verify that the begin and length parameters designate a subarray of positive length
     * and the weights are all non-negative, non-NaN, finite, and not all zero.
     * <p>
     * <ul>
     * <li>returns <code>true</code> iff the parameters designate a subarray of
     * non-negative length and the weights array contains legitimate values.</li>
     * <li><code></code> if any of the following are true:
     * <ul><li>the values array is NULL</li>
     *     <li>the weights array is NULL</li>
     *     <li>the weights array does not have the same length as the values array</li>
     *     <li>the weights array contains one or more infinite values</li>
     *     <li>the weights array contains one or more NaN values</li>
     *     <li>the weights array contains negative values</li>
     *     <li>the start and length arguments do not determine a valid array</li></ul>
     * </li>
     * <li>returns <code>false</li> if the array is non-null, but
     * <code>length</code> is 0 unless <code>allow_empty</code> is <code>true</code>.
     * </ul>
     *
     * @param values the input array.
     * @param weights the weights array.
     * @param begin index of the first array element to include.
     * @param length the number of elements to include.
     * @param allow_empty if {@code true} than allow zero length arrays to pass.
     * @return {@code true} if the parameters are valid.
     * @Null_Argument_Exception if either of the arrays are NULL
     * @ if the array indices are not valid, * the weights array contains NaN, infinite or negative elements, or there
     * are no positive weights.
     */
    static bool verify_values(const std::vector<double>& values, const std::vector<double>& weights, const int& begin, const int& length, const bool allow_empty)  
    {
        //Math_Utils::check_not_null(weights, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
        //Math_Utils::check_not_null(values, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);

        check_equal_length(weights, values);

        bool contains_positive_weight{};
        for (int i{ begin }; i < begin + length; i++)
        {
            const double weight = weights[i];
            if (std::isnan(weight)) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::NAN_ELEMENT_AT_INDEX, static_cast<int>(i));
            }
            if (std::isinfinite(weight)) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::INFINITE_ARRAY_ELEMENT, static_cast<double>(weight), static_cast<int>(i));
            }
            if (weight < 0) 
            {
                throw (hipparchus::exception::Localized_Core_Formats_Type::NEGATIVE_ELEMENT_AT_INDEX, static_cast<int>(i), static_cast<double>(weight));
            }
            if (!contains_positive_weight && weight > 0.0) 
            {
                contains_positive_weight = true;
            }
        }

        if (!contains_positive_weight) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::WEIGHT_AT_LEAST_ONE_NON_ZERO);
        }

        return verify_values(values, begin, length, allow_empty);
    }

    /**
     * Concatenates a sequence of arrays. The return array consists of the
     * entries of the input arrays concatenated in the order they appear in
     * the argument list.  Null arrays cause Null_Pointer_Exceptions; zero
     * length arrays are allowed (contributing nothing to the output array).
     *
     * @param x list of std::vector<double> arrays to concatenate
     * @return a array consisting of the entries of the argument arrays
     * @Null_Pointer_Exception if any of the arrays are NULL
     */
    static std::vector<double> concatenate(const std::vector<double>... x) 
    {
        int combined_length{};
        for (std::vector<double> a : x) 
        {
            combined_length += a.size();
        }
        int offset{};
        const std::vector<double> combined = std::vector<double>(combined_length];
        for (int i{}; i < x.size(); i++) 
        {
            const int cur_length = x[i].size();
            System.arraycopy(x[i], 0, combined, offset, cur_length);
            offset += cur_length;
        }
        return combined;
    }

    /**
     * Returns an array consisting of the unique values in {@code data}.
     * The return array is sorted in descending order.  Empty arrays
     * are allowed, but NULL arrays result in Null_Pointer_Exception.
     * Infinities are allowed.  NaN values are allowed with maximum
     * sort order - i.e., if there are NaN values in {@code data}, * {@codeNAN} will be the first element of the output array, * even if the array also contains {@code INFINITY}.
     *
     * @param data array to scan
     * @return descending list of values included in the input array
     * @Null_Pointer_Exception if data is NULL
     */
    static std::vector<double> unique(std::vector<double>& data) 
    {
        Tree_Set<Double> values = Tree_Set<>();
        for (int i{}; i < data.size(); i++) 
        {
            values.add(data[i]);
        }
        const int count = values.size();
        auto out = std::vector<double>(count];
        Iterator<Double> iterator = values.descending_iterator();
        int i{};
        while (iterator.has_next()) 
        {
            out[i++] = iterator.next();
        }
        return out;
    }
};