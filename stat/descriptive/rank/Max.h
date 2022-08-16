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

#include <limits>
#include <cmath>
#include <vector>
#include "../../descriptive/AbstractStorelessUnivariateStatistic.h"
#include "../../descriptive/AggregatableStatistic.hpp"
#include "../../../core/util/MathArrays.h"

/**
 * Returns the maximum of the available values.
 * <p>
 * <ul>
 * <li>The result is <code>NaN</code> iff all values are <code>NaN</code>
 * (i.e. <code>NaN</code> values have no impact on the value of the statistic).</li>
 * <li>If any of the values equals <code>INFINITY</code>, * the result is <code>INFINITY.</code></li>
 * </ul>
* <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
class Max : public Abstract_Storeless_Univariate_Statistic //, public Aggregatable_Statistic<Max>
{
private:


    /** Number of values that have been added */
    long my_n;

    /** Current value of the statistic */
    double my_value;

public:
    /**
     * Create a Max instance.
     */
    Max() : my_n{ 0 }, my_value{ std::numeric_limits<double>::quiet_NaN() } {};

    /**
     * Copy constructor, creates a {@code Max} identical
     * to the {@code original}.
     *
     * @param original the {@code Max} instance to copy
     * @ if original is NULL
     */
    Max(const Max& original) : my_n{ original.get_n() }, my_value{ original.get_result()}
    {
        //Math_Utils::check_not_null(original);
    }

    /** {@inherit_doc} */
    //override
    void increment(const double& d) 
    {
        if (d > my_value || std::isnan(my_value)) 
        {
            my_value = d;
        }
        my_n++;
    }

    /** {@inherit_doc} */
    void clear() override
    {
        my_value = std::numeric_limits<double>::quiet_NaN();
        my_n = 0;
    }

    /** {@inherit_doc} */
    //override
    double get_result() const
    {
        return my_value;
    }

    /** {@inherit_doc} */
    //override
    long get_n() const
    {
        return my_n;
    }

    /** {@inherit_doc} */
    //override
    void aggregate(Max other) 
    {
        //Math_Utils::check_not_null(other);
        if (other.get_n() > 0)
        {
            if (other.get_result() > my_value || std::isnan(my_value))
            {
                my_value = other.get_result();
            }
            my_n += other.get_n();
        }
    }

    /**
     * Returns the maximum of the entries in the specified portion of
     * the input array, or <code>Double.NaN</code> if the designated subarray
     * is empty.
     * <p>
     * Throws <code></code> if the array is NULL or
     * the array index parameters are not valid.
     * <p>
     * <ul>
     * <li>The result is <code>NaN</code> iff all values are <code>NaN</code>
     * (i.e. <code>NaN</code> values have no impact on the value of the statistic).</li>
     * <li>If any of the values equals <code>INFINITY</code>, * the result is <code>INFINITY.</code></li>
     * </ul>
     *
     * @param values the input array
     * @param begin index of the first array element to include
     * @param length the number of elements to include
     * @return the maximum of the values orNAN if length = 0
     * @ if the array is NULL or the array index
     *  parameters are not valid
     */
    //override
    double evaluate(const std::vector<double>& values, const int& begin, const int& length)
    {
        auto max = std::numeric_limits<double>::quiet_NaN();
        if (Math_Arrays::verify_values(values, begin, length)) 
        {
            max = values[begin];
            for (int i{ begin }; i < begin + length; i++)
            {
                if (!std::isnan(values[i])) 
                {
                    max = (max > values[i]) ? max : values[i];
                }
            }
        }
        return max;
    }

    /** {@inherit_doc} */
    //override
    //Max copy() 
    //{
    //    return Max(this);
    //}
};