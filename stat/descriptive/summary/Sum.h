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
//package org.hipparchus.stat.descriptive.summary;

//import java.io.Serializable;

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.stat.descriptive.Abstract_Storeless_Univariate_Statistic;
//import org.hipparchus.stat.descriptive.Aggregatable_Statistic;
//import org.hipparchus.stat.descriptive.Weighted_Evaluation;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <limits>
#include <cmath>
#include <vector>
#include "../../descriptive/AbstractStorelessUnivariateStatistic.h"
#include "../WeightedEvaluation.h"
#include "../../descriptive/AggregatableStatistic.hpp"
#include "../../../core/util/MathArrays.h"


/**
  * Returns the sum of the available values.
 * <p>
 * If there are no values in the dataset, then 0 is returned.
 * If any of the values are
 * <code>NaN</code>, then <code>NaN</code> is returned.
 * <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
class Sum : public Abstract_Storeless_Univariate_Statistic, public Aggregatable_Statistic<Sum>, public Weighted_Evaluation
{
private:
    /** The number of values that have been added */
    long my_n;

    /** The currently running sum */
    double my_value;

public:
    /**
     * Create a Sum instance.
     */
    Sum() : my_n{}, my_value{} {};

    /**
     * Copy constructor, creates a {@code Sum} identical
     * to the {@code original}.
     *
     * @param original the {@code Sum} instance to copy
     * @Null_Argument_Exception if original is NULL
     */
    Sum(const Sum& original) : my_n{ original.get_n() }, my_value{ original.get_result() }
    {
        //Math_Utils::check_not_null(original);
    };

    /** {@inherit_doc} */
    //override
     void increment(const double& d) 
    {
        my_value += d;
        my_n++;
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
    void clear() 
    {
        my_value = 0;
        my_n = 0;
    }

    /** {@inherit_doc} */
    //override
    void aggregate(const Sum& other) 
    {
        //Math_Utils::check_not_null(other);
        if (other.get_n() > 0)
        {
            my_n += other.get_n();
            my_value += other.get_result();
        }
    }

    /**
     * The sum of the entries in the specified portion of the input array, * or 0 if the designated subarray is empty.
     *
     * @param values the input array
     * @param begin index of the first array element to include
     * @param length the number of elements to include
     * @return the sum of the values or 0 if length = 0
     * @ if the array is NULL or the array index
     *  parameters are not valid
     */
    //override
    double evaluate(const std::vector<double>& values, const int& begin, const int& length)
         
        {

        double sum = std::numeric_limits<double>::quiet_NaN();
        if (Math_Arrays::verify_values(values, begin, length, true)) 
        {
            sum = 0.0;
            for (int i{ begin }; i < begin + length; i++) 
            {
                sum += values[i];
            }
        }
        return sum;
    }

    /**
     * The weighted sum of the entries in the specified portion of
     * the input array, or 0 if the designated subarray
     * is empty.
     * <p>
     * Throws <code></code> if any of the following are true:
     * <ul><li>the values array is NULL</li>
     *     <li>the weights array is NULL</li>
     *     <li>the weights array does not have the same length as the values array</li>
     *     <li>the weights array contains one or more infinite values</li>
     *     <li>the weights array contains one or more NaN values</li>
     *     <li>the weights array contains negative values</li>
     *     <li>the start and length arguments do not determine a valid array</li>
     * </ul></p>
     * <p>
     * Uses the formula, <pre>
     *    weighted sum = &Sigma;(values[i] * weights[i])
     * </pre></p>
     *
     * @param values the input array
     * @param weights the weights array
     * @param begin index of the first array element to include
     * @param length the number of elements to include
     * @return the sum of the values or 0 if length = 0
     * @ if the parameters are not valid
     */
    //override
    double evaluate(const std::vector<double>& values, const std::vector<double>& weights, const int& begin, const int& length)  
    {
        double sum = std::numeric_limits<double>::quiet_NaN();
        if (Math_Arrays::verify_values(values, weights, begin, length, true)) 
        {
            sum = 0.0;
            for (int i{ begin }; i < begin + length; i++) 
            {
                sum += values[i] * weights[i];
            }
        }
        return sum;
    }

    /** {@inherit_doc} */
    //override
    //Sum copy() 
    //{
    //    return Sum(this);
    //}

};