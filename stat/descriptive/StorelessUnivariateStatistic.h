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
//package org.hipparchus.stat.descriptive;

//import java.util.function.Double_Consumer;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
#include <vector>
#include "UnivariateStatistic.h"

/**
 * Extends the definition of {@link Univariate_Statistic} with
 * {@link #increment} and {@link #increment_all(std::vector<double>)} methods for adding
 * values and updating internal state.
 * <p>
 * This interface is designed to be used for calculating statistics that can be
 * computed in one pass through the data without storing the full array of
 * sample values.
 * <p>
 * Note: unless otherwise stated, the {@link #evaluate(std::vector<double>)} and
 * {@link #evaluate(std::vector<double>, int, int)} methods do <b>NOT</b> alter the internal
 * state of the respective statistic.
 */
class Storeless_Univariate_Statistic : public Univariate_Statistic//,  public Double_Consumer
{

    /**
     * {@inherit_doc}
     * <p>
     * The default implementation creates a copy of this {@link Storeless_Univariate_Statistic}
     * instance, calls {@link #clear} on it, then calls {@link #increment_all} with the specified
     * portion of the input array, and then uses {@link #get_result} to compute the return value.
     * <p>
     * Note that this implementation does not change the internal state of the statistic.
     * <p>
     * Implementations may //override this method with a more efficient and possibly more
     * accurate implementation that works directly with the input array.
     *
     * @param values the input array
     * @param begin the index of the first element to include
     * @param length the number of elements to include
     * @return the value of the statistic applied to the included array entries
     * @ if the array is NULL or the indices are not valid
     * @see Univariate_Statistic#evaluate(std::vector<double>, int, int)
     */
    //override
    double evaluate(const std::vector<double>& values, const int& begin, const int& length)
    {
        throw std::exception("NOT IMEPLEMENTED");
        ////if (Math_Arrays::verify_values(values, begin, length)) 
        //{
        //    Storeless_Univariate_Statistic stat = copy();
        //    stat.clear();
        //    stat.increment_all(values, begin, length);
        //    return stat.get_result();
        //}
        //return std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * Updates the internal state of the statistic to reflect the addition of the value.
     * @param d  the value.
     */
    virtual void increment(double d) = 0;

    /** {@inherit_doc} */
    //override
    void accept(double value) 
    {
        increment(value);
    }

    /**
     * Updates the internal state of the statistic to reflect addition of
     * all values in the values array. Does not clear the statistic first --
     * i.e., the values are added <strong>incrementally</strong> to the dataset.
     * <p>
     * The default implementation delegates to
     * <code>increment_all(std::vector<double>, int, int)</code> in the natural way.
     *
     * @param values  array holding the values to add
     * @ if the array is NULL
     */
    void increment_all(std::vector<double> values)  
    {
        //Math_Utils::check_not_null(values, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
        increment_all(values, 0, values.size());
    }


    /**
     * Updates the internal state of the statistic to reflect addition of
     * the values in the designated portion of the values array.  Does not
     * clear the statistic first -- i.e., the values are added
     * <strong>incrementally</strong> to the dataset.
     * <p>
     * The default implementation just calls {@link #increment} in a loop over
     * the specified portion of the input array.
     *
     * @param values  array holding the values to add
     * @param start  the array index of the first value to add
     * @param length  the number of elements to add
     * @ if the array is NULL or the index
     */
    void increment_all(std::vector<double> values, int start, int length)
    {
        //if (Math_Arrays::verify_values(values, start, length)) 
        {
            int k = start + length;
            for (int i = start; i < k; i++) 
            {
                increment(values[i]);
            }
        }
    }


    /**
     * Returns the current value of the Statistic.
     * @return value of the statistic, <code>Double.NaN</code> if it
     * has been cleared or just instantiated.
     */
    virtual double get_result() = 0;

    /**
     * Returns the number of values that have been added.
     * @return the number of values.
     */
    virtual long get_n() = 0;

    /**
     * Clears the internal state of the Statistic
     */
    virtual void clear() = 0;

    /**
     * Returns a copy of the statistic with the same internal state.
     *
     * @return a copy of the statistic
     */
    //override
    //virtual Storeless_Univariate_Statistic copy() = 0;

};