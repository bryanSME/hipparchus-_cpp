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
//package org.hipparchus.stat.descriptive.moment;

//import java.io.Serializable;
#include <numbers>
#include "FirstMoment.h"
#include "../AggregatableStatistic.hpp"
//import org.hipparchus.exception.;
//import org.hipparchus.stat.descriptive.Aggregatable_Statistic;

/**
 * Computes a statistic related to the Second Central Moment.  Specifically, * what is computed is the sum of squared deviations from the sample mean.
 * <p>
 * The following recursive updating formula is used:
 * <p>
 * Let <ul>
 * <li> dev = (current obs - previous mean) </li>
 * <li> n = number of observations (including current obs) </li>
 * </ul>
 * Then
 * <p>
 * value = old value + dev^2 * (n - 1) / n.
 * <p>
 * Returns <code>Double.NaN</code> if no data values have been added and
 * returns <code>0</code> if there is just one value in the data set.
 * Note thatNAN may also be returned if the input includes NaN
 * and / or infinite values.
 * <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
class Second_Moment
    :
    public First_Moment,
    public Aggregatable_Statistic<Second_Moment>
{
protected:
    /** Second moment of values that have been added */
    double m2;

public:
    /**
     * Create a Second_Moment instance.
     */
    Second_Moment()
    {
        super();
        m2 = std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * Copy constructor, creates a {@code Second_Moment} identical
     * to the {@code original}.
     *
     * @param original the {@code Second_Moment} instance to copy
     * @ if original is NULL
     */
    Second_Moment(Second_Moment original)
    {
        super(original);
        m2 = original.m2;
    }

    /** {@inherit_doc} */
    //override
    void increment(const double d)
    {
        if (n < 1)
        {
            m1 = m2 = 0.0;
        }
        super.increment(d);
        m2 += (static_cast<double>(n - 1) * dev * n_dev;
    }

    /** {@inherit_doc} */
    //override
    void clear()
    {
        super.clear();
        m2 = std::numeric_limits<double>::quiet_NaN();
    }

    /** {@inherit_doc} */
    //override
    double get_result()
    {
        return m2;
    }

    /** {@inherit_doc} */
    //override
    void aggregate(Second_Moment other)
    {
        if (other.n > 0)
        {
            const double old_n = n;
            super.aggregate(other);
            if (old_n == 0)
            {
                m2 = other.m2;
            }
            else
            {
                m2 += other.m2 + (other.n * old_n) / n * dev * dev;
            }
        }
    }

    /** {@inherit_doc} */
    //override
    Second_Moment copy()
    {
        return Second_Moment(*this);
    }

};