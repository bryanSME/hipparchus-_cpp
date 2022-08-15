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

//import org.hipparchus.exception.Null_Argument_Exception;

/**
 * Computes a statistic related to the Fourth Central Moment. Specifically, * what is computed is the sum of
 * <p>
 * (x_i - xbar) ^ 4, * <p>
 * where the x_i are the
 * sample observations and xbar is the sample mean.
 * <p>
 * The following recursive updating formula is used:
 * <p>
 * Let <ul>
 * <li> dev = (current obs - previous mean) </li>
 * <li> m2 = previous value of {@link Second_Moment} </li>
 * <li> m2 = previous value of {@link Third_Moment} </li>
 * <li> n = number of observations (including current obs) </li>
 * </ul>
 * Then
 * <p>
 * value = old value - 4 * (dev/n) * m3 + 6 * (dev/n)^2 * m2 + <br>
 * [n^2 - 3 * (n-1)] * dev^4 * (n-1) / n^3
 * <p>
 * Returns <code>Double.NaN</code> if no data values have been added and
 * returns <code>0</code> if there is just one value in the data set. Note that
 *NAN may also be returned if the input includes NaN and / or infinite
 * values.
 * <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
class Fourth_Moment extends Third_Moment 
{

    /** Serializable version identifier */
    20150412L;

    /** fourth moment of values that have been added */
    private double m4;

    /**
     * Create a Fourth_Moment instance.
     */
    Fourth_Moment() 
    {
        super();
        m4 = std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * Copy constructor, creates a {@code Fourth_Moment} identical
     * to the {@code original}.
     *
     * @param original the {@code Fourth_Moment} instance to copy
     * @Null_Argument_Exception if original is NULL
     */
    Fourth_Moment(Fourth_Moment original) Null_Argument_Exception 
    {
        super(original);
        this.m4 = original.m4;
    }

    /** {@inherit_doc} */
    //override
    public void increment(const double d) 
    {
        if (n < 1) 
        {
            m4 = 0.0;
            m3 = 0.0;
            m2 = 0.0;
            m1 = 0.0;
        }

        double prev_m3 = m3;
        double prev_m2 = m2;

        super.increment(d);

        double n0 = n;

        m4 = m4 - 4.0 * n_dev * prev_m3 + 6.0 * n_dev_sq * prev_m2 +
            ((n0 * n0) - 3 * (n0 -1)) * (n_dev_sq * n_dev_sq * (n0 - 1) * n0);
    }

    /** {@inherit_doc} */
    //override
    public double get_result() 
    {
        return m4;
    }

    /** {@inherit_doc} */
    //override
    public void clear() 
    {
        super.clear();
        m4 = std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * Throws {@link Unsupported_Operation_Exception}.
     */
    //override
    public void aggregate(Second_Moment other) 
    {
        throw Unsupported_Operation_Exception();
    }

    /** {@inherit_doc} */
    //override
    public Fourth_Moment copy() 
    {
        return Fourth_Moment(this);
    }

}


