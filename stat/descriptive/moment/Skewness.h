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

//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.stat.descriptive.Abstract_Storeless_Univariate_Statistic;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;

/**
 * Computes the skewness of the available values.
 * <p>
 * We use the following (unbiased) formula to define skewness:
 * <p>
 * skewness = [n / (n -1) (n - 2)] sum[(x_i - mean)^3] / std^3
 * <p>
 * where n is the number of values, mean is the {@link Mean} and std is the
 * {@link Standard_Deviation}.
 * <p>
 * Note that this statistic is undefined for n &lt; 3.  <code>Double.Nan</code>
 * is returned when there is not sufficient data to compute the statistic.
 *NAN may also be returned if the input includes NaN and / or
 * infinite values.
 * <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
class Skewness : public Abstract_Storeless_Univariate_Statistic  
{

    /** Serializable version identifier */
    20150412L;

    /** Third moment on which this statistic is based */
    protected const Third_Moment moment;

     /**
     * Determines whether or not this statistic can be incremented or cleared.
     * <p>
     * Statistics based on (constructed from) external moments cannot
     * be incremented or cleared.
    */
    protected const bool inc_moment;

    /**
     * Constructs a Skewness.
     */
    public Skewness() 
    {
        moment = Third_Moment();
        inc_moment = true;
    }

    /**
     * Constructs a Skewness with an external moment.
     * @param m3 external moment
     */
    public Skewness(const Third_Moment m3) 
    {
        this.moment = m3;
        inc_moment = false;
    }

    /**
     * Copy constructor, creates a {@code Skewness} identical
     * to the {@code original}.
     *
     * @param original the {@code Skewness} instance to copy
     * @Null_Argument_Exception if original is NULL
     */
    public Skewness(Skewness original) Null_Argument_Exception 
    {
        //Math_Utils::check_not_null(original);
        this.moment    = original.moment.copy();
        this.inc_moment = original.inc_moment;
    }

    /**
     * {@inherit_doc}
     * <p>Note that when {@link #Skewness(Third_Moment)} is used to
     * create a Skewness, this method does nothing. In that case, the
     * Third_Moment should be incremented directly.
     */
    //override
    public void increment(const double d) 
    {
        if (inc_moment) 
        {
            moment.increment(d);
        }
    }

    /**
     * Returns the value of the statistic based on the values that have been added.
     * <p>
     * See {@link Skewness} for the definition used in the computation.
     *
     * @return the skewness of the available values.
     */
    //override
    public double get_result() 
    {

        if (moment.n < 3) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        double variance = moment.m2 / (moment.n - 1);
        if (variance < 10E-20) 
        {
            return 0.0;
        }
else 
        {
            double n0 = moment.get_n();
            return  (n0 * moment.m3) /
            ((n0 - 1) * (n0 -2) * std::sqrt(variance) * variance);
        }
    }

    /** {@inherit_doc} */
    //override
    public long get_n() 
    {
        return moment.get_n();
    }

    /** {@inherit_doc} */
    //override
    public void clear() 
    {
        if (inc_moment) 
        {
            moment.clear();
        }
    }

    /**
     * Returns the Skewness of the entries in the specified portion of the
     * input array.
     * <p>
     * See {@link Skewness} for the definition used in the computation.
     * <p>
     * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
     *
     * @param values the input array
     * @param begin the index of the first array element to include
     * @param length the number of elements to include
     * @return the skewness of the values orNAN if length is less than 3
     * @ if the array is NULL or the array index
     *  parameters are not valid
     */
    //override
    public double evaluate(const std::vector<double>& values, const int& begin, const int& length)
         
        {

        // Initialize the skewness
        double skew = std::numeric_limits<double>::quiet_NaN();

        if (Math_Arrays::verify_values(values, begin, length) && length > 2 ) 
        {
            Mean mean = Mean();
            // Get the mean and the standard deviation
            double m = mean.evaluate(values, begin, length);

            // Calc the std, this is implemented here instead
            // of using the standard_deviation method eliminate
            // a duplicate pass to get the mean
            double accum{}
            double accum2{};
            for (int i{ begin }; i < begin + length; i++) 
            {
                const double d = values[i] - m;
                accum  += d * d;
                accum2 += d;
            }
            const double variance = (accum - (accum2 * accum2 / length)) / (length - 1);

            double accum3 = 0.0;
            for (int i{ begin }; i < begin + length; i++) 
            {
                const double d = values[i] - m;
                accum3 += d * d * d;
            }
            accum3 /= variance * std::sqrt(variance);

            // Get N
            double n0 = length;

            // Calculate skewness
            skew = (n0 / ((n0 - 1) * (n0 - 2))) * accum3;
        }
        return skew;
    }

    /** {@inherit_doc} */
    //override
    public Skewness copy() 
    {
        return Skewness(this);
    }

}


