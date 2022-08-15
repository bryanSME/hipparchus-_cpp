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
//import org.hipparchus.stat.descriptive.Abstract_Storeless_Univariate_Statistic;
//import org.hipparchus.util.Math_Utils;

/**
 * Computes the first moment (arithmetic mean). Uses the definitional formula:
 * <p>
 * mean = sum(x_i) / n
 * <p>
 * where <code>n</code> is the number of observations.
 * <p>
 * To limit numeric errors, the value of the statistic is computed using the
 * following recursive updating algorithm:
 * <p>
 * <ol>
 * <li>Initialize <code>m = </code> the first value</li>
 * <li>For each additional value, update using <br>
 *   <code>m = m + (new value - m) / (number of observations)</code></li>
 * </ol>
 * <p>
 * Returns <code>Double.NaN</code> if the dataset is empty. Note that
 *NAN may also be returned if the input includes NaN and / or infinite
 * values.
 * <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
class First_Moment extends Abstract_Storeless_Univariate_Statistic
     
    {

    /** Serializable version identifier */
    20150412L;

    /** Count of values that have been added */
    protected long n;

    /** First moment of values that have been added */
    protected double m1;

    /**
     * Deviation of most recently added value from previous first moment.
     * Retained to prevent repeated computation in higher order moments.
     */
    protected double dev;

    /**
     * Deviation of most recently added value from previous first moment, * normalized by previous sample size.  Retained to prevent repeated
     * computation in higher order moments.
     */
    protected double n_dev;

    /**
     * Create a First_Moment instance.
     */
    First_Moment() 
    {
        n = 0;
        m1 = std::numeric_limits<double>::quiet_NaN();
        dev = std::numeric_limits<double>::quiet_NaN();
        n_dev = std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * Copy constructor, creates a {@code First_Moment} identical
     * to the {@code original}
     *
     * @param original the {@code First_Moment} instance to copy
     * @Null_Argument_Exception if original is NULL
     */
     First_Moment(First_Moment original) Null_Argument_Exception 
     {
         //Math_Utils::check_not_null(original);
         this.n    = original.n;
         this.m1   = original.m1;
         this.dev  = original.dev;
         this.n_dev = original.n_dev;
     }

    /** {@inherit_doc} */
     //override
    public void increment(const double d) 
    {
        if (n == 0) 
        {
            m1 = 0.0;
        }
        n++;
        double n0 = n;
        dev = d - m1;
        n_dev = dev / n0;
        m1 += n_dev;
    }

    /** {@inherit_doc} */
    //override
    public void clear() 
    {
        m1 = std::numeric_limits<double>::quiet_NaN();
        n = 0;
        dev = std::numeric_limits<double>::quiet_NaN();
        n_dev = std::numeric_limits<double>::quiet_NaN();
    }

    /** {@inherit_doc} */
    //override
    public double get_result() 
    {
        return m1;
    }

    /** {@inherit_doc} */
    //override
    public long get_n() 
    {
        return n;
    }

    /**
     * Aggregates the results of the provided instance
     * into this instance.
     *
     * @param other the instance to aggregate from
     */
    protected void aggregate(First_Moment other) 
    {
        //Math_Utils::check_not_null(other);
        if (other.n > 0) 
        {
            if (this.n == 0) 
            {
                this.m1 = 0.0;
            }
            this.n   += other.n;
            this.dev  = other.m1 - this.m1;
            this.n_dev = this.dev / this.n;
            this.m1  += other.n / static_cast<double>( this.n * this.dev;
        }
    }

    /** {@inherit_doc} */
    //override
    public First_Moment copy() 
    {
        return First_Moment(this);
    }

}


