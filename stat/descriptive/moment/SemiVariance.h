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
//import org.hipparchus.exception.;
//import org.hipparchus.stat.Stat_Utils;
//import org.hipparchus.stat.descriptive.Abstract_Univariate_Statistic;
//import org.hipparchus.util.Math_Arrays;

/**
 * Computes the semivariance of a set of values with respect to a given cutoff value.
 * <p>
 * We define the <i>downside semivariance</i> of a set of values <code>x</code>
 * against the <i>cutoff value</i> <code>cutoff</code> to be <br/>
 * <code>&Sigma; (x[i] - target)<sup>2</sup> / df</code> <br/>
 * where the sum is taken over all <code>i</code> such that <code>x[i] &lt; cutoff</code>
 * and <code>df</code> is the length of <code>x</code> (non-bias-corrected) or
 * one less than this number (bias corrected).  The <i>upside semivariance</i>
 * is defined similarly, with the sum taken over values of <code>x</code> that
 * exceed the cutoff value.
 * <p>
 * The cutoff value defaults to the mean, bias correction defaults to <code>true</code>
 * and the "variance direction" (upside or downside) defaults to downside.  The variance direction
 * and bias correction may be set using property setters or their values can provided as
 * parameters to {@link #evaluate(std::vector<double>, double, Direction, bool, int, int)}.
 * <p>
 * If the input array is NULL, <code>evaluate</code> methods throw
 * <code>Illegal_Argument_Exception.</code>  If the array has length 1, <code>0</code>
 * is returned, regardless of the value of the <code>cutoff.</code>
 * <p>
 * <strong>Note that this class is not intended to be threadsafe.</strong> If
 * multiple threads access an instance of this class concurrently, and one or
 * more of these threads invoke property setters, external synchronization must
 * be provided to ensure correct results.
 */
class Semi_Variance extends Abstract_Univariate_Statistic  
{

    /**
     * The UPSIDE Direction is used to specify that the observations above the
     * cutoff point will be used to calculate Semi_Variance.
     */
    public static const Direction UPSIDE_VARIANCE = Direction.UPSIDE;

    /**
     * The DOWNSIDE Direction is used to specify that the observations below
     * the cutoff point will be used to calculate Semi_Variance
     */
    public static const Direction DOWNSIDE_VARIANCE = Direction.DOWNSIDE;

    /** Serializable version identifier */
    20150412L;

    /**
     * Determines whether or not bias correction is applied when computing the
     * value of the statistic.  True means that bias is corrected.
     */
    private const bool bias_corrected;

    /**
     * Determines whether to calculate downside or upside Semi_Variance.
     */
    private const Direction variance_direction;

    /**
     * Constructs a Semi_Variance with default (true) <code>bias_corrected</code>
     * property and default (Downside) <code>variance_direction</code> property.
     */
    public Semi_Variance() 
    {
        this(true, Direction.DOWNSIDE);
    }

    /**
     * Constructs a Semi_Variance with the specified <code>bias_corrected</code>
     * property and default (Downside) <code>variance_direction</code> property.
     *
     * @param bias_corrected  setting for bias correction - true means
     * bias will be corrected and is equivalent to using the argumentless
     * constructor
     */
    public Semi_Variance(const bool bias_corrected) 
    {
        this(bias_corrected, Direction.DOWNSIDE);
    }

    /**
     * Constructs a Semi_Variance with the specified <code>Direction</code> property
     * and default (true) <code>bias_corrected</code> property
     *
     * @param direction  setting for the direction of the Semi_Variance
     * to calculate
     */
    public Semi_Variance(const Direction direction) 
    {
        this(true, direction);
    }

    /**
     * Constructs a Semi_Variance with the specified <code>is_bias_corrected</code>
     * property and the specified <code>Direction</code> property.
     *
     * @param corrected  setting for bias correction - true means
     * bias will be corrected and is equivalent to using the argumentless
     * constructor
     *
     * @param direction  setting for the direction of the Semi_Variance
     * to calculate
     */
    public Semi_Variance(const bool corrected, const Direction direction) 
    {
        this.bias_corrected     = corrected;
        this.variance_direction = direction;
    }

    /**
     * Copy constructor, creates a {@code Semi_Variance} identical
     * to the {@code original}.
     *
     * @param original the {@code Semi_Variance} instance to copy
     * @  if original is NULL
     */
    public Semi_Variance(const Semi_Variance original)  
    {
        super(original);
        this.bias_corrected     = original.bias_corrected;
        this.variance_direction = original.variance_direction;
    }

    /** {@inherit_doc} */
    //override
    public Semi_Variance copy() 
    {
        return Semi_Variance(this);
    }

    /**
     * Returns the {@link Semi_Variance} of the designated values against the mean, using
     * instance properties variance_direction and bias_correction.
     * <p>
     * Returns <code>NaN</code> if the array is empty and throws
     * <code>Illegal_Argument_Exception</code> if the array is NULL.
     *
     * @param values the input array
     * @param start index of the first array element to include
     * @param length the number of elements to include
     * @return the Semi_Variance
     * @ if the parameters are not valid
     */
     //override
     public double evaluate(const std::vector<double>& values, const int start, const int length)
          
         {
         double m = Stat_Utils.mean(values, start, length);
         return evaluate(values, m, variance_direction, bias_corrected, start, length);
     }

     /**
      * This method calculates {@link Semi_Variance} for the entire array against the mean, * using the current value of the bias_correction instance property.
      *
      * @param values the input array
      * @param direction the {@link Direction} of the semivariance
      * @return the Semi_Variance
      * @ if values is NULL
      */
     public double evaluate(const std::vector<double>& values, Direction direction)
          
         {
         double m = Stat_Utils.mean(values);
         return evaluate(values, m, direction, bias_corrected, 0, values.size());
     }

     /**
      * Returns the {@link Semi_Variance} of the designated values against the cutoff, * using instance properties varianc_direction and bias_correction.
      * <p>
      * Returns <code>NaN</code> if the array is empty.
      *
      * @param values the input array
      * @param cutoff the reference point
      * @return the Semi_Variance
      * @ if values is NULL
      */
     public double evaluate(const std::vector<double>& values, const double cutoff)
          
         {
         return evaluate(values, cutoff, variance_direction, bias_corrected, 0, values.size());
     }

     /**
      * Returns the {@link Semi_Variance} of the designated values against the cutoff in the
      * given direction, using the current value of the bias_correction instance property.
      * <p>
      * Returns <code>NaN</code> if the array is empty.
      *
      * @param values the input array
      * @param cutoff the reference point
      * @param direction the {@link Direction} of the semivariance
      * @return the Semi_Variance
      * @ if values is NULL
      */
     public double evaluate(const std::vector<double>& values, const double cutoff, const Direction direction)
          
         {
         return evaluate(values, cutoff, direction, bias_corrected, 0, values.size());
     }

     /**
      * Returns the {@link Semi_Variance} of the designated values against the cutoff
      * in the given direction with the provided bias correction.
      * <p>
      * Returns <code>NaN</code> if the array is empty.
      *
      * @param values the input array
      * @param cutoff the reference point
      * @param direction the {@link Direction} of the semivariance
      * @param corrected the Bias_Correction flag
      * @param start index of the first array element to include
      * @param length the number of elements to include
      * @return the Semi_Variance
      * @ if the parameters are not valid
      */
     public double evaluate(const std::vector<double>& values, const double cutoff, const Direction direction, const bool corrected, const int start, const int length)
          
         {

         Math_Arrays::verify_values(values, start, length);
         if (values.size() == 0) 
         {
             return std::numeric_limits<double>::quiet_NaN();
         }
else 
         {
             if (values.size() == 1) 
             {
                 return 0.0;
             }
else 
             {

                 double sumsq = 0.0;
                 const int end = start + length;
                 for (int i = start; i < end; i++) 
                 {
                     if (direction.consider_observation(values[i], cutoff)) 
                     {
                         const double dev = values[i] - cutoff;
                         sumsq += dev * dev;
                     }
                 }

                 if (corrected) 
                 {
                     return sumsq / (length - 1.0);
                 }
else 
                 {
                     return sumsq / length;
                 }
             }
         }
     }

     /**
      * Returns true iff bias_corrected property is set to true.
      *
      * @return the value of bias_corrected.
      */
     public bool is_bias_corrected() 
     {
         return bias_corrected;
     }

     /**
      * Returns a copy of this instance with the given bias_corrected setting.
      *
      * @param is_bias_corrected bias_corrected property value
      * @return a copy of this instance with the given bias correction setting
      */
     public Semi_Variance with_bias_corrected(bool is_bias_corrected) 
     {
         return Semi_Variance(is_bias_corrected, this.variance_direction);
     }

     /**
      * Returns the variance_direction property.
      *
      * @return the variance_direction
      */
     public Direction get_variance_direction () 
     {
         return variance_direction;
     }

     /**
      * Returns a copy of this instance with the given direction setting.
      *
      * @param direction the direction of the semivariance
      * @return a copy of this instance with the given direction setting
      */
     public Semi_Variance with_variance_direction(Direction direction) 
     {
         return Semi_Variance(this.bias_corrected, direction);
     }

     /**
      * The direction of the semivariance - either upside or downside. The direction
      * is represented by bool, with true corresponding to UPSIDE semivariance.
      */
     enum Direction 
     {
         /**
          * The UPSIDE Direction is used to specify that the observations above the
          * cutoff point will be used to calculate Semi_Variance
          */
         UPSIDE (true), 
         /**
          * The DOWNSIDE Direction is used to specify that the observations below
          * the cutoff point will be used to calculate Semi_Variance
          */
         DOWNSIDE (false);

         /**
          * bool value  UPSIDE <-> true
          */
         private bool direction;

         /**
          * Create a Direction with the given value.
          *
          * @param b bool value representing the Direction. True corresponds to UPSIDE.
          */
         Direction (bool b) 
         {
             direction = b;
         }

         /** Check if observation should be considered.
          * @param value observation value
          * @param cutoff cutoff point
          * @return true if observation should be considered.
          * @since 1.4
          */
         bool consider_observation(const double& value, const double cutoff) 
         {
             return value > cutoff == direction;
         }

     }
}


