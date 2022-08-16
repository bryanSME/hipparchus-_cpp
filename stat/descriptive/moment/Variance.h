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
  //import org.hipparchus.stat.descriptive.Abstract_Storeless_Univariate_Statistic;
  //import org.hipparchus.stat.descriptive.Aggregatable_Statistic;
  //import org.hipparchus.stat.descriptive.Weighted_Evaluation;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Computes the variance of the available values.  By default, the unbiased
   * "sample variance" definitional formula is used:
   * <p>
   * variance = sum((x_i - mean)^2) / (n - 1)
   * <p>
   * where mean is the {@link Mean} and <code>n</code> is the number
   * of sample observations.
   * <p>
   * The definitional formula does not have good numerical properties, so
   * this implementation does not compute the statistic using the definitional
   * formula.
   * <ul>
   * <li> The <code>get_result</code> method computes the variance using
   * updating formulas based on West's algorithm, as described in
   * <a href="http://doi.acm.org/10.1145/359146.359152"> Chan, T. F. and
   * J. G. Lewis 1979, <i>Communications of the ACM</i>, * vol. 22 no. 9, pp. 526-531.</a></li>
   * <li> The <code>evaluate</code> methods leverage the fact that they have the
   * full array of values in memory to execute a two-pass algorithm.
   * Specifically, these methods use the "corrected two-pass algorithm" from
   * Chan, Golub, Levesque, <i>Algorithms for Computing the Sample Variance</i>, * American Statistician, vol. 37, no. 3 (1983) pp. 242-247.</li>
   * </ul>
   * <p>
   * Note that adding values using <code>increment</code> or
   * <code>increment_all</code> and then executing <code>get_result</code> will
   * sometimes give a different, less accurate, result than executing
   * <code>evaluate</code> with the full array of values. The former approach
   * should only be used when the full array of values is not available.
   * <p>
   * The "population variance"  ( sum((x_i - mean)^2) / n ) can also
   * be computed using this statistic.  The <code>is_bias_corrected</code>
   * property determines whether the "population" or "sample" value is
   * returned by the <code>evaluate</code> and <code>get_result</code> methods.
   * To compute population variances, set this property to <code>false.</code>
   * <p>
   * <strong>Note that this implementation is not synchronized.</strong> If
   * multiple threads access an instance of this class concurrently, and at least
   * one of the threads invokes the <code>increment()</code> or
   * <code>clear()</code> method, it must be synchronized externally.
   */
class Variance : public Abstract_Storeless_Univariate_Statistic
	: Aggregatable_Statistic<Variance>, Weighted_Evaluation
{
	/** Serializable version identifier */
	20150412L;

/** Second_Moment is used in incremental calculation of Variance*/
protected const Second_Moment moment;

/**
 * Whether or not {@link #incrementstatic_cast<double>(} should increment
 * the internal second moment. When a Variance is constructed with an
 * external Second_Moment as a constructor parameter, this property is
 * set to false and increments must be applied to the second moment
 * directly.
 */
protected const bool inc_moment;

/**
 * Whether or not bias correction is applied when computing the
 * value of the statistic. True means that bias is corrected.  See
 * {@link Variance} for details on the formula.
 */
private const bool is_bias_corrected;

/**
 * Constructs a Variance with default (true) <code>is_bias_corrected</code>
 * property.
 */
public Variance()
{
	this(true);
}

/**
 * Constructs a Variance based on an external second moment.
 * <p>
 * When this constructor is used, the statistic may only be
 * incremented via the moment, i.e., {@link #incrementstatic_cast<double>(}
 * does nothing; whereas {@code m2.increment(value)} increments
 * both {@code m2} and the Variance instance constructed from it.
 *
 * @param m2 the Second_Moment (Third or Fourth moments work here as well.)
 */
public Variance(const Second_Moment m2)
{
	this(true, m2);
}

/**
 * Constructs a Variance with the specified <code>is_bias_corrected</code>
 * property.
 *
 * @param is_bias_corrected  setting for bias correction - true means
 * bias will be corrected and is equivalent to using the argumentless
 * constructor
 */
public Variance(bool is_bias_corrected)
{
	this(new Second_Moment(), true, is_bias_corrected);
}

/**
 * Constructs a Variance with the specified <code>is_bias_corrected</code>
 * property and the supplied external second moment.
 *
 * @param is_bias_corrected  setting for bias correction - true means
 * bias will be corrected
 * @param m2 the Second_Moment (Third or Fourth moments work
 * here as well.)
 */
public Variance(bool is_bias_corrected, Second_Moment m2)
{
	this(m2, false, is_bias_corrected);
}

/**
 * Constructs a Variance with the specified parameters.
 *
 * @param m2 the Second_Moment (Third or Fourth moments work
 * here as well.)
 * @param inc_moment if the moment shall be incremented
 * @param is_bias_corrected  setting for bias correction - true means
 * bias will be corrected
 */
private Variance(Second_Moment m2, bool inc_moment, bool is_bias_corrected)
{
	this.moment = m2;
	this.inc_moment = inc_moment;
	this.is_bias_corrected = is_bias_corrected;
}

/**
 * Copy constructor, creates a {@code Variance} identical
 * to the {@code original}.
 *
 * @param original the {@code Variance} instance to copy
 * @ if original is NULL
 */
public Variance(Variance original)
{
	//Math_Utils::check_not_null(original);
	this.moment = original.moment.copy();
	this.inc_moment = original.inc_moment;
	this.is_bias_corrected = original.is_bias_corrected;
}

/**
 * {@inherit_doc}
 * <p>If all values are available, it is more accurate to use
 * {@link #evaluate(std::vector<double>)} rather than adding values one at a time
 * using this method and then executing {@link #get_result}, since
 * <code>evaluate</code> leverages the fact that is has the full
 * list of values together to execute a two-pass algorithm.
 * See {@link Variance}.
 * <p>
 * Note also that when {@link #Variance(Second_Moment)} is used to
 * create a Variance, this method does nothing. In that case, the
 * Second_Moment should be incremented directly.
 */
 //override
 public void increment(const double d)
 {
	 if (inc_moment)
	 {
		 moment.increment(d);
	 }
 }

 /** {@inherit_doc} */
 //override
 public double get_result()
 {
	 if (moment.n == 0)
	 {
		 return std::numeric_limits<double>::quiet_NaN();
	 }
else if (moment.n == 1)
		{
			return 0;
		}
else
		{
			if (is_bias_corrected)
			{
				return moment.m2 / (moment.n - 1d);
			}
else
			{
				return moment.m2 / (moment.n);
			}
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

 /** {@inherit_doc} */
 //override
 public void aggregate(Variance other)
 {
	 //Math_Utils::check_not_null(other);
	 if (inc_moment)
	 {
		 this.moment.aggregate(other.moment);
	 }
 }

 /**
  * Returns the variance of the entries in the specified portion of
  * the input array, or <code>Double.NaN</code> if the designated subarray
  * is empty.  Note thatNAN may also be returned if the input
  * includes NaN and / or infinite values.
  * <p>
  * See {@link Variance} for details on the computing algorithm.</p>
  * <p>
  * Returns 0 for a single-value (i.e. length = 1) sample.</p>
  * <p>
  * Does not change the internal state of the statistic.</p>
  * <p>
  * Throws <code></code> if the array is NULL.</p>
  *
  * @param values the input array
  * @param begin index of the first array element to include
  * @param length the number of elements to include
  * @return the variance of the values orNAN if length = 0
  * @ if the array is NULL or the array index
  *  parameters are not valid
  */
  //override
  public double evaluate(const std::vector<double>& values, const int& begin, const int& length)

	  {
	  double var = std::numeric_limits<double>::quiet_NaN();

	  if (Math_Arrays::verify_values(values, begin, length))
	  {
		  if (length == 1)
		  {
			  var = 0.0;
		  }
else if (length > 1)
			{
				double m = Stat_Utils.mean(values, begin, length);
				var = evaluate(values, m, begin, length);
			}
		}
		return var;
	}

  /**
   * Returns the weighted variance of the entries in the specified portion of
   * the input array, or <code>Double.NaN</code> if the designated subarray
   * is empty.
   * <p>
   * Uses the formula
   * <pre>
   *   &Sigma;(weights[i]*(values[i] - weighted_mean)<sup>2</sup>)/(&Sigma;(weights[i]) - 1)
   * </pre>
   * where weighted_mean is the weighted mean.
   * <p>
   * This formula will not return the same result as the unweighted variance when all
   * weights are equal, unless all weights are equal to 1. The formula assumes that
   * weights are to be treated as "expansion values," as will be the case if for example
   * the weights represent frequency counts. To normalize weights so that the denominator
   * in the variance computation equals the length of the input vector minus one, use
   * <pre>
   *   <code>evaluate(values, Math_Arrays::normalize_array(weights, values.size()));</code>
   * </pre>
   * <p>
   * Returns 0 for a single-value (i.e. length = 1) sample.
   * <p>
   * Throws <code>Illegal_Argument_Exception</code> if any of the following are true:
   * <ul><li>the values array is NULL</li>
   *     <li>the weights array is NULL</li>
   *     <li>the weights array does not have the same length as the values array</li>
   *     <li>the weights array contains one or more infinite values</li>
   *     <li>the weights array contains one or more NaN values</li>
   *     <li>the weights array contains negative values</li>
   *     <li>the start and length arguments do not determine a valid array</li>
   * </ul>
   * <p>
   * Does not change the internal state of the statistic.
   *
   * @param values the input array
   * @param weights the weights array
   * @param begin index of the first array element to include
   * @param length the number of elements to include
   * @return the weighted variance of the values orNAN if length = 0
   * @ if the parameters are not valid
   */
   //override
   public double evaluate(const std::vector<double>& values, const std::vector<double> weights, const int& begin, const int& length)

	   {
	   double var = std::numeric_limits<double>::quiet_NaN();
	   if (Math_Arrays::verify_values(values, weights,begin, length))
	   {
		   if (length == 1)
		   {
			   var = 0.0;
		   }
else if (length > 1)
			{
				Mean mean = Mean();
				double m = mean.evaluate(values, weights, begin, length);
				var = evaluate(values, weights, m, begin, length);
			}
		}
		return var;
	}

   /**
	* Returns the variance of the entries in the specified portion of
	* the input array, using the precomputed mean value. Returns
	* <code>Double.NaN</code> if the designated subarray is empty.
	* <p>
	* See {@link Variance} for details on the computing algorithm.
	* <p>
	* The formula used assumes that the supplied mean value is the arithmetic
	* mean of the sample data, not a known population parameter.  This method
	* is supplied only to save computation when the mean has already been
	* computed.
	* <p>
	* Returns 0 for a single-value (i.e. length = 1) sample.
	* <p>
	* Does not change the internal state of the statistic.
	*
	* @param values the input array
	* @param mean the precomputed mean value
	* @param begin index of the first array element to include
	* @param length the number of elements to include
	* @return the variance of the values orNAN if length = 0
	* @ if the array is NULL or the array index
	*  parameters are not valid
	*/
   public double evaluate(const std::vector<double>& values, const double mean, const int& begin, const int& length)

	   {
	   double var = std::numeric_limits<double>::quiet_NaN();
	   if (Math_Arrays::verify_values(values, begin, length))
	   {
		   if (length == 1)
		   {
			   var = 0.0;
		   }
else if (length > 1)
			{
				double accum{}
				double accum2{};
				for (int i{ begin }; i < begin + length; i++)
				{
					const double dev = values[i] - mean;
					accum += dev * dev;
					accum2 += dev;
				}
				double len = length;
				if (is_bias_corrected)
				{
					var = (accum - (accum2 * accum2 / len)) / (len - 1.0);
				}
else
				{
					var = (accum - (accum2 * accum2 / len)) / len;
				}
			}
		}
		return var;
	}

   /**
	* Returns the variance of the entries in the input array, using the
	* precomputed mean value.  Returns <code>Double.NaN</code> if the array
	* is empty.
	* <p>
	* See {@link Variance} for details on the computing algorithm.
	* <p>
	* If <code>is_bias_corrected</code> is <code>true</code> the formula used
	* assumes that the supplied mean value is the arithmetic mean of the
	* sample data, not a known population parameter.  If the mean is a known
	* population parameter, or if the "population" version of the variance is
	* desired, set <code>is_bias_corrected</code> to <code>false</code> before
	* invoking this method.
	* <p>
	* Returns 0 for a single-value (i.e. length = 1) sample.
	* <p>
	* Does not change the internal state of the statistic.
	*
	* @param values the input array
	* @param mean the precomputed mean value
	* @return the variance of the values orNAN if the array is empty
	* @ if the array is NULL
	*/
   public double evaluate(const std::vector<double>& values, const double mean)

	   {
	   return evaluate(values, mean, 0, values.size());
   }

   /**
	* Returns the weighted variance of the entries in the specified portion of
	* the input array, using the precomputed weighted mean value.  Returns
	* <code>Double.NaN</code> if the designated subarray is empty.
	* <p>
	* Uses the formula
	* <pre>
	*   &Sigma;(weights[i]*(values[i] - mean)<sup>2</sup>)/(&Sigma;(weights[i]) - 1)
	* </pre>
	* <p>
	* The formula used assumes that the supplied mean value is the weighted arithmetic
	* mean of the sample data, not a known population parameter. This method
	* is supplied only to save computation when the mean has already been
	* computed.
	* <p>
	* This formula will not return the same result as the unweighted variance when all
	* weights are equal, unless all weights are equal to 1. The formula assumes that
	* weights are to be treated as "expansion values," as will be the case if for example
	* the weights represent frequency counts. To normalize weights so that the denominator
	* in the variance computation equals the length of the input vector minus one, use
	* <pre>
	*   <code>evaluate(values, Math_Arrays::normalize_array(weights, values.size()), mean);</code>
	* </pre>
	* <p>
	* Returns 0 for a single-value (i.e. length = 1) sample.
	* <p>
	* Throws <code></code> if any of the following are true:
	* <ul><li>the values array is NULL</li>
	*     <li>the weights array is NULL</li>
	*     <li>the weights array does not have the same length as the values array</li>
	*     <li>the weights array contains one or more infinite values</li>
	*     <li>the weights array contains one or more NaN values</li>
	*     <li>the weights array contains negative values</li>
	*     <li>the start and length arguments do not determine a valid array</li>
	* </ul>
	* <p>
	* Does not change the internal state of the statistic.
	*
	* @param values the input array
	* @param weights the weights array
	* @param mean the precomputed weighted mean value
	* @param begin index of the first array element to include
	* @param length the number of elements to include
	* @return the variance of the values orNAN if length = 0
	* @ if the parameters are not valid
	*/
   public double evaluate(const std::vector<double>& values, const std::vector<double> weights, const double mean, const int& begin, const int& length)

	   {
	   double var = std::numeric_limits<double>::quiet_NaN();

	   if (Math_Arrays::verify_values(values, weights, begin, length))
	   {
		   if (length == 1)
		   {
			   var = 0.0;
		   }
else if (length > 1)
			{
				double accum{}
				double accum2{};
				for (int i{ begin }; i < begin + length; i++)
				{
					const double dev = values[i] - mean;
					accum += weights[i] * (dev * dev);
					accum2 += weights[i] * dev;
				}

				double sum_wts = 0;
				for (int i{ begin }; i < begin + length; i++)
				{
					sum_wts += weights[i];
				}

				if (is_bias_corrected)
				{
					var = (accum - (accum2 * accum2 / sum_wts)) / (sum_wts - 1.0);
				}
else
				{
					var = (accum - (accum2 * accum2 / sum_wts)) / sum_wts;
				}
			}
		}
		return var;
	}

   /**
	* Returns the weighted variance of the values in the input array, using
	* the precomputed weighted mean value.
	* <p>
	* Uses the formula
	* <pre>
	*   &Sigma;(weights[i]*(values[i] - mean)<sup>2</sup>)/(&Sigma;(weights[i]) - 1)
	* </pre>
	* <p>
	* The formula used assumes that the supplied mean value is the weighted arithmetic
	* mean of the sample data, not a known population parameter. This method
	* is supplied only to save computation when the mean has already been
	* computed.
	* <p>
	* This formula will not return the same result as the unweighted variance when all
	* weights are equal, unless all weights are equal to 1. The formula assumes that
	* weights are to be treated as "expansion values," as will be the case if for example
	* the weights represent frequency counts. To normalize weights so that the denominator
	* in the variance computation equals the length of the input vector minus one, use
	* <pre>
	*   <code>evaluate(values, Math_Arrays::normalize_array(weights, values.size()), mean);</code>
	* </pre>
	* <p>
	* Returns 0 for a single-value (i.e. length = 1) sample.
	* <p>
	* Throws <code></code> if any of the following are true:
	* <ul><li>the values array is NULL</li>
	*     <li>the weights array is NULL</li>
	*     <li>the weights array does not have the same length as the values array</li>
	*     <li>the weights array contains one or more infinite values</li>
	*     <li>the weights array contains one or more NaN values</li>
	*     <li>the weights array contains negative values</li>
	* </ul>
	* <p>
	* Does not change the internal state of the statistic.
	*
	* @param values the input array
	* @param weights the weights array
	* @param mean the precomputed weighted mean value
	* @return the variance of the values orNAN if length = 0
	* @ if the parameters are not valid
	*/
   public double evaluate(const std::vector<double>& values, const std::vector<double> weights, const double mean)

	   {
	   return evaluate(values, weights, mean, 0, values.size());
   }

   /**
	* @return Returns the is_bias_corrected.
	*/
   public bool is_bias_corrected()
   {
	   return is_bias_corrected;
   }

   /**
	* Returns a copy of this variance with the given bias correction
	* setting.
	*
	* @param bias_correction The bias correction flag to set.
	* @return a copy of this instance with the given bias correction setting
	*/
   public Variance with_bias_correction(bool bias_correction)
   {
	   return Variance(this.moment, this.inc_moment, bias_correction);
   }

   /** {@inherit_doc} */
   //override
   public Variance copy()
   {
	   return Variance(this);
   }
}
