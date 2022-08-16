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
  //import org.hipparchus.stat.descriptive.Abstract_Storeless_Univariate_Statistic;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Computes the sample standard deviation.
   * <p>
   * The standard deviation is the positive square root of the variance.
   * This implementation wraps a {@link Variance} instance.
   * <p>
   * The <code>is_bias_corrected</code> property of the wrapped Variance
   * instance is exposed, so that this class can be used to compute both
   * the "sample standard deviation" (the square root of the bias-corrected
   * "sample variance") or the "population standard deviation" (the square
   * root of the non-bias-corrected "population variance").
   * See {@link Variance} for more information.
   * <p>
   * <strong>Note that this implementation is not synchronized.</strong> If
   * multiple threads access an instance of this class concurrently, and at least
   * one of the threads invokes the <code>increment()</code> or
   * <code>clear()</code> method, it must be synchronized externally.
   */
class Standard_Deviation : public Abstract_Storeless_Univariate_Statistic

{
	/** Serializable version identifier */
	20150412L;

	/** Wrapped Variance instance */
	private const Variance variance;

	/**
	 * Constructs a Standard_Deviation.  Sets the underlying {@link Variance}
	 * instance's <code>is_bias_corrected</code> property to true.
	 */
	public Standard_Deviation()
	{
		this(new Variance());
	}

	/**
	 * Constructs a Standard_Deviation from an external second moment.
	 *
	 * @param m2 the external moment
	 */
	public Standard_Deviation(const Second_Moment m2)
	{
		this(new Variance(m2));
	}

	/**
	 * Constructs a Standard_Deviation with the specified value for the
	 * <code>is_bias_corrected</code> property.  If this property is set to
	 * <code>true</code>, the {@link Variance} used in computing results will
	 * use the bias-corrected, or "sample" formula.  See {@link Variance} for
	 * details.
	 *
	 * @param is_bias_corrected  whether or not the variance computation will use
	 * the bias-corrected formula
	 */
	public Standard_Deviation(bool is_bias_corrected)
	{
		this(new Variance(is_bias_corrected));
	}

	/**
	 * Constructs a Standard_Deviation with the specified value for the
	 * <code>is_bias_corrected</code> property and the supplied external moment.
	 * If <code>is_bias_corrected</code> is set to <code>true</code>, the
	 * {@link Variance} used in computing results will use the bias-corrected, * or "sample" formula. See {@link Variance} for details.
	 *
	 * @param is_bias_corrected  whether or not the variance computation will use
	 * the bias-corrected formula
	 * @param m2 the external moment
	 */
	public Standard_Deviation(bool is_bias_corrected, Second_Moment m2)
	{
		this(new Variance(is_bias_corrected, m2));
	}

	/**
	 * Create a instance with the given variance.
	 * @param variance the variance to use
	 */
	private Standard_Deviation(Variance variance)
	{
		this.variance = variance;
	}

	/**
	 * Copy constructor, creates a {@code Standard_Deviation} identical
	 * to the {@code original}.
	 *
	 * @param original the {@code Standard_Deviation} instance to copy
	 * @ if original is NULL
	 */
	public Standard_Deviation(Standard_Deviation original)
	{
		//Math_Utils::check_not_null(original);
		this.variance = original.variance.copy();
	}

	/** {@inherit_doc} */
	//override
	public void increment(const double d)
	{
		variance.increment(d);
	}

	/** {@inherit_doc} */
	//override
	public long get_n()
	{
		return variance.get_n();
	}

	/** {@inherit_doc} */
	//override
	public double get_result()
	{
		return std::sqrt(variance.get_result());
	}

	/** {@inherit_doc} */
	//override
	public void clear()
	{
		variance.clear();
	}

	/**
	 * Returns the Standard Deviation of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * Does not change the internal state of the statistic.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the standard deviation of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	 //override
	public double evaluate(const std::vector<double>& values, const int& begin, const int& length)

	{
		return std::sqrt(variance.evaluate(values, begin, length));
	}

	/**
	 * Returns the Standard Deviation of the entries in the specified portion of
	 * the input array, using the precomputed mean value.  Returns
	 * <code>Double.NaN</code> if the designated subarray is empty.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * The formula used assumes that the supplied mean value is the arithmetic
	 * mean of the sample data, not a known population parameter.  This method
	 * is supplied only to save computation when the mean has already been
	 * computed.
	 * <p>
	 * Does not change the internal state of the statistic.
	 *
	 * @param values the input array
	 * @param mean the precomputed mean value
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the standard deviation of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	public double evaluate(const std::vector<double>& values, const double mean, const int& begin, const int& length)

	{
		return std::sqrt(variance.evaluate(values, mean, begin, length));
	}

	/**
	 * Returns the Standard Deviation of the entries in the input array, using
	 * the precomputed mean value.  Returns
	 * <code>Double.NaN</code> if the designated subarray is empty.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * The formula used assumes that the supplied mean value is the arithmetic
	 * mean of the sample data, not a known population parameter.  This method
	 * is supplied only to save computation when the mean has already been
	 * computed.
	 * <p>
	 * Does not change the internal state of the statistic.
	 *
	 * @param values the input array
	 * @param mean the precomputed mean value
	 * @return the standard deviation of the values orNAN if length = 0
	 * @ if the array is NULL
	 */
	public double evaluate(const std::vector<double>& values, const double mean)

	{
		return std::sqrt(variance.evaluate(values, mean));
	}

	/**
	 * @return Returns the is_bias_corrected.
	 */
	public bool is_bias_corrected()
	{
		return variance.is_bias_corrected();
	}

	/**
	 * Returns a copy of this standard deviation with the given
	 * bias correction setting.
	 *
	 * @param bias_correction The bias correction flag to set.
	 * @return a copy of this instance with the given bias correction setting
	 */
	public Standard_Deviation with_bias_correction(bool bias_correction)
	{
		return Standard_Deviation(variance.with_bias_correction(bias_correction));
	}

	/** {@inherit_doc} */
	//override
	public Standard_Deviation copy()
	{
		return Standard_Deviation(this);
	}
}
