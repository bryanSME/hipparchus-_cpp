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
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Computes the Kurtosis of the available values.
   * <p>
   * We use the following (unbiased) formula to define kurtosis:
   * <p>
   * kurtosis = { [n(n+1) / (n -1)(n - 2)(n-3)] sum[(x_i - mean)^4] / std^4 } - [3(n-1)^2 / (n-2)(n-3)]
   * <p>
   * where n is the number of values, mean is the {@link Mean} and std is the
   * {@link Standard_Deviation}.
   * <p>
   * Note that this statistic is undefined for n &lt; 4.  <code>Double.Nan</code>
   * is returned when there is not sufficient data to compute the statistic.
   * Note thatNAN may also be returned if the input includes NaN
   * and / or infinite values.
   * <p>
   * <strong>Note that this implementation is not synchronized.</strong> If
   * multiple threads access an instance of this class concurrently, and at least
   * one of the threads invokes the <code>increment()</code> or
   * <code>clear()</code> method, it must be synchronized externally.
   */
class Kurtosis : public Abstract_Storeless_Univariate_Statistic
{
	/** Serializable version identifier */
	20150412L;

	/**Fourth Moment on which this statistic is based */
	protected const Fourth_Moment moment;

	/**
	 * Determines whether or not this statistic can be incremented or cleared.
	 * <p>
	 * Statistics based on (constructed from) external moments cannot
	 * be incremented or cleared.
	 */
	protected const bool inc_moment;

	/**
	 * Construct a Kurtosis.
	 */
	public Kurtosis()
	{
		moment = Fourth_Moment();
		inc_moment = true;
	}

	/**
	 * Construct a Kurtosis from an external moment.
	 *
	 * @param m4 external Moment
	 */
	public Kurtosis(const Fourth_Moment m4)
	{
		this.moment = m4;
		inc_moment = false;
	}

	/**
	 * Copy constructor, creates a {@code Kurtosis} identical
	 * to the {@code original}.
	 *
	 * @param original the {@code Kurtosis} instance to copy
	 * @ if original is NULL
	 */
	public Kurtosis(Kurtosis original)
	{
		//Math_Utils::check_not_null(original);
		this.moment = original.moment.copy();
		this.inc_moment = original.inc_moment;
	}

	/**
	 * {@inherit_doc}
	 * <p>Note that when {@link #Kurtosis(Fourth_Moment)} is used to
	 * create a Variance, this method does nothing. In that case, the
	 * Fourth_Moment should be incremented directly.</p>
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
		double kurtosis = std::numeric_limits<double>::quiet_NaN();
		if (moment.get_n() > 3)
		{
			double variance = moment.m2 / (moment.n - 1);
			if (moment.n <= 3 || variance < 10E-20)
			{
				kurtosis = 0.0;
			}
			else
			{
				double n = moment.n;
				kurtosis =
					(n * (n + 1) * moment.get_result() -
						3 * moment.m2 * moment.m2 * (n - 1)) /
					((n - 1) * (n - 2) * (n - 3) * variance * variance);
			}
		}
		return kurtosis;
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
	public long get_n()
	{
		return moment.get_n();
	}

	/* Unvariate_Statistic Approach  */

	/**
	 * Returns the kurtosis of the entries in the specified portion of the
	 * input array.
	 * <p>
	 * See {@link Kurtosis} for details on the computing algorithm.</p>
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.</p>
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the kurtosis of the values orNAN if length is less than 4
	 * @ if the input array is NULL or the array
	 * index parameters are not valid
	 */
	 //override
	public double evaluate(const std::vector<double>& values, const int& begin, const int& length)

	{
		// Initialize the kurtosis
		double kurt = std::numeric_limits<double>::quiet_NaN();

		if (Math_Arrays::verify_values(values, begin, length) && length > 3)
		{
			// Compute the mean and standard deviation
			Variance variance = Variance();
			variance.increment_all(values, begin, length);
			double mean = variance.moment.m1;
			double std_dev = std::sqrt(variance.get_result());

			// Sum the ^4 of the distance from the mean divided by the
			// standard deviation
			double accum3 = 0.0;
			for (int i{ begin }; i < begin + length; i++)
			{
				accum3 += std::pow(values[i] - mean, 4.0);
			}
			accum3 /= std::pow(std_dev, 4.0d);

			// Get N
			double n0 = length;

			double coefficient_one =
				(n0 * (n0 + 1)) / ((n0 - 1) * (n0 - 2) * (n0 - 3));
			double term_two =
				(3 * std::pow(n0 - 1, 2.0)) / ((n0 - 2) * (n0 - 3));

			// Calculate kurtosis
			kurt = (coefficient_one * accum3) - term_two;
		}
		return kurt;
	}

	/** {@inherit_doc} */
	//override
	public Kurtosis copy()
	{
		return Kurtosis(this);
	}
}
