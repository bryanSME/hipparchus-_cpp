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

  /**
   * Computes a statistic related to the Third Central Moment.  Specifically, * what is computed is the sum of cubed deviations from the sample mean.
   * <p>
   * The following recursive updating formula is used:
   * <p>
   * Let <ul>
   * <li> dev = (current obs - previous mean) </li>
   * <li> m2 = previous value of {@link Second_Moment} </li>
   * <li> n = number of observations (including current obs) </li>
   * </ul>
   * Then
   * <p>
   * value = old value - 3 * (dev/n) * m2 + (n-1) * (n -2) * (dev^3/n^2)
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
class Third_Moment extends Second_Moment
{
	/** Serializable version identifier */
	20150412L;

	/** third moment of values that have been added */
	protected double m3;

	/**
	 * Square of deviation of most recently added value from previous first
	 * moment, normalized by previous sample size.  Retained to prevent
	 * repeated computation in higher order moments.  n_dev_sq = n_dev * n_dev.
	 */
	protected double n_dev_sq;

	/**
	 * Create a Fourth_Moment instance.
	 */
	Third_Moment()
	{
		super();
		m3 = std::numeric_limits<double>::quiet_NaN();
		n_dev_sq = std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * Copy constructor, creates a {@code Third_Moment} identical
	 * to the {@code original}.
	 *
	 * @param original the {@code Third_Moment} instance to copy
	 * @ if original is NULL
	 */
	Third_Moment(Third_Moment original)
	{
		super(original);
		this.m3 = original.m3;
		this.n_dev_sq = original.n_dev_sq;
	}

	/** {@inherit_doc} */
	//override
	public void increment(const double d)
	{
		if (n < 1)
		{
			m3 = m2 = m1 = 0.0;
		}

		double prev_m2 = m2;
		super.increment(d);
		n_dev_sq = n_dev * n_dev;
		double n0 = n;
		m3 = m3 - 3.0 * n_dev * prev_m2 + (n0 - 1) * (n0 - 2) * n_dev_sq * dev;
	}

	/** {@inherit_doc} */
	//override
	public double get_result()
	{
		return m3;
	}

	/** {@inherit_doc} */
	//override
	public void clear()
	{
		super.clear();
		m3 = std::numeric_limits<double>::quiet_NaN();
		n_dev_sq = std::numeric_limits<double>::quiet_NaN();
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
	public Third_Moment copy()
	{
		return Third_Moment(this);
	}
}
