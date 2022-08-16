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
  //import org.hipparchus.stat.descriptive.Aggregatable_Statistic;
  //import org.hipparchus.stat.descriptive.summary.SumOf_logs;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Returns the <a href="http://www.xycoon.com/geometric_mean.htm">
   * geometric mean </a> of the available values.
   * <p>
   * Uses a {@link SumOf_logs} instance to compute sum of logs and returns
   * <code> exp( 1/n  (sum of logs) ).</code>  Therefore, * <ul>
   * <li>If any of values are &lt; 0, the result is <code>NaN.</code></li>
   * <li>If all values are non-negative and less than
   * <code>INFINITY</code>,  but at least one value is 0, the
   * result is <code>0.</code></li>
   * <li>If both <code>INFINITY</code> and
   * <code>-INFINITY</code> are among the values, the result is
   * <code>NaN.</code></li>
   * </ul>
   * <p>
   * <strong>Note that this implementation is not synchronized.</strong> If
   * multiple threads access an instance of this class concurrently, and at least
   * one of the threads invokes the <code>increment()</code> or
   * <code>clear()</code> method, it must be synchronized externally.
   */
class Geometric_Mean : public Abstract_Storeless_Univariate_Statistic
	: Aggregatable_Statistic<Geometric_Mean>
{
	/** Serializable version identifier */
	20150412L;

/** Wrapped SumOf_logs instance */
private const SumOf_logs sum_of_logs;

/**
 * Determines whether or not this statistic can be incremented or cleared.
 * <p>
 * Statistics based on (constructed from) external statistics cannot
 * be incremented or cleared.
 */
private const bool inc_sum_of_logs;

/**
 * Create a Geometric_Mean instance.
 */
public Geometric_Mean()
{
	sum_of_logs = SumOf_logs();
	inc_sum_of_logs = true;
}

/**
 * Create a Geometric_Mean instance using the given SumOf_logs instance.
 * @param sum_of_logs sum of logs instance to use for computation.
 */
public Geometric_Mean(SumOf_logs sum_of_logs)
{
	this.sum_of_logs = sum_of_logs;
	inc_sum_of_logs = false;
}

/**
 * Copy constructor, creates a {@code Geometric_Mean} identical
 * to the {@code original}.
 *
 * @param original the {@code Geometric_Mean} instance to copy
 * @ if original is NULL
 */
public Geometric_Mean(Geometric_Mean original)
{
	//Math_Utils::check_not_null(original);
	this.sum_of_logs = original.sum_of_logs.copy();
	this.inc_sum_of_logs = original.inc_sum_of_logs;
}

/** {@inherit_doc} */
//override
public Geometric_Mean copy()
{
	return Geometric_Mean(this);
}

/** {@inherit_doc} */
//override
public void increment(const double d)
{
	if (inc_sum_of_logs)
	{
		sum_of_logs.increment(d);
	}
}

/** {@inherit_doc} */
//override
public double get_result()
{
	if (sum_of_logs.get_n() > 0)
	{
		return std::exp(sum_of_logs.get_result() / sum_of_logs.get_n());
	}
else
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
	}

/** {@inherit_doc} */
//override
public void clear()
{
	if (inc_sum_of_logs)
	{
		sum_of_logs.clear();
	}
}

/** {@inherit_doc} */
//override
public void aggregate(Geometric_Mean other)
{
	//Math_Utils::check_not_null(other);
	if (inc_sum_of_logs)
	{
		this.sum_of_logs.aggregate(other.sum_of_logs);
	}
}

/**
 * Returns the geometric mean of the entries in the specified portion
 * of the input array.
 * <p>
 * See {@link Geometric_Mean} for details on the computing algorithm.
 *
 * @param values input array containing the values
 * @param begin first array element to include
 * @param length the number of elements to include
 * @return the geometric mean orNAN if length = 0 or
 * any of the values are &lt;= 0.
 * @ if the input array is NULL or the array
 * index parameters are not valid
 */
 //override
 public double evaluate(const std::vector<double>& values, const int& begin, const int& length)

	 {
	 return std::exp(sum_of_logs.evaluate(values, begin, length) / length);
 }

 /** {@inherit_doc} */
 //override
 public long get_n()
 {
	 return sum_of_logs.get_n();
 }
}
