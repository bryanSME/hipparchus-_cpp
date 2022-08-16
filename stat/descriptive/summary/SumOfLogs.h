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
  //package org.hipparchus.stat.descriptive.summary;

  //import java.io.Serializable;

  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.;
  //import org.hipparchus.stat.descriptive.Abstract_Storeless_Univariate_Statistic;
  //import org.hipparchus.stat.descriptive.Aggregatable_Statistic;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include "../../../core/util/MathArrays.h"
#include "../../descriptive/AbstractStorelessUnivariateStatistic.h"
#include "../../descriptive/AggregatableStatistic.hpp"

/**
 * Returns the sum of the natural logs for this collection of values.
 * <p>
 * Uses {@link org.hipparchus.util.FastMath#logstatic_cast<double>(} to compute the logs.
 * Therefore, * <ul>
 * <li>If any of values are &lt; 0, the result is <code>NaN.</code></li>
 * <li>If all values are non-negative and less than
 * <code>INFINITY</code>,  but at least one value is 0, the
 * result is <code>-INFINITY.</code></li>
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
class SumOf_logs : public Abstract_Storeless_Univariate_Statistic //, public Aggregatable_Statistic<SumOf_logs>
{
private:

	/** Number of values that have been added */
	int my_n;

	/** The currently running value */
	double my_value;

public:
	/**
	 * Create a SumOf_logs instance.
	 */
	SumOf_logs() : my_value{}, my_n{} {};

	/**
	 * Copy constructor, creates a {@code SumOf_logs} identical
	 * to the {@code original}.
	 *
	 * @param original the {@code SumOf_logs} instance to copy
	 * @ if original is NULL
	 */
	SumOf_logs(const SumOf_logs& original) : my_n{ original.get_result() }, my_value{ original.get_result() }
	{
		//Math_Utils::check_not_null(original);
	}

	/** {@inherit_doc} */
	//override
	void increment(const double d)
	{
		my_value += std::log(d);
		my_n++;
	}

	/** {@inherit_doc} */
	//override
	double get_result() const
	{
		return my_value;
	}

	/** {@inherit_doc} */
	//override
	long get_n() const
	{
		return my_n;
	}

	/** {@inherit_doc} */
	void clear() override
	{
		my_value = 0;
		my_n = 0;
	}

	/** {@inherit_doc} */
	//override
	void aggregate(const SumOf_logs& other)
	{
		//Math_Utils::check_not_null(other);
		if (other.get_n() > 0)
		{
			my_n += other.get_n();
			my_value += other.get_result();
		}
	}

	/**
	 * Returns the sum of the natural logs of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the sum of the natural logs of the values or 0 if
	 * length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	 //override
	double evaluate(const std::vector<double>& values, const int& begin, const int& length)
	{
		double sum_log = std::numeric_limits<double>::quiet_NaN();
		if (Math_Arrays::verify_values(values, begin, length, true))
		{
			sum_log = 0.0;
			for (int i{ begin }; i < begin + length; i++)
			{
				sum_log += std::log(values[i]);
			}
		}
		return sum_log;
	}

	/** {@inherit_doc} */
	//override
	/*SumOf_logs copy()
	{
		return SumOf_logs(this);
	}*/
};