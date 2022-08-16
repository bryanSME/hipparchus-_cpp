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
  //import org.hipparchus.stat.descriptive.summary.Sum;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <numbers>
#include<vector>
#include "../AbstractStorelessUnivariateStatistic.h"
#include "../AggregatableStatistic.hpp"
#include "../WeightedEvaluation.h"
#include "FirstMoment.h"
#include "../../../core/util/MathArrays.h"
#include "../moment/Mean.h"

/**
 * Computes the arithmetic mean of a set of values. Uses the definitional
 * formula:
 * <p>
 * mean = sum(x_i) / n
 * <p>
 * where <code>n</code> is the number of observations.
 * <p>
 * When {@link #incrementstatic_cast<double>(} is used to add data incrementally from a
 * stream of (unstored) values, the value of the statistic that
 * {@link #get_result()} returns is computed using the following recursive
 * updating algorithm:
 * <ol>
 * <li>Initialize <code>m = </code> the first value</li>
 * <li>For each additional value, update using <br>
 *   <code>m = m + (new value - m) / (number of observations)</code></li>
 * </ol>
 * <p>
 * If {@link #evaluate(std::vector<double>)} is used to compute the mean of an array
 * of stored values, a two-pass, corrected algorithm is used, starting with
 * the definitional formula computed using the array of stored values and then
 * correcting this by adding the mean deviation of the data values from the
 * arithmetic mean. See, e.g. "Comparison of Several Algorithms for Computing
 * Sample Means and Variances," Robert F. Ling, Journal of the American
 * Statistical Association, Vol. 69, No. 348 (Dec., 1974), pp. 859-866.
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
class Mean : public Abstract_Storeless_Univariate_Statistic, public Aggregatable_Statistic<Mean>, public Weighted_Evaluation
{
protected:
	/** First moment on which this statistic is based. */
	const First_Moment moment;

	/**
	 * Determines whether or not this statistic can be incremented or cleared.
	 * <p>
	 * Statistics based on (constructed from) external moments cannot
	 * be incremented or cleared.
	 */
	const bool my_inc_moment;

public:
	/** Constructs a Mean. */
	Mean()
		:
		moment{ First_Moment() },
		my_inc_moment{ true }
	{
	}

	/**
	 * Constructs a Mean with an External Moment.
	 *
	 * @param m1 the moment
	 */
	Mean(const First_Moment m1)
	{
		my_moment = m1;
		my_inc_moment = false;
	}

	/**
	 * Copy constructor, creates a {@code Mean} identical
	 * to the {@code original}.
	 *
	 * @param original the {@code Mean} instance to copy
	 * @ if original is NULL
	 */
	Mean(Mean original)
	{
		//Math_Utils::check_not_null(original);
		my_moment = original.moment.copy();
		my_inc_moment = original.inc_moment;
	}

	/**
	 * {@inherit_doc}
	 * <p>
	 * Note that when {@link #Mean(First_Moment)} is used to
	 * create a Mean, this method does nothing. In that case, the
	 * First_Moment should be incremented directly.
	 */
	 //override
	void increment(const double& d)
	{
		if (my_inc_moment)
		{
			my_moment.increment(d);
		}
	}

	/** {@inherit_doc} */
	//override
	void clear()
	{
		if (my_inc_moment)
		{
			my_moment.clear();
		}
	}

	/** {@inherit_doc} */
	//override
	double get_result()
	{
		return my_moment.m1;
	}

	/** {@inherit_doc} */
	//override
	long get_n()
	{
		return my_moment.get_n();
	}

	/** {@inherit_doc} */
	//override
	void aggregate(Mean other)
	{
		//Math_Utils::check_not_null(other);
		if (my_inc_moment)
		{
			my_moment.aggregate(other.moment);
		}
	}

	/**
	 * Returns the arithmetic mean of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the mean of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	 //override
	double evaluate(const std::vector<double>& values, const int& begin, const int& length)
	{
		if (Math_Arrays::verify_values(values, begin, length))
		{
			double sample_size = length;

			// Compute initial estimate using definitional formula
			double xbar = Stat_Utils.sum(values, begin, length) / sample_size;

			// Compute correction factor in second pass
			double correction = 0;
			for (int i{ begin }; i < begin + length; i++)
			{
				correction += values[i] - xbar;
			}
			return xbar + (correction / sample_size);
		}
		return std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * Returns the weighted arithmetic mean of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if either array is NULL.
	 * <p>
	 * See {@link Mean} for details on the computing algorithm. The two-pass algorithm
	 * described above is used here, with weights applied in computing both the original
	 * estimate and the correction factor.
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
	 *
	 * @param values the input array
	 * @param weights the weights array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the mean of the values orNAN if length = 0
	 * @ if the parameters are not valid
	 */
	 //override
	double evaluate(const std::vector<double>& values, const std::vector<double> weights, const int& begin, const int& length)
	{
		if (Math_Arrays::verify_values(values, weights, begin, length))
		{
			Sum sum = Sum();

			// Compute initial estimate using definitional formula
			double sumw = sum.evaluate(weights, begin, length);
			double xbarw = sum.evaluate(values, weights, begin, length) / sumw;

			// Compute correction factor in second pass
			double correction = 0;
			for (int i{ begin }; i < begin + length; i++)
			{
				correction += weights[i] * (values[i] - xbarw);
			}
			return xbarw + (correction / sumw);
		}
		return std::numeric_limits<double>::quiet_NaN();
	}

	/** {@inherit_doc} */
	//override
	Mean copy()
	{
		return Mean(*this);
	}
};