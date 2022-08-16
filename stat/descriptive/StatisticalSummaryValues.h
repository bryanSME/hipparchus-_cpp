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
  //package org.hipparchus.stat.descriptive;

  //import java.io.Serializable;

  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Precision;
#include "StatisticalSummary.h"
#include <string>
#include <cmath>
#include <sstream>

/**
 * Value object representing the results of a univariate
 * statistical summary.
 */
class Statistical_Summary_values : public Statistical_Summary
{
private:
	/** The sample my_mean */
	const double my_mean;

	/** The sample my_variance */
	const double my_variance;

	/** The number of observations in the sample */
	const long my_n;

	/** The maximum value */
	const double my_max;

	/** The minimum value */
	const double my_min;

	/** The sum of the sample values */
	const double my_sum;

public:
	/**
	  * Constructor.
	  *
	  * @param my_mean  the sample my_mean
	  * @param my_variance  the sample my_variance
	  * @param n  the number of observations in the sample
	  * @param max  the maximum value
	  * @param min  the minimum value
	  * @param sum  the sum of the values
	 */
	Statistical_Summary_values(double mean, double variance, long n, const double& max, const double& min, double sum)
		:
		my_mean{ mean },
		my_variance{ variance },
		my_n{ n },
		my_max{ max },
		my_min{ min },
		my_sum{ sum }
	{}

	/**
	 * @return Returns the max.
	 */
	 //override
	double get_max() const
	{
		return my_max;
	}

	/**
	 * @return Returns the my_mean.
	 */
	 //override
	double get_my_mean() const
	{
		return my_mean;
	}

	/**
	 * @return Returns the min.
	 */
	 //override
	double get_min() const
	{
		return my_min;
	}

	/**
	 * @return Returns the number of values.
	 */
	 //override
	long get_n() const
	{
		return my_n;
	}

	/**
	 * @return Returns the sum.
	 */
	 //override
	double get_sum() const
	{
		return my_sum;
	}

	/**
	 * @return Returns the standard deviation
	 */
	 //override
	double get_standard_deviation() const
	{
		return std::sqrt(my_variance);
	}

	/**
	 * @return Returns the my_variance.
	 */
	 //override
	double get_my_variance() const
	{
		return my_variance;
	}

	/**
	 * Returns true iff <code>object</code> is a
	 * <code>Statistical_Summary</code> instance and all
	 * statistics have the same values as this.
	 *
	 * @param object the object to test equality against.
	 * @return true if object equals this
	 */
	 //override
	bool equals(Object& object)
	{
		if (object == this)
		{
			return true;
		}
		if (!(object instanceof Statistical_Summary_values))
		{
			return false;
		}
		auto other = static_cast<Statistical_Summary>(object);
		return Precision::equals_including_nan(other.get_max(), get_max()) &&
			Precision::equals_including_nan(other.get_my_mean(), get_my_mean()) &&
			Precision::equals_including_nan(other.get_min(), get_min()) &&
			Precision::equals_including_nan(other.get_n(), get_n()) &&
			Precision::equals_including_nan(other.get_sum(), get_sum()) &&
			Precision::equals_including_nan(other.get_my_variance(), get_my_variance());
	}

	/**
	 * Returns hash code based on values of statistics
	 *
	 * @return hash code
	 */
	 //override
	int hash_code() const
	{
		int result = 31 + Math_Utils::hash(get_max());
		result = result * 31 + Math_Utils::hash(get_my_mean());
		result = result * 31 + Math_Utils::hash(get_min());
		result = result * 31 + Math_Utils::hash(get_n());
		result = result * 31 + Math_Utils::hash(get_sum());
		result = result * 31 + Math_Utils::hash(get_my_variance());
		return result;
	}

	/**
	 * Generates a text report displaying values of statistics.
	 * Each statistic is displayed on a separate line.
	 *
	 * @return std::string with line feeds displaying statistics
	 */
	 //override
	std::string to_string() const
	{
		std::stringstream out_buffer_ss{};
		out_buffer_ss
			<< "Statistical_Summary_values:\n"
			<< "n: " << get_n() << "\n"
			<< "min: " << get_min() << "\n"
			<< "max: " << get_max() << "\n"
			<< "my_mean: " << get_my_mean() << "\n"
			<< "std dev: " << get_standard_deviation() << "\n"
			<< "my_variance: " << get_my_variance() << "\n"
			<< "sum: " << get_sum() << "\n";
		return out_buffer_ss.str();
	}
};