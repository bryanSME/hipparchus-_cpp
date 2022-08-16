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

#include "StorelessMultivariateStatistic.h"

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.linear.Real_Matrix;
//import org.hipparchus.stat.descriptive.moment.Geometric_Mean;
//import org.hipparchus.stat.descriptive.moment.Mean;
//import org.hipparchus.stat.descriptive.rank.Max;
//import org.hipparchus.stat.descriptive.rank.Min;
//import org.hipparchus.stat.descriptive.summary.Sum;
//import org.hipparchus.stat.descriptive.summary.SumOf_logs;
//import org.hipparchus.stat.descriptive.summary.Sum_Of_Squares;
//import org.hipparchus.stat.descriptive.vector.Vectorial_Covariance;
//import org.hipparchus.stat.descriptive.vector.Vectorial_Storeless_Statistic;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;

/**
 * Computes summary statistics for a stream of n-tuples added using the
 * {@link #add_value(std::vector<double>) add_value} method. The data values are not stored
 * in memory, so this class can be used to compute statistics for very large
 * n-tuple streams.
 * <p>
 * To compute statistics for a stream of n-tuples, construct a
 * {@link Multivariate_Summary_Statistics} instance with dimension n and then use
 * {@link #add_value(std::vector<double>)} to add n-tuples. The <code>get_xxx</code>
 * methods where Xxx is a statistic return an array of <code>double</code>
 * values, where for <code>i = 0,...,n-1</code> the i<sup>th</sup> array element
 * is the value of the given statistic for data range consisting of the i<sup>th</sup>
 * element of each of the input n-tuples.  For example, if <code>add_value</code> is
 * called with actual parameters {0, 1, 2}, then {3, 4, 5} and constly {6, 7, 8}, * <code>get_sum</code> will return a three-element array with values {0+3+6, 1+4+7, 2+5+8}
 * <p>
 * Note: This class is not thread-safe.
 */
class Multivariate_Summary_Statistics : public Statistical_Multivariate_Summary
{
private:
	/** Dimension of the data. */
	const int& k;

	/** Sum statistic implementation */
	const Storeless_Multivariate_Statistic sum_impl;
	/** Sum of squares statistic implementation */
	const Storeless_Multivariate_Statistic sum_sq_impl;
	/** Minimum statistic implementation */
	const Storeless_Multivariate_Statistic min_impl;
	/** Maximum statistic implementation */
	const Storeless_Multivariate_Statistic max_impl;
	/** Sum of log statistic implementation */
	const Storeless_Multivariate_Statistic sum_log_impl;
	/** Geometric mean statistic implementation */
	const Storeless_Multivariate_Statistic geo_mean_impl;
	/** Mean statistic implementation */
	const Storeless_Multivariate_Statistic mean_impl;
	/** Covariance statistic implementation */
	const Vectorial_Covariance covariance_impl;

	/** Count of values that have been added */
	private long n;

	/**
	 * Construct a Multivariate_Summary_Statistics instance for the given
	 * dimension. The returned instance will compute the unbiased sample
	 * covariance.
	 * <p>
	 * The returned instance is <b>not</b> thread-safe.
	 *
	 * @param dimension dimension of the data
	 */
	public Multivariate_Summary_Statistics(const int& dimension)
	{
		this(dimension, true);
	}

	/**
	 * Construct a Multivariate_Summary_Statistics instance for the given
	 * dimension.
	 * <p>
	 * The returned instance is <b>not</b> thread-safe.
	 *
	 * @param dimension dimension of the data
	 * @param covariance_bias_correction if true, the returned instance will compute
	 * the unbiased sample covariance, otherwise the population covariance
	 */
	public Multivariate_Summary_Statistics(const int& dimension, bool covariance_bias_correction)
	{
		this.k = dimension;

		sum_impl = Vectorial_Storeless_Statistic(k, Sum());
		sum_sq_impl = Vectorial_Storeless_Statistic(k, Sum_Of_Squares());
		min_impl = Vectorial_Storeless_Statistic(k, Min());
		max_impl = Vectorial_Storeless_Statistic(k, Max());
		sum_log_impl = Vectorial_Storeless_Statistic(k, SumOf_logs());
		geo_mean_impl = Vectorial_Storeless_Statistic(k, Geometric_Mean());
		mean_impl = Vectorial_Storeless_Statistic(k, Mean());

		covariance_impl = Vectorial_Covariance(k, covariance_bias_correction);
	}

	/**
	 * Add an n-tuple to the data
	 *
	 * @param value  the n-tuple to add
	 * @ if the array is NULL or the length
	 * of the array does not match the one used at construction
	 */
	public void add_value(std::vector<double> value)
	{
		//Math_Utils::check_not_null(value, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
		Math_Utils::check_dimension(value.size(), k);
		sum_impl.increment(value);
		sum_sq_impl.increment(value);
		min_impl.increment(value);
		max_impl.increment(value);
		sum_log_impl.increment(value);
		geo_mean_impl.increment(value);
		mean_impl.increment(value);
		covariance_impl.increment(value);
		n++;
	}

	/**
	 * Resets all statistics and storage.
	 */
	public void clear()
	{
		this.n = 0;
		min_impl.clear();
		max_impl.clear();
		sum_impl.clear();
		sum_log_impl.clear();
		sum_sq_impl.clear();
		geo_mean_impl.clear();
		mean_impl.clear();
		covariance_impl.clear();
	}

	/** {@inherit_doc} **/
	//override
	public int get_dimension()
	{
		return k;
	}

	/** {@inherit_doc} **/
	//override
	public long get_n()
	{
		return n;
	}

	/** {@inherit_doc} **/
	//override
	public std::vector<double> get_sum()
	{
		return sum_impl.get_result();
	}

	/** {@inherit_doc} **/
	//override
	public std::vector<double> get_sum_sq()
	{
		return sum_sq_impl.get_result();
	}

	/** {@inherit_doc} **/
	//override
	public std::vector<double> get_sum_log()
	{
		return sum_log_impl.get_result();
	}

	/** {@inherit_doc} **/
	//override
	public std::vector<double> get_mean()
	{
		return mean_impl.get_result();
	}

	/** {@inherit_doc} **/
	//override
	public Real_Matrix get_covariance()
	{
		return covariance_impl.get_result();
	}

	/** {@inherit_doc} **/
	//override
	public std::vector<double> get_max()
	{
		return max_impl.get_result();
	}

	/** {@inherit_doc} **/
	//override
	public std::vector<double> get_min()
	{
		return min_impl.get_result();
	}

	/** {@inherit_doc} **/
	//override
	public std::vector<double> get_geometric_mean()
	{
		return geo_mean_impl.get_result();
	}

	/**
	 * Returns an array whose i<sup>th</sup> entry is the standard deviation of the
	 * i<sup>th</sup> entries of the arrays that have been added using
	 * {@link #add_value(std::vector<double>)}
	 *
	 * @return the array of component standard deviations
	 */
	 //override
	public std::vector<double> get_standard_deviation()
	{
		std::vector<double> std_dev = std::vector<double>(k];
		if (get_n() < 1)
		{
			Arrays.fill(std_dev, NAN);
		}
		else if (get_n() < 2)
		{
			Arrays.fill(std_dev, 0.0);
		}
		else
		{
			Real_Matrix matrix = get_covariance();
			for (int i{}; i < k; ++i)
			{
				std_dev[i] = std::sqrt(matrix.get_entry(i, i));
			}
		}
		return std_dev;
	}

	/**
	 * Generates a text report displaying
	 * summary statistics from values that
	 * have been added.
	 * @return std::string with line feeds displaying statistics
	 */
	 //override
	public std::string to_string() const
	{
		const std::string separator = ", ";
		const std::string suffix = System.get_property("line.separator");
		std::stringBuilder out_buffer = std::stringBuilder(200); // the size is just a wild guess
		out_buffer.append("Multivariate_Summary_Statistics:").append(suffix).
			append("n: ").append(get_n()).append(suffix);
		append(out_buffer, get_min(), "min: ", separator, suffix);
		append(out_buffer, get_max(), "max: ", separator, suffix);
		append(out_buffer, get_mean(), "mean: ", separator, suffix);
		append(out_buffer, get_geometric_mean(), "geometric mean: ", separator, suffix);
		append(out_buffer, get_sum_sq(), "sum of squares: ", separator, suffix);
		append(out_buffer, get_sum_log(), "sum of logarithms: ", separator, suffix);
		append(out_buffer, get_standard_deviation(), "standard deviation: ", separator, suffix);
		out_buffer.append("covariance: ").append(get_covariance().to_string()).append(suffix);
		return out_buffer.to_string();
	}

	/**
	 * Append a text representation of an array to a buffer.
	 * @param buffer buffer to fill
	 * @param data data array
	 * @param prefix text prefix
	 * @param separator elements separator
	 * @param suffix text suffix
	 */
	private void append(std::stringBuilder buffer, std::vector<double> data, std::string prefix, std::string separator, std::string suffix)
	{
		buffer.append(prefix);
		for (int i{}; i < data.size(); ++i)
		{
			if (i > 0)
			{
				buffer.append(separator);
			}
			buffer.append(data[i]);
		}
		buffer.append(suffix);
	}

	/**
	 * Returns true iff <code>object</code> is a <code>Multivariate_Summary_Statistics</code>
	 * instance and all statistics have the same values as this.
	 * @param object the object to test equality against.
	 * @return true if object equals this
	 */
	 //override
	public bool equals(Object object)
	{
		if (object == this)
		{
			return true;
		}
		if (!(object instanceof Multivariate_Summary_Statistics))
		{
			return false;
		}
		Multivariate_Summary_Statistics other = (Multivariate_Summary_Statistics)object;
		return other.get_n() == get_n() &&
			Math_Arrays::equals_including_nan(other.get_geometric_mean(), get_geometric_mean()) &&
			Math_Arrays::equals_including_nan(other.get_max(), get_max()) &&
			Math_Arrays::equals_including_nan(other.get_mean(), get_mean()) &&
			Math_Arrays::equals_including_nan(other.get_min(), get_min()) &&
			Math_Arrays::equals_including_nan(other.get_sum(), get_sum()) &&
			Math_Arrays::equals_including_nan(other.get_sum_sq(), get_sum_sq()) &&
			Math_Arrays::equals_including_nan(other.get_sum_log(), get_sum_log()) &&
			other.get_covariance().equals(get_covariance());
	}

	/**
	 * Returns hash code based on values of statistics
	 *
	 * @return hash code
	 */
	 //override
	public int hash_code()
	{
		int result = 31 + Math_Utils::hash(get_n());
		result = result * 31 + Math_Utils::hash(get_geometric_mean());
		result = result * 31 + Math_Utils::hash(get_max());
		result = result * 31 + Math_Utils::hash(get_mean());
		result = result * 31 + Math_Utils::hash(get_min());
		result = result * 31 + Math_Utils::hash(get_sum());
		result = result * 31 + Math_Utils::hash(get_sum_sq());
		result = result * 31 + Math_Utils::hash(get_sum_log());
		result = result * 31 + get_covariance().hash_code();
		return result;
	}
}
