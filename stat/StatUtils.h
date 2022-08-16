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
  //package org.hipparchus.stat;

  //import java.util.Arrays;
  //import java.util.List;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.stat.descriptive.Descriptive_Statistics;
  //import org.hipparchus.stat.descriptive.Univariate_Statistic;
  //import org.hipparchus.stat.descriptive.moment.Geometric_Mean;
  //import org.hipparchus.stat.descriptive.moment.Mean;
  //import org.hipparchus.stat.descriptive.moment.Variance;
  //import org.hipparchus.stat.descriptive.rank.Max;
  //import org.hipparchus.stat.descriptive.rank.Min;
  //import org.hipparchus.stat.descriptive.rank.Percentile;
  //import org.hipparchus.stat.descriptive.summary.Product;
  //import org.hipparchus.stat.descriptive.summary.Sum;
  //import org.hipparchus.stat.descriptive.summary.SumOf_logs;
  //import org.hipparchus.stat.descriptive.summary.Sum_Of_Squares;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <vector>
#include "../core/util/MathArrays.h"
#include "../core/util/MathUtils.h"
#include "../core/exception/LocalizedCoreFormats.h"
#include "descriptive/UnivariateStatistic.h"
#include "descriptive/rank/Max.h"
#include "descriptive/rank/Min.h"
//#include "descriptive/rank/"

/**
 * Stat_Utils provides static methods for computing statistics based on data
 * stored in std::vector<double> arrays.
 */
class Stat_Utils
{
private:
	/** sum */
	static const Univariate_Statistic SUM = Sum();

	/** sum_sq */
	static const Univariate_Statistic SUM_OF_SQUARES = Sum_Of_Squares();

	/** prod */
	static const Univariate_Statistic PRODUCT = Product();

	/** sum_log */
	static const Univariate_Statistic SUM_OF_LOGS = SumOf_logs();

	/** min */
	static const Univariate_Statistic MIN = Min();

	/** max */
	static const Univariate_Statistic MAX = Max();

	/** mean */
	static const Univariate_Statistic MEAN = Mean();

	/** variance */
	static const Variance VARIANCE = Variance();

	/** percentile */
	static const Percentile PERCENTILE = Percentile();

	/** geometric mean */
	static const Geometric_Mean GEOMETRIC_MEAN = Geometric_Mean();

	/**
	 * Private helper method.
	 * Assumes parameters have been validated.
	 * @param values input data
	 * @param begin index (0-based) of the first array element to include
	 * @param length the number of elements to include
	 * @return array of array of the most frequently occurring element(s) sorted in ascending order.
	 */
	static std::vector<double> get_mode(std::vector<double> values, const int& begin, const int& length)
	{
		// Add the values to the frequency table
		Frequency<Double> freq = Frequency<>();

		Arrays.stream(values, begin, begin + length)
			.filter(d -> !std::isnan(d))
			.for_each(freq::add_value);

		List<Double> list = freq.get_mode();
		// Convert the list to an array of primitive double
		return list.stream()
			.map_to_double(Double::double_value)
			.to_array();
	}

	/**
	 * Private Constructor
	 */
	Stat_Utils()
	{
	}

public:
	/**
	 * Returns the sum of the values in the input array, or
	 * <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the input array is NULL.
	 *
	 * @param values  array of values to sum
	 * @return the sum of the values or <code>Double.NaN</code> if the array is empty
	 * @ if the array is NULL
	 */
	static double sum(const double... values)
	{
		return SUM.evaluate(values);
	}

	/**
	 * Returns the sum of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the sum of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	static double sum(const std::vector<double>& values, const int& begin, const int& length)
	{
		return SUM.evaluate(values, begin, length);
	}

	/**
	 * Returns the sum of the squares of the entries in the input array, or
	 * <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 *
	 * @param values  input array
	 * @return the sum of the squared values or <code>Double.NaN</code> if the array is empty
	 * @ if the array is NULL
	 */
	static double sum_sq(const double... values)
	{
		return SUM_OF_SQUARES.evaluate(values);
	}

	/**
	 * Returns the sum of the squares of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the sum of the squares of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	static double sum_sq(const std::vector<double>& values, const int& begin, const int& length)
	{
		return SUM_OF_SQUARES.evaluate(values, begin, length);
	}

	/**
	 * Returns the product of the entries in the input array, or
	 * <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 *
	 * @param values the input array
	 * @return the product of the values orNAN if the array is empty
	 * @ if the array is NULL
	 */
	static double product(const double... values)
	{
		return PRODUCT.evaluate(values);
	}

	/**
	 * Returns the product of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the product of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	static double product(const std::vector<double>& values, const int& begin, const int& length)
	{
		return PRODUCT.evaluate(values, begin, length);
	}

	/**
	 * Returns the sum of the natural logs of the entries in the input array, or
	 * <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.summary.SumOf_logs}.
	 *
	 * @param values the input array
	 * @return the sum of the natural logs of the values orNAN if the array is empty
	 * @ if the array is NULL
	 */
	static double sum_log(const double... values)
	{
		return SUM_OF_LOGS.evaluate(values);
	}

	/**
	 * Returns the sum of the natural logs of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.summary.SumOf_logs}.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the sum of the natural logs of the values orNAN if
	 * length = 0
	 * @ if the array is NULL or the array index
	 * parameters are not valid
	 */
	static double sum_log(const std::vector<double>& values, const int& begin, const int& length)
	{
		return SUM_OF_LOGS.evaluate(values, begin, length);
	}

	/**
	 * Returns the arithmetic mean of the entries in the input array, or
	 * <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Mean} for
	 * details on the computing algorithm.
	 *
	 * @param values the input array
	 * @return the mean of the values orNAN if the array is empty
	 * @ if the array is NULL
	 */
	static double mean(const double... values)
	{
		return MEAN.evaluate(values);
	}

	/**
	 * Returns the arithmetic mean of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Mean Mean} for
	 * details on the computing algorithm.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the mean of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 * parameters are not valid
	 */
	static double mean(const std::vector<double>& values, const int& begin, const int& length)
	{
		return MEAN.evaluate(values, begin, length);
	}

	/**
	 * Returns the geometric mean of the entries in the input array, or
	 * <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Geometric_Mean Geometric_Mean}
	 * for details on the computing algorithm.
	 *
	 * @param values the input array
	 * @return the geometric mean of the values orNAN if the array is empty
	 * @ if the array is NULL
	 */
	static double geometric_mean(const double... values)
	{
		return GEOMETRIC_MEAN.evaluate(values);
	}

	/**
	 * Returns the geometric mean of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 * <p>
	 * Throws <code>Illegal_Argument_Exception</code> if the array is NULL.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Geometric_Mean Geometric_Mean}
	 * for details on the computing algorithm.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the geometric mean of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	static double geometric_mean(const std::vector<double>& values, const int& begin, const int& length)
	{
		return GEOMETRIC_MEAN.evaluate(values, begin, length);
	}

	/**
	 * Returns the variance of the entries in the input array, or
	 * <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * This method returns the bias-corrected sample variance (using {@code n - 1} in
	 * the denominator). Use {@link #population_variance(std::vector<double>)} for the non-bias-corrected
	 * population variance.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Variance Variance} for
	 * details on the computing algorithm.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * Throws <code></code> if the array is NULL.
	 *
	 * @param values the input array
	 * @return the variance of the values orNAN if the array is empty
	 * @ if the array is NULL
	 */
	static double variance(const double... values)
	{
		return VARIANCE.evaluate(values);
	}

	/**
	 * Returns the variance of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 * <p>
	 * This method returns the bias-corrected sample variance (using {@code n - 1} in
	 * the denominator). Use {@link #population_variance(std::vector<double>, int, int)} for the non-bias-corrected
	 * population variance.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Variance Variance} for
	 * details on the computing algorithm.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * Throws <code></code> if the array is NULL or the
	 * array index parameters are not valid.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the variance of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	static double variance(const std::vector<double>& values, const int& begin, const int& length)
	{
		return VARIANCE.evaluate(values, begin, length);
	}

	/**
	 * Returns the variance of the entries in the specified portion of
	 * the input array, using the precomputed mean value.  Returns
	 * <code>Double.NaN</code> if the designated subarray is empty.
	 * <p>
	 * This method returns the bias-corrected sample variance (using {@code n - 1} in
	 * the denominator). Use {@link #population_variance(std::vector<double>, double, int, int)} for
	 * the non-bias-corrected population variance.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Variance Variance} for
	 * details on the computing algorithm.
	 * <p>
	 * The formula used assumes that the supplied mean value is the arithmetic
	 * mean of the sample data, not a known population parameter.  This method
	 * is supplied only to save computation when the mean has already been
	 * computed.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * Throws <code></code> if the array is NULL or the
	 * array index parameters are not valid.
	 *
	 * @param values the input array
	 * @param mean the precomputed mean value
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the variance of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	static double variance(const std::vector<double>& values, const double mean, const int& begin, const int& length)
	{
		return VARIANCE.evaluate(values, mean, begin, length);
	}

	/**
	 * Returns the variance of the entries in the input array, using the
	 * precomputed mean value.  Returns <code>Double.NaN</code> if the array
	 * is empty.
	 * <p>
	 * This method returns the bias-corrected sample variance (using {@code n - 1} in
	 * the denominator).  Use {@link #population_variance(std::vector<double>, double)} for the
	 * non-bias-corrected population variance.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Variance Variance} for
	 * details on the computing algorithm.
	 * <p>
	 * The formula used assumes that the supplied mean value is the arithmetic
	 * mean of the sample data, not a known population parameter.  This method
	 * is supplied only to save computation when the mean has already been
	 * computed.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * Throws <code></code> if the array is NULL.
	 *
	 * @param values the input array
	 * @param mean the precomputed mean value
	 * @return the variance of the values orNAN if the array is empty
	 * @ if the array is NULL
	 */
	static double variance(const std::vector<double>& values, const double mean)
	{
		return VARIANCE.evaluate(values, mean);
	}

	/**
	 * Returns the <a href="http://en.wikibooks.org/wiki/Statistics/Summary/Variance">
	 * population variance</a> of the entries in the input array, or
	 * <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Variance Variance} for
	 * details on the formula and computing algorithm.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * Throws <code></code> if the array is NULL.
	 *
	 * @param values the input array
	 * @return the population variance of the values orNAN if the array is empty
	 * @ if the array is NULL
	 */
	static double population_variance(const double... values)
	{
		return Variance(false).evaluate(values);
	}

	/**
	 * Returns the <a href="http://en.wikibooks.org/wiki/Statistics/Summary/Variance">
	 * population variance</a> of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Variance Variance} for
	 * details on the computing algorithm.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * Throws <code></code> if the array is NULL or the
	 * array index parameters are not valid.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the population variance of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	static double population_variance(const std::vector<double>& values, const int& begin, const int& length)
	{
		return Variance(false).evaluate(values, begin, length);
	}

	/**
	 * Returns the <a href="http://en.wikibooks.org/wiki/Statistics/Summary/Variance">
	 * population variance</a> of the entries in the specified portion of
	 * the input array, using the precomputed mean value.  Returns
	 * <code>Double.NaN</code> if the designated subarray is empty.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Variance Variance} for
	 * details on the computing algorithm.
	 * <p>
	 * The formula used assumes that the supplied mean value is the arithmetic
	 * mean of the sample data, not a known population parameter.  This method
	 * is supplied only to save computation when the mean has already been
	 * computed.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * Throws <code></code> if the array is NULL or the
	 * array index parameters are not valid.
	 *
	 * @param values the input array
	 * @param mean the precomputed mean value
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the population variance of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	static double population_variance(const std::vector<double>& values, const double mean, const int& begin, const int& length)
	{
		return Variance(false).evaluate(values, mean, begin, length);
	}

	/**
	 * Returns the <a href="http://en.wikibooks.org/wiki/Statistics/Summary/Variance">
	 * population variance</a> of the entries in the input array, using the precomputed
	 * mean value. Returns <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.moment.Variance Variance} for
	 * details on the computing algorithm.
	 * <p>
	 * The formula used assumes that the supplied mean value is the arithmetic
	 * mean of the sample data, not a known population parameter. This method is
	 * supplied only to save computation when the mean has already been computed.
	 * <p>
	 * Returns 0 for a single-value (i.e. length = 1) sample.
	 * <p>
	 * Throws <code></code> if the array is NULL.
	 *
	 * @param values the input array
	 * @param mean the precomputed mean value
	 * @return the population variance of the values orNAN if the array is empty
	 * @ if the array is NULL
	 */
	static double population_variance(const std::vector<double>& values, const double mean)
	{
		return Variance(false).evaluate(values, mean);
	}

	/**
	 * Returns the maximum of the entries in the input array, or
	 * <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * Throws <code></code> if the array is NULL.
	 * <p>
	 * <ul>
	 * <li>The result is <code>NaN</code> iff all values are <code>NaN</code>
	 * (i.e. <code>NaN</code> values have no impact on the value of the statistic).</li>
	 * <li>If any of the values equals <code>INFINITY</code>, * the result is <code>INFINITY.</code></li>
	 * </ul>
	 *
	 * @param values the input array
	 * @return the maximum of the values orNAN if the array is empty
	 * @ if the array is NULL
	 */
	static double max(const double... values)
	{
		return MAX.evaluate(values);
	}

	/**
	 * Returns the maximum of the entries in the specified portion of the input array, * or <code>Double.NaN</code> if the designated subarray is empty.
	 * <p>
	 * Throws <code></code> if the array is NULL or
	 * the array index parameters are not valid.
	 * <p>
	 * <ul>
	 * <li>The result is <code>NaN</code> iff all values are <code>NaN</code>
	 * (i.e. <code>NaN</code> values have no impact on the value of the statistic).</li>
	 * <li>If any of the values equals <code>INFINITY</code>, * the result is <code>INFINITY.</code></li>
	 * </ul>
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the maximum of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	static double max(const std::vector<double>& values, const int& begin, const int& length)
	{
		return MAX.evaluate(values, begin, length);
	}

	/**
	 * Returns the minimum of the entries in the input array, or
	 * <code>Double.NaN</code> if the array is empty.
	 * <p>
	 * Throws <code></code> if the array is NULL.
	 * <p>
	 * <ul>
	 * <li>The result is <code>NaN</code> iff all values are <code>NaN</code>
	 * (i.e. <code>NaN</code> values have no impact on the value of the statistic).</li>
	 * <li>If any of the values equals <code>-INFINITY</code>, * the result is <code>-INFINITY.</code></li>
	 * </ul>
	 *
	 * @param values the input array
	 * @return the minimum of the values orNAN if the array is empty
	 * @ if the array is NULL
	 */
	static double min(const double... values)
	{
		return MIN.evaluate(values);
	}

	/**
	 * Returns the minimum of the entries in the specified portion of the input array, * or <code>Double.NaN</code> if the designated subarray is empty.
	 * <p>
	 * Throws <code></code> if the array is NULL or
	 * the array index parameters are not valid.
	 * <p>
	 * <ul>
	 * <li>The result is <code>NaN</code> iff all values are <code>NaN</code>
	 * (i.e. <code>NaN</code> values have no impact on the value of the statistic).</li>
	 * <li>If any of the values equals <code>-INFINITY</code>, * the result is <code>-INFINITY.</code></li>
	 * </ul>
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the minimum of the values orNAN if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	static double min(const std::vector<double>& values, const int& begin, const int& length)
	{
		return MIN.evaluate(values, begin, length);
	}

	/**
	 * Returns an estimate of the <code>p</code>th percentile of the values
	 * in the <code>values</code> array.
	 * <p>
	 * <ul>
	 * <li>Returns <code>Double.NaN</code> if <code>values</code> has length
	 *  <code>0</code></li>
	 * <li>Returns (for any value of <code>p</code>) <code>values[0]</code>
	 *  if <code>values</code> has length <code>1</code></li>
	 * <li>Throws <code>Illegal_Argument_Exception</code> if <code>values</code>
	 *  is NULL  or p is not a valid quantile value (p must be greater than 0
	 *  and less than or equal to 100)</li>
	 * </ul>
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.rank.Percentile Percentile}
	 * for a description of the percentile estimation algorithm used.
	 *
	 * @param values input array of values
	 * @param p the percentile value to compute
	 * @return the percentile value orNAN if the array is empty
	 * @ if <code>values</code> is NULL or p is invalid
	 */
	static double percentile(const std::vector<double>& values, const double p)
	{
		return PERCENTILE.evaluate(values, p);
	}

	/**
	 * Returns an estimate of the <code>p</code>th percentile of the values
	 * in the <code>values</code> array, starting with the element in (0-based)
	 * position <code>begin</code> in the array and including <code>length</code>
	 * values.
	 * <p>
	 * <ul>
	 * <li>Returns <code>Double.NaN</code> if <code>length = 0</code></li>
	 * <li>Returns (for any value of <code>p</code>) <code>values[begin]</code>
	 *  if <code>length = 1 </code></li>
	 * <li>Throws <code></code> if <code>values</code>
	 *  is NULL, <code>begin</code> or <code>length</code> is invalid, or
	 *  <code>p</code> is not a valid quantile value (p must be greater than 0
	 *  and less than or equal to 100)</li>
	 * </ul>
	 * <p>
	 * See {@link org.hipparchus.stat.descriptive.rank.Percentile Percentile}
	 * for a description of the percentile estimation algorithm used.
	 *
	 * @param values array of input values
	 * @param p the percentile to compute
	 * @param begin the first (0-based) element to include in the computation
	 * @param length the number of array elements to include
	 * @return the percentile value
	 * @ if the parameters are not valid or the input array is NULL
	 */
	static double percentile(const std::vector<double>& values, const int& begin, const int& length, const double p)
	{
		return PERCENTILE.evaluate(values, begin, length, p);
	}

	/**
	 * Returns the sum of the (signed) differences between corresponding elements of the
	 * input arrays -- i.e., sum(sample1[i] - sample2[i]).
	 *
	 * @param sample1  the first array
	 * @param sample2  the second array
	 * @return sum of paired differences
	 * @ if the arrays do not have the same (positive) length.
	 * @ if the sample arrays are empty.
	 */
	static double sum_difference(const std::vector<double>& sample1, const std::vector<double>& sample2)
	{
		int n = sample1.size();
		Math_Arrays::check_equal_length(sample1, sample2);
		if (n <= 0)
		{
			throw (Localized_Core_Formats_Types::INSUFFICIENT_DIMENSION);
		}
		double result{};
		for (int i{}; i < n; i++)
		{
			result += sample1[i] - sample2[i];
		}
		return result;
	}

	/**
	 * Returns the mean of the (signed) differences between corresponding elements of the
	 * input arrays -- i.e., sum(sample1[i] - sample2[i]) / sample1.size().
	 *
	 * @param sample1  the first array
	 * @param sample2  the second array
	 * @return mean of paired differences
	 * @ if the arrays do not have the same (positive) length.
	 * @ if the sample arrays are empty.
	 */
	static double mean_difference(const std::vector<double> sample1, const std::vector<double> sample2)
	{
		return sum_difference(sample1, sample2) / sample1.size();
	}

	/**
	 * Returns the variance of the (signed) differences between corresponding elements of the
	 * input arrays -- i.e., var(sample1[i] - sample2[i]).
	 *
	 * @param sample1  the first array
	 * @param sample2  the second array
	 * @param mean_difference   the mean difference between corresponding entries
	 * @return variance of paired differences
	 * @ if the arrays do not have the same length.
	 * @ if the arrays length is less than 2.
	 * @see #mean_difference(std::vector<double>,std::vector<double>)
	 */
	static double variance_difference(const std::vector<double> sample1, const std::vector<double> sample2, double mean_difference)
	{
		double sum1{};
		double sum2{};
		int n = sample1.size();
		Math_Arrays::check_equal_length(sample1, sample2);
		if (n < 2)
		{
			throw (Localized_Core_Format_Types::NUMBER_TOO_SMALL, n, 2);
		}
		for (int i{}; i < n; i++)
		{
			const double diff = sample1[i] - sample2[i];
			sum1 += (diff - mean_difference) * (diff - mean_difference);
			sum2 += diff - mean_difference;
		}
		return (sum1 - (sum2 * sum2 / n)) / (n - 1);
	}

	/**
	 * Normalize (standardize) the sample, so it is has a mean of 0 and a standard deviation of 1.
	 *
	 * @param sample Sample to normalize.
	 * @return normalized (standardized) sample.
	 */
	static std::vector<double> normalize(const double... sample)
	{
		Descriptive_Statistics stats = Descriptive_Statistics();

		// Add the data from the series to stats
		for (int i{}; i < sample.size(); i++)
		{
			stats.add_value(sample[i]);
		}

		// Compute mean and standard deviation
		double mean = stats.get_mean();
		double standard_deviation = stats.get_standard_deviation();

		// initialize the standardized_sample, which has the same length as the sample
		std::vector<double> standardized_sample = std::vector<double>(sample.size()];

		for (int i{}; i < sample.size(); i++)
		{
			// z = (x- mean)/standard_deviation
			standardized_sample[i] = (sample[i] - mean) / standard_deviation;
		}
		return standardized_sample;
	}

	/**
	 * Returns the sample mode(s).
	 * <p>
	 * The mode is the most frequently occurring value in the sample.
	 * If there is a unique value with maximum frequency, this value is returned
	 * as the only element of the output array. Otherwise, the returned array
	 * contains the maximum frequency elements in increasing order.
	 * <p>
	 * For example, if {@code sample} is {0, 12, 5, 6, 0, 13, 5, 17}, * the returned array will have length two, with 0 in the first element and
	 * 5 in the second.
	 * <p>
	 * NaN values are ignored when computing the mode - i.e., NaNs will never
	 * appear in the output array.  If the sample includes only NaNs or has
	 * length 0, an empty array is returned.
	 *
	 * @param sample input data
	 * @return array of array of the most frequently occurring element(s) sorted in ascending order.
	 * @ if the indices are invalid or the array is NULL
	 */
	static std::vector<double> mode(double... sample)
	{
		//Math_Utils::check_not_null(sample, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
		return get_mode(sample, 0, sample.size());
	}

	/**
	 * Returns the sample mode(s).
	 * <p>
	 * The mode is the most frequently occurring value in the sample.
	 * If there is a unique value with maximum frequency, this value is returned
	 * as the only element of the output array. Otherwise, the returned array
	 * contains the maximum frequency elements in increasing order.
	 * <p>
	 * For example, if {@code sample} is {0, 12, 5, 6, 0, 13, 5, 17}, * the returned array will have length two, with 0 in the first element and
	 * 5 in the second.
	 * <p>
	 * NaN values are ignored when computing the mode - i.e., NaNs will never
	 * appear in the output array.  If the sample includes only NaNs or has
	 * length 0, an empty array is returned.
	 *
	 * @param sample input data
	 * @param begin index (0-based) of the first array element to include
	 * @param length the number of elements to include
	 * @return array of array of the most frequently occurring element(s) sorted in ascending order.
	 * @ if the indices are invalid or the array is NULL
	 */
	static std::vector<double> mode(std::vector<double> sample, const int& begin, const int& length)
	{
		//Math_Utils::check_not_null(sample, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);

		if (begin < 0)
		{
			throw std::exception("not implemented");
			// throw (hipparchus::exception::Localized_Core_Formats_Type::START_POSITION, static_cast<int>(begin));
		}

		if (length < 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::hipparchus::exception::Localized_Core_Formats_Type::size(), static_cast<int>(length));
		}

		return get_mode(sample, begin, length);
	}
};