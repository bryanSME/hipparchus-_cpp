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
  //import java.util.Arrays;
  //import java.util.function.Double_Consumer;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.stat.descriptive.moment.Geometric_Mean;
  //import org.hipparchus.stat.descriptive.moment.Kurtosis;
  //import org.hipparchus.stat.descriptive.moment.Mean;
  //import org.hipparchus.stat.descriptive.moment.Skewness;
  //import org.hipparchus.stat.descriptive.moment.Variance;
  //import org.hipparchus.stat.descriptive.rank.Max;
  //import org.hipparchus.stat.descriptive.rank.Min;
  //import org.hipparchus.stat.descriptive.rank.Percentile;
  //import org.hipparchus.stat.descriptive.summary.Sum;
  //import org.hipparchus.stat.descriptive.summary.Sum_Of_Squares;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Resizable_Double_Array;
#include <sstream>
#include "StatisticalSummary.h"
#include <vector>
#include <string>
#include <cmath>
#include "UnivariateStatistic.h"
#include "DescriptiveStatistics.h"
#include "../../stat/descriptive/rank/Percentile.h"

/**
 * Maintains a dataset of values of a single variable and computes descriptive
 * statistics based on stored data.
 * <p>
 * The {@link #get_window_size() window_size} property sets a limit on the number
 * of values that can be stored in the dataset. The default value, INFINITE_WINDOW, * puts no limit on the size of the dataset. This value should be used with
 * caution, as the backing store will grow without bound in this case.
 * <p>
 * For very large datasets, {@link Streaming_Statistics}, which does not store
 * the dataset, should be used instead of this class. If <code>window_size</code>
 * is not INFINITE_WINDOW and more values are added than can be stored in the
 * dataset, values are added in a "rolling" manner, with values replacing
 * the "oldest" values in the dataset.
 * <p>
 * Note: this class is not threadsafe.
 */
class Descriptive_Statistics : public Statistical_Summary //, public Double_Consumer
{
protected:
	/**
	 * Represents an infinite window size.  When the {@link #get_window_size()}
	 * returns this value, there is no limit to the number of data values
	 * that can be stored in the dataset.
	 */
	static constexpr int INFINITE_WINDOW{ -1 };

	/**
	 * Copy constructor.
	 * <p>
	 * Construct a Descriptive_Statistics instance that
	 * is a copy of original.
	 *
	 * @param original Descriptive_Statistics instance to copy
	 * @org.hipparchus.exception. if original is NULL
	 */
	Descriptive_Statistics(const Descriptive_Statistics& original)
	{
		//Math_Utils::check_not_null(original);

		// Copy data and window size
		my_window_size = original.get_window_size();
		my_eDA = original.get_eDA();

		// Copy implementations
		my_max_impl = original.get_max_impl();
		my_min_impl = original.get_min_impl();
		my_mean_impl = original.get_mean_impl();
		my_sum_impl = original.get_sum_impl();
		my_sum_of_squares_impl = original.get_sum_of_squares_impl();
		my_variance_impl = original.get_variance_impl();
		my_geometric_mean_impl = original.get_geometric_mean_impl();
		my_kurtosis_impl = original.get_kurtosis_impl();
		my_skewness_impl = original.get_skewness_impl();
		my_percentile_impl = original.get_percentile_impl();
	}

private:
	/** The statistic used to calculate the population variance - fixed. */
	static const Univariate_Statistic POPULATION_VARIANCE = Variance(false);

	/** Maximum statistic implementation. */
	Univariate_Statistic          my_max_impl;
	/** Minimum statistic implementation. */
	Univariate_Statistic          my_min_impl;
	/** Sum statistic implementation. */
	Univariate_Statistic          my_sum_impl;
	/** Sum of squares statistic implementation. */
	Univariate_Statistic          my_sum_of_squares_impl;
	/** Mean statistic implementation. */
	Univariate_Statistic          my_mean_impl;
	/** Variance statistic implementation. */
	Univariate_Statistic          my_variance_impl;
	/** Geometric mean statistic implementation. */
	Univariate_Statistic          my_geometric_mean_impl;
	/** Kurtosis statistic implementation. */
	Univariate_Statistic          my_kurtosis_impl;
	/** Skewness statistic implementation. */
	Univariate_Statistic          my_skewness_impl;
	/** Percentile statistic implementation. */
	Percentile                    my_percentile_impl;

	/** holds the window size. */
	int my_window_size;

	/** Stored data values. */
	Resizable_Double_Array my_eDA;

public:
	/**
	 * Construct a Descriptive_Statistics instance with an infinite window.
	 */
	Descriptive_Statistics()
	{
		Descriptive_Statistics(INFINITE_WINDOW);
	}

	/**
	 * Construct a Descriptive_Statistics instance with the specified window.
	 *
	 * @param size the window size.
	 * @ if window size is less than 1 but
	 * not equal to {@link #INFINITE_WINDOW}
	 */
	Descriptive_Statistics(const int& size)
	{
		Descriptive_Statistics(size, false, NULL);
	}

	/**
	 * Construct a Descriptive_Statistics instance with an infinite window
	 * and the initial data values in std::vector<double> initial_double_array.
	 *
	 * @param initial_double_array the initial std::vector<double>.
	 * @org.hipparchus.exception. if the input array is NULL
	 */
	Descriptive_Statistics(std::vector<double> initial_double_array)
	{
		Descriptive_Statistics(INFINITE_WINDOW, true, initial_double_array);
	}

	/**
	 * Construct a Descriptive_Statistics instance with the specified window.
	 *
	 * @param window_size the window size
	 * @param has_initial_values if initial values have been provided
	 * @param initial_values the initial values
	 * @org.hipparchus.exception. if initial_values is NULL
	 * @ if window size is less than 1 but
	 * not equal to {@link #INFINITE_WINDOW}
	 */
	Descriptive_Statistics(const int& window_size, bool has_initial_values, const std::vector<double>& initial_values)
	{
		if (window_size < 1 && window_size != INFINITE_WINDOW)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_WINDOW_SIZE, window_size);
		}

		if (has_initial_values)
		{
			//Math_Utils::check_not_null(initial_values, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
		}

		my_window_size = window_size;
		int initial_capacity = my_window_size < 0
			? 100
			: my_window_size;
		my_eDA = has_initial_values
			? Resizable_Double_Array(initial_values)
			: Resizable_Double_Array(initial_capacity);

		max_impl = Max();
		my_min_impl = Min();
		my_sum_impl = Sum();
		my_sum_of_squares_impl = Sum_Of_Squares();
		my_mean_impl = Mean();
		my_variance_impl = Variance();
		my_geometric_mean_impl = Geometric_Mean();
		my_kurtosis_impl = Kurtosis();
		my_skewness_impl = Skewness();
		my_percentile_impl = Percentile();
	}

	/**
	 * Returns a copy of this Descriptive_Statistics instance with the same internal state.
	 *
	 * @return a copy of this
	 */
	Descriptive_Statistics copy()
	{
		return Descriptive_Statistics(*this);
	}

	/**
	 * Adds the value to the dataset. If the dataset is at the maximum size
	 * (i.e., the number of stored elements equals the currently configured
	 * window_size), the first (oldest) element in the dataset is discarded
	 * to make room for the value.
	 *
	 * @param v the value to be added
	 */
	void add_value(double v)
	{
		if (my_window_size != INFINITE_WINDOW)
		{
			if (get_n() == my_window_size)
			{
				my_eDA.add_element_rolling(v);
			}
			else if (get_n() < my_window_size)
			{
				my_eDA.add_element(v);
			}
		}
		else
		{
			my_eDA.add_element(v);
		}
	}

	/** {@inherit_doc} */
	//override
	void accept(double v)
	{
		add_value(v);
	}

	/**
	 * Resets all statistics and storage.
	 */
	void clear()
	{
		my_eDA.clear();
	}

	/**
	 * Removes the most recent value from the dataset.
	 *
	 * @Math_Illegal_State_Exception if there are no elements stored
	 */
	void remove_most_recent_value()
	{
		throw std::exception("not implemented");
		/*try
		{
			eDA.discard_most_recent_elements(1);
		}
		catch ( ex)
		{
			throw Math_Illegal_State_Exception(ex, hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
		}*/
	}

	/**
	 * Replaces the most recently stored value with the given value.
	 * There must be at least one element stored to call this method.
	 *
	 * @param v the value to replace the most recent stored value
	 * @return replaced value
	 * @Math_Illegal_State_Exception if there are no elements stored
	 */
	double replace_most_recent_value(const double& v)
	{
		return my_eDA.substitute_most_recent_element(v);
	}

	/**
	 * Apply the given statistic to the data associated with this set of statistics.
	 * @param stat the statistic to apply
	 * @return the computed value of the statistic.
	 */
	double apply(const Univariate_Statistic& stat)
	{
		// No try-catch or advertised exception here because arguments
		// are guaranteed valid.
		return my_eDA.compute(stat);
	}

	/** {@inherit_doc} */
	//override
	double get_mean()
	{
		return apply(my_mean_impl);
	}

	/**
	 * Returns the geometric mean of the available values.
	 * <p>
	 * See {@link Geometric_Mean} for details on the computing algorithm.
	 *
	 * @see <a href="http://www.xycoon.com/geometric_mean.htm">
	 * Geometric mean</a>
	 *
	 * @return The geometric_mean,NAN if no values have been added, * or if any negative values have been added.
	 */
	double get_geometric_mean()
	{
		return apply(my_geometric_mean_impl);
	}

	/**
	 * Returns the standard deviation of the available values.
	 * @return The standard deviation,NAN if no values have been added
	 * or 0.0 for a single value set.
	 */
	 //override
	double get_standard_deviation()
	{
		double std_dev = std::numeric_limits<double>::quiet_NaN();
		if (get_n() > 0)
		{
			if (get_n() > 1)
			{
				std_dev = std::sqrt(get_variance());
			}
			else
			{
				std_dev = 0.0;
			}
		}
		return std_dev;
	}

	/**
	 * Returns the quadratic mean of the available values.
	 *
	 * @see <a href="http://mathworld.wolfram.com/Root-Mean-Square.html">
	 * Root Mean Square</a>
	 *
	 * @return The quadratic mean or {@codeNAN} if no values
	 * have been added.
	 */
	double get_quadratic_mean()
	{
		const long n = get_n();
		return n > 0
			? std::sqrt(get_sum_of_squares() / n)
			: std::numeric_limits<double>::quiet_NaN();
	}

	/** {@inherit_doc} */
	//override
	double get_variance()
	{
		return apply(my_variance_impl);
	}

	/**
	 * Returns the population variance of the available values.
	 *
	 * @see <a href="http://en.wikibooks.org/wiki/Statistics/Summary/Variance">
	 * Population variance</a>
	 *
	 * @return The population variance,NAN if no values have been added, * or 0.0 for a single value set.
	 */
	double get_population_variance()
	{
		return apply(POPULATION_VARIANCE);
	}

	/**
	 * Returns the skewness of the available values. Skewness is a
	 * measure of the asymmetry of a given distribution.
	 *
	 * @return The skewness,NAN if less than 3 values have been added.
	 */
	double get_skewness()
	{
		return apply(my_skewness_impl);
	}

	/**
	 * Returns the Kurtosis of the available values. Kurtosis is a
	 * measure of the "peakedness" of a distribution.
	 *
	 * @return The kurtosis,NAN if less than 4 values have been added.
	 */
	double get_kurtosis()
	{
		return apply(my_kurtosis_impl);
	}

	/** {@inherit_doc} */
	//override
	double get_max()
	{
		return apply(my_max_impl);
	}

	/** {@inherit_doc} */
	//override
	double get_min()
	{
		return apply(my_min_impl);
	}

	/** {@inherit_doc} */
	//override
	double get_sum()
	{
		return apply(my_sum_impl);
	}

	/**
	 * Returns the sum of the squares of the available values.
	 * @return The sum of the squares orNAN if no
	 * values have been added.
	 */
	double get_sum_of_squares()
	{
		return apply(my_sum_of_squares_impl);
	}

	/**
	 * Returns an estimate for the pth percentile of the stored values.
	 * <p>
	 * The implementation provided here follows the first estimation procedure presented
	 * <a href="http://www.itl.nist.gov/div898/handbook/prc/section2/prc252.htm">here.</a>
	 * </p><p>
	 * <strong>Preconditions</strong>:<ul>
	 * <li><code>0 &lt; p &le; 100</code> (otherwise an
	 * <code></code> is thrown)</li>
	 * <li>at least one value must be stored (returns <code>Double.NaN
	 *     </code> otherwise)</li>
	 * </ul>
	 *
	 * @param p the requested percentile (scaled from 0 - 100)
	 * @return An estimate for the pth percentile of the stored data
	 * @ if p is not a valid quantile
	 */
	double get_percentile(const double& p)
	{
		my_percentile_impl.set_quantile(p);
		return apply(my_percentile_impl);
	}

	/** {@inherit_doc} */
	//override
	long get_n()
	{
		return my_eDA.get_num_elements();
	}

	/**
	 * Returns the maximum number of values that can be stored in the
	 * dataset, or INFINITE_WINDOW (-1) if there is no limit.
	 *
	 * @return The current window size or -1 if its Infinite.
	 */
	int get_window_size() const
	{
		return my_window_size;
	}

	/**
	 * Window_Size controls the number of values that contribute to the
	 * reported statistics.  For example, if window_size is set to 3 and the
	 * values {1,2,3,4,5} have been added <strong> in that order</strong> then
	 * the <i>available values</i> are {3,4,5} and all reported statistics will
	 * be based on these values. If {@code window_size} is decreased as a result
	 * of this call and there are more than the value of elements in the
	 * current dataset, values from the front of the array are discarded to
	 * reduce the dataset to {@code window_size} elements.
	 *
	 * @param window_size sets the size of the window.
	 * @ if window size is less than 1 but
	 * not equal to {@link #INFINITE_WINDOW}
	 */
	void set_window_size(const int& window_size)
	{
		if (window_size < 1 && window_size != INFINITE_WINDOW)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NOT_POSITIVE_WINDOW_SIZE, window_size);
		}

		my_window_size = window_size;

		// We need to check to see if we need to discard elements
		// from the front of the array.  If the window_size is less than
		// the current number of elements.
		if (window_size != INFINITE_WINDOW && window_size < my_eDA.get_num_elements())
		{
			my_eDA.discard_front_elements(my_eDA.get_num_elements() - window_size);
		}
	}

	/**
	 * Returns the current set of values in an array of double primitives.
	 * The order of addition is preserved.  The returned array is a fresh
	 * copy of the underlying data -- i.e., it is not a reference to the
	 * stored data.
	 *
	 * @return the current set of numbers in the order in which they
	 * were added to this set
	 */
	std::vector<double> get_values() const
	{
		return my_eDA.get_elements();
	}

	/**
	 * Returns the current set of values in an array of double primitives, * sorted in ascending order.  The returned array is a fresh
	 * copy of the underlying data -- i.e., it is not a reference to the
	 * stored data.
	 * @return returns the current set of
	 * numbers sorted in ascending order
	 */
	std::vector<double> get_sorted_values()
	{
		auto sort = get_values();
		Arrays.sort(sort);
		return sort;
	}

	/**
	 * Returns the element at the specified index
	 * @param index The Index of the element
	 * @return return the element at the specified index
	 */
	double get_element(const int& index)
	{
		return my_eDA.get_element(index);
	}

	/**
	 * Generates a text report displaying univariate statistics from values
	 * that have been added.  Each statistic is displayed on a separate line.
	 *
	 * @return std::string with line feeds displaying statistics
	 */
	 //override
	std::string to_string() const
	{
		std::stringstream out_buffer{};
		out_buffer
			<< "Descriptive_Statistics:\n"
			<< "n: " << get_n() << "\n"
			<< "min: " << get_min() << "\n"
			<< "max: " << get_max() << "\n"
			<< "mean: " << get_mean() << "\n"
			<< "std dev: " << get_standard_deviation() << "\n";
		//try
		//{
			// No catch for MIAE because actual parameter is valid below
		out_buffer << "median: " << get_percentile(50) << "\n";
		//}
		//catch (Math_Illegal_State_Exception ex)
		//{
		//    out_buffer.append("median: unavailable").append(endl);
		//}
		out_buffer
			<< "skewness: " << get_skewness() << "\n"
			<< "kurtosis: " << get_kurtosis() << "\n";
		return out_buffer.str();
	}

	int get_window_size() const
	{
		return my_window_size;
	}

	Univariate_Statistic get_eDA() const
	{
		return my_eDA;
	}

	Univariate_Statistic get_max_impl() const
	{
		return my_max_impl;
	}

	Univariate_Statistic get_min_impl() const
	{
		return my_min_impl;
	}

	Univariate_Statistic get_mean_impl() const
	{
		return my_mean_impl;
	}

	Univariate_Statistic get_sum_impl() const
	{
		return my_sum_impl;
	}

	Univariate_Statistic get_sum_of_squares_impl() const
	{
		return my_sum_of_squares_impl
	}

	Univariate_Statistic get_variance_impl() const
	{
		return my_variance_impl;
	}

	Univariate_Statistic get_geometric_mean_impl() const
	{
		return my_geometric_mean_impl;
	}

	Univariate_Statistic get_kurtosis_impl() const
	{
		return my_kurtosis_impl;
	}

	Univariate_Statistic get_skewness_impl() const
	{
		return my_skewness_impl;
	}

	Univariate_Statistic get_percentile_impl() const
	{
		return my_percentile_impl;
	}
};