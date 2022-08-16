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
  //import java.util.function.Double_Consumer;

  //import org.hipparchus.exception.;
  //import org.hipparchus.stat.descriptive.moment.Geometric_Mean;
  //import org.hipparchus.stat.descriptive.moment.Mean;
  //import org.hipparchus.stat.descriptive.moment.Second_Moment;
  //import org.hipparchus.stat.descriptive.moment.Variance;
  //import org.hipparchus.stat.descriptive.rank.Max;
  //import org.hipparchus.stat.descriptive.rank.Min;
  //import org.hipparchus.stat.descriptive.rank.Random_Percentile;
  //import org.hipparchus.stat.descriptive.summary.Sum;
  //import org.hipparchus.stat.descriptive.summary.SumOf_logs;
  //import org.hipparchus.stat.descriptive.summary.Sum_Of_Squares;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Precision;
#include "summary/Sum.h"
#include "summary/SumOfLogs.h"
#include "summary/SumOfSquares.h"
#include "rank/"

/**
 * Computes summary statistics for a stream of data values added using the
 * {@link #add_valuestatic_cast<double>( add_value} method. The data values are not stored in
 * memory, so this class can be used to compute statistics for very large data
 * streams.
 * <p>
 * By default, all statistics other than percentiles are maintained.  Percentile
 * calculations use an embedded {@link Random_Percentile} which carries more memory
 * and compute overhead than the other statistics, so it is disabled by default.
 * To enable percentiles, either pass {@code true} to the constructor or use a
 * {@link Streaming_StatisticsBuilder} to configure an instance with percentiles turned
 * on. Other stats can also be selectively disabled using
 * {@code Streaming_StatisticsBulder}.
 * <p>
 * Note: This class is not thread-safe.
 */
class Streaming_Statistics : Statistical_Summary, Aggregatable_Statistic<Streaming_Statistics>, Double_Consumer
{
private:
	/** count of values that have been added */
	long n;

	/** Second_Moment is used to compute the mean and variance */
	const Second_Moment second_moment;
	/** min of values that have been added */
	const Min min_impl;
	/** max of values that have been added */
	const Max max_impl;
	/** sum of values that have been added */
	const Sum sum_impl;
	/** sum of the square of each value that has been added */
	const Sum_Of_Squares sum_of_squares_impl;
	/** sum_log of values that have been added */
	const SumOf_logs sum_of_logs_impl;
	/** mean of values that have been added */
	const Mean mean_impl;
	/** variance of values that have been added */
	const Variance variance_impl;
	/** geo_mean of values that have been added */
	const Geometric_Mean geo_mean_impl;
	/** population variance of values that have been added */
	const Variance population_variance;
	/** source of percentiles */
	const Random_Percentile random_percentile;

	/** whether or not moment stats (sum, mean, variance) are maintained */
	const bool compute_moments;
	/** whether or not sum of squares and quadratic mean are maintained */
	const bool compute_sum_of_squares;
	/** whether or not sum of logs and geometric mean are maintained */
	const bool compute_sum_of_logs;
	/** whether or not percentiles are maintained */
	const bool compute_percentiles;
	/** whether or not min and max are maintained */
	const bool compute_extrema;

	/**
	 * Private constructor used by {@link Streaming_StatisticsBuilder}.
	 *
	 * @param compute_percentiles whether or not percentiles are maintained
	 * @param compute_moments whether or not moment stats (mean, sum, variance) are maintained
	 * @param compute_sum_of_logs whether or not sum of logs and geometric mean are maintained
	 * @param compute_sum_of_squares whether or not sum of squares and quadratic mean are maintained
	 * @param compute_extrema whether or not min and max are maintained
	 */
	Streaming_Statistics(bool compute_percentiles, bool compute_moments, bool compute_sum_of_logs, bool compute_sum_of_squares, bool compute_extrema)
	{
		compute_moments = compute_moments;
		compute_sum_of_logs = compute_sum_of_logs;
		compute_sum_of_squares = compute_sum_of_squares;
		compute_percentiles = compute_percentiles;
		compute_extrema = compute_extrema;

		second_moment = compute_moments ? Second_Moment() : NULL;
		max_impl = compute_extrema ? Max() : NULL;
		min_impl = compute_extrema ? Min() : NULL;
		sum_impl = compute_moments ? Sum() : NULL;
		sum_of_squares_impl = compute_sum_of_squares ? Sum_Of_Squares() : NULL;
		sum_of_logs_impl = compute_sum_of_logs ? SumOf_logs() : NULL;
		mean_impl = compute_moments ? Mean(this.second_moment) : NULL;
		variance_impl = compute_moments ? Variance(this.second_moment) : NULL;
		geo_mean_impl = compute_sum_of_logs ? Geometric_Mean(this.sum_of_logs_impl) : NULL;
		population_variance = compute_moments ? Variance(false, this.second_moment) : NULL;
		random_percentile = compute_percentiles ? Random_Percentile() : NULL;
	}

public:
	/**
	 * Construct a Streaming_Statistics instance, maintaining all statistics
	 * other than percentiles.
	 */
	Streaming_Statistics()
	{
		Streaming_Statistics(false);
	}

	/**
	 * Construct a Streaming_Statistics instance, maintaining all statistics
	 * other than percentiles and with/without percentiles per the argument.
	 *
	 * @param compute_percentiles whether or not percentiles are maintained
	 */
	Streaming_Statistics(bool compute_percentiles)
	{
		Streaming_Statistics(compute_percentiles, true, true, true, true);
	}

	/**
	 * A copy constructor. Creates a deep-copy of the {@code original}.
	 *
	 * @param original the {@code Streaming_Statistics} instance to copy
	 * @ if original is NULL
	 */
	Streaming_Statistics(Streaming_Statistics original)
	{
		//Math_Utils::check_not_null(original);

		this.n = original.n;
		this.second_moment = original.compute_moments ? original.second_moment.copy() : NULL;
		this.max_impl = original.compute_extrema ? original.max_impl.copy() : NULL;
		this.min_impl = original.compute_extrema ? original.min_impl.copy() : NULL;
		this.sum_impl = original.compute_moments ? original.sum_impl.copy() : NULL;
		this.sum_of_logs_impl = original.compute_sum_of_logs ? original.sum_of_logs_impl.copy() : NULL;
		this.sum_of_squares_impl = original.compute_sum_of_squares ? original.sum_of_squares_impl.copy() : NULL;

		// Keep statistics with embedded moments in synch
		this.mean_impl = original.compute_moments ? Mean(this.second_moment) : NULL;
		this.variance_impl = original.compute_moments ? Variance(this.second_moment) : NULL;
		this.geo_mean_impl = original.compute_sum_of_logs ? Geometric_Mean(this.sum_of_logs_impl) : NULL;
		this.population_variance = original.compute_moments ? Variance(false, this.second_moment) : NULL;
		this.random_percentile = original.compute_percentiles ? original.random_percentile.copy() : NULL;

		this.compute_moments = original.compute_moments;
		this.compute_sum_of_logs = original.compute_sum_of_logs;
		this.compute_sum_of_squares = original.compute_sum_of_squares;
		this.compute_percentiles = original.compute_percentiles;
		this.compute_extrema = original.compute_extrema;
	}

	/**
	 * Returns a copy of this Streaming_Statistics instance with the same internal state.
	 *
	 * @return a copy of this
	 */
	public Streaming_Statistics copy()
	{
		return Streaming_Statistics(this);
	}

	/**
	 * Return a {@link Statistical_Summary_values} instance reporting current
	 * statistics.
	 * @return Current values of statistics
	 */
	public Statistical_Summary get_summary()
	{
		return Statistical_Summary_values(get_mean(), get_variance(), get_n(), get_max(), get_min(), get_sum());
	}

	/**
	 * Add a value to the data
	 * @param value the value to add
	 */
	public void add_value(double value)
	{
		if (compute_moments)
		{
			second_moment.increment(value);
			sum_impl.increment(value);
		}
		if (compute_extrema)
		{
			min_impl.increment(value);
			max_impl.increment(value);
		}
		if (compute_sum_of_squares)
		{
			sum_of_squares_impl.increment(value);
		}
		if (compute_sum_of_logs)
		{
			sum_of_logs_impl.increment(value);
		}
		if (compute_percentiles)
		{
			random_percentile.increment(value);
		}
		n++;
	}

	/** {@inherit_doc} */
	//override
	public void accept(double value)
	{
		add_value(value);
	}

	/**
	 * Resets all statistics and storage.
	 */
	public void clear()
	{
		this.n = 0;
		if (compute_extrema)
		{
			min_impl.clear();
			max_impl.clear();
		}
		if (compute_moments)
		{
			sum_impl.clear();
			second_moment.clear();
		}
		if (compute_sum_of_logs)
		{
			sum_of_logs_impl.clear();
		}
		if (compute_sum_of_squares)
		{
			sum_of_squares_impl.clear();
		}
		if (compute_percentiles)
		{
			random_percentile.clear();
		}
	}

	/** {@inherit_doc} */
	//override
	public long get_n()
	{
		return n;
	}

	/** {@inherit_doc} */
	//override
	public double get_max()
	{
		return compute_extrema ? max_impl.get_result() : std::numeric_limits<double>::quiet_NaN();
	}

	/** {@inherit_doc} */
	//override
	public double get_min()
	{
		return compute_extrema ? min_impl.get_result() : std::numeric_limits<double>::quiet_NaN();
	}

	/** {@inherit_doc} */
	//override
	public double get_sum()
	{
		return compute_moments ? sum_impl.get_result() : std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * Returns the sum of the squares of the values that have been added.
	 * <p>
	 *NAN is returned if no values have been added.
	 *
	 * @return The sum of squares
	 */
	public double get_sum_of_squares()
	{
		return compute_sum_of_squares ? sum_of_squares_impl.get_result() : std::numeric_limits<double>::quiet_NaN();
	}

	/** {@inherit_doc} */
	//override
	public double get_mean()
	{
		return compute_moments ? mean_impl.get_result() : std::numeric_limits<double>::quiet_NaN();
	}

	/** {@inherit_doc} */
	//override
	public double get_variance()
	{
		return compute_moments ? variance_impl.get_result() : std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * Returns the <a href="http://en.wikibooks.org/wiki/Statistics/Summary/Variance">
	 * population variance</a> of the values that have been added.
	 * <p>
	 *NAN is returned if no values have been added.
	 *
	 * @return the population variance
	 */
	public double get_population_variance()
	{
		return compute_moments ? population_variance.get_result() : std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * Returns the geometric mean of the values that have been added.
	 * <p>
	 *NAN is returned if no values have been added.
	 *
	 * @return the geometric mean
	 */
	public double get_geometric_mean()
	{
		return compute_sum_of_logs ? geo_mean_impl.get_result() : std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * Returns the sum of the logs of the values that have been added.
	 * <p>
	 *NAN is returned if no values have been added.
	 *
	 * @return the sum of logs
	 */
	public double get_sum_of_logs()
	{
		return compute_sum_of_logs ? sum_of_logs_impl.get_result() : std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * Returns a statistic related to the Second Central Moment. Specifically, * what is returned is the sum of squared deviations from the sample mean
	 * among the values that have been added.
	 * <p>
	 * Returns <code>Double.NaN</code> if no data values have been added and
	 * returns <code>0</code> if there is just one value in the data set.
	 *
	 * @return second central moment statistic
	 */
	public double get_second_moment()
	{
		return compute_moments ? second_moment.get_result() : std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * Returns the quadratic mean, a.k.a.
	 * <a href="http://mathworld.wolfram.com/Root-Mean-Square.html">
	 * root-mean-square</a> of the available values
	 *
	 * @return The quadratic mean or {@codeNAN} if no values
	 * have been added.
	 */
	public double get_quadratic_mean()
	{
		if (compute_sum_of_squares)
		{
			long size = get_n();
			return size > 0 ? std::sqrt(get_sum_of_squares() / size) : std::numeric_limits<double>::quiet_NaN();
		}
		else
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
	}

	/**
	 * Returns the standard deviation of the values that have been added.
	 * <p>
	 *NAN is returned if no values have been added.
	 *
	 * @return the standard deviation
	 */
	 //override
	public double get_standard_deviation()
	{
		long size = get_n();
		if (compute_moments)
		{
			if (size > 0)
			{
				return size > 1 ? std::sqrt(get_variance()) : 0.0;
			}
			else
			{
				return std::numeric_limits<double>::quiet_NaN();
			}
		}
		else
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
	}

	/**
	 * Returns an estimate of the median of the values that have been entered.
	 * See {@link Random_Percentile} for a description of the algorithm used for large
	 * data streams.
	 *
	 * @return the median
	 */
	public double get_median()
	{
		return random_percentile != NULL ? random_percentile.get_result(50d) : std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * Returns an estimate of the given percentile of the values that have been entered.
	 * See {@link Random_Percentile} for a description of the algorithm used for large
	 * data streams.
	 *
	 * @param percentile the desired percentile (must be between 0 and 100)
	 * @return estimated percentile
	 */
	public double get_percentile(double percentile)
	{
		return random_percentile != NULL ? random_percentile.get_result(percentile) : std::numeric_limits<double>::quiet_NaN();
	}

	/**
	 * {@inherit_doc}
	 * Statistics are aggregated only when both this and other are maintaining them.  For example, * if this.compute_moments is false, but other.compute_moments is true, the moment data in other
	 * will be lost.
	 */
	 //override
	public void aggregate(Streaming_Statistics other)
	{
		//Math_Utils::check_not_null(other);

		if (other.n > 0)
		{
			this.n += other.n;
			if (compute_moments && other.compute_moments)
			{
				this.second_moment.aggregate(other.second_moment);
				this.sum_impl.aggregate(other.sum_impl);
			}
			if (compute_extrema && other.compute_extrema)
			{
				this.min_impl.aggregate(other.min_impl);
				this.max_impl.aggregate(other.max_impl);
			}
			if (compute_sum_of_logs && other.compute_sum_of_logs)
			{
				this.sum_of_logs_impl.aggregate(other.sum_of_logs_impl);
			}
			if (compute_sum_of_squares && other.compute_sum_of_squares)
			{
				this.sum_of_squares_impl.aggregate(other.sum_of_squares_impl);
			}
			if (compute_percentiles && other.compute_percentiles)
			{
				this.random_percentile.aggregate(other.random_percentile);
			}
		}
	}

	/**
	 * Generates a text report displaying summary statistics from values that
	 * have been added.
	 *
	 * @return std::string with line feeds displaying statistics
	 */
	 //override
	public std::string to_string() const
	{
		std::stringBuilder out_buffer = std::stringBuilder(200); // the size is just a wild guess
		std::string endl = "\n";
		out_buffer.append("Streaming_Statistics:").append(endl).
			append("n: ").append(get_n()).append(endl).
			append("min: ").append(get_min()).append(endl).
			append("max: ").append(get_max()).append(endl).
			append("sum: ").append(get_sum()).append(endl).
			append("mean: ").append(get_mean()).append(endl).
			append("variance: ").append(get_variance()).append(endl).
			append("population variance: ").append(get_population_variance()).append(endl).
			append("standard deviation: ").append(get_standard_deviation()).append(endl).
			append("geometric mean: ").append(get_geometric_mean()).append(endl).
			append("second moment: ").append(get_second_moment()).append(endl).
			append("sum of squares: ").append(get_sum_of_squares()).append(endl).
			append("sum of logs: ").append(get_sum_of_logs()).append(endl);
		return out_buffer.to_string();
	}

	/**
	 * Returns true iff <code>object</code> is a <code>Streaming_Statistics</code>
	 * instance and all statistics have the same values as this.
	 *
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
		if (!(object instanceof Streaming_Statistics))
		{
			return false;
		}
		Streaming_Statistics other = (Streaming_Statistics)object;
		return other.get_n() == get_n() &&
			Precision::equals_including_nan(other.get_max(), get_max()) &&
			Precision::equals_including_nan(other.get_min(), get_min()) &&
			Precision::equals_including_nan(other.get_sum(), get_sum()) &&
			Precision::equals_including_nan(other.get_geometric_mean(), get_geometric_mean()) &&
			Precision::equals_including_nan(other.get_mean(), get_mean()) &&
			Precision::equals_including_nan(other.get_sum_of_squares(), get_sum_of_squares()) &&
			Precision::equals_including_nan(other.get_sum_of_logs(), get_sum_of_logs()) &&
			Precision::equals_including_nan(other.get_variance(), get_variance()) &&
			Precision::equals_including_nan(other.get_median(), get_median());
	}

	/**
	 * Returns hash code based on values of statistics.
	 * @return hash code
	 */
	 //override
	public int hash_code()
	{
		int result = 31 + Math_Utils::hash(get_n());
		result = result * 31 + Math_Utils::hash(get_max());
		result = result * 31 + Math_Utils::hash(get_min());
		result = result * 31 + Math_Utils::hash(get_sum());
		result = result * 31 + Math_Utils::hash(get_geometric_mean());
		result = result * 31 + Math_Utils::hash(get_mean());
		result = result * 31 + Math_Utils::hash(get_sum_of_squares());
		result = result * 31 + Math_Utils::hash(get_sum_of_logs());
		result = result * 31 + Math_Utils::hash(get_variance());
		result = result * 31 + Math_Utils::hash(get_median());
		return result;
	}

	/**
	 * Returns a {@link Streaming_StatisticsBuilder} to source configured
	 * {@code Streaming_Statistics} instances.
	 *
	 * @return a Streaming_StatisticsBuilder instance
	 */
	public static Streaming_StatisticsBuilder builder()
	{
		return Streaming_StatisticsBuilder();
	}

	/**
	 * Builder for Streaming_Statistics instances.
	 */
	public static class Streaming_StatisticsBuilder
	{
		/** whether or not moment statistics are maintained by instances created by this factory */
		private bool compute_moments;
		/** whether or not sum of squares and quadratic mean are maintained by instances created by this factory */
		private bool compute_sum_of_squares;
		/** whether or not sum of logs and geometric mean are maintained by instances created by this factory */
		private bool compute_sum_of_logs;
		/** whether or not percentiles are maintained by instances created by this factory */
		private bool compute_percentiles;
		/** whether or not min and max are maintained by instances created by this factory */
		private bool compute_extrema;

		/** Simple constructor.
		 */
		public Streaming_StatisticsBuilder()
		{
			compute_moments = true;
			compute_sum_of_squares = true;
			compute_sum_of_logs = true;
			compute_percentiles = false;
			compute_extrema = true;
		}

		/**
		 * Sets the compute_moments setting of the factory
		 *
		 * @param arg whether or not instances created using {@link #build()} will
		 * maintain moment statistics
		 * @return a factory with the given compute_moments property set
		 */
		public Streaming_StatisticsBuilder moments(bool arg)
		{
			this.compute_moments = arg;
			return this;
		}

		/**
		 * Sets the compute_sum_of_logs setting of the factory
		 *
		 * @param arg whether or not instances created using {@link #build()} will
		 * maintain log sums
		 * @return a factory with the given compute_sum_of_logs property set
		 */
		public Streaming_StatisticsBuilder sum_of_logs(bool arg)
		{
			this.compute_sum_of_logs = arg;
			return this;
		}

		/**
		 * Sets the compute_sum_of_squares setting of the factory.
		 *
		 * @param arg whether or not instances created using {@link #build()} will
		 * maintain sums of squares
		 * @return a factory with the given compute_sum_of_squares property set
		 */
		public Streaming_StatisticsBuilder sum_of_squares(bool arg)
		{
			this.compute_sum_of_squares = arg;
			return this;
		}

		/**
		 * Sets the compute_percentiles setting of the factory.
		 *
		 * @param arg whether or not instances created using {@link #build()} will
		 * compute percentiles
		 * @return a factory with the given compute_percentiles property set
		 */
		public Streaming_StatisticsBuilder percentiles(bool arg)
		{
			this.compute_percentiles = arg;
			return this;
		}

		/**
		 * Sets the compute_extrema setting of the factory.
		 *
		 * @param arg whether or not instances created using {@link #build()} will
		 * compute min and max
		 * @return a factory with the given compute_extrema property set
		 */
		public Streaming_StatisticsBuilder extrema(bool arg)
		{
			this.compute_extrema = arg;
			return this;
		}

		/**
		 * Builds a Streaming_Statistics instance with currently defined properties.
		 *
		 * @return newly configured Streaming_Statistics instance
		 */
		public Streaming_Statistics build()
		{
			return Streaming_Statistics(compute_percentiles, compute_moments, compute_sum_of_logs, compute_sum_of_squares, compute_extrema);
		}
	}
};