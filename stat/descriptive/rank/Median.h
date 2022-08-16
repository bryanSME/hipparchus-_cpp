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
  //package org.hipparchus.stat.descriptive.rank;

  //import java.io.Serializable;

  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.;
  //import org.hipparchus.stat.descriptive.Abstract_Univariate_Statistic;
  //import org.hipparchus.stat.descriptive.rank.Percentile.Estimation_Type;
  //import org.hipparchus.stat.ranking.NaN_Strategy;
  //import org.hipparchus.util.Kth_Selector;
#include <vector>
#include "../../../core/util/KthSelector.h"
#include "Percentile.h"
#include "../../ranking/NaNStrategy.h"

/**
 * Returns the median of the available values.  This is the same as the 50th percentile.
 * See {@link Percentile} for a description of the algorithm used.
 * <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
class Median : Abstract_Univariate_Statistic
{
private:

	/** Fixed quantile. */
	static constexpr double FIXED_QUANTILE_50{ 50.0 };

	/** The percentile impl to calculate the median. */
	const Percentile my_percentile;

	/**
	 * Constructs a Median with the specific {@link Estimation_Type}, * {@link NaN_Strategy} and {@link Kth_Selector}.
	 *
	 * @param estimation_type one of the percentile {@link Estimation_Type estimation types}
	 * @param nan_strategy one of {@link NaN_Strategy} to handle with NaNs
	 * @param kth_selector {@link Kth_Selector} to use for pivoting during search
	 * @ if p is not within (0,100]
	 * @ if type or NaN_Strategy passed is NULL
	 */
	Median(const Estimation_Type& estimation_type, const NaN_Strategy& nan_strategy, const Kth_Selector& kth_selector) : my_percentile{ Percentile(FIXED_QUANTILE_50, estimation_type, nan_strategy, kth_selector) } {};

public:
	/**
	 * Default constructor.
	 */
	Median() : my_percentile{ Percentile(FIXED_QUANTILE_50) } {};

	Percentile get_percentile() const
	{
		return my_percentile;
	}

	/**
	 * Copy constructor, creates a {@code Median} identical
	 * to the {@code original}
	 *
	 * @param original the {@code Median} instance to copy
	 * @ if original is NULL
	 */
	Median(const Median& original)
	{
		Abstract_Univariate_Statistic(original);
		//super(original);
		my_percentile = original.get_percentile();
	}

	/** {@inherit_doc} */
	//override
	double evaluate(std::vector<double> values, int begin, int length)

	{
		return my_percentile.evaluate(values, begin, length);
	}

	/** {@inherit_doc} */
	//override
	Median copy()
	{
		return Median(*this);
	}

	/**
	 * Get the estimation {@link Estimation_Type type} used for computation.
	 *
	 * @return the {@code estimation_type} set
	 */
	Estimation_Type get_estimation_type()
	{
		return percentile.get_estimation_type();
	};

	/**
	 * Build a instance similar to the current one except for the
	 * {@link Estimation_Type estimation type}.
	 *
	 * @param new_estimation_type estimation type for the instance
	 * @return a instance, with changed estimation type
	 * @ when new_estimation_type is NULL
	 */
	Median with_estimation_type(const Estimation_Type new_estimation_type)
	{
		return Median(new_estimation_type, my_percentile.get_nan_strategy(), my_percentile.get_kth_selector());
	};

	/**
	 * Get the {@link NaN_Strategy NaN Handling} strategy used for computation.
	 * @return {@code NaN Handling} strategy set during construction
	 */
	NaN_Strategy get_nan_strategy()
	{
		return my_percentile.get_nan_strategy();
	};

	/**
	 * Build a instance similar to the current one except for the
	 * {@link NaN_Strategy NaN handling} strategy.
	 *
	 * @param new_nan_strategy NaN strategy for the instance
	 * @return a instance, with changed NaN handling strategy
	 * @ when new_nan_strategy is NULL
	 */
	Median with_na_n_strategy(const NaN_Strategy new_nan_strategy)
	{
		return Median(my_percentile.get_estimation_type(), new_nan_strategy, my_percentile.get_kth_selector());
	};

	/**
	 * Get the {@link Kth_Selector kth_selector} used for computation.
	 * @return the {@code kth_selector} set
	 */
	Kth_Selector get_kth_selector()
	{
		return my_percentile.get_kth_selector();
	};

	/**
	 * Build a instance similar to the current one except for the
	 * {@link Kth_Selector kth_selector} instance specifically set.
	 *
	 * @param new_kth_selector Kth_Selector for the instance
	 * @return a instance, with changed Kth_Selector
	 * @ when new_kth_selector is NULL
	 */
	Median with_kth_selector(const Kth_Selector& new_kth_selector)
	{
		return Median(my_percentile.get_estimation_type(), my_percentile.get_nan_strategy(), new_kth_selector);
	};
};