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
  //import org.hipparchus.stat.descriptive.Weighted_Evaluation;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <limits>
#include <cmath>
#include <vector>
#include "../../descriptive/AbstractStorelessUnivariateStatistic.h"
#include "../../descriptive/AggregatableStatistic.hpp"
#include "../../../core/util/MathArrays.h"
#include "../../../core/util/MathUtils.h"
#include "../../descriptive/WeightedEvaluation.h"

/**
 * Returns the product of the available values.
 * <p>
 * If there are no values in the dataset, then 1 is returned.
 * If any of the values are
 * <code>NaN</code>, then <code>NaN</code> is returned.
 * <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
class Product : public Abstract_Storeless_Univariate_Statistic, public Aggregatable_Statistic<Product>, public Weighted_Evaluation
{
private:

	/** The number of values that have been added */
	long my_n;

	/** The current Running Product */
	double my_value;

public:
	/**
	 * Create a Product instance.
	 */
	Product() : my_n{}, my_value{ 1 }{};

	/**
	 * Copy constructor, creates a {@code Product} identical
	 * to the {@code original}.
	 *
	 * @param original the {@code Product} instance to copy
	 * @  if original is NULL
	 */
	Product(const Product& original) : my_n{ original.get_n() }, my_value{ original.get_status() }
	{
		//Math_Utils::check_not_null(original);
	}

	/** {@inherit_doc} */
	//override
	void increment(const double& d)
	{
		my_value *= d;
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
		my_value = 1;
		my_n = 0;
	}

	/** {@inherit_doc} */
	//override
	void aggregate(const Product& other)
	{
		//Math_Utils::check_not_null(other);
		if (other.get_n() > 0)
		{
			my_n += other.get_n();
			my_value *= other.get_result();
		}
	}

	/**
	 * Returns the product of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 *
	 * @param values the input array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the product of the values or 1 if length = 0
	 * @ if the array is NULL or the array index
	 *  parameters are not valid
	 */
	 //override
	double evaluate(const std::vector<double>& values, const int& begin, const int& length)
	{
		auto product = std::numeric_limits<double>::quiet_NaN();
		if (Math_Arrays::verify_values(values, begin, length, true))
		{
			product = 1.0;
			for (int i{ begin }; i < begin + length; i++)
			{
				product *= values[i];
			}
		}
		return product;
	}

	/**
	 * Returns the weighted product of the entries in the specified portion of
	 * the input array, or <code>Double.NaN</code> if the designated subarray
	 * is empty.
	 * <p>
	 * Throws <code></code> if any of the following are true:
	 * <ul><li>the values array is NULL</li>
	 *     <li>the weights array is NULL</li>
	 *     <li>the weights array does not have the same length as the values array</li>
	 *     <li>the weights array contains one or more infinite values</li>
	 *     <li>the weights array contains one or more NaN values</li>
	 *     <li>the weights array contains negative values</li>
	 *     <li>the start and length arguments do not determine a valid array</li>
	 * </ul>
	 * <p>
	 * Uses the formula, * <pre>
	 *    weighted product = &prod;values[i]<sup>weights[i]</sup>
	 * </pre>
	 * <p>
	 * that is, the weights are applied as exponents when computing the weighted product.
	 *
	 * @param values the input array
	 * @param weights the weights array
	 * @param begin index of the first array element to include
	 * @param length the number of elements to include
	 * @return the product of the values or 1 if length = 0
	 * @ if the parameters are not valid
	 */
	 //override
	double evaluate(const std::vector<double>& values, const std::vector<double>& weights, const int& begin, const int& length)
	{
		auto product = std::numeric_limits<double>::quiet_NaN();
		if (Math_Arrays::verify_values(values, weights, begin, length, true))
		{
			product = 1.0;
			for (int i{ begin }; i < begin + length; i++)
			{
				product *= std::pow(values[i], weights[i]);
			}
		}
		return product;
	}

	/** {@inherit_doc} */
	//override
	//Product copy()
	//{
	//    return Product(this);
	//}
}
