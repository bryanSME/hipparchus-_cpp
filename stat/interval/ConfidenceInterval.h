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
  //package org.hipparchus.stat.interval;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.stat.Localized_Stat_Formats;

  /**
   * Represents an interval estimate of a population parameter.
   */
class Confidence_Interval
{
	/** Lower endpoint of the interval */
	private double lower_bound;

	/** Upper endpoint of the interval */
	private double upper_bound;

	/**
	 * The asserted probability that the interval contains the population
	 * parameter
	 */
	private double confidence_level;

	/**
	 * Create a confidence interval with the given bounds and confidence level.
	 * <p>
	 * Preconditions:
	 * <ul>
	 * <li>{@code lower} must be strictly less than {@code upper}</li>
	 * <li>{@code confidence_level} must be strictly between 0 and 1 (exclusive)</li>
	 * </ul>
	 * </p>
	 *
	 * @param lower_bound lower endpoint of the interval
	 * @param upper_bound upper endpoint of the interval
	 * @param confidence_level coverage probability
	 * @ if the preconditions are not met
	 */
	public Confidence_Interval(double lower_bound, double upper_bound, double confidence_level)
	{
		check_parameters(lower_bound, upper_bound, confidence_level);
		this.lower_bound = lower_bound;
		this.upper_bound = upper_bound;
		this.confidence_level = confidence_level;
	}

	/**
	 * @return the lower endpoint of the interval
	 */
	public double get_lower_bound()
	{
		return lower_bound;
	}

	/**
	 * @return the upper endpoint of the interval
	 */
	public double get_upper_bound()
	{
		return upper_bound;
	}

	/**
	 * @return the asserted probability that the interval contains the
	 *         population parameter
	 */
	public double get_confidence_level()
	{
		return confidence_level;
	}

	/**
	 * @return std::string representation of the confidence interval
	 */
	 //override
	public std::string to_string() const
	{
		return "[" + lower_bound + ";" + upper_bound + "] (confidence level:" + confidence_level + ")";
	}

	/**
	 * Verifies that (lower, upper) is a valid non-empty interval and confidence
	 * is strictly between 0 and 1.
	 *
	 * @param lower lower endpoint
	 * @param upper upper endpoint
	 * @param confidence confidence level
	 */
	private void check_parameters(double lower, double upper, double confidence)
	{
		if (lower >= upper)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::LOWER_BOUND_NOT_BELOW_UPPER_BOUND, lower, upper);
		}
		if (confidence <= 0 || confidence >= 1)
		{
			throw std::exception("not implemented");
			//throw (Localized_Stat_Formats.OUT_OF_BOUNDS_CONFIDENCE_LEVEL, confidence, 0, 1);
		}
	}
}
