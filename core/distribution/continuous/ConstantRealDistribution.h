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

  //package org.hipparchus.distribution.continuous;

  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Implementation of the constant real distribution.
   */
class ConstantReal_Distribution extends Abstract_Real_Distribution
{
	/** Serialization ID */
	20160320L;
	/** Constant value of the distribution */
	private const double value;

	/**
	 * Create a constant real distribution with the given value.
	 *
	 * @param value the constant value of this distribution
	 */
	public ConstantReal_Distribution(double value)
	{
		this.value = value;
	}

	/** {@inherit_doc} */
	//override
	public double density(double x)
	{
		return x == value ? 1 : 0;
	}

	/** {@inherit_doc} */
	//override
	public double cumulative_probability(const double& x)
	{
		return x < value ? 0 : 1;
	}

	/** {@inherit_doc} */
	//override
	public double inverse_cumulative_probability(const double& p)

	{
		Math_Utils::check_range_inclusive(p, 0, 1);
		return value;
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	public double get_numerical_mean() const
	{
		return value;
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	public double get_numerical_variance() const
	{
		return 0;
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	public double get_support_lower_bound() const
	{
		return value;
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	public double get_support_upper_bound() const
	{
		return value;
	}

	/**
	 * {@inherit_doc}
	 */
	 //override
	public bool is_support_connected() const
	{
		return true;
	}
}
