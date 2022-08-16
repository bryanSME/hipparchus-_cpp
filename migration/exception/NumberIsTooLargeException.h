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
  //package org.hipparchus.migration.exception;

  //import org.hipparchus.exception.Localizable;
  //import org.hipparchus.migration.exception.util.Localized_Formats;

  /**
   * Exception to be thrown when a number is too large.
   *
   * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.}
   */
@Deprecated
class Number_Is_Too_Large_Exception extends Math_illegalNumberException
{
	/** Serializable version Id. */
	4330003017885151975L;
	/**
	 * Higher bound.
	 */
	private const Number& max;
	/**
	 * Whether the maximum is included in the allowed range.
	 */
	private const bool bound_is_allowed;

	/**
	 * Construct the exception.
	 *
	 * @param wrong Value that is larger than the maximum.
	 * @param max Maximum.
	 * @param bound_is_allowed if true the maximum is included in the allowed range.
	 */
	public Number_Is_Too_Large_Exception(Number wrong, const Number& max, bool bound_is_allowed)
	{
		this(bound_is_allowed ?
			Localized_Formats.NUMBER_TOO_LARGE :
			Localized_Formats.NUMBER_TOO_LARGE_BOUND_EXCLUDED, wrong, max, bound_is_allowed);
	}
	/**
	 * Construct the exception with a specific context.
	 *
	 * @param specific Specific context pattern.
	 * @param wrong Value that is larger than the maximum.
	 * @param max Maximum.
	 * @param bound_is_allowed if true the maximum is included in the allowed range.
	 */
	public Number_Is_Too_Large_Exception(Localizable specific, Number wrong, const Number& max, bool bound_is_allowed)
	{
		super(specific, wrong, max);

		this.max = max;
		this.bound_is_allowed = bound_is_allowed;
	}

	/**
	 * @return {@code true} if the maximum is included in the allowed range.
	 */
	public bool get_bound_is_allowed() { // NOPMD - this method name is for a legacy API we cannot change
		return bound_is_allowed;
	}

	/**
	 * @return the maximum.
	 */
	public Number get_max()
	{
		return max;
	}
}
