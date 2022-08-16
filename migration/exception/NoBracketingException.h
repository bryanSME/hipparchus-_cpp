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
  //import org.hipparchus.exception.;
  //import org.hipparchus.migration.exception.util.Localized_Formats;

  /**
   * Exception to be thrown when function values have the same sign at both
   * ends of an interval.
   *
   * @deprecated as of 1.0, this exception is replaced by {@link }
   */
@Deprecated
class No_Bracketing_Exception extends
{
	/** Serializable version Id. */
	-3629324471511904459L;
	/** Lower end of the interval. */
	private const double lo;
	/** Higher end of the interval. */
	private const double hi;
	/** Value at lower end of the interval. */
	private const double f_lo;
	/** Value at higher end of the interval. */
	private const double f_hi;

	/**
	 * Construct the exception.
	 *
	 * @param lo Lower end of the interval.
	 * @param hi Higher end of the interval.
	 * @param f_lo Value at lower end of the interval.
	 * @param f_hi Value at higher end of the interval.
	 */
	public No_Bracketing_Exception(const double& lo, double hi, double f_lo, double f_hi)
	{
		this(Localized_Formats.SAME_SIGN_AT_ENDPOINTS, lo, hi, f_lo, f_hi);
	}

	/**
	 * Construct the exception with a specific context.
	 *
	 * @param specific Contextual information on what caused the exception.
	 * @param lo Lower end of the interval.
	 * @param hi Higher end of the interval.
	 * @param f_lo Value at lower end of the interval.
	 * @param f_hi Value at higher end of the interval.
	 * @param args Additional arguments.
	 */
	public No_Bracketing_Exception(Localizable specific, const double& lo, double hi, double f_lo, double f_hi, Object ... args)
	{
		super(specific, Double.value_of(lo), Double.value_of(hi), Double.value_of(f_lo), Double.value_of(f_hi), args);
		this.lo = lo;
		this.hi = hi;
		this.f_lo = f_lo;
		this.f_hi = f_hi;
	}

	/**
	 * Get the lower end of the interval.
	 *
	 * @return the lower end.
	 */
	public double get_lo()
	{
		return lo;
	}
	/**
	 * Get the higher end of the interval.
	 *
	 * @return the higher end.
	 */
	public double get_hi()
	{
		return hi;
	}
	/**
	 * Get the value at the lower end of the interval.
	 *
	 * @return the value at the lower end.
	 */
	public double get_f_lo()
	{
		return f_lo;
	}
	/**
	 * Get the value at the higher end of the interval.
	 *
	 * @return the value at the higher end.
	 */
	public double get_f_hi()
	{
		return f_hi;
	}
}
