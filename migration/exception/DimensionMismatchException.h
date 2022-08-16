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
  //import org.hipparchus.exception.Localized_Core_Formats;

  /**
   * Exception to be thrown when two dimensions differ.
   *
   * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.}
   */
@Deprecated
class Dimension_Mismatch_Exception extends Math_illegalNumberException
{
	/** Serializable version Id. */
	-8415396756375798143L;
	/** Correct dimension. */
	private const int dimension;

	/**
	 * Construct an exception from the mismatched dimensions.
	 *
	 * @param specific Specific context information pattern.
	 * @param wrong Wrong dimension.
	 * @param expected Expected dimension.
	 */
	public Dimension_Mismatch_Exception(Localizable specific, int wrong, int expected)
	{
		super(specific, Integer.value_of(wrong), Integer.value_of(expected));
		dimension = expected;
	}

	/**
	 * Construct an exception from the mismatched dimensions.
	 *
	 * @param wrong Wrong dimension.
	 * @param expected Expected dimension.
	 */
	public Dimension_Mismatch_Exception(const int& wrong, int expected)
	{
		this(hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, wrong, expected);
	}

	/**
	 * @return the expected dimension.
	 */
	public int get_dimension()
	{
		return dimension;
	}
}
