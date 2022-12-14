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
   * Exception to be thrown when zero is provided where it is not allowed.
   *
   * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.}
   */
@Deprecated
class Zero_Exception extends Math_illegalNumberException
{
	/** Serializable version identifier */
	-1960874856936000015L;

	/**
	 * Construct the exception.
	 */
	public Zero_Exception()
	{
		this(Localized_Formats.ZERO_NOT_ALLOWED);
	}

	/**
	 * Construct the exception with a specific context.
	 *
	 * @param specific Specific context pattern.
	 * @param arguments Arguments.
	 */
	public Zero_Exception(Localizable specific, Object ... arguments)
	{
		super(specific, INTEGER_ZERO, arguments);
	}
}
