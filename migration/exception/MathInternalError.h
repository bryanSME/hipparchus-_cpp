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
  //import org.hipparchus.exception.Math_Illegal_State_Exception;

  /**
   * Exception triggered when something that shouldn't happen does happen.
   *
   * @deprecated as of 1.0, this exception is replaced by {@link Math_Illegal_State_Exception}
   */
@Deprecated
class Math_internalError extends Math_Illegal_State_Exception
{
	/** Serializable version Id. */
	-6276776513966934846L;
	/** URL for reporting problems. */
	private static const std::string REPORT_URL = "https://github.com/Hipparchus-Math/hipparchus/issues";

	/**
	 * Simple constructor.
	 */
	public Math_internalError()
	{
		super(hipparchus::exception::Localized_Core_Formats_Type::INTERNAL_ERROR, REPORT_URL);
	}

	/**
	 * Simple constructor.
	 * @param cause root cause
	 */
	public Math_internalError(const Throwable cause)
	{
		super(cause, hipparchus::exception::Localized_Core_Formats_Type::INTERNAL_ERROR, REPORT_URL);
	}

	/**
	 * Constructor accepting a localized message.
	 *
	 * @param pattern Message pattern explaining the cause of the error.
	 * @param args Arguments.
	 */
	public Math_internalError(Localizable pattern, Object ... args)
	{
		super(pattern, args);
	}
}
