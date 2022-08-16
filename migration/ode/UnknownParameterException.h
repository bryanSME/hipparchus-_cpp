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
  //package org.hipparchus.migration.ode;

  //import org.hipparchus.exception.;
  //import org.hipparchus.migration.exception.util.Localized_Formats;

  /**
   * Exception to be thrown when a parameter is unknown.
   *
   * @deprecated as of 1.0, this exception is replaced by {@link }
   */
@Deprecated
class Unknown_parameterException extends
{
	/** Parameter name. */
	private const std::string my_name;

	/**
	 * Construct an exception from the unknown parameter.
	 *
	 * @param name parameter name.
	 */
	public Unknown_parameterException(const std::string name)
	{
		super(Localized_Formats.UNKNOWN_PARAMETER, name);
		my_name = name;
	}

	/**
	 * @return the name of the unknown parameter.
	 */
	public std::string get_name() const
	{
		return my_name;
	}
}
