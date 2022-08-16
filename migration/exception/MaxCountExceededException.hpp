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

#include <type_traits>
#include "../LocalizedMigrationFormats.h"
#include "../../core/exception/Localizable.h"
//import org.hipparchus.exception.Localizable;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.migration.exception.util.Localized_Formats;

/**
 * Exception to be thrown when some counter maximum value is exceeded.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link Math_Illegal_State_Exception}
 */
 //@Deprecated
template<
	typename T, //real type
	typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type
>
class Max_Count_Exceeded_Exception //: public Math_Illegal_State_Exception
{
private:
	/**
	 * Maximum number of evaluations.
	 */
	const T max;

public:
	/**
	 * Construct the exception.
	 *
	 * @param max Maximum.
	 */
	Max_Count_Exceeded_Exception(const T& max)
	{
		Max_Count_Exceeded_Exception(Localized_Formats::MAX_COUNT_EXCEEDED, max);
	}
	/**
	 * Construct the exception with a specific context.
	 *
	 * @param specific Specific context pattern.
	 * @param max Maximum.
	 * @param args Additional arguments.
	 */
	Max_Count_Exceeded_Exception(const Localizable& specific, const T& max, Object ... args)
		:
		my_max{ max }
	{
		super(specific, max, args);
	}

	/**
	 * @return the maximum number of evaluations.
	 */
	T get_max() const
	{
		return my_max;
	}
};