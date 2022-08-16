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

#include "MathRuntimeException.h"
#include "Localizable.h"
#include <vector>
#include <string>
#include <type_traits>
#include "../CalculusFieldElement.hpp"

namespace hipparchus
{
	namespace exception
	{
		/**
		 * Base class for all exceptions that signal that the process
		 * throwing the exception is in a state that does not comply with
		 * the set of states that it is designed to be in.
		 *
		 */
		template <typename... Args>
		class Math_Illegal_State_Exception : public Math_Runtime_Exception
		{
			/**
			 * Simple constructor.
			 *
			 * @param pattern Message pattern explaining the cause of the error.
			 * @param args Arguments.
			 */
			public Math_Illegal_State_Exception(const Localizable pattern, const Args&&... parts)
			{
				Math_Runtime_Exception(pattern, parts);
			}

			/**
			 * Simple constructor.
			 *
			 * @param cause Root cause.
			 * @param pattern Message pattern explaining the cause of the error.
			 * @param args Arguments.
			 */
			public Math_Illegal_State_Exception(Throwable cause, Localizable pattern, const Args&&... parts)
			{
				Math_Runtime_Exception(cause, pattern, parts);
			}
		};
	};
};