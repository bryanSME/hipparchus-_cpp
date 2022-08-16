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
  //package org.hipparchus.exception;

  //import java.text.Message_Format;
  //import java.util.Locale;
#include <vector>
#include <string>
#include <format>
#include "LocalizedException.h"
#include "LocalizedCoreFormats.h"
#include <type_traits>
#include "../CalculusFieldElement.hpp"

namespace hipparchus
{
	namespace exception
	{
		/**
		 * All conditions checks that fail due to a {@code NULL} argument must throw
		 * this exception.
		 * This class is meant to signal a precondition violation ("null is an illegal
		 * argument") and so does not extend the standard {@code Null_Pointer_Exception}.
		 * Propagation of {@code Null_Pointer_Exception} from within Hipparchus is
		 * construed to be a bug.
		 * <p>
		 * Note: from 1.0 onwards, this class extends {@link Null_Pointer_Exception} instead
		 * of {@link }.
		 *
		 */
		template <typename... Args>
		class : public Localized_Exception //, Null_Pointer_Exception
		{
		public:
			/**
			 * Default constructor.
			 */
			()
		{
			(Localized_Core_Formats::NULL_NOT_ALLOWED);
		}

			/** Simple constructor.
			 * @param specifier format specifier (to be translated).
			 * @param parts parts to insert in the format (no translation).
			 */
			(const Localizable specifier, const Args&&... parts)
			{
				my_specifier = specifier;
				my_parts = (std::is_null_pointer(parts)) ? Object[0] : parts.clone();
			}

			/** {@inherit_doc} */
			//override
			std::string get_message(const Locale& locale)
			{
				return build_message(locale, specifier, parts);
			}

			/** {@inherit_doc} */
			//override
			std::string get_message()
			{
				return get_message(Locale.US);
			}

			/** {@inherit_doc} */
			//override
			std::string get_localized_message()
			{
				return get_message(Locale.get_default());
			}

			/** {@inherit_doc} */
			//override
			Localizable get_specifier()
			{
				return specifier;
			}

			/** {@inherit_doc} */
			//override
			std::vector<Object> get_parts() const
			{
				return my_parts;
			}

		private:
			/** Format specifier (to be translated). */
			private const Localizable my_specifier;

			/** Parts to insert in the format (no translation). */
			private const std::vector<Object> my_parts;

			/**
			 * Builds a message string by from a pattern and its arguments.
			 * @param locale Locale in which the message should be translated
			 * @param specifier format specifier (to be translated)
			 * @param parts parts to insert in the format (no translation)
			 * @return a message string
			 */
			static std::string build_message(const Locale& locale, const Localizable specifier, const Args&&... parts)
			{
				return std::vformat(rt_fmt_str, std::make_format_args(args...));
				return (specifier == NULL) ? "" : Message_Format(specifier.get_localized_string(locale), locale).format(parts);
			}
		};
	};
};