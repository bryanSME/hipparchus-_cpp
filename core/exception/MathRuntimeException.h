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
#include <string>
#include <vector>
#include "Localizable.h"
#include <locale>
#include "LocalizedException.h"
#include "LocalizedCoreFormats.h"
#include <format>
#include <type_traits>
#include "../CalculusFieldElement.hpp"

namespace hipparchus
{
    namespace exception
    {
        /**
         * All exceptions thrown by the Hipparchus code inherit from this class.
         */

        template <typename... Args>
        class Math_Runtime_Exception : public Localized_Exception //, public Runtime_Exception,
        {
        private:
            /** URL for reporting problems for internal errors. */
            inline const static std::string REPORT_URL = ""; // https://github.com/Hipparchus-Math/hipparchus/issues

            /** Format specifier (to be translated). */
            Localizable my_specifier;

            /** Parts to insert in the format (no translation). */
            inline static Args&&... my_parts;

            /**
             * Builds a message string by from a pattern and its arguments.
             * @param locale Locale in which the message should be translated
             * @return a message string
             */
             //@Suppress_Warnings("PMD.Avoid_Catching_Generic_Exception") // catching Exception is intentional here
            template <typename... Args>
            std::string build_message(const std::locale& locale)
            {
                /*if (my_specifier == NULL)
                {
                    return "";
                }*/
                // CHECKSTYLE: stop Illegal_Catch
                try
                {
                    return std::vformat(rt_fmt_str, std::make_format_args(args...));
                }
                catch (std::exception& e)
                {
                    
                    this.add_suppressed(e);
                    return my_specifier.get_source_string();
                }
                // CHECKSTYLE: resume Illegal_Catch
            }

        public:
            /** Simple constructor.
             * @param specifier format specifier (to be translated).
             * @param parts parts to insert in the format (no translation).
             */
            template <typename... Args>
            Math_Runtime_Exception(const Localizable& specifier, const Args&&... parts) : my_specifier{ specifier }, my_parts{ parts } {};

            /** Simple constructor.
             * @param cause root cause.
             * @param specifier format specifier (to be translated).
             * @param parts parts to insert in the format (no translation).
             */
            template <typename... Args>
            Math_Runtime_Exception(const char* cause, const Localizable& specifier, const Args&&... parts) : my_specifier{ specifier }, my_parts{ parts }
            {
                Localized_Exception(cause);
            }

            /** Create an exception for an internal error.
             * @return a runtime exception indicating an internal error
             */
            static Math_Runtime_Exception create_internal_error()
            {
                return Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::INTERNAL_ERROR, REPORT_URL);
            }

            /** Create an exception for an internal error.
             * @param cause root cause
             * @return a runtime exception, indicating an internal error and wrapping the
             * given throwable
             */
            static Math_Runtime_Exception create_internal_error(const char* cause)
            {
                return Math_Runtime_Exception(cause, Localized_Core_Formats_Type::INTERNAL_ERROR, REPORT_URL);
            }

            /** {@inherit_doc} */
            //override
            std::string get_message(const std::locale& locale)
            {
                return build_message(locale);
            }

            /** {@inherit_doc} */
            //override
            std::string get_message()
            {
                return get_message(std::locale("en_US"));
            }

            /** {@inherit_doc} */
            //override
            std::string get_localized_message()
            {
                return get_message(std::locale());
            }

            /** {@inherit_doc} */
            //override
            Localizable get_specifier() const
            {
                return my_specifier;
            }

            /** {@inherit_doc} */
            //override
            std::vector<std::string> get_parts() const
            {
                return my_parts;
            }
        };
    };
};