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
#include <type_traits>
#include "../CalculusFieldElement.hpp"

namespace hipparchus
{
    namespace exception
    {
        /**
         * Base class for all preconditions violation exceptions.
         * In most cases, this class should not be instantiated directly: it should
         * serve as a base class to create all the exceptions that have the semantics
         * of the standard {@link Illegal_Argument_Exception}.
         *
         */
        template<typename T>
        class Math_Illegal_Argument_Exception : public Math_Runtime_Exception, public std::exception
        {
        private:
            char* my_message;

        public:
            /**
             * @param pattern Message pattern explaining the cause of the error.
             * @param args Arguments.
             */
            char* what(Localizable pattern, T args)
            {
                super(pattern, args);
            }

            Math_Illegal_Argument_Exception(const char* msg) : my_message{ msg } {};

            char* what()
            {
                return my_message;
            }

            /** Simple constructor.
             * @param cause root cause.
             * @param specifier format specifier (to be translated).
             * @param parts parts to insert in the format (no translation).
             * @since 1.4
             */
            public (const Throwable cause, const Localizable& specifier, const Object ... parts)
            {
                super(cause, specifier, parts);
            }

        };
    };
};