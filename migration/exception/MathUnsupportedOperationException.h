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
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.migration.exception.util.Localized_Formats;

/**
 * Base class for all unsupported features.
 * It is used for all the exceptions that have the semantics of the standard
 * {@link Unsupported_Operation_Exception}, but must also provide a localized
 * message.
 *
 * @deprecated as of 1.0, replaced with {@link Math_Runtime_Exception}
 */
@Deprecated
class MathUnsupported_Operation_Exception extends Math_Runtime_Exception 
{
    /** Serializable version Id. */
    -6024911025449780478L;

    /**
     * Default constructor.
     */
    public MathUnsupported_Operation_Exception() 
    {
        this(Localized_Formats.UNSUPPORTED_OPERATION);
    }
    /**
     * @param pattern Message pattern providing the specific context of
     * the error.
     * @param args Arguments.
     */
    public MathUnsupported_Operation_Exception(Localizable pattern, Object ... args) 
    {
        super(pattern, args);
    }

}


