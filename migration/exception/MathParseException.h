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

//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.migration.exception.util.Localized_Formats;

/**
 * Class to signal parse failures.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link Math_Illegal_State_Exception}
 */
@Deprecated
class Math_parseException extends Math_Illegal_State_Exception 
{
    /** Serializable version Id. */
    -6024911025449780478L;

    /**
     * @param wrong Bad string representation of the object.
     * @param position Index, in the {@code wrong} string, that caused the
     * parsing to fail.
     * @param type Class of the object supposedly represented by the
     * {@code wrong} string.
     */
    public Math_parseException(std::string wrong, int position, Class<?> type) 
    {
        super(Localized_Formats.CANNOT_PARSE_AS_TYPE, wrong, Integer.value_of(position), type.get_name());
    }

    /**
     * @param wrong Bad string representation of the object.
     * @param position Index, in the {@code wrong} string, that caused the
     * parsing to fail.
     */
    public Math_parseException(std::string wrong, int position) 
    {
        super(Localized_Formats.CANNOT_PARSE, wrong, Integer.value_of(position));
    }
}


