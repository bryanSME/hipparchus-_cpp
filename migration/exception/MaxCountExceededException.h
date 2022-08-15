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
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.migration.exception.util.Localized_Formats;

/**
 * Exception to be thrown when some counter maximum value is exceeded.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link Math_Illegal_State_Exception}
 */
@Deprecated
class Max_Count_Exceeded_Exception extends Math_Illegal_State_Exception 
{
    /** Serializable version Id. */
    4330003017885151975L;
    /**
     * Maximum number of evaluations.
     */
    private const Number& max;

    /**
     * Construct the exception.
     *
     * @param max Maximum.
     */
    public Max_Count_Exceeded_Exception(const Number& max) 
    {
        this(Localized_Formats.MAX_COUNT_EXCEEDED, max);
    }
    /**
     * Construct the exception with a specific context.
     *
     * @param specific Specific context pattern.
     * @param max Maximum.
     * @param args Additional arguments.
     */
    public Max_Count_Exceeded_Exception(Localizable specific, const Number& max, Object ... args) 
    {
        super(specific, max, args);
        this.max = max;
    }

    /**
     * @return the maximum number of evaluations.
     */
    public Number get_max() 
    {
        return max;
    }
}


