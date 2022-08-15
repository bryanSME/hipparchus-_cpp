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
//import org.hipparchus.migration.exception.util.Localized_Formats;

/**
 * Exception to be thrown when a number is not finite.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.}
 */
@Deprecated
class Not_Finite_Number_Exception extends Math_illegalNumberException 
{
    /** Serializable version Id. */
    -6100997100383932834L;

    /**
     * Construct the exception.
     *
     * @param wrong Value that is infinite or NaN.
     * @param args Optional arguments.
     */
    public Not_Finite_Number_Exception(Number wrong, Object ... args) 
    {
        this(Localized_Formats.NOT_FINITE_NUMBER, wrong, args);
    }

    /**
     * Construct the exception with a specific context.
     *
     * @param specific Specific context pattern.
     * @param wrong Value that is infinite or NaN.
     * @param args Optional arguments.
     */
    public Not_Finite_Number_Exception(Localizable specific, Number wrong, Object ... args) 
    {
        super(specific, wrong, args);
    }
}


