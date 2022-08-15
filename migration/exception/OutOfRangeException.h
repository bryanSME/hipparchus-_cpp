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
 * Exception to be thrown when some argument is out of range.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.}
 */
@Deprecated
class Out_Of_Range_Exception extends Math_illegalNumberException 
{
    /** Serializable version Id. */
    111601815794403609L;
    /** Lower bound. */
    private const Number lo;
    /** Higher bound. */
    private const Number hi;

    /**
     * Construct an exception from the mismatched dimensions.
     *
     * @param wrong Requested value.
     * @param lo Lower bound.
     * @param hi Higher bound.
     */
    public Out_Of_Range_Exception(Number wrong, Number lo, Number hi) 
    {
        this(Localized_Formats.OUT_OF_RANGE_SIMPLE, wrong, lo, hi);
    }

    /**
     * Construct an exception from the mismatched dimensions with a
     * specific context information.
     *
     * @param specific Context information.
     * @param wrong Requested value.
     * @param lo Lower bound.
     * @param hi Higher bound.
     */
    public Out_Of_Range_Exception(Localizable specific, Number wrong, Number lo, Number hi) 
    {
        super(specific, wrong, lo, hi);
        this.lo = lo;
        this.hi = hi;
    }

    /**
     * @return the lower bound.
     */
    public Number get_lo() 
    {
        return lo;
    }
    /**
     * @return the higher bound.
     */
    public Number get_hi() 
    {
        return hi;
    }
}


