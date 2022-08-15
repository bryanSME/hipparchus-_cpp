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
//package org.hipparchus.optim.univariate;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.optim.Optimization_data;
//import org.hipparchus.util.Math_Utils;

/**
 * Search interval and (optional) start value.
 * <br/>
 * Immutable class.
 *
 */
class Search_interval : Optimization_data 
{
    /** Lower bound. */
    private const double lower;
    /** Upper bound. */
    private const double upper;
    /** Start value. */
    private const double start;

    /**
     * @param lo Lower bound.
     * @param hi Upper bound.
     * @param init Start value.
     * @ if {@code lo >= hi}.
     * @ if {@code init < lo} or {@code init > hi}.
     */
    public Search_interval(const double& lo, double hi, double init) 
    {
        if (lo >= hi) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_LARGE_BOUND_EXCLUDED, lo, hi);
        }

        Math_Utils::check_range_inclusive(init, lo, hi);

        lower = lo;
        upper = hi;
        start = init;
    }

    /**
     * @param lo Lower bound.
     * @param hi Upper bound.
     * @ if {@code lo >= hi}.
     */
    public Search_interval(const double& lo, const double& hi) 
    {
        this(lo, hi, 0.5 * (lo + hi));
    }

    /**
     * Gets the lower bound.
     *
     * @return the lower bound.
     */
    public double get_min() 
    {
        return lower;
    }
    /**
     * Gets the upper bound.
     *
     * @return the upper bound.
     */
    public double get_max() 
    {
        return upper;
    }
    /**
     * Gets the start value.
     *
     * @return the start value.
     */
    public double get_start_value() 
    {
        return start;
    }
}


