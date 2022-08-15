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
//package org.hipparchus.optim;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;

/**
 * Maximum number of evaluations of the function to be optimized.
 *
 */
class Max_Eval : Optimization_data 
{
    /** Allowed number of evalutations. */
    private const int max;

    /**
     * @param max Allowed number of evalutations.
     * @ if {@code max <= 0}.
     */
    public Max_Eval(const int& max) 
    {
        if (max <= 0) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL_BOUND_EXCLUDED, max, 0);
        }

        this.max = max;
    }

    /**
     * Gets the maximum number of evaluations.
     *
     * @return the allowed number of evaluations.
     */
    public int get_max_eval() 
    {
        return max;
    }

    /**
     * Factory method that creates instance of this class that represents
     * a virtually unlimited number of evaluations.
     *
     * @return a instance suitable for allowing {@link Integer#MAX_VALUE}
     * evaluations.
     */
    public static Max_Eval unlimited() 
    {
        return Max_Eval(Integer.MAX_VALUE);
    }
}


