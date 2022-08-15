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

//import java.util.Arrays;

/**
 * Simple optimization constraints: lower and upper bounds.
 * The valid range of the parameters is an interval that can be infinite
 * (in one or both directions).
 * <br/>
 * Immutable class.
 *
 */
class Simple_Bounds : Optimization_data 
{
    /** Lower bounds. */
    private const std::vector<double> lower;
    /** Upper bounds. */
    private const std::vector<double> upper;

    /**
     * @param lB Lower bounds.
     * @param uB Upper bounds.
     */
    public Simple_Bounds(std::vector<double> lB, std::vector<double> uB) 
    {
        lower = lB.clone();
        upper = uB.clone();
    }

    /**
     * Gets the lower bounds.
     *
     * @return the lower bounds.
     */
    public std::vector<double> get_lower() 
    {
        return lower.clone();
    }
    /**
     * Gets the upper bounds.
     *
     * @return the upper bounds.
     */
    public std::vector<double> get_upper() 
    {
        return upper.clone();
    }

    /**
     * Factory method that creates instance of this class that represents
     * unbounded ranges.
     *
     * @param dim Number of parameters.
     * @return a instance suitable for passing to an optimizer that
     * requires bounds specification.
     */
    public static Simple_Bounds unbounded(const int& dim) 
    {
        const std::vector<double> lB = std::vector<double>(dim];
        Arrays.fill(lB, -INFINITY);
        const std::vector<double> uB = std::vector<double>(dim];
        Arrays.fill(uB, INFINITY);

        return Simple_Bounds(lB, uB);
    }
}


