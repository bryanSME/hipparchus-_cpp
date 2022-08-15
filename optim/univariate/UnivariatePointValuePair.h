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

//import java.io.Serializable;

/**
 * This class holds a point and the value of an objective function at this
 * point.
 * This is a simple immutable container.
 *
 */
class UnivariatePoint_valuePair  
{
    /** Serializable version identifier. */
    private static const long serial_version_uid = 1003888396256744753L;
    /** Point. */
    private const double point;
    /** Value of the objective function at the point. */
    private const double value;

    /**
     * Build a point/objective function value pair.
     *
     * @param point Point.
     * @param value Value of an objective function at the point
     */
    public UnivariatePoint_valuePair(const double point, const double value) 
    {
        this.point = point;
        this.value = value;
    }

    /**
     * Get the point.
     *
     * @return the point.
     */
    public double get_point() 
    {
        return point;
    }

    /**
     * Get the value of the objective function.
     *
     * @return the stored value of the objective function.
     */
    public double get_value() 
    {
        return value;
    }
}


