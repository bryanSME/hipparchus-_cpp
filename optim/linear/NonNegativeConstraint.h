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
//package org.hipparchus.optim.linear;

//import org.hipparchus.optim.Optimization_data;

/**
 * A constraint for a linear optimization problem indicating whether all
 * variables must be restricted to non-negative values.
 *
 */
class Non_Negative_Constraint : Optimization_data 
{
    /** Whether the variables are all positive. */
    private const bool is_restricted;

    /**
     * @param restricted If {@code true}, all the variables must be positive.
     */
    public Non_Negative_Constraint(bool restricted) 
    {
        is_restricted = restricted;
    }

    /**
     * Indicates whether all the variables must be restricted to non-negative
     * values.
     *
     * @return {@code true} if all the variables must be positive.
     */
    public bool is_restricted_to_non_negative() 
    {
        return is_restricted;
    }
}


