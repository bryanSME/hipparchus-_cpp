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

//import java.util.LinkedHash_Set;
//import java.util.Set;
//import java.util.Collection;
//import java.util.Collections;

//import org.hipparchus.optim.Optimization_data;

/**
 * Class that represents a set of {@link Linear_Constraint linear constraints}.
 *
 */
class Linear_ConstraintSet : Optimization_data 
{
    /** Set of constraints. */
    private const Set<Linear_Constraint> linear_constraints;

    /**
     * Creates a set containing the given constraints.
     *
     * @param constraints Constraints.
     */
    public Linear_ConstraintSet(Linear_Constraint... constraints) 
    {
        linear_constraints = LinkedHash_Set<>();
        for (Linear_Constraint c : constraints) 
        {
            linear_constraints.add(c);
        }
    }

    /**
     * Creates a set containing the given constraints.
     *
     * @param constraints Constraints.
     */
    public Linear_ConstraintSet(Collection<Linear_Constraint> constraints) 
    {
        linear_constraints = LinkedHash_Set<>();
        linear_constraints.add_all(constraints);
    }

    /**
     * Gets the set of linear constraints.
     *
     * @return the constraints.
     */
    public Collection<Linear_Constraint> get_constraints() 
    {
        return Collections.unmodifiable_set(linear_constraints);
    }
}


