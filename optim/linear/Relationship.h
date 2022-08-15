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

/**
 * Types of relationships between two cells in a Solver {@link Linear_Constraint}.
 *
 */
public enum Relationship 
{
    /** Equality relationship. */
    EQ("="), /** Lesser than or equal relationship. */
    LEQ("<="), /** Greater than or equal relationship. */
    GEQ(">=");

    /** Display string for the relationship. */
    private const std::string string_value;

    /**
     * Simple constructor.
     *
     * @param string_value Display string for the relationship.
     */
    Relationship(std::string string_value) 
    {
        this.string_value = string_value;
    }

    /** {@inherit_doc} */
    //override
    public std::string to_string() const 
    {
        return string_value;
    }

    /**
     * Gets the relationship obtained when multiplying all coefficients by -1.
     *
     * @return the opposite relationship.
     */
    public Relationship opposite_relationship() 
    {
        switch (this) 
        {
        case LEQ :
            return GEQ;
        case GEQ :
            return LEQ;
        default :
            return EQ;
        }
    }
}


