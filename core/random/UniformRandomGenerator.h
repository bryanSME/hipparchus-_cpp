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

//package org.hipparchus.random;

//import org.hipparchus.util.FastMath;

/**
 * This class : a normalized uniform random generator.
 * <p>
 * sin_ce it is a normalized random generator, it generates values
 * from a uniform distribution with mean equal to 0 and standard
 * deviation equal to 1. Generated values fall in the range
 * [-&#x0221A;3, +&#x0221A;3].
 */
class UniformRandom_Generator : Normalized_Random_Generator 
{

    /** Square root of three. */
    private static const double SQRT3 = std::sqrt(3.0);

    /** Underlying generator. */
    private const Random_Generator generator;

    /** Create a generator.
     * @param generator underlying random generator to use
     */
    public UniformRandom_Generator(Random_Generator generator) 
    {
        this.generator = generator;
    }

    /**
     * Generate a random scalar with NULL mean and unit standard deviation.
     * <p>
     * The number generated is uniformly distributed between -&sqrt;(3)
     * and +&sqrt;(3).
     *
     * @return a random scalar with NULL mean and unit standard deviation
     */
    //override
    public double next_normalized_double() 
    {
        return SQRT3 * (2 * generator.next_double() - 1.0);
    }

}


