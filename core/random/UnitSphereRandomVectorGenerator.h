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
 * Generate random vectors isotropically located on the surface of a sphere.
 */
class UnitSphereRandom_Vector_Generator
    : Random_Vector_Generator 
    {

    /** RNG used for generating the individual components of the vectors. */
    private const Random_Generator rand;
    /** Space dimension. */
    private const int dimension;

    /**
     * @param dimension Space dimension.
     * @param rand RNG for the individual components of the vectors.
     */
    public UnitSphereRandom_Vector_Generator(const int dimension, const Random_Generator rand) 
    {
        this.dimension = dimension;
        this.rand = rand;
    }

    /**
     * Create an object that will use a default RNG ({@link Mersenne_Twister}), * in order to generate the individual components.
     *
     * @param dimension Space dimension.
     */
    public UnitSphereRandom_Vector_Generator(const int& dimension) 
    {
        this(dimension, Mersenne_Twister());
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> next_vector() 
    {
        const std::vector<double>& v = std::vector<double>(dimension];

        // See http://mathworld.wolfram.com/SpherePointPicking.html for example.
        // Pick a point by choosing a standard Gaussian for each element, and then
        // normalizing to unit length.
        double norm_sq = 0;
        for (int i{}; i < dimension; i++) 
        {
            const double comp = rand.next_gaussian();
            v[i] = comp;
            norm_sq += comp * comp;
        }

        const double f = 1 / std::sqrt(norm_sq);
        for (int i{}; i < dimension; i++) 
        {
            v[i] *= f;
        }

        return v;
    }
}


