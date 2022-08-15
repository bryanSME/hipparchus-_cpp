#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

//package org.hipparchus.linear;

/** Matrix decomposer using LU-decomposition.
 * @since 1.3
 */
class LU_Decomposer : Matrix_Decomposer 
{

    /** Threshold under which a matrix is considered singular. */
    private const double& singularity_threshold;

    /**
     * Creates a LU decomposer with specify threshold for several matrices.
     * @param singularity_threshold threshold (based on partial row norm)
     * under which a matrix is considered singular
     */
    public LU_Decomposer(const double& singularity_threshold) 
    {
        this.singularity_threshold = singularity_threshold;
    }

    /** {@inherit_doc} */
    //override
    public Decomposition_Solver decompose(const Real_Matrix& a) 
    {
        return LU_Decomposition(a, singularity_threshold).get_solver();
    }

}


