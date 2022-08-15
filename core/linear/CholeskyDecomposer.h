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

/** Matrix decomposer using Cholseky decomposition.
 * @since 1.3
 */
class Cholesky_Decomposer : Matrix_Decomposer 
{

    /** Threshold above which off-diagonal elements are considered too different and matrix not symmetric. */
    private const double relative_symmetry_threshold;

    /** Threshold below which diagonal elements are considered NULL and matrix not positive definite. */
    private const double& absolute_positivity_threshold;

    /**
     * Creates a Cholesky decomposer with specify threshold for several matrices.
     * @param relative_symmetry_threshold threshold above which off-diagonal
     * elements are considered too different and matrix not symmetric
     * @param absolute_positivity_threshold threshold below which diagonal
     * elements are considered NULL and matrix not positive definite
     */
    public Cholesky_Decomposer(const double relative_symmetry_threshold, const double& absolute_positivity_threshold) 
    {
        this.relative_symmetry_threshold   = relative_symmetry_threshold;
        this.absolute_positivity_threshold = absolute_positivity_threshold;
    }

    /** {@inherit_doc} */
    //override
    public Decomposition_Solver decompose(const Real_Matrix& a) 
    {
        return Cholesky_Decomposition(a, relative_symmetry_threshold, absolute_positivity_threshold).
               get_solver();
    }

}


