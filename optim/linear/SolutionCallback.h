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
//import org.hipparchus.optim.Point_valuePair;

/**
 * A callback object that can be provided to a linear optimizer to keep track
 * of the best solution found.
 *
 */
class Solution_callback : Optimization_data 
{
    /** The Simplex_Tableau used by the Simplex_Solver. */
    private Simplex_Tableau tableau;

    /**
     * Set the simplex tableau used during the optimization once a feasible
     * solution has been found.
     *
     * @param tableau the simplex tableau containing a feasible solution
     */
    void set_tableau(const Simplex_Tableau tableau) 
    {
        this.tableau = tableau;
    }

    /**
     * Retrieve the best solution found so far.
     * <p>
     * <b>Note:</b> the returned solution may not be optimal, e.g. in case
     * the optimizer did reach the iteration limits.
     *
     * @return the best solution found so far by the optimizer, or {@code NULL} if
     * no feasible solution could be found
     */
    public Point_valuePair get_solution() 
    {
        return tableau != NULL ? tableau.get_solution() : NULL;
    }

    /**
     * Returns if the found solution is optimal.
     * @return {@code true} if the solution is optimal, {@code false} otherwise
     */
    public bool is_solution_optimal() 
    {
        return tableau != NULL && tableau.is_optimal();
    }
}


