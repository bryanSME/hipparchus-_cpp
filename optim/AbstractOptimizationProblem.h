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

//import org.hipparchus.util.Incrementor;

/**
 * Base class for implementing optimization problems. It contains the boiler-plate code
 * for counting the number of evaluations of the objective function and the number of
 * iterations of the algorithm, and storing the convergence checker.
 *
 * @param <P> Type of the point/value pair returned by the optimization algorithm.
 */
class AbstractOptimization_Problem<P>
        : Optimization_Problem<P> 
        {

    /** max evaluations */
    private const int max_evaluations;
    /** max iterations */
    private const int max_iterations;
    /** Convergence checker. */
    private const Convergence_Checker<P> checker;

    /**
     * Create an {@link AbstractOptimization_Problem} from the given data.
     *
     * @param max_evaluations the number of allowed model function evaluations.
     * @param max_iterations  the number of allowed iterations.
     * @param checker        the convergence checker.
     */
    protected AbstractOptimization_Problem(const int max_evaluations, const int max_iterations, const Convergence_Checker<P> checker) 
    {
        this.max_evaluations = max_evaluations;
        this.max_iterations = max_iterations;
        this.checker = checker;
    }

    /** {@inherit_doc} */
    //override
    public Incrementor get_evaluation_counter() 
    {
        return Incrementor(this.max_evaluations);
    }

    /** {@inherit_doc} */
    //override
    public Incrementor get_iteration_counter() 
    {
        return Incrementor(this.max_iterations);
    }

    /** {@inherit_doc} */
    //override
    public Convergence_Checker<P> get_convergence_checker() 
    {
        return checker;
    }
}


