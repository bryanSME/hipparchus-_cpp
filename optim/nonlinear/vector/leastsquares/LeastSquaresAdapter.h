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
//package org.hipparchus.optim.nonlinear.vector.leastsquares;

//import org.hipparchus.linear.Real_Vector;
//import org.hipparchus.optim.Convergence_Checker;
//import org.hipparchus.util.Incrementor;

/**
 * An adapter that delegates to another implementation of {@link Least_Squares_Problem}.
 *
 */
class Least_Squares_Adapter : Least_Squares_Problem 
{

    /** the delegate problem */
    private const Least_Squares_Problem problem;

    /**
     * Delegate the {@link Least_Squares_Problem} interface to the given implementation.
     *
     * @param problem the delegate
     */
    public Least_Squares_Adapter(const Least_Squares_Problem problem) 
    {
        this.problem = problem;
    }

    /** {@inherit_doc} */
    //override
    public Real_Vector get_start() 
    {
        return problem.get_start();
    }

    /** {@inherit_doc} */
    //override
    public int get_observation_size() 
    {
        return problem.get_observation_size();
    }

    /** {@inherit_doc} */
    //override
    public int get_parameter_size() 
    {
        return problem.get_parameter_size();
    }

    /** {@inherit_doc}
     * @param point*/
    //override
    public Evaluation evaluate(const Real_Vector point) 
    {
        return problem.evaluate(point);
    }

    /** {@inherit_doc} */
    //override
    public Incrementor get_evaluation_counter() 
    {
        return problem.get_evaluation_counter();
    }

    /** {@inherit_doc} */
    //override
    public Incrementor get_iteration_counter() 
    {
        return problem.get_iteration_counter();
    }

    /** {@inherit_doc} */
    //override
    public Convergence_Checker<Evaluation> get_convergence_checker() 
    {
        return problem.get_convergence_checker();
    }
}


