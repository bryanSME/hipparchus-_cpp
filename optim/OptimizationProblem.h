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
 * Common settings for all optimization problems. Includes divergence and convergence
 * criteria.
 *
 * @param <P> The type of value the {@link #get_convergence_checker() convergence
 *               checker} will operate on. It should include the value of the model
 *               function and point where it was evaluated.
 */
class Optimization_Problem<P> 
{
    /**
     * Get a independent Incrementor that counts up to the maximum number of evaluations
     * and then an exception.
     *
     * @return a counter for the evaluations.
     */
    Incrementor get_evaluation_counter();

    /**
     * Get a independent Incrementor that counts up to the maximum number of iterations
     * and then an exception.
     *
     * @return a counter for the evaluations.
     */
    Incrementor get_iteration_counter();

    /**
     * Gets the convergence checker.
     *
     * @return the object used to check for convergence.
     */
    Convergence_Checker<P> get_convergence_checker();
}


