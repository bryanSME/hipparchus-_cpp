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

/**
 * An algorithm that can be applied to a non-linear least squares problem.
 */
class Least_Squares_Optimizer 
{

    /**
     * Solve the non-linear least squares problem.
     *
     *
     * @param least_squares_problem the problem definition, including model function and
     *                            convergence criteria.
     * @return The optimum.
     */
    Optimum optimize(Least_Squares_Problem least_squares_problem);

    /**
     * The optimum found by the optimizer. This object contains the point, its value, and
     * some metadata.
     */
    interface Optimum extends Least_Squares_Problem.Evaluation 
    {

        /**
         * Get the number of times the model was evaluated in order to produce this
         * optimum.
         *
         * @return the number of model (objective) function evaluations
         */
        int get_evaluations();

        /**
         * Get the number of times the algorithm iterated in order to produce this
         * optimum. In general least squares it is common to have one {@link
         * #get_evaluations() evaluation} per iterations.
         *
         * @return the number of iterations
         */
        int get_iterations();

        /**
         * Create a optimum from an evaluation and the values of the counters.
         *
         * @param value       the function value
         * @param evaluations number of times the function was evaluated
         * @param iterations  number of iterations of the algorithm
         * @return a optimum based on the given data.
         */
        static Optimum of(const Least_Squares_Problem.Evaluation value, const int evaluations, const int iterations) 
        {
            return Optimum_Impl(value, evaluations, iterations);
        }

    }

}


