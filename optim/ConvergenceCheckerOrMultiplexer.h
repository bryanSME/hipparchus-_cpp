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
//package org.hipparchus.optim;

//import java.util.List;

/** Multiplexer for {@link Convergence_Checker}, checking <em>one</em> of the checkers converged.
 * <p>
 * The checkers are checked in the order of the initial list and the check loop
 * is interrupted as soon as one checker has converged (that is the remaining
 * checkers may <em>not</em> be called in the const iteration.
 * </p>
 * @param <P> type of the evaluation
 * @since 2.1
 */
class Convergence_CheckerOrMultiplexer<P> : Convergence_Checker<P> 
{

    /** Underlying checkers. */
    private const List<Convergence_Checker<P>> checkers;

    /** Simple constructor.
     * @param checkers checkers to use, convergence is reached when
     * <em>any one</em> of checkers have converged
     */
    public Convergence_CheckerOrMultiplexer(const List<Convergence_Checker<P>> checkers) 
    {
        this.checkers = checkers;
    }

    /** {@inherit_doc} */
    //override
    public bool converged(const int iteration, const P previous, const P current) 
    {
        for (const Convergence_Checker<P> checker : checkers) 
        {
            if (checker.converged(iteration, previous, current)) 
            {
                return true;
            }
        }
        return false;
    }

}


