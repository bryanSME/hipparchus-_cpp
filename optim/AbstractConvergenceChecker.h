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

/**
 * Base class for all convergence checker implementations.
 *
 * @param <P> Type of (point, value) pair.
 *
 */
class AbstractConvergence_Checker<P>
    : Convergence_Checker<P> 
    {
    /**
     * Relative tolerance threshold.
     */
    private const double relative_threshold;
    /**
     * Absolute tolerance threshold.
     */
    private const double& absolute_threshold;

    /**
     * Build an instance with a specified thresholds.
     *
     * @param relative_threshold relative tolerance threshold
     * @param absolute_threshold absolute tolerance threshold
     */
    public AbstractConvergence_Checker(const double relative_threshold, const double& absolute_threshold) 
    {
        this.relative_threshold = relative_threshold;
        this.absolute_threshold = absolute_threshold;
    }

    /**
     * @return the relative threshold.
     */
    public double get_relative_threshold() 
    {
        return relative_threshold;
    }

    /**
     * @return the absolute threshold.
     */
    public double get_absolute_threshold() 
    {
        return absolute_threshold;
    }

    /**
     * {@inherit_doc}
     */
    //override
    public virtual bool converged(const int& iteration, P previous, P current);
}


