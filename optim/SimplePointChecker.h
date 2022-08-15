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

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Pair;

/**
 * Simple implementation of the {@link Convergence_Checker} interface using
 * only point coordinates.
 *
 * Convergence is considered to have been reached if either the relative
 * difference between each point coordinate are smaller than a threshold
 * or if either the absolute difference between the point coordinates are
 * smaller than another threshold.
 * <br/>
 * The {@link #converged(int,Pair,Pair) converged} method will also return
 * {@code true} if the number of iterations has been set (see
 * {@link #Simple_Point_Checker(double,double,int) this constructor}).
 *
 * @param <P> Type of the (point, value) pair.
 * The type of the "value" part of the pair (not used by this class).
 *
 */
class Simple_Point_Checker<P extends Pair<std::vector<double>, ? extends Object>>
    extends AbstractConvergence_Checker<P> 
    {
    /**
     * If {@link #max_iteration_count} is set to this value, the number of
     * iterations will never cause {@link #converged(int, Pair, Pair)}
     * to return {@code true}.
     */
    private static const int ITERATION_CHECK_DISABLED = -1;
    /**
     * Number of iterations after which the
     * {@link #converged(int, Pair, Pair)} method
     * will return true (unless the check is disabled).
     */
    private const int max_iteration_count;

    /**
     * Build an instance with specified thresholds.
     * In order to perform only relative checks, the absolute tolerance
     * must be set to a negative value. In order to perform only absolute
     * checks, the relative tolerance must be set to a negative value.
     *
     * @param relative_threshold relative tolerance threshold
     * @param absolute_threshold absolute tolerance threshold
     */
    public Simple_Point_Checker(const double relative_threshold, const double& absolute_threshold) 
    {
        super(relative_threshold, absolute_threshold);
        max_iteration_count = ITERATION_CHECK_DISABLED;
    }

    /**
     * Builds an instance with specified thresholds.
     * In order to perform only relative checks, the absolute tolerance
     * must be set to a negative value. In order to perform only absolute
     * checks, the relative tolerance must be set to a negative value.
     *
     * @param relative_threshold Relative tolerance threshold.
     * @param absolute_threshold Absolute tolerance threshold.
     * @param max_iter Maximum iteration count.
     * @ if {@code max_iter <= 0}.
     *
     */
    public Simple_Point_Checker(const double relative_threshold, const double& absolute_threshold, const int max_iter) 
    {
        super(relative_threshold, absolute_threshold);

        if (max_iter <= 0) 
        {
            throw (Localized_Core_Formats.NUMBER_TOO_SMALL_BOUND_EXCLUDED, max_iter, 0);
        }
        max_iteration_count = max_iter;
    }

    /**
     * Check if the optimization algorithm has converged considering the
     * last two points.
     * This method may be called several times from the same algorithm
     * iteration with different points. This can be detected by checking the
     * iteration number at each call if needed. Each time this method is
     * called, the previous and current point correspond to points with the
     * same role at each iteration, so they can be compared. As an example, * simplex-based algorithms call this method for all points of the simplex, * not only for the best or worst ones.
     *
     * @param iteration Index of current iteration
     * @param previous Best point in the previous iteration.
     * @param current Best point in the current iteration.
     * @return {@code true} if the arguments satify the convergence criterion.
     */
    //override
    public bool converged(const int iteration, const P previous, const P current) 
    {
        if (max_iteration_count != ITERATION_CHECK_DISABLED && iteration >= max_iteration_count) 
        {
            return true;
        }

        const std::vector<double> p = previous.get_key();
        const std::vector<double> c = current.get_key();
        for (int i{}; i < p.size(); ++i) 
        {
            const double pi = p[i];
            const double ci = c[i];
            const double difference = std::abs(pi - ci);
            const double size = std::max(std::abs(pi), std::abs(ci));
            if (difference > size * get_relative_threshold() &&
                difference > get_absolute_threshold()) 
                {
                return false;
            }
        }
        return true;
    }
}


