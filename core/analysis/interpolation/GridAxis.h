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
//package org.hipparchus.analysis.interpolation;
#include <vector>
#include <atomic>
#include "../../util/MathArrays.h"
//import java.io.Serializable;
//import java.util.concurrent.atomic.Atomic_Integer;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;

/**
 * Helper for finding interpolation nodes along one axis of grid data.
 * <p>
 * This class is intended to be used for interpolating inside grids.
 * It works on any sorted data without duplication and size at least
 * {@code n} where {@code n} is the number of points required for
 * interpolation (i.e. 2 for linear interpolation, 3 for quadratic...)
 * <p>
 * </p>
 * The method uses linear interpolation to select the nodes indices.
 * It should be O(1) for sufficiently regular data, therefore much faster
 * than bisection. It also features caching, which improves speed when
 * interpolating several points in raw in the close locations, i.e. when
 * successive calls have a high probability to return the same interpolation
 * nodes. This occurs for example when scanning with small steps a loose
 * grid. The method also works on non-regular grids, but may be slower in
 * this case.
 * </p>
 * <p>
 * This class is thread-safe.
 * </p>
 * @since 1.4
 */
class Grid_Axis  
{
private:
    /** All the coordinates of the interpolation points, sorted in increasing order. */
    std::vector<double> my_grid;

    /** Number of points required for interpolation. */
    int my_n;

    /** Cached value of last x index. */
    std::atomic<int> my_cache;

public:
    /** Simple constructor.
     * @param grid coordinates of the interpolation points, sorted in increasing order
     * @param n number of points required for interpolation, i.e. 2 for linear, 3
     * for quadratic...
     * @exception  if grid size is smaller than {@code n}
     * or if the grid is not sorted in strict increasing order
     */
    Grid_Axis(const std::vector<double>& grid, const int& n)
        : my_grid{ grid }, my_n{ n }, my_cache{ std::atomic<int>(0) }
    {
        //// safety checks
        //if (my_grid.size() < n) 
        //{
        //    throw (hipparchus::exception::Localized_Core_Formats_Type::INSUFFICIENT_DIMENSION, grid.size(), n);
        //}
        //Math_Arrays::check_order(grid);
    }

    /** Get the number of points of the grid.
     * @return number of points of the grid
     */
    int size() const 
    {
        return my_grid.size();
    }

    /** Get the number of points required for interpolation.
     * @return number of points required for interpolation
     */
    int get_n() const 
    {
        return my_n;
    }

    /** Get the interpolation node at specified index.
     * @param index node index
     * @return coordinate of the node at specified index
     */
    double node(const int& index) const
    {
        return my_grid[index];
    }

    /** Get the index of the first interpolation node for some coordinate along the grid.
     * <p>
     * The index return is the one for the lowest interpolation node suitable for
     * {@code t}. This means that if {@code i} is returned the nodes to use for
     * interpolation at coordinate {@code t} are at indices {@code i}, {@code i+1}, * ..., {@code i+n-1}, where {@code n} is the number of points required for
     * interpolation passed at construction.
     * </p>
     * <p>
     * The index is selected in order to have the subset of nodes from {@code i} to
     * {@code i+n-1} as balanced as possible around {@code t}:
     * </p>
     * <ul>
     *   <li>
     *     if {@code t} is inside the grid and sufficiently far from the endpoints
     *     <ul>
     *       <li>
     *         if {@code n} is even, the returned nodes will be perfectly balanced:
     *         there will be {@code n/2} nodes smaller than {@code t} and {@code n/2}
     *         nodes larger than {@code t}
     *       </li>
     *       <li>
     *         if {@code n} is odd, the returned nodes will be slightly unbalanced by
     *         one point: there will be {@code (n+1)/2} nodes smaller than {@code t}
     *         and {@code (n-1)/2} nodes larger than {@code t}
     *       </li>
     *     </ul>
     *   </li>
     *   <li>
     *     if {@code t} is inside the grid and close to endpoints, the returned nodes
     *     will be unbalanced: there will be less nodes on the endpoints side and
     *     more nodes on the interior side
     *   </li>
     *   <li>
     *     if {@code t} is outside of the grid, the returned nodes will completely
     *     off balance: all nodes will be on the same side with respect to {@code t}
     *   </li>
     * </ul>
     * <p>
     * It is <em>not</em> an error to call this method with {@code t} outside of the grid, * it simply implies that the interpolation will become an extrapolation and accuracy
     * will decrease as {@code t} goes farther from the grid points. This is intended so
     * interpolation does not fail near the end of the grid.
     * </p>
     * @param t coordinate of the point to interpolate
     * @return index {@code i} such {@link #nodestatic_cast<int>( node(i)}, {@link #nodestatic_cast<int>( node(i+1)}, * ... {@link #nodestatic_cast<int>( node(i+n-1)} can be used for interpolating a value at
     * coordinate {@code t}
     * @since 1.4
     */
    int interpolation_index(const double& t) 
    {

        const int middle_offset = (my_n - 1) / 2;
        int i_inf = middle_offset;
        int i_sup = my_grid.size() - (my_n - 1) + middle_offset;

        // first try to simply reuse the cached index, // for faster return in a common case
        auto    cached = my_cache.load(); // .get();
        const auto    middle = cached + middle_offset;
        const double& a_mid0  = my_grid[middle];
        const double& a_mid1  = my_grid[middle + 1];
        if (t < a_mid0) 
        {
            if (middle == i_inf) 
            {
                // we are in the unbalanced low area
                return cached;
            }
        }
        else if (t < a_mid1) 
        {
            // we are in the balanced middle area
            return cached;
        }
        else 
        {
            if (middle == i_sup - 1) 
            {
                // we are in the unbalanced high area
                return cached;
            }
        }

        // we need to find a index
        auto a_inf = my_grid[i_inf];
        auto a_sup = my_grid[i_sup];
        while (i_sup - i_inf > 1) 
        {
            const int i_interp = static_cast<int>( ((i_inf * (a_sup - t) + i_sup * (t - a_inf)) / (a_sup - a_inf));
            const int i_med    = std::max(i_inf + 1, std::min(i_interp, i_sup - 1));
            if (t < my_grid[i_med]) 
            {
                // keeps looking in the lower part of the grid
                i_sup = i_med;
                a_sup = my_grid[i_sup];
            }
            else 
            {
                // keeps looking in the upper part of the grid
                i_inf = i_med;
                a_inf = my_grid[i_inf];
            }
        }

       const auto new_cached = i_inf - middle_offset;
       my_cache.compare_exchange_strong(cached, new_cached);
       return new_cached;

    }

};