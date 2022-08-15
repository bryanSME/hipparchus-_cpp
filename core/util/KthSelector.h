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
//package org.hipparchus.util;

//import java.io.Serializable;
//import java.util.Arrays;
#include<vector>
#include<cmath>
#include<algorithm>
#include "PivotingStrategy.h"
//import org.hipparchus.exception.Null_Argument_Exception;


/**
 * A Simple K<sup>th</sup> selector implementation to pick up the
 * K<sup>th</sup> ordered element from a work array containing the
 * input numbers.
 */
class Kth_Selector  
{
private:

    /** Minimum selection size for insertion sort rather than selection. */
    static constexpr int MIN_SELECT_SIZE{ 15 };

    /** A {@link Pivoting_Strategy} used for pivoting.  */
    Pivoting_Strategy my_pivoting_strategy;
   
    /**
     * Partition an array slice around a pivot.Partitioning exchanges array
     * elements such that all elements smaller than pivot are before it and
     * all elements larger than pivot are after it.
     *
     * @param work work array
     * @param begin index of the first element of the slice of work array
     * @param end index after the last element of the slice of work array
     * @param pivot initial index of the pivot
     * @return index of the pivot after partition
     */
    int partition(std::vector<double>& work, const int& begin, const int& end, const int& pivot)
    {
        auto value = work[pivot];
        work[pivot] = work[begin];

        int i = begin + 1;
        int j = end - 1;
        while (i < j)
        {
            while (i < j && (work[j] - value) > 0)
            {
                --j;
            }
            while (i < j && (work[i] - value) < 0)
            {
                ++i;
            }

            if (i < j)
            {
                const double tmp = work[i];
                work[i++] = work[j];
                work[j--] = tmp;
            }
        }

        if (i >= end || (work[i] - value) > 0)
        {
            --i;
        }
        work[begin] = work[i];
        work[i] = value;
        return i;
    }

public:
    /**
     * Constructor with default {@link Pivoting_Strategy#MEDIAN_OF_3 median of 3}
     * pivoting strategy.
     */
    Kth_Selector() : my_pivoting_strategy{ /*Pivoting_Strategy::MEDIAN_OF_3 */} {};

    /**
     * Constructor with specified pivoting strategy
     *
     * @param pivoting_strategy pivoting strategy to use
     * @Null_Argument_Exception when pivoting_strategy is NULL
     */
    Kth_Selector(const Pivoting_Strategy& pivoting_strategy)
    {
        //Math_Utils::check_not_null(pivoting_strategy);
        my_pivoting_strategy = pivoting_strategy;
    }

    /** Get the pivoting strategy.
     * @return pivoting strategy
     */
    Pivoting_Strategy get_pivoting_strategy() const
    {
        return my_pivoting_strategy;
    }

    /**
     * Select K<sup>th</sup> value in the array.
     *
     * @param work work array to use to find out the K<sup>th</sup> value
     * @param pivots_heap cached pivots heap that can be used for efficient estimation
     * @param k the index whose value in the array is of interest
     * @return K<sup>th</sup> value
     */
    double select(std::vector<double>& work, std::vector<int>& pivots_heap, const int& k) 
    {
        int begin{};
        int end = work.size();
        int node{};
        while (end - begin > Kth_Selector::MIN_SELECT_SIZE)
        {
            int pivot;

            if (node < pivots_heap.size() && pivots_heap[node] >= 0) 
            {
                // the pivot has already been found in a previous call
                // and the array has already been partitioned around it
                pivot = pivots_heap[node];
            }
            else 
            {
                // select a pivot and partition work array around it
                pivot = partition(work, begin, end, my_pivoting_strategy.pivot_index(work, begin, end));
                if (node < pivots_heap.size()) 
                {
                    pivots_heap[node] = pivot;
                }
            }

            if (k == pivot) 
            {
                // the pivot was exactly the element we wanted
                return work[k];
            }
            else if (k < pivot) 
            {
                // the element is in the left partition
                end  = pivot;
                node = std::min(static_cast<std::size_t>(2) * node + 1, pivots_heap.size());
            }
            else 
            {
                // the element is in the right partition
                begin = pivot + 1;
                node  = std::min(static_cast<std::size_t>(2) * node + 2, pivots_heap.size());
            }
        }
        std::sort(work.begin(), work.end());
        return work[k];
    }
};