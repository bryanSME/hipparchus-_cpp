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
//package org.hipparchus.stat;

//import java.util.Comparator;
#include "Frequency.hpp"
#include "../core/util/Comparable.h"

/**
 * Maintains a frequency distribution of long values.
 * <p>
 * Accepts byte, short, int, long primitive or Integer and long values.
 * <p>
 * Integer values (byte, short, int, long, Integer, long) are not
 * distinguished by type, i.e. {@code add_value(static_cast<long>(2)), * add_value(2), add_value(2L)} all have the same effect (similarly
 * for arguments to {@code get_count()} etc.).
 * <p>
 * NOTE: std::byte and short values will be implicitly converted to int values
 * by the compiler, thus there are no explicit overloaded methods for these
 * primitive types.
 * <p>
 * The values are ordered using the default (natural order), unless a
 * {@code Comparator} is supplied in the constructor.
 */
class Long_Frequency : Frequency<long>
{
public:
    /**
     * Default constructor.
     */
    Long_Frequency()
    {
        // This constructor is intentionally empty. Nothing special is needed here.
    }

    /**
     * Constructor allowing values Comparator to be specified.
     *
     * @param comparator Comparator used to order values
     */
    Long_Frequency(const Comparable<long>& comparator)
    {
        Frequency(comparator);
    }

    /**
     * Adds 1 to the frequency count for v.
     *
     * @param v the value to add.
     */
    void add_value(const int& v)
    {
        increment_value(static_cast<long>(v), 1);
    }

    /**
     * Increments the frequency count for v.
     *
     * @param v the value to add.
     * @param increment the amount by which the value should be incremented
     */
    void increment_value(const int& v, const long& increment)
    {
        increment_value(static_cast<long>(v), increment);
    }

    //-------------------------------------------------------------------------

    /**
     * Returns the number of values equal to v.
     *
     * @param v the value to lookup.
     * @return the frequency of v.
     */
    long get_count(const int& v)
    {
        return get_count(static_cast<long>(v));
    }

    /**
     * Returns the percentage of values that are equal to v
     * (as a proportion between 0 and 1).
     * <p>
     * Returns {@codeNAN} if no values have been added.
     *
     * @param v the value to lookup
     * @return the proportion of values equal to v
     */
    double get_pct(const int& v)
    {
        return get_pct(static_cast<long>(v));
    }

    //-----------------------------------------------------------------------------------------

    /**
     * Returns the cumulative frequency of values less than or equal to v.
     *
     * @param v the value to lookup.
     * @return the proportion of values equal to v
     */
    long get_cum_freq(const int& v)
    {
        return get_cum_freq(static_cast<long>(v));
    }

    //----------------------------------------------------------------------------------------------

    /**
     * Returns the cumulative percentage of values less than or equal to v
     * (as a proportion between 0 and 1).
     * <p>
     * Returns {@codeNAN} if no values have been added.
     *
     * @param v the value to lookup
     * @return the proportion of values less than or equal to v
     */
    double get_cum_pct(const int& v)
    {
        return get_cum_pct(static_cast<long>(v));
    }

};