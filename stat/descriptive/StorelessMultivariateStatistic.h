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

#include <vector>

/**
 * Base interface implemented by storeless multivariate statistics.
 */
class Storeless_Multivariate_Statistic
{

    /**
     * Updates the internal state of the statistic to reflect
     * the addition of the value.
     *
     * @param d the value
     */
    virtual void increment(std::vector<double> d) = 0;

    /**
     * Returns the current value of the Statistic.
     * @return value of the statistic, <code>Double.NaN</code> if it
     * has been cleared or just instantiated.
     */
    virtual std::vector<double> get_result() = 0;

    /**
     * Returns the number of values that have been added.
     * @return the number of values.
     */
    virtual long get_n() = 0;

    /**
     * Clears the internal state of the statistic.
     */
    virtual void clear() = 0;

    /**
     * Returns the dimension of the statistic.
     * @return the dimension of the statistic
     */
    virtual int get_dimension() = 0;
};