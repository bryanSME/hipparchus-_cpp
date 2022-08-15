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

//import org.hipparchus.exception.Localized_Core_Formats;
#include <vector>

/**
 * Base interface implemented by all statistics.
 */
class Univariate_Statistic //: public Math_Arrays::Function 
{
    /**
     * Returns the result of evaluating the statistic over the input array.
     * <p>
     * The default implementation delegates to
     * <code>evaluate(std::vector<double>, int, int)</code> in the natural way.
     *
     * @param values input array
     * @return the value of the statistic applied to the input array
     * @  if values is NULL
     */
    //override
    double evaluate(std::vector<double> values)  
    {
        //Math_Utils::check_not_null(values, hipparchus::exception::Localized_Core_Formats_Type::INPUT_ARRAY);
        return evaluate(values, 0, values.size());
    }

    /**
     * Returns the result of evaluating the statistic over the specified entries
     * in the input array.
     *
     * @param values the input array
     * @param begin the index of the first element to include
     * @param length the number of elements to include
     * @return the value of the statistic applied to the included array entries
     * @ if values is NULL or the indices are invalid
     */
    //override
    virtual double evaluate(std::vector<double> values, int begin, int length) = 0;

    /**
     * Returns a copy of the statistic with the same internal state.
     *
     * @return a copy of the statistic
     */
    virtual Univariate_Statistic copy() = 0;
};


