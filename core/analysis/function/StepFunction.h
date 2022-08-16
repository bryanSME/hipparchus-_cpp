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
#include "../UnivariateFunction.h"
#include "../../util/MathArrays.h"

/**
 * <a href="http://en.wikipedia.org/wiki/Step_function">
 *  Step function</a>.
 *
 */
class Step_Function : Univariate_Function
{
private:
    /** Abscissae. */
    const std::vector<double> my_abscissa;
    /** Ordinates. */
    const std::vector<double> my_ordinate;

public:
    /**
     * Builds a step function from a list of arguments and the corresponding
     * values. Specifically, returns the function h(x) defined by <pre><code>
     * h(x) = y[0] for all x &lt; x[1]
     *        y[1] for x[1] &le; x &lt; x[2]
     *        ...
     *        y[y.size() - 1] for x &ge; x[x.size() - 1]
     * </code></pre>
     * The value of {@code x[0]} is ignored, but it must be strictly less than
     * {@code x[1]}.
     *
     * @param x Domain values where the function changes value.
     * @param y Values of the function.
     * @
     * if the {@code x} array is not sorted in strictly increasing order.
     * @Null_Argument_Exception if {@code x} or {@code y} are {@code NULL}.
     * @ if {@code x} or {@code y} are zero-length.
     * @ if {@code x} and {@code y} do not
     * have the same length.
     */
    Step_Function(const std::vector<double>& x, const std::vector<double>& y)
        :
        my_abscissa{ x },
        my_ordinate{ y }
    {
        if (x.size() == 0 || y.size() == 0)
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        Math_Arrays::check_equal_length(y, x);
        Math_Arrays::check_order(x);
    }

    /** {@inherit_doc} */
    //override
    double value(const double& x)
    {
        const int index = Arrays.binary_search(my_abscissa, x);

        if (index < -1)
        {
            // "x" is between "abscissa[-index-2]" and "abscissa[-index-1]".
            return my_ordinate[-index - 2];
        }
        if (index >= 0)
        {
            // "x" is exactly "abscissa[index]".
            return my_ordinate[index];
        }
        // Otherwise, "x" is smaller than the first value in "abscissa"
        // (hence the returned value should be "ordinate[0]").
        return my_ordinate[0];
    }
};