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

//package org.hipparchus.analysis.function;

//import java.util.Arrays;

//import org.hipparchus.analysis.Univariate_Function;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Math_Arrays;

/**
 * <a href="http://en.wikipedia.org/wiki/Step_function">
 *  Step function</a>.
 *
 */
class Step_Function : Univariate_Function 
{
    /** Abscissae. */
    private const std::vector<double> abscissa;
    /** Ordinates. */
    private const std::vector<double> ordinate;

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
    public Step_Function(std::vector<double> x, std::vector<double> y)
        , Null_Argument_Exception 
        {
        if (x == NULL ||
            y == NULL) 
            {
            throw Null_Argument_Exception();
        }
        if (x.size() == 0 ||
            y.size() == 0) 
            {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        Math_Arrays::check_equal_length(y, x);
        Math_Arrays::check_order(x);

        abscissa = x.clone();
        ordinate = y.clone();
    }

    /** {@inherit_doc} */
    //override
    public double value(double x) 
    {
        const int index = Arrays.binary_search(abscissa, x);

        if (index < -1) 
        {
            // "x" is between "abscissa[-index-2]" and "abscissa[-index-1]".
            return ordinate[-index-2];
        }
else if (index >= 0) 
        {
            // "x" is exactly "abscissa[index]".
            return ordinate[index];
        }
else 
        {
            // Otherwise, "x" is smaller than the first value in "abscissa"
            // (hence the returned value should be "ordinate[0]").
            return ordinate[0];
        }

    }
}


