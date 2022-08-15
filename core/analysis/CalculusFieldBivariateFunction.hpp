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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

//import org.hipparchus.CalculusFieldElement;
#include <type_traits>
#include "../CalculusFieldElement.hpp"

/**
 * An interface representing a bivariate field function.
 * @param <T> type of the field elements
 * @since 1.5
 */
template<typename T, typename std::enable_if<std::is_base_of<CalculusFieldElement<T>, T>::value>::type* = nullptr>
class CalculusFieldBivariateFunction
{

    /**
     * Compute the value for the function.
     *
     * @param x Abscissa for which the function value should be computed.
     * @param y Ordinate for which the function value should be computed.
     * @return the value.
     */
    virtual T value(T x, T y) = 0;

};