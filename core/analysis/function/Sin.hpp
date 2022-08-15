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

#include <type_traits>
#include <cmath>
#include "../differentiation/UnivariateDifferentiableFunction.h"
#include "../differentiation/Derivative.h"

/**
 * Sine function.
 *
 */
class Sin : public Univariate_Differentiable_Function
{
public:
    /** {@inherit_doc} */
    //override
    double value(const double& x) const
    {
        return std::sin(x);
    }

    /** {@inherit_doc} */
    //override
    template<typename T, typename std::enable_if<std::is_base_of<Derivative<T>, T>::value>::type* = nullptr>
    T value(const T& x) const
    {
        return x.sin();
    }

};