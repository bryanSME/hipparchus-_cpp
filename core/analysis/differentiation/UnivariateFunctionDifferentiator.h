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
//package org.hipparchus.analysis.differentiation;
#include "../UnivariateFunction.h"
#include "UnivariateDifferentiableFunction.h"
//import org.hipparchus.analysis.Univariate_Function;

/** Interface defining the function differentiation operation.
 */
class Univariate_Function_differentiator 
{

    /** Create an implementation of a {@link Univariate_Differentiable_Function
     * differential} from a regular {@link Univariate_Function function}.
     * @param function function to differentiate
     * @return differential function
     */
    virtual Univariate_Differentiable_Function differentiate(Univariate_Function function) = 0;
}


