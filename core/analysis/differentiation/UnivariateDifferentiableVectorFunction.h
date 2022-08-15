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

//import org.hipparchus.analysis.Univariate_Vector_Function;
//import org.hipparchus.exception.;
#include "../UnivariateVectorFunction.h"
#include <vector>

/**
 * Extension of {@link Univariate_Vector_Function} representing a univariate differentiable vectorial function.
 *
 */
class Univariate_Differentiable_Vector_Function : public Univariate_Vector_Function
{

    /**
     * Compute the value for the function.
     * @param x the point for which the function value should be computed
     * @param <T> the type of the field elements
     * @return the value
     * @exception  if {@code x} does not
     * satisfy the function's constraints (argument out of bound, or unsupported
     * derivative order for example)
     */
    template<typename T>
    <T extends Derivative<T>> std::vector<T> value(T x);

};