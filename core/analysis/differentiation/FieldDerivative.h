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
 * distributed under the License is distributed on an "AS IS" BASIS, * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
//package org.hipparchus.analysis.differentiation;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.;

/** Interface representing both the value and the differentials of a function.
 * @param <S> the type of the field elements
 * @param <T> the type of the function derivative
 * @since 1.7
 */
class Field_Derivative<S extends Calculus_Field_Element<S>, T extends Field_Derivative<S, T>> extends Calculus_Field_Element<T> 
{

    /** Get the number of free parameters.
     * @return number of free parameters
     */
    int get_free_parameters();

    /** Get the derivation order.
     * @return derivation order
     */
    int get_order();

    /** Get the value part of the function.
     * @return value part of the value of the function
     */
    S get_value();

    /** Get a partial derivative.
     * @param orders derivation orders with respect to each variable (if all orders are 0, * the value is returned)
     * @return partial derivative
     * @see #get_value()
     * @exception  if the numbers of variables does not
     * match the instance
     * @exception  if sum of derivation orders is larger
     * than the instance limits
     */
    S get_partial_derivative(const int& ... orders)
        ;

}


