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

//package org.hipparchus.ode.nonstiff;

//import org.hipparchus.Calculus_Field_Element;
#include <type_traits>
#include "../../core/CalculusFieldElement.h"

/** This interface represents an integrator  based on Butcher arrays.
 * @see Runge_Kutta_Field_Integrator
 * @see EmbeddedRunge_Kutta_Field_Integrator
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = NULLptr>
class FieldButcher_Array_Provider
{

    /** Get the time steps from Butcher array (without the first zero).
     * @return time steps from Butcher array (without the first zero
     */
    std::vector<T> get_c();

    /** Get the internal weights from Butcher array (without the first empty row).
     * @return internal weights from Butcher array (without the first empty row)
     */
    std::vector<std::vector<T>> get_a();

    /** Get the external weights for the high order method from Butcher array.
     * @return external weights for the high order method from Butcher array
     */
    std::vector<T> get_b();

};