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
//package org.hipparchus.special.elliptic.jacobi;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"
//import org.hipparchus.Calculus_Field_Element;

/** Values of {@link Field_Jacobi_Theta Jacobi theta} functions.
 * <p>
 * This is a container for the four Jacobi theta functions
 * \xce\xb8\xe2\x82\x81(z|\xcf\x84), \xce\xb8\xe2\x82\x82(z|\xcf\x84), \xce\xb8\xe2\x82\x83(z|\xcf\x84), and \xce\xb8\xe2\x82\x84(z|\xcf\x84).
 * </p>
 * @param <T> the type of the field elements
 * @see Field_Jacobi_Theta
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Theta
{

    /** Value of the \xce\xb8\xe2\x82\x81(z|\xcf\x84) function. */
    private const T theta1;

    /** Value of the \xce\xb8\xe2\x82\x82(z|\xcf\x84) function. */
    private const T theta2;

    /** Value of the \xce\xb8\xe2\x82\x83(z|\xcf\x84) function. */
    private const T theta3;

    /** Value of the \xce\xb8\xe2\x82\x84(z|\xcf\x84) function. */
    private const T theta4;

    /** Simple constructor.
     * @param theta1 value of the \xce\xb8\xe2\x82\x81(z|\xcf\x84) function
     * @param theta2 value of the \xce\xb8\xe2\x82\x82(z|\xcf\x84) function
     * @param theta3 value of the \xce\xb8\xe2\x82\x83(z|\xcf\x84) function
     * @param theta4 value of the \xce\xb8\xe2\x82\x84(z|\xcf\x84) function
     */
    Field_Theta(const T theta1, const T theta2, const T theta3, const T theta4)
    {
        this.theta1 = theta1;
        this.theta2 = theta2;
        this.theta3 = theta3;
        this.theta4 = theta4;
    }

    /** Get the value of the \xce\xb8\xe2\x82\x81(z|\xcf\x84) function.
     * @return \xce\xb8\xe2\x82\x81(z|\xcf\x84)
     */
    public T theta1()
    {
        return theta1;
    }

    /** Get the value of the \xce\xb8\xe2\x82\x82(z|\xcf\x84) function.
     * @return \xce\xb8\xe2\x82\x82(z|\xcf\x84)
     */
    public T theta2()
    {
        return theta2;
    }

    /** Get the value of the \xce\xb8\xe2\x82\x83(z|\xcf\x84) function.
     * @return \xce\xb8\xe2\x82\x83(z|\xcf\x84)
     */
    public T theta3()
    {
        return theta3;
    }

    /** Get the value of the \xce\xb8\xe2\x82\x84(z|\xcf\x84) function.
     * @return \xce\xb8\xe2\x82\x84(z|\xcf\x84)
     */
    public T theta4()
    {
        return theta4;
    }

};