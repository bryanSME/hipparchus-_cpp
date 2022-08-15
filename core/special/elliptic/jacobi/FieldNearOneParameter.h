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

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sinh_Cosh;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/** Algorithm for computing the principal Jacobi functions for parameters slightly below one.
 * <p>
 * The algorithm for evaluating the functions is based on approximation
 * in terms of hyperbolic functions. It is given in Abramowitz and Stegun, * sections 16.15.
 * </p>
 * @param <T> the type of the field elements
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Near_One_Parameter extends Field_Jacobi_Elliptic<T>
{

    /** Complementary parameter of the Jacobi elliptic function. */
    private const T m1_fourth;

    /** Simple constructor.
     * @param m parameter of the Jacobi elliptic function (must be one or slightly below one here)
     */
    Field_Near_One_Parameter(const T m)
    {
        super(m);
        this.m1_fourth = m.get_field().get_one().subtract(m).multiply(0.25);
    }

    /** {@inherit_doc} */
    //override
    public Field_Copolar_N<T> values_n(const T u)
    {
        const Field_Sinh_Cosh<T> sch = std::sinh_cosh(u);
        const T                sech = sch.cosh().reciprocal();
        const T                t = sch.sinh().multiply(sech);
        const T                factor = sch.sinh().multiply(sch.cosh()).subtract(u).multiply(sech).multiply(m1_fourth);
        return Field_Copolar_N<>(t.add(factor.multiply(sech)),  // equation 16.15.1
            sech.subtract(factor.multiply(t)),        // equation 16.15.2
            sech.add(factor.multiply(t)));            // equation 16.15.3
    }

};