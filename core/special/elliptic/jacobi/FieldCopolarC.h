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
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/** Copolar trio with pole at point c in Glaisher\xe2\x80\x99s Notation.
 * <p>
 * This is a container for the three subsidiary Jacobi elliptic functions
 * {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
 * </p>
 * @param <T> the type of the field elements
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldCopolar_C
{

    /** Value of the dc function. */
    private const T dc;

    /** Value of the nc function. */
    private const T nc;

    /** Value of the sc function. */
    private const T sc;

    /** Simple constructor.
     * @param trio_n copolar trio with pole at point n in Glaisher\xe2\x80\x99s Notation
     */
    FieldCopolar_C(const Field_Copolar_N<T> trio_n)
    {
        this.nc = trio_n.cn().reciprocal();
        this.sc = nc.multiply(trio_n.sn());
        this.dc = nc.multiply(trio_n.dn());
    }

    /** Get the value of the dc function.
     * @return dc(u|m)
     */
    public T dc()
    {
        return dc;
    }

    /** Get the value of the nc function.
     * @return nc(u|m)
     */
    public T nc()
    {
        return nc;
    }

    /** Get the value of the sc function.
     * @return sc(u|m)
     */
    public T sc()
    {
        return sc;
    }

};