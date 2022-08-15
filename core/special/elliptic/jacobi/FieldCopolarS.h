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

/** Copolar trio with pole at point s in Glaisher\xe2\x80\x99s Notation.
 * <p>
 * This is a container for the three subsidiary Jacobi elliptic functions
 * {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
 * </p>
 * @param <T> the type of the field elements
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class FieldCopolar_S
{

    /** Value of the cs function. */
    private const T cs;

    /** Value of the dn function. */
    private const T ds;

    /** Value of the ns function. */
    private const T ns;

    /** Simple constructor.
     * @param trio_n copolar trio with pole at point n in Glaisher\xe2\x80\x99s Notation
     */
    FieldCopolar_S(const Field_Copolar_N<T> trio_n)
    {
        this.ns = trio_n.sn().reciprocal();
        this.cs = ns.multiply(trio_n.cn());
        this.ds = ns.multiply(trio_n.dn());
    }

    /** Get the value of the cs function.
     * @return cs(u|m)
     */
    public T cs()
    {
        return cs;
    }

    /** Get the value of the ds function.
     * @return ds(u|m)
     */
    public T ds()
    {
        return ds;
    }

    /** Get the value of the ns function.
     * @return ns(u|m)
     */
    public T ns()
    {
        return ns;
    }

};