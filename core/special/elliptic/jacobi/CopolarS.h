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
#include "CopolarN.h"

/** Copolar trio with pole at point s in Glaisher’s Notation.
 * <p>
 * This is a container for the three subsidiary Jacobi elliptic functions
 * {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
 * </p>
 * @since 2.0
 */
class Copolar_S 
{
private:
    /** Value of the cs function. */
    const double my_cs;

    /** Value of the dn function. */
    const double my_ds;

    /** Value of the ns function. */
    const double my_ns;

public:
    /** Simple constructor.
     * @param trio_n copolar trio with pole at point n in Glaisher’s Notation
     */
    Copolar_S(const Copolar_N& trio_n) 
        :
        my_ns{ 1.0 / trio_n.sn() },
        my_cs{ my_ns * trio_n.cn() },
        my_ds{ my_ns * trio_n.dn() }
    {};

    /** Get the value of the cs function.
     * @return cs(u|m)
     */
    double cs() const
    {
        return my_cs;
    }

    /** Get the value of the ds function.
     * @return ds(u|m)
     */
    double ds() const
    {
        return my_ds;
    }

    /** Get the value of the ns function.
     * @return ns(u|m)
     */
    double ns() const
    {
        return my_ns;
    }

}


